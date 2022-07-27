from __future__ import annotations
import numpy as np
import dask
import statsmodels.api as sm
import math
import pandas as pd

import matplotlib.pyplot as plt

import gps

import pdb

## !! To-do - move over ell2xyz etc here then remove Kinematic/Postprocess refs

def create_time_index(
    data
    ) -> pd.DataFrame:
    """
    Requires TRACK-originated dataframe with YY, DOY and Seconds columns.
    Can be used with parquet files.

    Could move this to along-side/within read_track_file?
    """
    ix = pd.to_datetime(data['YY'].astype(str)) + \
            pd.to_timedelta(data['DOY'].round().astype(int), unit='days') + \
            pd.to_timedelta(data['Seconds'].round().astype(int), unit='sec')
    return ix


def calculate_local_neu(
    data : pd.DataFrame,
    x0 : float,
    y0 : float, 
    z0 : float, 
    lat0 : float, 
    lon0 : float
    ) -> pd.DataFrame:
    """
    Given local NEU origin coordinates (x, y, z), transform all lat/lon/height 
    into this coordinate reference system.

    Usually, origin coordinates have already been calculated by calculate_local_origin.py.

    :param data:
    :param x0: ECEF XYZ X origin
    :param y0: ECEF XYZ Y origin
    :param z0: ECEF XYZ Z origin
    :param lat0: WGS84 latitude origin (corresponding to y0)
    :param lon0: WGS84 longitude origin (corresponding to x0)
    """

    pp = gps.PostProcess()

    ## Transform ellipsoidal coordinates to cartesian
    # First convert to radians
    tmp_y = data['Latitude'] * (math.pi/180)
    tmp_x = data['Longitude'] * (math.pi/180)
    # Convert to XYZ
    xyz = pp.ell2xyz(tmp_y, tmp_x, data['Height'])

    ## Transform absolute cartesian coordinates to local origin cartesian
    dX = xyz['x'] - x0
    dY = xyz['y'] - y0
    dZ = xyz['z'] - z0

    ## Transform to Local Geodetic: North-East-Up.
    # Matrix reshaping is to satisfy the requirements of the function.
    neu = pp.ct2lg(dX, dY, dZ,
                   lat0 * (math.pi / 180),
                   lon0 * (math.pi / 180))
    
    # Format to dataframe
    neu = pd.DataFrame(neu, index=data.index)
    # Keep z named as z.
    neu = neu.rename({'y':'North', 'x':'East'}, axis='columns')
    return neu
    

def apply_pole_corrections():
    # not needed yet
    pass
    

def apply_exclusions(
    data : pd.DataFrame,
    exclusion_file : str
    ) -> pd.DataFrame:
    
    excl = pd.read_csv(exclusion_file, parse_dates=['excl_start', 'excl_end'])
    for ix, row in excl.iterrows():
        print('Excluding %s - %s' %(row.excl_start, row.excl_end))
        data.at[row.excl_start:row.excl_end, :] = np.nan

    return data
    

def filter_positions(
    data : dask.DataFrame,
    thresh_rms : float=50.0,
    thresh_h : float=9.0
    ) -> dask.DataFrame:
    """
    :param thresh_rms: Threshold RMS value in mm to retain.
    :param thresh_h: Threshold height std. deviation in cm to retain.
    """
    # Am I streaming the filtered data back to disk?
    
    # Filter by RMS
    data = data[data['RMS'] <= thresh_rms]

    # Filter by height std. dev.
    data = data[data['SigH'] <= thresh_h]

    # Filter by removing TRACK-interpolated values
    data = data[data['N'] == 0]

    return data


def calculate_displacement_trajectory(
    data: pd.DataFrame | dask.DataFrame
    ) -> tuple:
    """
    Compute 2-D velocity magnitude and co-variance - North and East.
    Export to disk.

    Computes on the contents of the whole dataframe.

    :param data: dataframe. Can pass in a subset.
    """

    # East
    # A.T. 2022-07-25 - not sure why the divide by length of year is needed.
    X = data['Fract_DOY'] / 365.25
    X = sm.add_constant(X)
    y = data['East']
    r = sm.OLS(y, X).fit()
    # Take the slope parameter
    print(r.summary())
    vel_e = r.params[0]

    # North
    X = data['Fract_DOY'] / 365.25
    X = sm.add_constant(X)
    y = data['North']
    r = sm.OLS(y, X).fit()
    # Take the slope parameter
    print(r.summary())
    vel_n = r.params[0]

    return (vel_n, vel_e)


def create_rot_matrix(
    directions: tuple | list
    ) -> np.array:
    """
    Create rotation matrix, to be used by rotate_to_displacements.
    """

    direc = math.atan2(*directions)
    print(direc)
    
    #%rotate --> xy(:,1)=along track, xy(:,2)=across track.
    R1 = np.array([
        [math.cos(-direc), -math.sin(-direc)],
        [math.sin(-direc), math.cos(-direc)]
        ])

    return R1


def rotate_to_displacements(
    east : pd.Series | np.array | list,
    north : pd.Series | np.array | list,
    R1 : np.array
    ):

    # Note the transposes below!!
    # (R1*[smap(:,19) smap(:,18)]')';
    xy = np.dot(R1, np.array([east, north])).T
    return pd.DataFrame(xy, index=east.index, columns=('x', 'y'))


def regularise(
    data : pd.DataFrame,
    interval : str,
    add_flag : str = "interpolated"
    ) -> pd.DataFrame:
    """
    Sample the x,y,z input onto the desired frequency and fill gaps with linear interpolation.

    :param data: DataFrame of X,Y,Z, indexed by time.
    """
    if add_flag is not None:
        data[add_flag] = False

    data = data.resample(interval).asfreq()

    # matlab script checked for identical time values ... should this ever happen though?
    
    data_iterp = data.filter(items=('x','y','z'), axis='columns').interpolate()
    data = pd.concat((data_iterp, data[add_flag]), axis='columns')
    #data[add_flag] = data[add_flag].where(data[add_flag] != False, True)

    return data


def remove_displacement_outliers(
    data,
    interval : str='10s',
    iterations : int=2,
    mt : dict={'x':0.08, 'y':0.04, 'z':0.15},
    median_win : str='2H',
    sigma_mult : float=2
    ):
    """
    Remove outliers by median filtering.

    Usually run with at least two iterations - the first iteration removes large
    absolute differences, the subsequent iterations remove outliers via sigma.

    :param interval: interval of data as pandas string
    :param iterations: number of iterations to run
    :param mt: dict of absolute difference thresholds for each dimension
    :param median_win: pandas window string for nth stage filtering
    :param sigma_mult: std.dev. multiplier for nth stage filtering
    """
    for n in range(0,iterations):
        print('Median Filter: Iteration %s' %n)
        if n == 0:
            smoothed = data.filter(items=('x','y','z'), axis=1) \
                        .resample('15Min').first() \
                        .rolling('24H', center=True).median() \
                        .resample(interval).asfreq() \
                        .interpolate()
            diffs = (data - smoothed).apply(np.abs)
            data = data[~(
                (diffs['x'] >= mt['x']) | 
                (diffs['y'] >= mt['y']) | 
                (diffs['z'] >= mt['z'])
                )]
        else:
            median = data.filter(items=('x','y','z'), axis=1) \
                        .rolling(median_win, center=True).median() 
            diffs = (data - median).apply(np.abs) 
            sigma = median.rolling(median_win, center=True).std()
            # Apply user-prescribed tolerance
            sigma_ = sigma * sigma_mult
            # Remove differences from the median that are larger than the 
            # median's sigma.
            data = data[~(
                (diffs['x'] >= sigma_['x']) | 
                (diffs['y'] >= sigma_['y']) | 
                (diffs['z'] >= sigma_['z'])
                )]

    # Fill gaps
    #data = regularise(data, interval)
    return data


def smooth_displacement(
    data
    ):
    pass

    # Smooth with Gaussian low-pass filter

    # Remove interpolated values??


def calculate_daily_velocities():
    pass


def calculate_short_velocities(
    window : int
    ):
    pass