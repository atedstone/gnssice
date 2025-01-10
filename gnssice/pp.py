""" Kinematic GNSS Post-Processing

HISTORY
Elements of this module were originally contained in gps.py.
Andrew Tedstone (andrew.tedstone@unifr.ch), July 2022

"""
from __future__ import annotations
import numpy as np
import statsmodels.api as sm
import math
import pandas as pd
from scipy.signal import filtfilt
from pandas.tseries.frequencies import to_offset
import matplotlib.pyplot as plt

from gnssice import gnss
from gnssice.gaussfiltcoef import gaussfiltcoef

# Length of year in fractional days
YEAR_LENGTH_DAYS = 365.26
# Length of Earth's semi-major axis in metres 
ELLPS_A = 6378137.0
# Length of Earth's semi-minor axis in metres 
ELLPS_B = 6356752.3142
# First numerical eccentricity
ELLPS_E2 = 1.0 - math.pow((ELLPS_B/ELLPS_A), 2)


def create_time_index(
    data : pd.DataFrame
    ) -> pd.DatetimeIndex:
    """
    Requires TRACK-originated dataframe with YY, DOY and Seconds columns.
    Can be used with parquet files.

    We need to subtract 1 from the supplied DOY:
    pd.to_datetime(2021)               -> datetime(2021, 1, 1)
    if DOY=1, then without subtraction -> datetime(2021, 1, 2)
    
    """
    ix = pd.to_datetime(data['YY'].astype(str)) + \
            pd.to_timedelta(data['DOY'].round().astype(int) - 1, unit='days') + \
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

    ## Transform ellipsoidal coordinates to cartesian
    # First convert to radians
    tmp_y = data['Latitude'] * (math.pi/180)
    tmp_x = data['Longitude'] * (math.pi/180)
    # Convert to XYZ
    xyz = ell2xyz(tmp_y, tmp_x, data['Height'])

    ## Transform absolute cartesian coordinates to local origin cartesian
    dX = xyz['x'] - x0
    dY = xyz['y'] - y0
    dZ = xyz['z'] - z0

    ## Transform to Local Geodetic: North-East-Up.
    # Matrix reshaping is to satisfy the requirements of the function.
    neu = ct2lg(dX, dY, dZ,
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
    """
    Apply temporal exclusions supplied by user, setting data in these bounds to NaN.
    
    The CSV file has the format: excl_start,excl_end,comment.
    The format of the dates/times in the file should be of the form
    yyyy-mm-ddThh:mm:ss (i.e. ISO standard, T is the separator between date and time).

    :param data: data to which to apply exclusions
    :param exclusion_file: path to the CSV file listing the exclusions.
    """
    
    excl = pd.read_csv(exclusion_file, parse_dates=['excl_start', 'excl_end'])
    for ix, row in excl.iterrows():
        print('Excluding %s - %s' %(row.excl_start, row.excl_end))
        data.loc[row.excl_start:row.excl_end, :] = np.nan

    return data
    

def filter_positions(
    data : pd.DataFrame,
    thresh_rms : float=50.0,
    thresh_h : float=10.0,
    thresh_N : int=3,
    thresh_NotF : int=0
    ) -> pd.DataFrame:
    """
    Filter (remove) bad positions based on their RMS, height standard deviation
    and if flagged as being interpolated by TRACK.

    :param thresh_rms: Threshold RMS value in mm to retain. (<=)
    :param thresh_h: Threshold height std. deviation in cm to retain. (<=)
    :param thresh_N: Threshold to apply to column N, retain above this value. (>)
    :param thresh_NotF : Threshold to apply to NotF, retain below this value (<=).

    Default parameters are based on Bartholomew/Sole/Tedstone (thresh_rm, thresh_h) 
    and Doyle 2014 (thresh_N, thresh_NotF).
    """
    
    # Filter by RMS
    data = data[data['RMS'] <= thresh_rms]

    # Filter by height std. dev.
    data = data[data['SigH'] <= thresh_h]

    # Filter by removing TRACK-interpolated values
    data = data[data['N'] > thresh_N]

    # Only retain epochs in which all ambiguities have been fixed (i.e. no unfixed ambiguities left).
    data = data[data['NotF'] <= thresh_NotF]

    return data


def calculate_displacement_trajectory(
    data: pd.DataFrame,
    verbose : bool=False
    ) -> tuple:
    """
    Compute 2-D velocity magnitude and co-variance - North and East.
    Export to disk.

    Computes on the contents of the whole dataframe.

    :param data: dataframe. Can pass in a subset.
    :param verbose: if True, print statsmodels linreg results.
    """

    # East
    # A.T. 2022-07-25 - not sure why the divide by length of year is needed.
    X = data['Fract_DOY'] / YEAR_LENGTH_DAYS
    X = sm.add_constant(X)
    y = data['East']
    rx = sm.OLS(y, X).fit()
    # Take the slope parameter
    vel_e = rx.params[0]

    # North
    X = data['Fract_DOY'] / YEAR_LENGTH_DAYS
    X = sm.add_constant(X)
    y = data['North']
    ry = sm.OLS(y, X).fit()
    # Take the slope parameter
    vel_n = ry.params[0]
    
    if verbose:
        print(rx.summary())
        print(ry.summary())

    return (vel_n, vel_e)


def create_rot_matrix(
    directions: tuple | list
    ) -> np.array:
    """
    Create rotation matrix, to be used by rotate_to_displacements.

    Designed to be called with the outputs of calculate_displacement_trajectory.

    :param directions: tuple with (north_direction, east_direction)
    """

    direc = math.atan2(*directions)
    
    R1 = np.array([
        [math.cos(-direc), -math.sin(-direc)],
        [math.sin(-direc), math.cos(-direc)]
        ])

    return R1


def rotate_to_displacements(
    east : pd.Series | np.array | list,
    north : pd.Series | np.array | list,
    R1 : np.array
    ) -> pd.DataFrame:
    """
    Rotate North-East data to along-track and across-track (x, y) data.

    :param east: East coordinates to rotate
    :param north: North coordinates to rotate
    :param R1: 2x2 rotation matrix (e.g. output of create_rot_matrix())
    """
    xy = np.dot(R1, np.array([east, north])).T
    return pd.DataFrame(xy, index=east.index, columns=('x', 'y'))


def regularise(
    data : pd.DataFrame,
    interval : str,
    add_flag : str = "interpolated"
    ) -> pd.DataFrame:
    """
    Sample the x,y,z input onto the desired frequency and fill gaps with linear 
    interpolation.

    :param data: DataFrame of X,Y,Z, indexed by time.
    :param interval: pandas Offset string.
    :param add_flag: None, or name of column to create containing 0 if original 
        data, 1 if interpolated.
    """
    if add_flag is not None:
        flag = pd.Series(0, index=data.index, name=add_flag, dtype=np.int32)
        flag[data.x.isna()] = np.nan
        flag = flag.resample(interval).asfreq()

    data = data.resample(interval).asfreq()
    
    data_iterp = data.filter(items=('x','y','z'), axis='columns').interpolate()
    data_iterp = pd.concat((data_iterp, flag), axis='columns')
    data_iterp[add_flag][data_iterp[add_flag].isna()] = 1

    return data_iterp


def remove_displacement_outliers(
    data,
    interval : str,
    iterations : int=2,
    mt : dict={'x':0.08, 'y':0.04, 'z':0.15},
    median_win : str='2H',
    sigma_mult : float=2
    ):
    """
    Remove outliers by median filtering.

    Usually run with at least two iterations - the first iteration removes large
    absolute differences, the subsequent iterations remove outliers via sigma.

    :param data: dataframe with x,y,z. Does not need to have even frequency.
    :param interval: sampling interval of data as pandas string
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
                        .rolling('24h', center=True).median() \
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

    return data


def smooth_displacement(
    data : pd.DataFrame,
    gauss_win_secs : int,
    gauss_win_z_mult : int=4
    ) -> pd.DataFrame:
    """ 
    Smooth x, y and z with a Gaussian forward-backward filter. Returns
    new dataframe with x,y,z columns.

    :param data: dataframe with x,y,z columns and regular time frequency
    :param gauss_win_secs: Size of the Gaussian filter window in seconds
    :param gauss_win_z_mult: Multiplier for size of z filter window
    """
    interval = data.index.freq.delta.seconds
    # Design filters
    gaus_coef = gaussfiltcoef((1/interval),(1/gauss_win_secs))
    gaus_coef_z = gaussfiltcoef((1/interval),(1/gauss_win_secs * gauss_win_z_mult))
    # Apply filters
    xs = pd.Series(filtfilt(gaus_coef, 1, data.x), index=data.index, name='x')
    ys = pd.Series(filtfilt(gaus_coef, 1, data.y), index=data.index, name='y')
    zs = pd.Series(filtfilt(gaus_coef_z, 1, data.z), index=data.index, name='z')
    return pd.concat((xs,ys,zs), axis='columns')


def position_by_regression(
    ser : pd.Series, 
    center: bool=False, 
    return_stats: bool=False
    ) -> float | pd.Series:
    """
    Estimate value at end of timestamped Series (default) or center (option) using 
    ordinary least squares regression through the Series.

    Only works with Series containing at least 10 values.

    Intended to be used in a resampling operation, e.g.:

    >>> df.resample('10D').apply(position_by_regression)

    >>> df.resample('10D').apply(position_by_regression, return_stats=True)

    :param center: default is to compute value for last index of Series. If True, compute
        value for center timestamp instead.
    :param return_stats: if True, return pd.Series containing the computed value, the
    OLS m and c parameters, the r-squared value.

    """
    # OLS regression does not work with NaNs.
    ser = ser.dropna()
    if len(ser) < 10:
        return np.nan

    # Fit model, using Julian calendar as exog
    juld = ser.index.to_julian_date()
    X = sm.add_constant(juld)
    y = ser.values
    m = sm.OLS(y, X)
    f = m.fit()

    # Compute value
    if center:
        t = juld[0] + ((juld[-1] - juld[0]) / 2)
    else:
        t = juld[-1]
    v = f.params[1] * t + f.params[0]

    if return_stats:
        return pd.Series([v, f.params[1], f.params[0], f.rsquared],
            index=['p', 'param_m', 'param_c', 'r2'])
    else:
        return v


def detrend_z(
    data : pd.DataFrame,
    slope : float=None
    ) -> (pd.Series, float):
    """
    Detrend z (height data) as a function of x. Optionally, can be used with a 
    pre-determined slope. This is useful if a slope has been calculated previously
    from an earlier block of data - have the option to use this earlier one instead.

    :param data: dataframe containing at least x,z.
    :param slope: float of slope.
    :returns: (detrended z, slope)
    """
    if slope is None:
        X = data['x']
        X = sm.add_constant(X)
        y = data['z']
        rz = sm.OLS(y, X).fit()
        slope = rz.params[2]

    z_detr = data['z'] - data['x'] * slope
    return (z_detr, slope)


def calculate_daily_velocities(
    x : pd.Series,
    tz : str=None,
    flag_iterp : pd.Series=None
    ) -> pd.Series:
    """
    Calculate 24-h along-track velocity using average position in first hour of each day.
    Units metres/year.

    if tz is provided, first localise the (assumed UTC) measurements to local
    timezone before calculating the velocities.

    :param x: Series of along-track displacement from which to calculate velocity.
    :param tz: A 'tz-database' time zone string, e.g. America/Nuuk.
    """
    if tz is not None and tz != '':
        x = x.tz_localize(tz)

    v_24h = x.resample('24h').apply(lambda x: x[x.index.hour == 0].mean())
    v_24h = (v_24h.shift(1) - v_24h) * YEAR_LENGTH_DAYS

    if flag_iterp is not None:
        if tz is not None and tz != '':
            flag_iterp = flag_iterp.tz_localize(tz)
            f = flag_iterp.resample('D').first()
            return pd.concat([v_24h, f], axis='columns')

    return v_24h
   

def v_mult(
    window : str
    ) -> float:
    """
    Given a pandas window/offset string, Calculate the multipler needed to 
    convert displacement to m/yr.

    Info on to_offset:
    https://stackoverflow.com/questions/40223470/how-do-i-get-at-the-pandas-offsets-object-given-an-offset-string
    """
    offset = to_offset(window)
    year = pd.Timedelta(days=YEAR_LENGTH_DAYS)
    multiplier = year / offset
    return multiplier


def calculate_short_velocities(
    x : pd.Series,
    window : str
    ) -> pd.Series:
    """
    Calculate 'instantaneous' along-track velocity across window. Units metres/year.

    :param data: Series of along-track displacement to convert to velocity.
    :param window: pandas Offset str
    """
    multiplier = v_mult(window)
    vs = (x.shift(freq=window) - x) * multiplier
    return vs 


def ell2xyz(lat,lon,h,a=ELLPS_A,e2=ELLPS_E2):
    """Convert lat lon height to local north-east-up.
    
    lat, lon and h must be numpy 1-d arrays.
    
    If a and e2 are None then their values will be obtained from the 
    class variables self.a and self.e2.
    
    Based on the matlab exchange function...:
    ELL2XYZ  Converts ellipsoidal coordinates to cartesian.
    Vectorized.
     Version: 2011-02-19
     Useage:  [x,y,z]=ell2xyz(lat,lon,h,a,e2)
              [x,y,z]=ell2xyz(lat,lon,h)
     Input:   lat - vector of ellipsoidal latitudes (radians)
              lon - vector of ellipsoidal E longitudes (radians)
              h   - vector of ellipsoidal heights (m)
              a   - ref. ellipsoid major semi-axis (m); default GRS80
              e2  - ref. ellipsoid eccentricity squared; default GRS80
     Output:  x \
              y  > vectors of cartesian coordinates in CT system (m)
              z /
    
     Original copyright (c) 2011, Michael R. Craymer Email: mike@craymer.com
                
    """
        
    v = a / np.sqrt(1 - e2 * np.sin(lat) * np.sin(lat))
    x = (v + h) * np.cos(lat) * np.cos(lon)
    y = (v + h) * np.cos(lat) * np.sin(lon)
    z = (v * (1 - e2) + h) * np.sin(lat)   
            
    toret = {}
    toret['x'] = x
    toret['y'] = y
    toret['z'] = z
    return toret
    

def ct2lg(dX,dY,dZ,lat,lon):
    """Converts CT coordinate differences to local geodetic.
    
    All inputs must be numpy arrays.
    
    Local origin at lat,lon,h. If lat,lon are vectors, dx,dy,dz
    are referenced to origin at lat,lon of same index. If
    astronomic lat,lon input, output is in local astronomic
    system. Vectorized in both dx,dy,dz and lat,lon. See also
    LG2CT.
    Version: 2011-02-19
    Useage:  [dx,dy,dz]=ct2lg(dX,dY,dZ,lat,lon)
    Input:   dX  - vector of X coordinate differences in CT
     dY  - vector of Y coordinate differences in CT
     dZ  - vector of Z coordinate differences in CT
     lat - lat(s) of local system origin (rad); may be vector
     lon - lon(s) of local system origin (rad); may be vector
    Output:  dx  - vector of x coordinates in local system (north)
     dy  - vector of y coordinates in local system (east)
     dz  - vector of z coordinates in local system (ht)
            
    Ported to Python from original Matlab exchange version...
    Original copyright (c) 2011, Michael R. Craymer 
    All rights reserved.
    Email: mike@craymer.com 

    """

    n = len(dX)
    if isinstance(lat, float):
        lat = np.ones((n,1)) * lat
        lon = np.ones((n,1)) * lon
    R = np.zeros((3,3,n))

    if isinstance(dX, pd.Series):
        dX = dX.to_numpy()
        dY = dY.to_numpy()
        dZ = dZ.to_numpy()
    
    R[0,0,:] = -np.sin(lat.T) * np.cos(lon.T)
    R[0,1,:] = -np.sin(lat.T) * np.sin(lon.T)
    R[0,2,:] = np.cos(lat.T)
    
    R[1,0,:] = -np.sin(lon.T)
    R[1,1,:] = np.cos(lon.T)
    R[1,2,:] = np.zeros((1,n))
    
    R[2,0,:] = np.cos(lat.T) * np.cos(lon.T)
    R[2,1,:] = np.cos(lat.T) * np.sin(lon.T)
    R[2,2,:] = np.sin(lat.T)
    
    RR = np.reshape(R[0,:,:],(3,n))
    dx = np.sum(RR.T * np.stack((dX,dY,dZ),axis=1),axis=1)
    RR = np.reshape(R[1,:,:],(3,n))
    dy = np.sum(RR.T * np.stack((dX,dY,dZ),axis=1),axis=1)
    RR = np.reshape(R[2,:,:],(3,n))
    dz = np.sum(RR.T * np.stack((dX,dY,dZ),axis=1),axis=1)
    
    toret = {}
    toret['x'] = dx
    toret['y'] = dy
    toret['z'] = dz
    return toret
    
    