from __future__ import annotations
import numpy as np
import dask
import statsmodels.api as sm

def calculate_local_neu():
    dX = xyz['x'] - np.median(xyz['x'][0:100]) 
    dY = xyz['y'] - np.median(xyz['y'][0:100])
    dZ = xyz['z'] - np.median(xyz['z'][0:100])
    # The calculated neu coordinates reference lat-lon origin point.
    # All the nasty matrix reshaping is to satisfy the requirements of
    # the function!
    neu = self.ct2lg(np.array([dX]).T,np.array([dY]).T,np.array([dZ]).T,
                   np.array([smap_all[0,3] * (math.pi / 180)]),
                   np.array([smap_all[0,4] * (math.pi / 180)]))

    # make displacements absolute - apply this directly to NEU values

    # Save back to columns named North, East, Up


def apply_pole_corrections():
    # not needed yet
    pass
    

def apply_exclusions():
    # not needed yet
    pass
    # This needs some TOML file trickery or equivalent


def filter_positions(
    data : dask.DataFrame,
    thresh_rms : float=50.0,
    thresh_h : float=9.0
    ) -> dask.DataFrame:
    """
    :param thresh_rms: Threshold RMS value in mm to retain.
    :param thresh_h: Threshold height std. deviation in cm to retain.
    """
    # Am I streaming the corrected data back to disk?
    
    # Filter by RMS
    data = data[data['RMS'] <= thresh_rms]

    # Filter by height std. dev.
    data = data[data['SigH'] <= thresh_h]

    # Filter by removing interpolated values
    data = data[data['#'] == 0]

    return data


def calculate_displacement_trajectory(
    data: dask.DataFrame
    ) -> None:
    """
    Compute 2-D velocity magnitude and co-variance - North and East.
    Export to disk.

    Computes on the contents of the whole dataframe.

    :param data: dataframe. Can pass in a subset.
    """

    # East
    # A.T. 2022-07-25 - not sure why the divide by length of year is needed.
    X = data['Fract DOY'] / 365.25
    X = sm.add_constant(X)
    y = data['East']
    r = sm.OLS(y, X).fit()
    # Take the slope parameter
    print(r)
    vel_e = r.params[0]

    # North
    X = data['Fract DOY'] / 365.25
    X = sm.add_constant(X)
    y = data['North']
    r = sm.OLS(y, X).fit()
    # Take the slope parameter
    print(r)
    vel_n = r.params[0]

    return (vel_n, vel_e)


def create_rot_matrix(
    directions: tuple | list
    ) -> np.array:
    direc = math.atan2(directions)
    
    #%rotate --> xy(:,1)=along track, xy(:,2)=across track.
    R1 = np.array([
        [math.cos(-direc), -math.sin(-direc)]
        [math.sin(-direc), math.cos(-direc)]
        ])

def rotate_to_displacements(
    east : np.array | list,
    north : np.array | list,
    directions : None=None | tuple | list,
    ):

    if directions is None:
        R1 = create_rot_matrix(directions)
    else:
        R1 = directions    

    # Note the transposes below!!
    # (R1*[smap(:,19) smap(:,18)]')';
    xy = R1 * np.array([east, north])
    return xy
    # After this, rename columns in df to x,y


def add_time_index(
    data
    ) -> pd.DataFrame:
    pass


def smooth_displacements(
    data,
    interval : str,
    iterations : int=2,
    median_win : str='2H',
    median_thresholds : dict={'x':0.08, 'y':0.04, 'z':0.15}
    ):

    data = data.resample(interval).asfreq()

    # matlab script checked for identical time values ... should this ever happen though?

    # This should apply across all of NEU?
    # Flag interpolated values ??
    data = data.interpolate()

    # Remove outliers
    # Median filter
    for n in range(0,iterations)
        print('Median Filter: Iteration %s' %n)
        smoothed = data.rolling(median_win).median()
        for dim, t in median_thresholds.values():
            diffs = np.abs(data[dim] - smoothed[dim])
            if n == 0:
                data[dim] = data[dim][diffs <= t]
            else:
                sigma = data.rolling(median_win).std()
                data[dim] = data[dim][diffs <= (sigma * sigma_mult)]
                # Actually, in MATLAB code, all data (X,Y,Z) is culled irrespectively of dimension which 'failed'.




    # Smooth with Gaussian low-pass filter

    # Remove interpolated values??


def calculate_daily_velocities():
    pass


def calculate_short_velocities(
    window : int
    )
    pass