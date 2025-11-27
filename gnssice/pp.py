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
import collections.abc
import pdb
from copy import deepcopy

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

REGEX_GEOD_BALE_FILE = r'[a-z0-9]{4}\_[a-z0-9]{4}\_[0-9]{4}\_[0-9]{3}\_[0-9]{3}_GEOD.parquet'
REGEX_L1_FILE = r'[a-z0-9]{4}\_[0-9]{4}\_[0-9]{3}\_[0-9]{4}\_[0-9]{3}_geod.parquet'


def update_legacy_geod_col_names(df):
    return df.rename(columns={
            'Latitude':'Latitude_deg', 
            'Longitude':'Longitude_deg', 
            'Height':'Height_m', 
            'SigN':'SigN_cm', 
            'SigE':'SigE_cm', 
            'SigH':'SigH_cm', 
            'RMS':'RMS_mm', 
            'Atm':'Atm_mm', 
            'plus_minus':'plus_minus_mm' 
            })   


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
    # NEU file
    neu_cols = ['YY', 'MM', 'DD', 'HR', 'MIN', 'Sec']
    if set(neu_cols).issubset(data.columns):
        _tmp = data.filter(items=neu_cols, axis='columns')
        _tmp = _tmp.rename(columns={'YY':'year', 'MM':'month', 'DD':'day', 'HR':'hour', 'MIN':'minute', 'Sec':'second'})
        ix = pd.to_datetime(_tmp)

    # GEOD file
    elif set(['YY', 'DOY', 'Seconds']).issubset(data.columns):
        ix = pd.to_datetime(data['YY'].astype(str)) + \
                pd.to_timedelta(data['DOY'].round().astype(int) - 1, unit='days') + \
                pd.to_timedelta(data['Seconds'].round().astype(int), unit='sec')

    else:
        raise ValueError('Could not find columns required to create time index.')

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
    tmp_y = data['Latitude_deg'] * (math.pi/180)
    tmp_x = data['Longitude_deg'] * (math.pi/180)
    # Convert to XYZ
    xyz = ell2xyz(tmp_y, tmp_x, data['Height_m'])

    ## Transform absolute cartesian coordinates to local origin cartesian
    dX = xyz['x_m'] - x0
    dY = xyz['y_m'] - y0
    dZ = xyz['z_m'] - z0

    ## Transform to Local Geodetic: North-East-Up.
    # Matrix reshaping is to satisfy the requirements of the function.
    neu = ct2lg(dX, dY, dZ,
                lat0 * (math.pi / 180),
                lon0 * (math.pi / 180))
    
    # Format to dataframe
    neu = pd.DataFrame(neu, index=data.index)
    # Keep z named as z.
    neu = neu.rename({'y_m':'North_m', 'x_m':'East_m'}, axis='columns')
    return neu
    

def correct_pole_changes(
    data : pd.DataFrame,
    correction_file : str,
    ) -> pd.DataFrame:
    """     
    Correct the XYZ displacement record for pole changes resulting from maintenance.

    The CSV file has the columns:
    - maintenance_start : timestamps in ISO format
    - maintenance_end : timestamps in ISO format
    - dx_m : auto or float in metres
    - dy_m : auto or float in metres
    - dz_m : auto or float in metres
    - comment : text note, not used in this analysis.

    For any change listed as 'auto', this function will calculate the difference in 
    the position of that dimension between the hour before the maintenance_start and the
    hour after the maintenance_end.

    :param data: data to which to apply pole changes
    :param correction_file: path to the CSV file containing the corrections.
    """

    def auto_difference(corr_info, col, window='1h'):
        win_start = row.maintenance_start - pd.Timedelta(window)
        win_end = row.maintenance_end + pd.Timedelta(window)
        pos_before = d.loc[win_start:row.maintenance_start, col].mean()
        pos_after = d.loc[row.maintenance_end:win_end, col].mean()
        diff = pos_after - pos_before
        return diff

    corrs = pd.read_csv(correction_file, parse_dates=['maintenance_start', 'maintenance_end'])

    d = deepcopy(data)

    for ix, row in corrs.iterrows():
        print('Correcting %s - %s' %(row.maintenance_start, row.maintenance_end))
        dx = row.dx_m
        if dx == 'auto':
            dx = auto_difference(row, 'x_m')
        d.loc[row.maintenance_end:, 'x_m'] -= float(dx)
        
        dy = row.dy_m
        if dy == 'auto':
            dy = auto_difference(row, 'y_m')
        d.loc[row.maintenance_end:, 'y_m'] -= float(dy)
        
        dz = row.dz_m
        if dz == 'auto':
            dz = auto_difference(row, 'z_m')
        d.loc[row.maintenance_end:, 'z_m'] -= float(dz)

        print(dx, dy, dz)
        # Remove any data acquired during the maintenance window
        d.loc[row.maintenance_start:row.maintenance_end] = np.nan

    return d

            

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
    data['exclude'] = False
    for ix, row in excl.iterrows():
        print('Excluding %s - %s' %(row.excl_start, row.excl_end))
        data.loc[row.excl_start:row.excl_end, 'exclude'] = True

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
    data = data[data['RMS_mm'] <= thresh_rms]

    # Filter by height std. dev.
    data = data[data['SigH_cm'] <= thresh_h]

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
    y = data['East_m']
    rx = sm.OLS(y, X).fit()
    # Take the slope parameter
    vel_e = rx.params[0]

    # North
    X = data['Fract_DOY'] / YEAR_LENGTH_DAYS
    X = sm.add_constant(X)
    y = data['North_m']
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
    return pd.DataFrame(xy, index=east.index, columns=('x_m', 'y_m'))


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
        flag[data.x_m.isna()] = np.nan
        flag = flag.resample(interval).asfreq()

    data = data.resample(interval).asfreq()
    
    data_iterp = data.filter(items=('x_m','y_m','z_m'), axis='columns').interpolate()
    data_iterp = pd.concat((data_iterp, flag), axis='columns')
    data_iterp.loc[data_iterp[add_flag].isna(), add_flag] = 1

    return data_iterp


def remove_displacement_outliers(
    data,
    interval : str,
    iterations : int=2,
    mt : dict={'x_m':0.08, 'y_m':0.04, 'z_m':0.15},
    median_win : str='2h',
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
            smoothed = data.filter(items=('x_m','y_m','z_m'), axis=1) \
                        .resample('15Min').first() \
                        .rolling('24h', center=True).median() \
                        .resample(interval).asfreq() \
                        .interpolate()
            diffs = (data - smoothed).apply(np.abs)
            # Only look at differences for uninterpolated data points
            diffs = diffs[diffs.index.isin(data.index)]
            #pdb.set_trace()
            data = data[~(
                (diffs['x_m'] >= mt['x_m']) | 
                (diffs['y_m'] >= mt['y_m']) | 
                (diffs['z_m'] >= mt['z_m'])
                )]
        else:
            median = data.filter(items=('x_m','y_m','z_m'), axis=1) \
                        .rolling(median_win, center=True).median() 
            diffs = (data - median).apply(np.abs) 
            sigma = median.rolling(median_win, center=True).std()
            # Apply user-prescribed tolerance
            sigma_ = sigma * sigma_mult
            # Remove differences from the median that are larger than the 
            # median's sigma.
            data = data[~(
                (diffs['x_m'] >= sigma_['x_m']) | 
                (diffs['y_m'] >= sigma_['y_m']) | 
                (diffs['z_m'] >= sigma_['z_m'])
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
    # Make sure data have seconds frequency
    assert data.index.freq.freqstr[-1] == 's'
    # Get seconds-based sampling frequency
    interval = data.index.freq.n
    # Design filters
    gaus_coef = gaussfiltcoef((1/interval),(1/gauss_win_secs))
    gaus_coef_z = gaussfiltcoef((1/interval),(1/gauss_win_secs * gauss_win_z_mult))
    # Apply filters
    xs = pd.Series(filtfilt(gaus_coef, 1, data.x_m), index=data.index, name='x_m')
    ys = pd.Series(filtfilt(gaus_coef, 1, data.y_m), index=data.index, name='y_m')
    zs = pd.Series(filtfilt(gaus_coef_z, 1, data.z_m), index=data.index, name='z_m')
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
        X = data['x_m']
        X = sm.add_constant(X)
        y = data['z_m']
        rz = sm.OLS(y, X).fit()
        slope = rz.params[2]

    z_detr = data['z_m'] - data['x_m'] * slope
    return (z_detr, slope)


def resample_with_window(ts, freq, half_window, agg='mean',
                                    offset=None, dropna=True):
    """
    Resample an irregular time series to a fixed grid (freq, optional offset),
    aggregating only points within +/- half_window of each grid boundary.

    Method:
    1) Mask points that are within +/- half_window (modulo freq) of the grid.
    2) Shift the masked timestamps forward by +half_window.
    3) Resample on the shifted timestamps (this puts the +/-half window inside each non-overlapping bin).
    4) Optionally drop NaN bins.

    Notes:
    - Works for fixed-length frequencies (timedeltas): seconds, minutes, hours, days, etc.
      It does NOT support variable-length freq like 'M', 'Q', 'A' without extra logic.
    - Requires half_window <= freq / 2 (logical constraint for symmetric half-windows).

    
    Generated by ChatGPT following (multiple!) prompts by AT, Nov 2025.
    """
    if isinstance(half_window, str):
        half_window = pd.Timedelta(half_window)
    if isinstance(freq, str):
        # try to convert freq to Timedelta for validation and arithmetic
        try:
            freq_td = pd.Timedelta(freq)
        except Exception as e:
            raise ValueError("freq must be a fixed-length offset convertible to Timedelta (e.g. '1D','1H').") from e
    else:
        freq_td = pd.Timedelta(freq)

    if half_window > freq_td / 2:
        raise ValueError("half_window must be <= freq/2 to avoid ambiguous windows.")

    ts = ts.sort_index()
    if ts.empty:
        # return an empty resampled index
        start = pd.Timestamp.now().floor(freq)  # arbitrary empty
        return pd.Series(dtype=float)

    # Build target grid (same anchor resample would use)
    start = ts.index.min().floor(freq)
    end = ts.index.max().ceil(freq)
    targets = pd.date_range(start, end, freq=freq)
    if offset is not None:
        offset_td = pd.Timedelta(offset)
        targets = targets + offset_td
    else:
        offset_td = pd.Timedelta(0)

    # Convert to int64 ns for vectorised modulo arithmetic
    ix_ns = ts.index.view('int64')
    first_target_ns = targets[0].value
    freq_ns = freq_td.value
    half_ns = half_window.value

    # relative position modulo freq
    rel = (ix_ns - first_target_ns) % freq_ns

    # mask points within +/- half_window of a boundary (either near 0 or near freq_ns)
    mask = (rel <= half_ns) | (rel >= (freq_ns - half_ns))
    ts_masked = ts.iloc[mask]

    if ts_masked.empty:
        # no points fall in any half-window
        empty_idx = pd.date_range(start, end, freq=freq)
        if offset is not None:
            empty_idx = empty_idx + offset_td
        if isinstance(ts, pd.DataFrame):
            return pd.DataFrame(index=empty_idx, columns=ts.columns).astype(float)
        else:
            return pd.Series(index=empty_idx, dtype=float)

    # shift masked timestamps forward by +half_window so window [t-half, t+half) -> [t, t+2*half)
    shifted_index = ts_masked.index + half_window
    ts_shifted = ts_masked.copy()
    ts_shifted.index = shifted_index

    # now resample on the shifted series; bins [t, t+freq) will contain points that were in [t-half, t+half)
    result = ts_shifted.resample(freq, offset=offset_td).agg(agg)

    # drop rows where all aggregated values are NaN (optional)
    if dropna:
        result = result.dropna(how='all')

    return result


def calculate_epoch_hoz_sigma(
    sigma_e, 
    sigma_n, 
    freq, 
    window, 
    **kwargs
    ) -> pd.Series:
    """
    Calculate the horizontal sigma from the East and North components at each epoch.
    
    freq : the desired epoch frequency, e.g. '1D'
    window : the half-window size either side of the frequency interval to examine for data.

    By default, epochs without any data will be dropped from the returned df.

    Passes kwargs directly to resample_with_fixed_window().
    """
    # Calculating uncertainties at midnight +/- 3 h every day, using a time series with irregularly spaced timestamps.
    uncs_e = resample_with_window(sigma_e, freq, window, **kwargs) 
    uncs_n  = resample_with_window(sigma_n, freq, window, **kwargs) 
    epoch_sigmas = np.sqrt(uncs_e**2 + uncs_n**2)    
    return epoch_sigmas


def calculate_vel_uncertainties(
    x : float | pd.Series,
    velocity : float | pd.Series, 
    epoch_uncertainty: pd.Series | tuple=(0.02, 0.02) 
    ) -> pd.Series:
    """ Calculate the velocity uncertainties per given time period.

    :param x: Series of displacements, or float displacement
    :param velocity: Series of velocities, with timestamps matching x, or float velocity
    :param epoch_uncertainty: the positional uncertainty at each epoch. 
        either tuple (t0, t1) or Pandas series of time points, which get shifted appropriately.
    :returns: pd.Series

    For reference, see Lokkegaard et al. (2024, Nat.Comms. Earth & Env.).
    """
    if isinstance(x, (collections.abc.Sequence, pd.Series, pd.DataFrame)):
        x = np.abs(x.shift(-1) - x).iloc[0:-1]

    if isinstance(epoch_uncertainty, pd.Series):
        eu1 = epoch_uncertainty
        # Negative shifting brings the next time point into this time point
        # Which produces a time series where the index denotes the start of the time period
        eu2 = epoch_uncertainty.shift(-1).iloc[:-1]
    else:
        eu1, eu2 = epoch_uncertainty

    return (np.sqrt(eu1**2 + eu2**2) / np.abs(x)) * velocity


def calculate_velocities_from_epochs(
    x : pd.Series, 
    sigmas : pd.Series
    ) -> tuple(pd.Series, pd.DataFrame):
    """ Calculate epoch-to-epoch velocities, where those epochs may be IRREGULARLY spaced. 

    The epochs used for the calculation are determined by the index of the `sigmas` df.
    
    """
    # Use the retained sigmas to find the corresponding X value.
    # We look in the smoothed data `filtd_disp`, NOT `filtd`, because the 
    # latter positions are really noisy and generate massive transient velocity spikes.
    epoch_x_pos_smoothed = x[x.index.isin(sigmas.index)]
    
    # Find the duration of the interval between each epoch - needed as the data are not evenly spaced in time. 
    w = epoch_x_pos_smoothed.index[1:]-epoch_x_pos_smoothed.index[:-1]
    # Calculate the corresponding multiplier for velocity m yr calc
    mults = v_mult(w)
    # Calculate the epoch-to-epoch displacement
    # For some reason this has to be a shift operation, not an iloc-offset based one.
    epoch_disp_x = np.abs(epoch_x_pos_smoothed.shift(-1) - epoch_x_pos_smoothed).iloc[:-1]

    # Calculate the corresponding velocity
    epoch_vel = epoch_disp_x * mults
    
    # Put a dummy velocity at the end of the series.
    # This denotes the end Timestamp of the epoch's velocity.
    epoch_vel.loc[epoch_x_pos_smoothed.index[-1]] = np.nan

    return (epoch_x_pos_smoothed, epoch_vel)


def calculate_epoch_velocities(
    x : pd.Series,
    sigma_e : pd.Series,
    sigma_n : pd.Series,
    freq: str,
    window='3h'
    ) -> pd.DataFrame:
    """
    Wrapper function to calculate epoch-to-epoch velocities and their corresponding uncertainties.

    """
    # Calculating uncertainties at midnight +/- 3 h every day, using a time series with irregularly spaced timestamps.
    epoch_sigmas = calculate_epoch_hoz_sigma(sigma_e, sigma_n, freq, window) 
    # Observed x (along-track displacement) per epoch and corresponding velocities from epoch to epoch.
    epoch_x, epoch_vel = calculate_velocities_from_epochs(x, epoch_sigmas)
    # Velocity uncertainites
    epoch_vel_unc = calculate_vel_uncertainties(epoch_x, epoch_vel, epoch_uncertainty=epoch_sigmas) 
    # Combine
    df = pd.concat([epoch_x, epoch_sigmas, epoch_vel, epoch_vel_unc], axis=1)
    df.columns = ['x_m', 'x_sigma_m', 'v_myr', 'v_uncertainty_myr']
    return df


def calculate_regular_velocities(
    x : pd.Series,
    period: str='24h',
    method : str='epoch',
    tz : str=None,
    flag_iterp : pd.Series=None,
    calc_uncertainty : bool=True, 
    sigmas : (pd.Series, pd.Series),
    **kwargs
    ) -> pd.DataFrame:
    """
    LEGACY FUNCTION AS OF NOV 2025, MAINTAINED FOR BACKWARD COMPATIBILITY WITH PROCESSING 2009-2013.
    USE THE EPOCH-TO-EPOCH FUNCTIONALITY INSTEAD!

    Calculate along-track velocity using either (a) first epoch or (b) average position in the 
    first hour of each period. By default, produces daily velocity time series. The time series are always 
    regularly spaced, so if a period has no observations available then the value will be filled
    by interpolation.
    
    Units metres/year.

    if tz is provided, first localise the (assumed UTC) measurements to local
    timezone before calculating the velocities.

    :param x: Series of along-track displacement from which to calculate velocity.
    :param method: str, either 'epoch' to difference between epochs separated by 24 h, 
    or 'average' to difference the mean position in the first hour of each day.
    :param tz: A 'tz-database' time zone string, e.g. America/Nuuk.
    :param flag_iterp:
    :param calc_uncertainty: if True, calculate uncertainty.
    :param sigmas: when calc_uncertainty is True, provide a tuple of sigma North and
    sigmaEast to  return uncertainties based on epoch-by-epoch
        sigma. Provide kwarg `half_window` to override defualt. 
        If None, 'standard' fixed epoch uncertainty will be used; override 
        default by supplying the kwarg `epoch_uncertainty` (see 
        .calculate_period_vel_uncertainties).

    """
    if tz is not None and tz != '':
        x = x.tz_localize(tz)

    if method == 'epoch':
        x = x.resample(period).first()
    elif method == 'average':
        x = x.resample(period).apply(lambda x: x[x.index.hour == 0].mean())
    else:
        raise ValueError('Unknown option provided for `method`.')
    
    # Shifting by -1 brings the next time point into this time point
    # Which produces a time series where the index denotes the start of the time period
    # and the stated velocity is for the period index:index+1
    v = np.abs(x.shift(-1) - x) * v_mult(period) 

    if flag_iterp is not None:
        if tz is not None and tz != '':
            flag_iterp = flag_iterp.tz_localize(tz)
            f = flag_iterp.resample(period).first()
            return pd.concat([v, f], axis='columns')

    if calc_uncertainty: 
        v.name = 'v_myr'
        if isinstance(sigmas, tuple):
            epoch_uncs = calculate_epoch_hoz_sigma(sigmas[0], sigmas[1], period, '3h', **kwargs)
            unc = calculate_vel_uncertainties(x, v, epoch_uncertainty=epoch_uncs)   
        else:
            unc = calculate_vel_uncertainties(x, v, **kwargs)
        v = v.to_frame()
        v = np.round(v, 2)
        v['unc_myr'] = unc
        v['disp_m'] = x
        if isinstance(sigmas, pd.DataFrame):
            v['sigma_epoch_m'] = epoch_uncs # = #np.round(

    return v
   

# def v_mult(
#     window : str
#     ) -> float:
#     """
#     Given a pandas window/offset string, Calculate the multipler needed to 
#     convert displacement to m/yr.

#     Info on to_offset:
#     https://stackoverflow.com/questions/40223470/how-do-i-get-at-the-pandas-offsets-object-given-an-offset-string
#     """
#     offset = to_offset(window)
#     year = pd.Timedelta(days=YEAR_LENGTH_DAYS)
#     multiplier = year / offset
#     return multiplier

def v_mult(windows):
    """
    Given a pandas window/offset string, Calculate the multipler needed to 
    convert displacement to m/yr.

    Info on to_offset:
    https://stackoverflow.com/questions/40223470/how-do-i-get-at-the-pandas-offsets-object-given-an-offset-string

    Parameters
    ----------
    windows : str | list[str] | np.ndarray | pd.Series
        One or more pandas window/offset strings (e.g., '1M', '2W', '3D').

    Returns
    -------
    float | np.ndarray
        The multiplier(s) needed to convert displacement to m/yr.
        Returns a scalar if input is scalar, otherwise an array.
    """
    scalar_input = False
    if isinstance(windows, str):
        # Handle scalar input
        windows = [windows]
        scalar_input = True

    # Convert to pandas Series for vectorized apply
    windows = pd.Series(windows)

    year = pd.Timedelta(days=YEAR_LENGTH_DAYS)
    multipliers = windows.apply(lambda w: year / to_offset(w))

    if scalar_input:
        return multipliers.iloc[0]
    return multipliers.to_numpy()


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
    toret['x_m'] = x
    toret['y_m'] = y
    toret['z_m'] = z
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
    toret['x_m'] = dx
    toret['y_m'] = dy
    toret['z_m'] = dz
    return toret
    
    