"""
Calculate seasonal and annual displacements from GNSS data.

As at June 2023, not ready to be used as a command line tool. Instead use magic
%run functionality in a Notebook or ipython env.

Andrew Tedstone (andrew.tedstone@unifr.ch), July 2022.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import numpy as np
from copy import deepcopy
import calendar
import seaborn as sns
sns.set_style('whitegrid')
import argparse
import os
from glob import glob


def load_multiple_xyz(files):
    store = []
    for f in files:
        store.append(pd.read_hdf(f, key='xyz'))
    xyz = pd.concat(store, axis=0)
    return xyz

def make_contiguous(df):
    # Coarsen the time series and interpolate over any gaps.
    xyzi = xyz.resample('1h').median()
    original = xyzi.x.notna() 
    xyzi = xyzi.interpolate()
    xyzi['original'] = original
    return xyzi

# +
def backdate(df, start='2021-05-01', freq='1h', sample=1000):
    # To estimate 1 May to 1 May displacement, for some sites we need to extend
    # the time range a bit.
    if df.index[0] > pd.Timestamp(start):
        # Make a synthetic time series between 1 May and the start of the series
        synth_ts = pd.date_range(start, df.index[0], freq=freq)

        # Use the first n. sampled hours.
        X = df.index[0:sample].to_julian_date()
        X = sm.add_constant(X)
        y = df.x.iloc[0:sample]
        m = sm.OLS(y, X)
        f = m.fit()

        # Predict the hourly displacements
        pred = f.predict(sm.add_constant(synth_ts.to_julian_date()))

        newx = pd.Series(pred, index=synth_ts, name='x').to_frame()

        toret = pd.concat((newx, df), axis=0)
        return toret
    else:
        return df
    
def fwddate(df, end='2024-05-01', freq='1h', sample=1000):
    # To estimate 1 May to 1 May displacement, for some sites we need to extend
    # the time range a bit.
    if df.index[-1] < pd.Timestamp(end):
        # Make a synthetic time series between 1 May and the start of the series
        synth_ts = pd.date_range(df.index[-1], end, freq=freq)

        # Use the first n. sampled hours.
        X = df.index[-sample:-1].to_julian_date()
        X = sm.add_constant(X)
        y = df.x.iloc[-sample:-1]
        m = sm.OLS(y, X)
        f = m.fit()

        # Predict the hourly displacements
        pred = f.predict(sm.add_constant(synth_ts.to_julian_date()))

        newx = pd.Series(pred, index=synth_ts, name='x').to_frame()
        toret = pd.concat((df, newx), axis=0)
        return toret
    else:
        return df


# -

def find_nearest_occupation(df, datetime):
    df = df[df.original == True]
    # returns absolute index into df e.g. array([5])
    iloc_idx = df.index.get_indexer([datetime], method='nearest')  
    # Get the named index
    loc_idx = df.index[iloc_idx]                    
    return loc_idx


def calculate_disps(df, years, periods):
    """
    In case of annual displacement, it corresponds to year->year+1.

    :param df: dataframe containing x column.
    :param years: list of years to estimate for.
    :param periods: { 'period_name': [(start_month, start_day), (end_month, end_day), tolerance_doys] }
    if tolerance_doys = -1 then no tolerance criteria are applied.
    """
    dates_store = []
    for year in years:
        dates_store_year = {}
        for period_name, bounds in periods.items():

            # Create start date and find nearest observation
            st1 = pd.Timestamp(year, bounds[0][0], bounds[0][1], 0, 0)
            st2 = pd.Timestamp(year, bounds[0][0], bounds[0][1], 23, 59)
            
            nearest_start = find_nearest_occupation(df, st1)
            diff1 = np.abs((nearest_start - st1).days).values
            print('START', period_name, bounds[0], nearest_start[0], diff1)

            # Create end date and find nearest observation
            if bounds[1][0] < bounds[0][0]:
                year_here = year + 1
            else:
                year_here = year
            en1 = pd.Timestamp(year_here, bounds[1][0], bounds[1][1], 0, 0)
            en2 = pd.Timestamp(year_here, bounds[1][0], bounds[1][1], 23, 59)
            nearest_end = find_nearest_occupation(df, en1)
            diff2 = np.abs((nearest_end - en1).days).values
            print('END  ', period_name, bounds[1], nearest_end[0], diff2)

            #!!! COME BACK TO THIS!!!
            last_year_obs = df[df.original == True].iloc[-1].name.year
            if period_name == 'Annual' and year_here > last_year_obs:
                continue

            if (~np.isnan(bounds[2])) and ((diff1 > bounds[2]) or (diff2 > bounds[2])):
                # Observations are outside permissible bounds
                #store_year[period_name] = np.nan
                continue
            else:
                # Observations are within permissible bounds

                # Calculate displacement through period
                disp = xyzi.x.loc[en1:en2].mean() - xyzi.x.loc[st1:st2].mean()
                disp = np.abs(np.round(disp, 2))

                # Calculate velocity of period in metres per year
                if calendar.isleap(year_here):
                    year_length = 366
                else:
                    year_length = 365
                period_length = (en2 - st1).days + 1
                vel = np.abs(np.round(disp / period_length * year_length, 2))

                if np.isnan(vel):
                    continue

                # Store results
                dates_store_year = dict(
                    period=period_name,
                    year=year,
                    disp_m=disp,
                    vel_m_yr=vel,
                    period_start_month_day=bounds[0], 
                    period_end_month_day=bounds[1], 
                    nearest_obs_start=nearest_start[0], 
                    nearest_obs_end=nearest_end[0]
                )
                dates_store.append(dates_store_year)        
        
    return dates_store


# +
#sites = ['camp', 'f003', 'f004', 'fs05', 'kanu', 'lev5', 'lev6']
sites = ['lev6']    
    
# Syntax:: period_name: [(start_month, start_day), (end_month, end_day), tolerance_days]
periods = {
    'Annual':[(5,1), (4,30), np.nan],
    'Summer':[(5,1), (8,30), 20],
    'Winter':[(9,1), (4,30), 20],
    'ES':[(5,1), (6,30), 5],
    'LS':[(7,1), (8,30), 5]
    }
# -
for site in sites:

    print(site)

    search_path = os.path.join(os.environ['GNSS_L2DIR'], site, '*_disp.h5')
    fn = glob(search_path)
    if len(fn) > 1:
        raise ValueError('More than one Level-2 displacement/velocity file (*_disp.h5) found.')
    elif len(fn) == 0:
        raise ValueError('No Level-2 displacement/velocity file (*_disp.h5) found.')
    fn = fn[0]    
    
    #xyz = load_multiple_xyz([fn])
    xyz = pd.read_hdf(fn, key='xyz')
    xyzi = make_contiguous(deepcopy(xyz))
    xyzi = backdate(xyzi)
    xyzi = fwddate(xyzi)

    year_start = xyzi.iloc[0].name.year
    year_end = xyzi.iloc[-1].name.year
    years = np.arange(year_start, year_end+1)
    
    disps = calculate_disps(xyzi, years, periods)
    
    # Convert to DataFrame
    disps_pd = pd.DataFrame(disps)
    # Add site ID
    disps_pd['site'] = site
    # Move site ID column to left-hand-side
    cols = disps_pd.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    disps_pd = disps_pd[cols]
    
    disps_pd.to_csv(os.path.join(os.environ['GNSS_L2DIR'], site, f'{site}_seasonal_annual_{year_start}_{year_end}.csv'), 
                    index=False)
    display(disps_pd)
# ---
# ## Reloading example

xyz.loc['2022-04']

store = []
for site in sites:
    # os.environ['GNSS_L2DIR']
    fglob = os.path.join('/scratch/hislide-level2/', f'{site}/{site}_seasonal_annual_*.csv')
    fs = glob(fglob)
    if fs is None:
        raise FileNotFoundError(f'No file for {site}.')
    elif len(fs) > 1:
        raise ValueError(f'More than one file found for {site}.')
    else:
        fs = fs[0]
    d = pd.read_csv(fs)
    store.append(d)
all_data = pd.concat(store, axis=0)

all_data

# To reduce to single...
vel_m_yr = all_data[all_data.period == 'Annual'].pivot(index='year', values='vel_m_yr', columns='site')
vel_m_yr

vel_m_yr.plot(marker='o')


