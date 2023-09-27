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

import argparse
import os


def load_xyz(files):
    store = []
    for f in files:
        store.append(pd.read_hdf(f, key='xyz'))
    xyz = pd.concat(store, axis=0)
    return xyz

def make_contiguous(df):
    # Coarsen the time series and interpolate over any gaps.
    xyzi = xyz.resample('1H').median()
    original = xyzi.x.notna() 
    xyzi = xyzi.interpolate()
    xyzi['original'] = original
    return xyzi

# +
def backdate(df, start='2021-05-01', freq='1H', sample=1000):
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
    
def fwddate(df, end='2023-05-01', freq='1H', sample=1000):
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

def calculate_disps(df, years, periods):
    """
    In case of annual displacement, it corresponds to year->year+1.

    :param df: dataframe containing x column.
    :param years: list of years to estimate for.
    :param periods: { 'period_name': [(start_month, start_day), (end_month, end_day)] }
    """
    store = {}
    for year in years:
        store_year = {}
        for period_name, bounds in periods.items():
            st = pd.Timestamp(year, bounds[0][0], bounds[0][1], 0, 0)
            if bounds[1][0] < bounds[0][0]:
                year_here = year + 1
            else:
                year_here = year
            en = pd.Timestamp(year_here, bounds[1][0], bounds[1][1], 23, 0)
            disp = xyzi.x.loc[en] - xyzi.x.loc[st]
            disp = np.abs(np.round(disp, 2))
            store_year[period_name] = disp
        store[year] = store_year
    return store


if __name__ == '__main__':        

    import seaborn as sns
    sns.set_style('whitegrid')

    # Eventually could wrap this logic up further, so that the script can process
    # every site then export to a spreadsheet.

    #files = ['lev5_rusb_2021_129_242_GEOD.h5', 'lev5_rusb_2022_137_138_GEOD.h5']
    files = ['lev6_rusb_2021_126_265_GEOD.h5', 'lev6_rusb_2022_134_135_GEOD.h5']
    #files = ['f003_rusb_2021_126_214_GEOD.h5', 'f003_rusb_2022_131_131_GEOD.h5']
    #files = ['kanu_2021_122_2022_135_geod.h5']
    #files = ['fs05_2021_126_2022_129_geod.h5']
    files = ['camp_2020_206_2023_126_disp.h5']
    files = ['/scratch/gnss/lev5/lev5_2021_130_2023_129_disp.h5']
    files = ['/scratch/gnss/lev6/lev6_2021_127_2023_123_disp.h5']
    files = ['/scratch/gnss/camp/camp_2021_126_2023_126_disp.h5']
    files = ['/scratch/gnss/f004/f004_2021_122_2023_123_disp.h5']
    files = ['/scratch/gnss/f003/f003_2021_127_2023_124_disp.h5']
    files = ['/scratch/gnss/kanu/kanu_2021_122_2023_119_disp.h5']

    periods = {
        'Annual':[(5,1), (4,30)],
        'Summer':[(5,1), (8,30)],
        'Winter':[(9,1), (4,30)],
        'ES':[(5,1), (6,30)],
        'LS':[(7,1), (8,30)]
        }

    years = [2021, 2022]

    xyz = load_xyz(files)
    xyzi = make_contiguous(xyz)
    xyzi = backdate(xyzi)
    xyzi = fwddate(xyzi)

    disps = calculate_disps(xyzi, years, periods)

    # To do: calculate average velocities for each period. 
    # Individual 24-h velocities are already available, should be OK to interpolate between the data blocks?
    # Velocities calculated here could otherwise differ slightly as they are not being computed from the full-temporal resolution gauss-smoothed xyz.

disps



