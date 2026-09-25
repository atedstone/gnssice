""" 
Helper functions to work with Canadian PPP outputs

Andrew Tedstone, Aug/Sep 2026
"""

import pandas as pd
import dask.dataframe as dd

def read_csv(p, round_to='1s'):
    """
    p : path to CSV file(s), can contain wildcards
    """
    # Load the CSV files, which contain decimal degrees coordinates
    ddf = dd.read_csv(p)
    df = ddf.compute()
    df.index = (pd.to_datetime(df['year'].astype(int).astype(str), format='%Y')
      + pd.to_timedelta(df['day_of_year'] - 1, unit='D')
      + pd.to_timedelta(df['decimal_hour'], unit='h'))
    if round_to is not None:
        df.index = df.index.round(round_to)
    return df

def read_pos(p):
    """
    p : path to CSV file(s), can contain wildcards
    """

    # Load the POS files, which don't have decimal degrees but do have sigmas info
    # To fix!: FutureWarning: Support for nested sequences for 'parse_dates' in pd.read_csv is deprecated. Combine the desired columns with pd.to_datetime after parsing instead.
    ddf = dd.read_csv(p, sep=r'\s+', skiprows=3, on_bad_lines='warn')
    df_pos = ddf.compute()
    # Concatenating with T in the middle makes these isoformat datetime strings
    df_pos.index = pd.to_datetime(ddf['YEAR-MM-DD'] + 'T' + ddf['HR:MN:SS.SS']) 
    return df_pos

def to_track_format(df_csv, df_pos):
    # Check that the dataframes match precisely in length
    assert len(df_csv) == len(df_pos)

    # Now sort both by time
    df_csv = df_csv.sort_index()
    df_pos = df_pos.sort_index()

    # The PPP outputs contain a single row for the next day, which results in duplicated time stamps
    # ...remove them.
    df_csv = df_csv[~df_pos.index.duplicated()]
    df_pos = df_pos[~df_pos.index.duplicated()]

    # Renaming 
    mapping = {
        'latitude_decimal_degree':'Latitude_deg',
        'longitude_decimal_degree':'Longitude_deg',
        'ellipsoidal_height_m':'Height_m',
        'day_of_year':'DOY',
        'year':'YY'
    }

    # Do the renaming
    df_csv = df_csv.rename(mapping, axis=1)

    # Remove the following columns
    drop = ['rcvr_clk_ns', 'decimal_hour']
    df_csv = df_csv.drop(columns=drop)

    # Calculate fractional DOY (for correspondence with TRACK outputs)
    fdoy = (df_csv.index.day_of_year) + ((df_csv.index - df_csv.index.floor('D')) / pd.Timedelta(24, 'h'))
    df_csv['Fract_DOY'] = fdoy

    # Calculate daily seconds elapsed (for correspondence with TRACK outputs)
    daily_seconds_elapsed = (df_csv.index - pd.to_datetime(df_csv.index.date)).total_seconds()
    df_csv['Seconds'] = daily_seconds_elapsed

    # Add the needed POS columns into the CSV frame, using standard Track names
    df_csv['SigN_cm'] = df_pos['SDLAT(95%)'] * 100
    df_csv['SigE_cm'] = df_pos['SDLON(95%)'] * 100
    df_csv['SigH_cm'] = df_pos['SDHGT(95%)'] * 100

    return df_csv