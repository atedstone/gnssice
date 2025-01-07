# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Convert Canadian PPP files to "Track" Parquet format
#
# This script is used to convert CSV position files originating from the Canadian PPP service into Parquet files with the same column names as those originating from the TRACK processing.
#
# AT, 18.12.2024

# %%
import pandas as pd

# %%
f = pd.read_csv('/Users/atedston/scratch/gnss_2023_2024/f003/f003_1230.csv')

# %%
mapping = {
    'latitude_decimal_degree':'Latitude',
    'longitude_decimal_degree':'Longitude',
    'ellipsoidal_height_m':'Height',
    'day_of_year':'DOY',
    'year':'YY'
}

drop = ['rcvr_clk_ns', 'decimal_hour']

# %%
timestamps = pd.to_datetime(f.year.astype(str)) + \
            (pd.to_timedelta(f.day_of_year.astype(int), unit='days') - pd.Timedelta(days=1)) + \
            pd.to_timedelta(f.decimal_hour, unit='h')

# %%
f

# %%
timestamps

# %%
fdoy = (timestamps.dt.day_of_year) + ((timestamps - timestamps.dt.floor('D')) / pd.Timedelta(24, 'h'))

# %%
daily_seconds_elapsed = (timestamps - pd.to_datetime(timestamps.dt.date)).dt.total_seconds()

# %%
ff = f.rename(mapping, axis=1)
ff['Fract_DOY'] = fdoy
ff['Seconds'] = daily_seconds_elapsed
ff = ff.drop(columns=drop)

# %%
ff

# %%
ff.to_parquet('/Users/atedston/scratch/gnss_2023_2024/f003/f003_1230_geod.parquet')

# %%
geod

# %%
plt.plot(geod.Fract_DOY, geod.Seconds, '.')

# %%
geod

# %%
