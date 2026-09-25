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
# Most functionality moved to `ppp.py` module 22.09.2026, AT

# %%
import pandas as pd
import matplotlib.pyplot as plt
from gnssice import ppp

# %%
csv = ppp.read_csv('/Users/atedston/scratch/flowstate-gnss-processing/ilhw/full_output/Reach-Base_raw_20260425171558.csv')
pos = ppp.read_pos('/Users/atedston/scratch/flowstate-gnss-processing/ilhw/full_output/Reach-Base_raw_20260425171558.pos')
df = ppp.to_track_format(csv, pos)

# %%
df = df.resample('10s').first()

# %%
plt.figure()
plt.plot(timestamps, df.Latitude, '.', alpha=0.2)

# %%
df.to_parquet('/scratch/flowstate-gnss-level1/ilhw/ilhw_pppp_2026_115_115_GEOD.parquet')
