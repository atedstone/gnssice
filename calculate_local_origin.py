#!/usr/bin/env python
"""
Calculate and save origin of local cartesian grid.

"""

import pandas as pd
from dask import dataframe as dd
import argparse
import math
import numpy as np

import gps

p = argparse.ArgumentParser('Calculate and save origin coordinates of local cartesian grid.')
p.add_argument('site', type=str, help='Name of site')
p.add_argument('geod_file', type=str, help='Filepath/name of GEOD parquet file from which to take origin coordinates')
args = p.parse_args()

# Load a concatenated GEOD file.
# Don't need much ...
# Let's assume we're working with Parquet only.
# Could load with Dask actualy

#geod = dd.read_parquet(args.geod_file)
geod = pd.read_parquet(args.geod_file)

# We have to recalculate local cartesian for pre_res first 100 points as these are not saved in pre_res.
pp = gps.PostProcess()
pre_res_xyz = pp.ell2xyz(geod['Latitude'].iloc[0:100] * (math.pi/180),
                           geod['Longitude'].iloc[0:100] * (math.pi/180),
                           geod['Height'].iloc[0:100]) 

x = np.median(pre_res_xyz['x'])
y = np.median(pre_res_xyz['y'])
z = np.median(pre_res_xyz['z'])

lat0 = geod['Latitude'].iloc[0:100].median()
lon0 = geod['Longitude'].iloc[0:100].median()

# save x,y,z.
df = pd.DataFrame({'x0':[x], 'y0':[y], 'z0':[z], 'lat0':lat0, 'lon0':lon0})
df.to_csv('%s_origin.csv' %args.site, index=False)
