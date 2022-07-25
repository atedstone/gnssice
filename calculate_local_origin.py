"""
Calculate and save origin of local cartesian grid.

"""

import pandas as pd
import dask
import argparse

import gps

p = argparse.ArgumentParser('Calculate and save origin coordinates of local cartesian grid.')
p.add_argument('site', type=str, help='Name of site')
p.add_argument('geod_file', type=str, help='Filepath/name of GEOD parquet file from which to take origin coordinates')
args = p.parse_args()

# Load a concatenated GEOD file.
# Don't need much ...
# Let's assume we're working with Parquet only.
# Could load with Dask actualy

geod = dask.dataframe.read_parquet(p.geod_file)

# We have to recalculate local cartesian for pre_res first 100 points as these are not saved in pre_res.
pre_res_xyz = self.ell2xyz(geod.iloc[0:100,'Latitude'] * (math.pi/180),
                           geod.iloc[0:100,'Longitude'] * (math.pi/180),
                           geod.iloc[0:100,'Height']) 

x = np.median(pre_res_xyz['x'])
y = np.median(pre_res_xyz['y'])
z = np.median(pre_res_xyz['z'])

# save x,y,z.
df = pd.DataFrame({'x':[x], 'y':[y], 'z':[z]})
df.to_csv('%s_origin.csv' %p.site, index=False)
