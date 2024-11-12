"""
Concatenative daily windowed GEOD files together from a given period into 
a single Parquet/CSV/TSV file.

Originally tried to do this via appending to a parquet file, but this does not
work - parquet files must be written in one go.  This is not a problem given that
this script only works over a maximum of a year and is intended for 10-sec data.

Andrew Tedstone, July 2022.
"""
from __future__ import annotations
import pandas as pd
import argparse
import os

import warnings
warnings.filterwarnings('ignore')

from gnssice import gnss as gps

def concatenate_daily_geod(
    base : str, 
    rover : str, 
    year : int, 
    start_doy : int, 
    finish_doy : int, 
    outformat : str, 
    exclude_doys : list
    ) -> None:

    output_fname = '%s_%s_%s_%s_%s_GEOD' %(rover, base, year, start_doy, finish_doy)

    if outformat == 'tsv':
        sep = '\t'
        output_fname += '.dat'
    elif outformat == 'csv':
        sep = ','
        output_fname += '.csv'
    elif outformat == 'parquet':
        output_fname += '.parquet'
    else:
        raise ValueError('Unknown file format')

    if exclude_doys is None:
        exclude_doys = []

    all_data = []        
    if start_doy == finish_doy:
        finish_it = start_doy + 1
    else:
        finish_it = finish_doy
    print ("Daily key: E=Excluded, N=No file, C=Concatenated.")
    for doy in range(start_doy, finish_it):
        if doy in exclude_doys:
            print ('E'+ str(doy) + ', ',)
            continue
        
        # Try to read in the daily GEOD file 
        try: 
            # First try new filename format
            saved_opts = dict(r=rover, b=base, y=year, d=str(doy).zfill(3))     
            fname = '{r}_{b}_{y}_{d}_{ftype}.dat'.format(ftype='GEOD', **saved_opts)
            if os.path.exists(fname):
                data = gps.read_track_geod_file(fname)
            else:
                # Try old filename format (deprecated)
                data = gps.read_track_geod_file(rover + '_' + base + '_' + str(doy).zfill(3) + 'GEOD.dat')
        except IOError:
            print ('S' + str(doy) + ', ',)
            continue
        print ('C' + str(doy) + ', ')
        # Remove the overlapping parts of the window
        data = data[(data['Fract_DOY'] >= doy) & (data['Fract_DOY'] < doy+1)]
        if len(data) == 0:
            continue

        all_data.append(data)

    if len(all_data) > 0:
        all_data = pd.concat(all_data, axis=0)

        if outformat in ['csv', 'tsv']:
            all_data.to_csv(output_fname, index=False, header=write_header, sep=sep)
        elif outformat == 'parquet':
            all_data.to_parquet(output_fname, index=False)

        print('Concatenated to %s.'%output_fname)
        return all_data

    else:
        print('No data to concatenate.')



def cli():

    p = argparse.ArgumentParser('Concatenate daily TRACK GEOD files. Run in the folder containing the files.')

    p.add_argument('base', type=str)
    p.add_argument('rover', type=str)
    p.add_argument('year', type=int)
    p.add_argument('start_doy', type=int)
    p.add_argument('finish_doy', type=int)

    p.add_argument('-outformat', choices=['parquet', 'csv', 'tsv'], default='parquet')
    p.add_argument('-exclude_doys', nargs='*', type=int)

    args = p.parse_args()
    
    data = concatenate_daily_geod(args.base, args.rover, args.year, args.start_doy, args.finish_doy,
        args.outformat, args.exclude_doys)