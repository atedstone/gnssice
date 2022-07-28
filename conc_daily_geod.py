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

import pdb

import gps

def concatenate_daily_geod(
    base : str, 
    rover : str, 
    year : int, 
    start_doy : int, 
    finish_doy : int, 
    outformat : str, 
    exclude_doys : list
    ) -> None:

    fname = 'geod_%s_%s_%s_%s_%s' %(base, rover, year, start_doy, end_doy)

    if outformat == 'tsv':
        sep = '\t'
        fname += '.dat'
    elif outformat == 'csv':
        sep = ','
        fname += '.csv'
    elif outformat == 'parquet':
        fname += '.parquet'
    else:
        raise ValueError('Unknown file format')

    if exclude_doys is None:
        exclude_doys = []

    all_data = []        
    print ("Daily key: E=Excluded, N=No file, C=Concatenated.")
    for doy in range(start_doy, finish_doy):
        if doy in exclude_doys:
            print ('E'+ str(doy) + ', ',)
            continue
        
        # Try to read in the daily GEOD file 
        try:      
            data = gps.read_track_file(rover + '_' + base + '_' + str(doy).zfill(3) + 'GEOD.dat')
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
        print(all_data.dtypes)

        if outformat in ['csv', 'tsv']:
            all_data.to_csv(fname, index=False, header=write_header, sep=sep)
        elif outformat == 'parquet':
            all_data.to_parquet(fname, index=False)

        print('Concatenated to %s.'%fname)
        return all_data

    else:
        print('No data to concatenate.')




if __name__ == '__main__':

    p = argparse.ArgumentParser('Concatenate daily TRACK GEOD files.')

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