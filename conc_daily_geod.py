"""
Concatenative daily windowed GEOD files together from a given period into 
a single Parquet/CSV/TSV file.

Andrew Tedstone, July 2022.
"""

from __future__ import annotations
import pandas as pd
import argparse

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

    fname = base + "_" + rover + "_" + str(year) + "_geod"

    if outformat == 'tsv':
        sep = '\t'
        open_as = 'a+'
        fname += '.dat'
    elif outformat == 'csv':
        sep = ','
        open_as = 'a+'
        fname += '.csv'
    elif outformat == 'parquet':
        open_as = 'ab+'
        fname += '.parquet'
    else:
        raise ValueError('Unknown file format')

    written = False
    with open(fname, open_as) as f:
        
        print ("Daily key: E=Excluded, N=No file, C=Concatenated.")
        for doy in range(start_doy, end_doy):
            if doy in exclude_doy:
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
            data = data[(data['Fract DOY'] >= doy) & (data['Fract DOY'] < doy+1)]
            if len(data) == 0:
                continue

            # Append to file
            if outformat in ['csv', 'tsv']:
                data.to_csv(f, index=False, header=write_header, sep=sep)
            elif outformat == 'parquet':
                data.to_parquet(f, index=False)

            written = True

    if written:
        print('Concatenated to %s.'%fname)
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

    concatenate_daily_geod(**args)