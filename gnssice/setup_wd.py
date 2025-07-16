#!/usr/bin/env python3
"""
Set up the basic working environment according to the shell environment variables set.

Andrew Tedstone, July 2025.
"""

from gnss import shellcmd
import os
from pathlib import Path
import datetime as dt
import argparse

def cli():
    parser = argparse.ArgumentParser(
        description="""
        Set up working directory and download VMF3, SP3 and IONEX files. 

        This script can also be run to extend the time series of these files in the existing working directory.
        Note: no checks for whether the files are already available are undertaken by this script.
        """
    )
    parser.add_argument("year", type=int, help="Year for downloads")
    parser.add_argument("doy_start", type=int, help="Day of year on which to start downloads")
    parser.add_argument("doy_end", type=int, help="Day of year on which to end downloads")
    parser.add_argument("-overlap", action='store_true', type=bool, help='Overlap the SP3 and IONEX files')
    
    args = parser.parse_args()
    run_setup(args.year, args.doy_start, args.doy_end, args.overlap)

def print_header(text):
    print('-------------------------------------------------------------------')
    print(f'{text}')

def run_setup(year, doy_start, doy_end, overlap):

    if os.environ['GNSS_WORK'] is None:
        raise RuntimeError('$GNSS_WORK not found in env variables.')
        
    # Derive time period timestamps from year/DOY
    date_start = datetime.strptime(f"{year}-{doy_start:03d}", "%Y-%j")
    date_end = datetime.strptime(f"{year}-{doy_end:03d}", "%Y-%j")

    # Create working directory
    Path(os.environ['GNSS_WORK']).mkdir(exist_ok=True)

    # Create directory for user to place raw data into
    Path(os.environ['GNSS_PATH_RAWDATA']).mkdir(exist_ok=True)

    # Download VMF grids
    print_header('Downloading VMF grids')
    sout, serr = shellcmd(f'get_vmf3 {date_start} {date_end}')

    # SP3 files
    print_header('Downloading SP3 orbit files, overlapping if requested')
    if overlap:
        _add = ' --overlap'
    else:
        _add = ''
    sout, serr = shellcmd(f'get_orbits {year} {doy_start} {doy_end}{_add}')

    # IONEX files
    print_header('Downloading IONEX files')
    Path(os.environ['GNSS_PATH_IONEX_DAILY']).mkdir(exist_ok=True)
    ndays = doy_end - doy_start
    sout, serr = shellcmd('cd ' + os.environ['GNSS_PATH_IONEX_DAILY'])
    sout, serr = shellcmd(f'sh_get_ion -yr {year} -doy {doy_start} -ndays {ndays}')
    sout, serr = shellcmd('cd ' + os.environ['GNSS_WORK'])
    if overlap:
        print_header('Overlapping IONEX files')
        for day in date_range:
            sout, serr = shellcmd('splice_ionex {day} 2')

    print_header('Finished initial setup.')
    print('Working directory is ready for RINEX processing.')

if __name__ == '__main__':
    cli()