#!/usr/bin/env python3
"""
Set up the basic working environment according to the shell environment variables set.

Andrew Tedstone, July 2025.
"""

from gnssice.gnss import shellcmd
import os
from pathlib import Path
import datetime as dt
import argparse
import pandas as pd

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
    parser.add_argument("-overlap", action='store_true', help='Overlap the SP3 and IONEX files')
    
    args = parser.parse_args()
    run_setup(args.year, args.doy_start, args.doy_end, args.overlap)

def print_header(text):
    print('-------------------------------------------------------------------')
    print(f'{text}')

def run_setup(year, doy_start, doy_end, overlap):

    if os.environ['GNSS_WORK'] is None:
        raise RuntimeError('$GNSS_WORK not found in env variables.')

    print_header('Set up of basic working environment')
    print('$GNSS_WORK: ' + os.environ['GNSS_WORK'])

    # Derive time period timestamps from year/DOY
    date_start = dt.datetime.strptime(f"{year}-{doy_start:03d}", "%Y-%j").strftime('%Y-%m-%d')
    date_end = dt.datetime.strptime(f"{year}-{doy_end:03d}", "%Y-%j").strftime('%Y-%m-%d')
    print(f'Corresponding start and end dates: {date_start} - {date_end}')

    # Create working directory
    Path(os.environ['GNSS_WORK']).mkdir(exist_ok=True)

    # Create directory for user to place raw data into
    Path(os.environ['GNSS_PATH_RAWDATA']).mkdir(exist_ok=True)

    # Download VMF grids
    print_header('Downloading VMF grids')
    d = os.path.join(os.environ['GNSS_WORK'], 'vmf')
    sout, serr = shellcmd(f'get_vmf3 --start-date {date_start} --end-date {date_end} --output-dir {d}')
    if serr != '':
        raise RuntimeError(f'get_vmf3 failed: {serr}')

    # SP3 files
    print_header('Downloading SP3 orbit files, overlapping if requested')
    if overlap:
        _add = ' --overlap'
        sout, err = shellcmd('ln -s ')
    else:
        _add = ''
    sout, serr = shellcmd(f'get_orbits {year} {doy_start} {doy_end}{_add}')
    if serr != '':
        raise RuntimeError(f'get_orbits failed: {serr}')

    # IONEX files
    print_header('Downloading IONEX files, overlapping if requested')
    sout, serr = shellcmd(f'get_ionex {year} {doy_start} {doy_end}{_add}')

    # Sym-links depend on whether the workflow is setup for overlapped processing or not.
    print_header('Setting symlinks')
    wd = os.environ['GNSS_WORK']
    if overlap:
        shellcmd('ln -s ' + Path(os.environ['GNSS_PATH_SP3_OVERLAP']).stem + ' sp3', cwd=wd)
        shellcmd('ln -s ' + Path(os.environ['GNSS_PATH_IONEX_OVERLAP']).stem + ' ionex', cwd=wd)
        shellcmd('ln -s ' + Path(os.environ['GNSS_PATH_RINEX_OVERLAP']).stem + ' rinex', cwd=wd)
        print('Symlinks to sp3, ionex, rinex set from $GNSS_PATH_*_OVERLAP paths')
    else:
        shellcmd('ln -s ' + Path(os.environ['GNSS_PATH_SP3_DAILY']).stem + ' sp3', cwd=wd)
        shellcmd('ln -s ' + Path(os.environ['GNSS_PATH_IONEX_DAILY']).stem + ' ionex', cwd=wd)
        shellcmd('ln -s ' + Path(os.environ['GNSS_PATH_RINEX_DAILY']).stem + ' rinex', cwd=wd)
        print('Symlinks to sp3, ionex, rinex set from $GNSS_PATH_*_DAILY paths')


    print_header('Finished initial setup.')
    print('Working directory is ready for RINEX processing.')

if __name__ == '__main__':
    cli()