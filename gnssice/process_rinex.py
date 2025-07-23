"""
Convert CryoLogger GVT and Leica mdb files to daily windowed rinex files, in preparation for
kinematic processing using track.

Created on Thu Feb 16 11:49:30 2012
Updated July 2025 for GVT/u-blox files.

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import sys
import argparse
import datetime as dt
from glob import glob
import os

from gnssice import gnss
rinex = gnss.RinexConvert()

def print_header(text):
    print('-------------------------------------------------------------------')
    print(f'\t {text}')
    print('-------------------------------------------------------------------')

def cli():

    if os.environ['GNSS_PATH_RAWDATA'] is None:
        raise ValueError('Environment variable GNSS_PATH_RAWDATA not found. Have you set up your environment?')

    parser = argparse.ArgumentParser("""
        Windowed RINEX file creation.
        -----------------------------
        Generate daily RINEX files from binary receiver formats.
        Optionally, produce overlapping RINEX files.

        When called with type F, LS, LD, this script always outputs daily RINEX files.
        """
    )

    parser.add_argument('site', type=str, help='Site to process')
    parser.add_argument('file_type', type=str,
        help=("""
            Input file type: one of 
            - G  : CryoLogger GVT daily UBX files
            - R  : daily RINEX files (located in $GNSS_PATH_RINEX_DAILY directory; only overlapping is possible)
            - LS : single Leica MDB (M0x) file 
            - LD : daily Leica MDB files'
        """)
    )

    parser.add_argument('-overlap', action='store_true', 
        help='Produce overlapping RINEX files beginning at -22:00:00 and extending 28 hours.')

    parser.add_argument('-filename', type=str, help='For Leica single files (LS) only: filename to process in DATA_PATH_BINARY/<site>/.')
    parser.add_argument('-start', 
        type=lambda s: dt.datetime.strptime(s, '%Y-%m-%d'),
        help='For Leica single files (LS) only: day to begin processing, yyyy-mm-dd'
    )
    parser.add_argument('-finish', 
        type=lambda s: dt.datetime.strptime(s, '%Y-%m-%d'),
        help='For Leica single files (LS) only: day to finish processing, yyyy-mm-dd'
    )

    args = parser.parse_args()
    rinex.observer = os.environ['GNSS_RINEX_OBSERVER']
    rinex.institution = os.environ['GNSS_RINEX_INSTITUTION']

    print_header('Windowed RINEX file creation')
    print ("Observer set as '" + rinex.observer + "'. Change in your env variables if this is not correct!")

    # CryoLogger GVT / u-blox files UBX to RINEX:
    if args.file_type == 'G':
        print_header('Converting GVT --> RINEX')
        # First convert to RINEX
        scan_path = os.path.join(os.environ['GNSS_PATH_RAWDATA'], args.site, '*.ubx')
        print(f'Scanning {scan_path}')
        ubx_files = glob(scan_path)
        print('Found {n} files...'.format(n=len(ubx_files)))
        for f in ubx_files:
            sout, serr = rinex.gvt_to_rinex(f, args.site, os.path.join(os.environ['GNSS_PATH_RINEX_DAILY'], args.site))
            if serr is not None:
                print(serr)
                raise IOError
        print('Finished.')

    # Single Leica File MDB to RINEX (e.g. L1200 systems):                     
    elif args.file_type == 'LS':
        print_header('Converting Leica MDB --> RINEX')
        pth = os.path.join(os.environ['GNSS_PATH_RAWDATA'], args.site, args.filename)
        rinex.leica_to_daily_rinex(pth,
                             site=args.site,
                             start_date=args.start, end_date=args.finish)

    # Daily Leica files MDB to RINEX (e.g. L500 systems):
    elif args.file_type == 'LD':

        print_header('Leica daily files: This option is out-of-date and has been DISABLED.')
        sys.exit()

        prefix = input("Enter the daily file series prefix (the bit before the .):")
        out = rinex.leica_joindaily(args.filename, args.site)
        rinex_files = out['joined_rinex']
        
        print ("Problem files: ") 
        print (out['problems'])

        for fn in rinex_files:
            print ("Window overlap on " + fn)
            rinex.window_overlap(fn, start_at, window) 
        print("Finished.")

    elif args.file_type == 'R':
        print_header('RINEX files option')

    # Invalid option
    else: 
        print ("Invalid option specified, exiting.")
        sys.exit()

    # Now do overlapping, always from daily RINEX files (never from any other type of file)
    if args.overlap:
        print_header('Commence window overlaps')
   
        rinex_files = glob(os.path.join(os.environ['GNSS_PATH_RINEX_DAILY'], args.site, '*o'))
        print('Found {n} files...'.format(n=len(rinex_files)))
        for f in rinex_files:
            rinex.window_overlap(f, st_timestart='220000', dh=28)
        print('Finished.')

    if (not args.overlap) and (args.file_type == 'R'):
        print('You specified RINEX files but did not ask them to be overlapped. There is no work to do.')