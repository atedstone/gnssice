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

def cli():

    parser = argparse.ArgumentParser('Windowed RINEX file creation')

    parser.add_argument('site', type=str, help='Site to process')
    parser.add_argument('-filename', type=str, help='Filename to process, needed if single Leica file.')
    parser.add_argument('file_type', type=str, help='File type: one of [g] CryoLogger GVT daily UBX files, [s]ingle Leica file, [r]inex daily files e.g. Kellyville from SOPAC, [d]aily Leica files')

    parser.add_argument('-start', 
        type=lambda s: dt.datetime.strptime(s, '%Y-%m-%d'),
        help='For Leica single files, day to begin processing, yyyy-mm-dd'
    )
    parser.add_argument('-finish', 
        type=lambda s: dt.datetime.strptime(s, '%Y-%m-%d'),
        help='For Leica single files, day to finish processing, yyyy-mm-dd'
    )

    parser.add_argument('-overlap', action='store_true', 
        help='Produce overlapping RINEX files beginning at -22:00:00 and extending 28 hours.')

    args = parser.parse_args()
    rinex.observer = os.environ['GNSS_RINEX_OBSERVER']
    rinex.institution = os.environ['GNSS_RINEX_INSTITUTION']

    print ("----------------------------\nWindowed rinex file creation\n----------------------------")
    print ("Observer set as '" + rinex.observer + "'. Change in your env variables if this is not correct!")
     
    if args.overlap:
        start_at = '-22:00:00'
        window = 28
    else:
        start_at = '00:00:00'
        window = 24

    if args.overlap:
        rinex_files = glob(os.path.join(os.environ['DATA_PATH_RINEX_DAILY'], args.site, '*'))
        for f in rinex_files:
            pass



    # CryoLogger GVT / u-blox files:
    if args.file_type == 'g':
        # First convert to RINEX
        ubx_files = glob(os.path.join(os.environ['DATA_PATH_BINARY'], args.site, '*.ubx'))
        print('Found {n} files...'.format(n=len(ubx_files)))
        for f in ubx_files:
            sout, serr = rinex.gvt_to_rinex(f, args.site, os.path.join(os.environ['DATA_PATH_RINEX_DAILY'], args.site))
            if serr is not None:
                print(serr)
                raise IOError  
    
    # Single Leica File (e.g. L1200 systems):                     
    elif args.file_type == 's':
        rinex.leica_to_daily_rinex(args.filename,
                             site=args.site,
                             start_date=args.start, end_date=args.finish)


    # Overlapping of daily Rinex files:
    elif args.file_type == 'r':
        # First need to concatenate files
        y = input("Enter the rinex year suffix (e.g. 11 for 2011): ")
        fn = "all_" + args.site + "." + str(y) + "o"
        print ("File will be saved to " + fn)
        cmd = "teqc " + args.site + "*." + str(y) + "o > " + fn
        status = gnss.shellcmd(cmd)
        print (status['stdout'])
        print (status['stderr'])
        
        print ("Starting window overlap.")
        rinex.window_overlap(fn, start_at, window) 
        print ("Finished.")


    # Daily Leica files (old L500 systems only):
    elif args.file_type == 'd':

        print('This option is out-of-date and has been disabled.')
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

    else: 
        print ("Invalid option specified, exiting.")
        sys.exit()