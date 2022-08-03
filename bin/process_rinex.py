#!/usr/bin/env python
"""
Convert leica mdb files to daily windowed rinex files, in preparation for
kinematic processing using track.

Created on Thu Feb 16 11:49:30 2012

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import sys
import argparse
import datetime as dt

import gps
rinex = gps.RinexConvert()

if __name__ == '__main__':

    parser = argparse.ArgumentParser('Windowed RINEX file creation')

    parser.add_argument('site', type=str, help='Site to process')
    parser.add_argument('filename', type=str, help='Filename/prefix to process')
    parser.add_argument('file_type', type=str, help='File type: one of [d]aily Leica files, [s]ingle Leica file, [r]inex daily files e.g. Kellyville from SOPAC')

    parser.add_argument('-start', 
        type=lambda s: dt.datetime.strptime(s, '%Y-%m-%d'),
        help='For single files, day to begin processing, yyyy-mm-dd'
    )
    parser.add_argument('-finish', 
        type=lambda s: dt.datetime.strptime(s, '%Y-%m-%d'),
        help='For single files, day to finish processing, yyyy-mm-dd'
    )

    args = parser.parse_args()
    rinex.observer = "Andrew Tedstone"
    rinex.institution = "University of Fribourg"
           

    print ("----------------------------\nWindowed rinex file creation\n----------------------------")
    print ("Observer set as '" + rinex.observer + "'. Change in process_rinex if this is not correct!")
     

    # Daily Leica Files:
    if args.file_type == 'd':
        prefix = input("Enter the daily file series prefix (the bit before the .):")
        out = rinex.leica_joindaily(args.filename, args.site)
        rinex_files = out['joined_rinex']
        
        print ("Problem files: ") 
        print (out['problems'])

        for fn in rinex_files:
            print ("Window overlap on " + fn)
            rinex.window_overlap(fn,"-22:00:00",28) 
        print("Finished.")

    # Single Leica File:                      
    elif args.file_type == 's':
        rinex.window_overlap(args.filename, "-22:00:00", 28, 
                             leica=True,
                             site=args.site,
                             start_date=args.start, end_date=args.finish)
        print ("Finished.")

    # Daily Rinex Files:
    elif args.file_type == 'r':
        # First need to concatenate files
        y = input("Enter the rinex year suffix (e.g. 11 for 2011): ")
        fn = "all_" + args.site + "." + str(y) + "o"
        print ("File will be saved to " + fn)
        cmd = "teqc " + args.site + "*." + str(y) + "o > " + fn
        status = gps.shellcmd(cmd)
        print (status['stdout'])
        print (status['stderr'])
        
        print ("Starting window overlap.")
        rinex.window_overlap(fn,"-22:00:00",28) 
        print ("Finished.")

    else: 
        print ("Invalid option specified, exiting.")
        sys.exit()