#!/usr/bin/python
"""
Convert leica mdb files to daily windowed rinex files, in preparation for
kinematic processing using track.

Created on Thu Feb 16 11:49:30 2012

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import gps
rinex = gps.RinexConvert()

rinex.observer = "Andrew Tedstone"
rinex.institution = "University of Edinburgh"
       
print "----------------------------\nWindowed rinex file creation\n----------------------------"
print "Observer set as '" + rinex.observer + "'. Change in process_rinex if this is not correct!"
 
print "Options:"
print "     d : Process from daily leica files."
print "     s : Process from single big leica file."
print "     r : Process from daily rinex files (e.g. kellyville files from SOPAC)." 
opt = raw_input("Choose an option (d,s,r): ")
s_site = raw_input("Enter the 4-character site identifier: ")

# Daily Leica Files:
if opt == 'd':
    prefix = raw_input("Enter the daily file series prefix (the bit before the .):")
    out = rinex.leica_joindaily(prefix,s_site)
    rinex_files = out['joined_rinex']
    
    print "Problem files: " 
    print out['problems']

    for fn in rinex_files:
        print "Doing window overlap on " + fn
        rinex.window_overlap(fn,"-22:00:00",28) 
    print "Finished."

# Single Leica File:                      
elif opt == 's':
    s_fn = raw_input("Enter the filename: ")
    s_start_date = input("Enter start date, in format [yyyy,m,d]: ")
    s_end_date = input("Enter end date, in format [yyyy,m,d]: ")
    rinex.window_overlap(s_fn,"-22:00:00",28, 
                         leica=True,
                         site=s_site,
                         start_date=s_start_date,end_date=s_end_date)
    print "Finished."

# Daily Rinex Files:
elif opt == 'r':
    # First need to concatenate files
    y = raw_input("Enter the rinex year suffix (e.g. 11 for 2011): ")
    fn = "all_" + s_site + "." + str(y) + "o"
    print "File will be saved to " + fn
    cmd = "teqc " + s_site + "*." + str(y) + "o > " + fn
    status = gps.shellcmd(cmd)
    print status['stdout']
    print status['stderr']
    
    print "Starting window overlap."
    rinex.window_overlap(fn,"-22:00:00",28) 
    print "Finished."

else: 
    print "Invalid option specified, exiting."
    exit