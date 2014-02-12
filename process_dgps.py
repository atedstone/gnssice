#! /usr/bin/env python
"""
Set up session for Leverett GPS track processing.

Created on Thu Feb 16 10:52:10 2012

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import gps
import ajtuseful as u
k = gps.Kinematic()

# If you want to set this interactively, set value to None.
# Otherwise provide in the form [x, y, z].
k.apriori = None

base = raw_input("Enter base for this processing session (levb, kely, lev4, cref available): ")

if base == "levb":
    k.ion_stats = 1
    k.MW_WL = 1
    k.LG = 0.15
elif base == "kely":
    k.ion_stats = 2
    k.MW_WL = 0.1
    k.LG = 0.15
elif base == "lev4":
    k.ion_stats = 1
    k.MW_WL = 0.1
    k.LG = 0.15
elif base == "cref":
    k.ion_stats = 1
    k.MW_WL = 1
    k.LG = 0.15
else:
    print "Unknown base station - cannot set processing parameters. Exiting."
    exit

rover = raw_input("Enter rover for this processing session: ")
doy_start = int(raw_input("Start on day: "))
doy_end = int(raw_input("End on day: ")) 

print "Enter apriori coordinates only if track has not already been run for \
the day previous to the one you want to start on."
if k.apriori == None:
    if u.confirm("Enter apriori coordinates?",resp=True):
        apr_x = float(raw_input("APR X: ")) 
         # track can't cope with -ve sign being provided on command line, hence
         # it has to be provided in cmd file.
        apr_y = float(raw_input("APR Y (-ve sign already in track cmd file): "))
        apr_z = float(raw_input("APR Z: "))
        k.apriori = [apr_x,apr_y,apr_z]

k.track(base,rover,doy_start,doy_end) 

    
