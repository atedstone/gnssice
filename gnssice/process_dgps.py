
"""
Set up session for Leverett GPS track processing.

Created on Thu Feb 16 10:52:10 2012

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import argparse
from gnssice import gnss

def cli():

    parser = argparse.ArgumentParser('Kinematic GNSS processing with TRACK')

    parser.add_argument('base', type=str, help='Name of base site')
    parser.add_argument('rover', type=str, help='Name of rover site')
    parser.add_argument('year', type=int, help='Year')
    parser.add_argument('doy_start', type=int, help='Start DOY')
    parser.add_argument('doy_end', type=int, help='End DOY')

    parser.add_argument('-ap', type=float, default=None, nargs=3, help='Rover a-priori coordinates (X Y Z')

    args = parser.parse_args()

    k = gnss.Kinematic()

    # Otherwise provide in the form [x, y, z].
    k.apriori = None

    if args.base == "levb":
        k.ion_stats = 1
        k.MW_WL = 1
        k.LG = 0.15
    elif args.base == "kely":
        k.ion_stats = 2
        k.MW_WL = 0.1
        k.LG = 0.15
    elif args.base == "rusb":
        k.ion_stats = 1 
        k.MW_WL = 0.1
        k.LG = 0.15
    elif args.base == "klsq":
        k.ion_stats = 2
        k.MW_WL = 0.1
        k.LG = 0.15
    else:
        print ("Unknown base station - cannot set processing parameters. Exiting.")
        exit

    k.apriori = args.ap

    k.track(args.base, args.rover, args.year, args.doy_start, args.doy_end) 

    
