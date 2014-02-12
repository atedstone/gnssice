# -*- coding: utf-8 -*-
"""
Created on Fri May 24 12:34:26 2013

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import numpy as np
import matplotlib.pyplot as plt

start_doy = 242. + (365*3)
end_doy = 244. + (365*3)
f = open("M:/scratch/gps_workspace_spr2013/lev3_2009_2013_geod.dat")
lines = []
for line in f.readlines():
    if float(line[13]) > start_doy and float(line[13]) < end_doy:
        lines.append(line)
f.close()