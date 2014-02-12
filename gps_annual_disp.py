# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 12:02:37 2014

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gps_api
import math

keys = gps_api.Geod().keys


print 'p12'
p12 = pd.read_table("/scratch/s1144267/gps_final_velocities_spring2013/lev6_2009_2013_geod.dat",sep=" ",names=keys,usecols=['North','East'],skiprows=1740764,nrows=245)

print 'p13'                    
p13 = pd.read_table("M:/scratch/gps_final_velocities_spring2013/lev6_2009_2013_geod.dat",sep=" ",names=keys,usecols=['North','East'],skiprows=1896143)
 
print 'pythag'                   
p12e = p12['East'].mean()
p12n = p12['North'].mean()

p13e = p13['East'].mean()
p13n = p13['North'].mean()

diff_n = p13n - p12n
diff_e = p13e - p12e

diff_vector = math.sqrt(diff_n**2 + diff_e**2)

print diff_vector
