# -*- coding: utf-8 -*-
"""
Modify this script on a per-need basis to pull subsets of data from the 
leverett transect h5 file out and save into mat files for other folk to use.

Created on Fri Mar 15 13:23:14 2013

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import numpy as np
import matplotlib.pyplot as plt

from tables import *
import numpy as np
import datetime as dt
import doy as doyf
import time

import scipy.io as sio


gpsf = openFile("gps_final_velocities_spring2013/lev_transect_2009-spr2013.h5", mode="r")


for site in ['lev0','lev1','lev2','lev3','lev4','lev5','lev6']:
    date_start = dt.datetime(2009,5,1)
    date_end = dt.datetime(2009,9,1)
    slice_start = time.mktime(date_start.timetuple())
    slice_end = time.mktime(date_end.timetuple())
    
    v_6h = []
    v_6h = [[doyf.datetime2doy(dt.datetime.utcfromtimestamp(row['timestamp'])),row['velocity']] for row in gpsf.getNode('/' + site + '/v_6h').where("(timestamp >= " + str(slice_start) + ") & (timestamp <= " + str(slice_end) + ")")]
    v_6h = np.array(zip(*v_6h)).T
    
    v_24h = []
    v_24h = [[doyf.datetime2doy(dt.datetime.utcfromtimestamp(row['timestamp'])),row['velocity']] for row in gpsf.getNode('/' + site + '/v_24h').where("(timestamp >= " + str(slice_start) + ") & (timestamp <= " + str(slice_end) + ")")]
    v_24h = np.array(zip(*v_24h)).T
    
    t = []
    t = [[doyf.datetime2doy(dt.datetime.utcfromtimestamp(row['timestamp'])),row['temperature']] for row in gpsf.getNode('/' + site + '/temperature').where("(timestamp >= " + str(slice_start) + ") & (timestamp <= " + str(slice_end) + ")")]
    t = np.array(zip(*t)).T
    
    xyzt = []
    xyzt = [[row['x'],row['y'],row['z'],doyf.datetime2doy(dt.datetime.utcfromtimestamp(row['timestamp']))] for row in gpsf.getNode('/' + site + '/xyz').where("(timestamp >= " + str(slice_start) + ") & (timestamp <= " + str(slice_end) + ")")]
    xyzt = np.array(zip(*xyzt)).T
    
    sio.savemat('M:/Data/werder_package_2009/' + site + '_20090501-20090901_spr13h5.mat',{'v_24h':v_24h,'v_6h':v_6h, 'xyzt':xyzt, 't':t})