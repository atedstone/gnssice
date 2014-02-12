# -*- coding: utf-8 -*-
"""

TO DO : add GEOD dat files to h5 file, which will compress them ready for archiving.
Created on Thu Feb 21 12:15:23 2013

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

from tables import *
import gc
import scipy.io

import datetime as dt
import calendar # fpr unix timestamps
import doy as doyf

import gps_api

sites = ['lev0','lev1','lev2','lev3','lev4','lev5','lev6']             

vel_f = ['lev0vel_2009_spr2013_UTC-2_20130607.mat',
         'lev1vel_2009_spr2013_UTC-2_20130607.mat',
         'lev2vel_2009_spr2013_UTC-2_20130607.mat',
         'lev3vel_2009_spr2013_UTC-2_20130607.mat',
         'lev4vel_2009_spr2013_UTC-2_20130607.mat',
         'lev5vel_2009_spr2013_UTC-2_20130607.mat',
         'lev6vel_2009_spr2013_UTC-2_20130607.mat']

for f in vel_f:
    mat = scipy.io.loadmat(f)
    if mat.has_key('v_24h_opt'):
        print f + ' - True'
    else:
        print f + ' - FALSE'
         
      

# Set up classes to store each type of data.
class v_24h(IsDescription):
    timestamp = Int32Col()
    velocity = Float32Col()
       
api = gps_api.gpsm()

h5file = openFile("lev_transect_2009-spr2013.h5", mode="a", title="Leverett transect processed final data")
root = h5file.root

# Set up a group for each transect site
counter = 0
for s in sites:
    
    group = h5file.getNode('/' + s)
    print group
        
    try:
        table = h5file.removeNode('/' + s,'v_24h')
        print 'v_24h node removed'
    except:
        print 'v_24h node does not yet exist'
    
          
            
    # Load the site mat file
    print ' Loading velocities mat file'
    mat = scipy.io.loadmat(vel_f[counter])
    
    
    # Set up 24h
    print ' Adding 24h velocities'
    t_24h = h5file.createTable(group,"v_24h",v_24h,"24h velocities for transect site ('opt' dataset)")
    t_24h.attrs.timezone = "UTC-2"
    t_24h_r = t_24h.row  
    
    # First convert cumulative DOY to proper DOY, then to python datetime
    reindexed_doy = api.reindex_year_doy(doy=mat['v_24h_opt'][:,0],start_year=2009,save=False)
    as_datetime = doyf.doy2datetime(reindexed_doy[0],reindexed_doy[1])
    
    i = 0
    row_c = 0
    for row in mat['v_24h_opt']:  
        if isinstance(as_datetime[i],dt.datetime):
            # Convert python datetime to unix format.
            timestamp = calendar.timegm(as_datetime[i].timetuple())
        else:
            timestamp = -9999

        t_24h_r['timestamp'] = timestamp
        t_24h_r['velocity'] = row[1]
        t_24h_r.append()
        row_c += 1
        if row_c > 1000:
            t_24h.flush()
            row_c = 0
        i += 1
    t_24h.flush()
    
   
    # Increment counter which indexes mat file list locations
    counter += 1
       
    
print h5file
h5file.close()