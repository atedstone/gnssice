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

sites = ['lev0','lev1','lev2','lev3','lev4','lev5','lev6','lev7','lev8']
site_descriptions = ['Site 1','Site 2','Site 3','Site 4','Site 5','Site 6','Site 7','Site 3N','Site 3M']
site_approx_elevation = [445,601,788,1051,1213,1473,1712,831,790]
site_approx_lat =  [67.06940169,
                    67.09262000,
                    67.10414350,
                    67.11504380,
                    67.12607527,
                    67.15336110,
                    67.16266243,
                    67.14634979,
                    67.12240000
                    ]
site_approx_lon =  [-50.13070662,
                    -50.03462000,
                    -49.81476781,
                    -49.40900344,
                    -49.01827330,
                    -48.37715363,
                    -47.55346856,
                    -49.81250749,
                    -49.80770000
                    ] 
site_approx_pos_from = '2012 May'                    

vel_f = ['lev0vel_2009_spr2013_UTC-2_20130528.mat',
         'lev1vel_2009_spr2013_UTC-2_20130528.mat',
         'lev2vel_2009_spr2013_UTC-2_20130528.mat',
         'lev3vel_2009_spr2013_UTC-2_20130528.mat',
         'lev4vel_2009_spr2013_UTC-2_20130528.mat',
         'lev5vel_2009_spr2013_UTC-2_20130528.mat',
         'lev6vel_2009_spr2013_UTC-2_20130528.mat',
         'lev7vel_2009_spr2013_UTC-2_20130528.mat',
         'lev8vel_2009_spr2013_UTC-2_20130528.mat']
         
xyzt_f = ['lev0xyzt_2009_spr2013_UTC-2_20130528.mat',
          'lev1xyzt_2009_spr2013_UTC-2_20130528.mat',
          'lev2xyzt_2009_spr2013_UTC-2_20130528.mat',
          'lev3xyzt_2009_spr2013_UTC-2_20130528.mat',
          'lev4xyzt_2009_spr2013_UTC-2_20130528.mat',
          'lev5xyzt_2009_spr2013_UTC-2_20130528.mat',
          'lev6xyzt_2009_spr2013_UTC-2_20130528.mat',
          'lev7xyzt_2009_spr2013_UTC-2_20130528.mat',
          'lev8xyzt_2009_spr2013_UTC-2_20130528.mat']         

# Set up classes to store each type of data.
class v_24h(IsDescription):
    timestamp = Int32Col()
    velocity = Float32Col()
    
class v_6h(IsDescription):
    timestamp = Int32Col()
    velocity = Float32Col()
    
class xyz(IsDescription):
    timestamp = Int32Col()
    x = Float32Col()
    y = Float32Col()
    z = Float32Col()
    
    
api = gps_api.gpsm()

h5file = openFile("lev_transect_2009-spr2013.h5", mode="w", title="Leverett transect processed final data")
root = h5file.root
root._v_attrs.created_by = 'Andrew Tedstone'
root._v_attrs.created_on = '2013 June 4'

# Set up a group for each transect site
groups = []
counter = 0
for s in sites:
    print "group for " + s
    group = h5file.createGroup(root,s,'Data for transect site with identifier ' + s)
    group._v_attrs.alt_name = site_descriptions[counter]
    group._v_attrs.approx_elevation_m = site_approx_elevation[counter]
    group._v_attrs.approx_lat = site_approx_lat[counter]
    group._v_attrs.approx_lon = site_approx_lon[counter]
    group._v_attrs.approx_pos_date = site_approx_pos_from
    groups.append(group)
        
        
    # PART 1 : VELOCITIES
    
    # Load the site mat file
    print ' Loading velocities mat file'
    mat = scipy.io.loadmat(vel_f[counter])
    
    
    # Set up 24h
    print ' Adding 24h velocities'
    t_24h = h5file.createTable(group,"v_24h",v_24h,"24h velocities for transect site")
    t_24h.attrs.timezone = "UTC-2"
    t_24h_r = t_24h.row  
    
    # First convert cumulative DOY to proper DOY, then to python datetime
    reindexed_doy = api.reindex_year_doy(doy=mat['v_24h'][:,0],start_year=2009,save=False)
    as_datetime = doyf.doy2datetime(reindexed_doy[0],reindexed_doy[1])
    
    i = 0
    row_c = 0
    for row in mat['v_24h']:  
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
    

    # Set up 6h
    print ' Adding 6h velocities'  
    tot_rows = mat['v_6h'][:,0].size
    t_6h = h5file.createTable(group,"v_6h",v_6h,"velcities smoothed with a 6-hour sliding window for transect site")
    t_6h.attrs.timezone = "UTC-2"
    t_6h_r = t_6h.row
     
    # First convert cumulative DOY to proper DOY, then to python datetime
    reindexed_doy = api.reindex_year_doy(doy=mat['v_6h'][:,0],start_year=2009,save=False)
    as_datetime = doyf.doy2datetime(reindexed_doy[0],reindexed_doy[1])

    i = 0
    row_c = 0
    flush_c = 0
    for row in mat['v_6h']:
        if isinstance(as_datetime[i],dt.datetime):
            # Convert python datetime to unix format.
            timestamp = calendar.timegm(as_datetime[i].timetuple())
        else:
            timestamp = 'NaN'

        t_6h_r['timestamp'] = timestamp
        t_6h_r['velocity'] = row[1]
        t_6h_r.append()
        row_c += 1
        if row_c > 1000:
            t_6h.flush()
            row_c = 0
            flush_c += 1
        i += 1
     
    # Attempt to clear the memory. 
    t_6h.flush()      
    mat = []
    gc.collect()
    
    
    
    # PART 2 : XYZ
    
    # Set up table
    t_xyz = h5file.createTable(group,"xyz",xyz,"X, Y, Z values for transect site")
    t_xyz.attrs.timezone = "UTC-2"
    t_xyz_r = t_xyz.row    
    
    # Load site xyzt
    print ' Loading xyzt mat file'
    mat = scipy.io.loadmat(xyzt_f[counter])
    tot_rows = mat['xyzt'][:,0].size
    
    # First convert cumulative DOY to proper DOY, then to python datetime
    reindexed_doy = api.reindex_year_doy(doy=mat['xyzt'][:,3],start_year=2009,save=False)
    as_datetime = doyf.doy2datetime(reindexed_doy[0],reindexed_doy[1])
    
    print ' Adding xyzt data'
    row_c = 0
    i = 0
    flush_c = 0
    for row in mat['xyzt']:
        if isinstance(as_datetime[i],dt.datetime):
            # Convert python datetime to unix format.
            timestamp = calendar.timegm(as_datetime[i].timetuple())
        else:
            timestamp = 'NaN'

        t_xyz_r['timestamp'] = timestamp
        t_xyz_r['x'] = row[0]
        t_xyz_r['y'] = row[1]
        t_xyz_r['z'] = row[2]
        t_xyz_r.append()
        row_c += 1
        if row_c > 1000:
            t_xyz.flush()
            row_c = 0
            flush_c += 1
        i += 1
    
    # Try to clear up    
    t_xyz.flush()
    mat = []
    gc.collect()
    
    
    # Increment counter which indexes mat file list locations
    counter += 1
    

    
    

    
    
print h5file
h5file.close()