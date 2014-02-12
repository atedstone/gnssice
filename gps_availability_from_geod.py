# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 10:45:05 2014

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gps_api
api = gps_api.gpsm()
import doy as doyf



keys = ['YY','DOY','Seconds','Latitude','Longitude','Height','SigN',
                    'SigE','SigH','RMS','#','Atm','Atm+-','Fract DOY','Epoch',
                    '#BF','NotF','NotFK','North','East','Set']

files = ['lev0_2009_2013_geod.dat',
         'lev1_2009_2013_geod.dat',
         'lev2_2009_2013_geod.dat',
         'lev3_2009_2013_geod.dat',
         'lev4_2009_2013_geod.dat',
         'lev5_2009_2013_geod.dat',
         'lev6_2009_2013_geod.dat',
         'lev7_2011_2013_geod.dat',
         'lev8_2012_2013_geod.dat']
         
siteids = ['lev0','lev1','lev2','lev3','lev4','lev5','lev6','lev7','lev8']         
   
all_avail = {}
master_ts = pd.date_range('2009-5-5','2013-5-12')
   
for siteid,f in zip(siteids,files):    
    print f             
    doys = pd.read_csv('M:/scratch/gps_final_velocities_spring2013/' + f,sep=' ',header=None,names=keys,usecols=[13])
    
    doys = doys.apply(np.floor)
    doys.drop_duplicates(inplace=True)
    
    # First convert cumulative DOY to proper DOY, then to python datetime
    if siteid == 'lev7':
        start = 2011
    elif siteid == 'lev8':
        start = 2012
    else:
        start = 2009
    reindexed_doy = api.reindex_year_doy(doy=doys['Fract DOY'].tolist(),start_year=start,save=False)
    as_datetime = doyf.doy2datetime(reindexed_doy[0],reindexed_doy[1])
    
    all_avail[siteid] = pd.Series(True,index=as_datetime)
    #flag = [1] * len(as_datetime)
    #avail = pd.Series(flag,index=as_datetime)

master = pd.DataFrame(index=master_ts,data=all_avail)
master.to_excel('M:/Data/lev_gps_availability.xls')

#k = 1
#fig = plt.figure(figsize=(10,5))
#ax = plt.subplot(111)
#for siteid in siteids:
#    ax.plot(master.index,master[siteid] * k,label=siteid,marker='o')
#    k += 1
#ax.set_yticks([0,1,2,3,4,5,6,7,8,9,10])
#ax.set_yticklabels(siteids)
#plt.ylim(0,10)
#ax.grid(color='lightgray', alpha=0.7)
#show_d3()
    


#doys = []
#for piece in chunker:
#    for row in piece:
#        doys.append(row['Fract DOY'])
#        
#obs_avail = pd.Series(doys)
#obs_avail.drop_duplicates(inplace=True)