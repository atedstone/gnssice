# -*- coding: utf-8 -*-
"""
http://stackoverflow.com/questions/13590484/calculating-euclidean-distance-between-consecutive-points-of-an-array-with-numpy
Created on Wed Jan 15 10:32:26 2014

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mpl_toolkits.basemap.pyproj as pyproj
import datetime as dt
from scipy import signal
import math

import gps_api
geod = gps_api.Geod()
api = gps_api.gpsm()

def gaussfiltcoef(SR,fco):
    """
    %GAUSSFILTCOEF  Return coefficients of Gaussian lowpass filter.
    % SR=sampling rate, fco=cutoff (-3dB) freq, both in Hz. 
    % Coeffs for FIR filter of length L (L always odd) are computed.
    % This symmetric FIR filter of length L=2N+1 has delay N/SR seconds.
    % Examples of use
    %    Compute Gaussian filter frequency response for SR=1000, fco=50 Hz:
    %    freqz(gaussfiltcoef(1000,50),1,256,1000);
    %    Filter signal X sampled at 5kHz with Gaussian filter with fco=500:
    %    y=filter(gaussfiltcoef(5000,500),1,X);
    % SR, fco are not sanity-checked.  WCR 2006-10-11.
    
    PORTED FROM THE MATLAB VERSION!
    """

    
    a = 3.011 * fco
    N = math.ceil(0.398*SR/fco) #filter half-width, excluding midpoint
    b = np.zeros((N*2+1,))
    #Width N corresponds to at least +-3 sigma which captures at least 99.75%
    #of area under Normal density function. sigma=1/(a*sqrt(2pi)).
    L=2*N+1            #full length of FIR filter
    for k in range(-N,N):
        b[k+N+1]=3.011*(fco/SR)*np.exp(-np.pi*(a*k/SR)**2)
    
    #b(k) coeffs computed above will add to almost exactly unity, but not 
    #quite exact due to finite sampling and truncation at +- 3 sigma.  
    #Next line adjusts to make coeffs b(k) sum to exactly unity.
    b=b/np.sum(b)
    
    return b
    
    

utm22n = pyproj.Proj("+init=EPSG:32622")

files = ['lev8_2012_2013_geod.dat']
         
siteids = ['lev8'] 
start_years = [2012]     

loy = 366

start = dt.datetime(2012,5,1)   
   
nsites = len(siteids)
n = 1 

all_data = {}
 
#for siteid,f,start_year in zip(siteids,files,start_years):    
f = files[0]
siteid = siteids[0]
start_year = start_years[0]

print 'Loading ' + f + '...'
        
# Load data            
data = geod.load_geod('/scratch/s1144267/gps_final_velocities_spring2013/' + f,
                   start_year=start_year,columns=['Latitude','Longitude','RMS','SigH','#'])
                                    

print 'Deleting bad data'
# Remove bad data    
data[data.index < start] = np.nan
data[data['RMS'] > 60] = np.nan
data[data['SigH'] > 10] = np.nan
data[data['#'] == 0] = np.nan 
data.dropna()

data.fillna(method='ffill',inplace=True)
   
# Convert lat-lon to UTM
print 'Converting to lat lon'
xy = utm22n(data['Longitude'].tolist(),data['Latitude'].tolist())
data['utmx'] = xy[0]
data['utmy'] = xy[1]

# Remove outliers
print 'Removing outliers...'
its = 3
for n in range(0,its):
    print ' Iteration ' + str(n+1) + ' of ' + str(its)
    data['utmxs'] = signal.medfilt(data['utmx'],kernel_size=49)
    data['utmys'] = signal.medfilt(data['utmy'],kernel_size=49)
    

# Smooth the data
print 'Smoothing...'
data['utmxs'] = signal.filtfilt((gaussfiltcoef((1./10),(1./(7200)))),[1],data['utmxs'])
data['utmys'] = signal.filtfilt((gaussfiltcoef((1./10),(1./(7200)))),[1],data['utmys'])

# Daily differences
print 'Calculating daily velocities...'
daily = data.tz_localize('UTC').tz_convert('America/Godthab').resample('24h',how='first')
pts = np.array([daily['utmxs'],daily['utmys']]).T
differences = np.diff(pts,axis=0)
segdists = np.hypot(differences[:,0],differences[:,1])
segdists[segdists > 1000] = np.nan
v_24h = pd.Series(segdists * loy,index=daily.index[:-1])

# Epoch differences
print 'Calculating epoch velocities...'
data = data.tz_localize('UTC').tz_convert('America/Godthab')
differences = np.diff(np.array([data['utmxs'],data['utmys']]).T,axis=0)
segdists = np.hypot(differences[:,0],differences[:,1])
segdists[segdists > 100] = np.nan
v = pd.Series(segdists * (6*60*24*loy),index=data.index[:-1])


    
    # Convert to local coord system
    #data['utmx'] -= data['utmx'][start].mean()
    #data['utmy'] -= data['utmy'][start].mean()
    
    
    
    
    

    
    
    
    
    