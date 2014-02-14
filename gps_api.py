# -*- coding: utf-8 -*-
"""
gpsm -- Helper functions/API to final Leverett GPS velocity MAT files, utilising 
pandas library.

Created on Wed Nov 28 12:20:22 2012

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import numpy as np

try:
    from pandas import DataFrame
    import pandas as pd
    import matplotlib.pyplot as plt
except:
    print 'Could not load pandas or matplotlib'
import scipy.io
import calendar as cal
import math
from copy import deepcopy
from ajtuseful import *
from doy import *

class gpsm:
    """ GPS manipulation functionality - load data, re-index timestamps ...
    
    load_vel : load single receivers velocity into pandas Series
    load_multi_vel : load multiple receivers velocity into pandas DataFrame
    reindex_year_doy : convert incremental DOY to year-doy
    reindex_datetime : convert year-doy to python datetimes
    
    """
    
    def __init__(self):
        self.d = None
        self.tmp = None
        self.subset_year = None
        

    def load_vel(self,filename,name,ds='v_24h',reindex=True,year=False):
        """ Load specific dataset (ds) from GPS velocity MAT file with filename.
        
        Parameters
        ----------
        filename : path and name of MAT file from which to load data.
        name : Name of Series.
        ds : name of dataset, i.e. v_24h or v_6h.
        reindex : if True, reindex to hierarchical year,doy index.
        year : if not false, return just data for year specified. Ignored if 
        reindex=False.
        
        Returns
        -------
        Velocity data as pandas Series object in self.d.
        
        """
        data = scipy.io.loadmat(filename,variable_names=[ds])
        
        if ds == 'xyzt':
            ts_idx = 3
            ds_idx = 2
        else:
            ts_idx = 0
            ds_idx = 1
            
        if reindex == True:
            index = self.reindex_year_doy(data[ds][:,ts_idx])
            asSeries = pd.Series(data=data[ds][:,ds_idx],index=[index[0],index[1]],name=name)
        else:
            asSeries = pd.Series(data=data[ds][:,ds_idx],index=data[ds][:,ts_idx])
            
        if year != False and reindex == True:
            print 'Extracting ' + str(year)
            asSeries = asSeries[year]
            
            years = [int(year)] * len(asSeries)
            new_idx = pd.Index([years,list(asSeries.index)])
            self.tmp = new_idx
            self.subset_year = year
        
            #asSeriesSubNewIdx = pd.Series(list(asSeriesSub),index=new_idx)
            #asSeries = asSeriesSubNewIdx
      
        self.d = asSeries
        
        
    def load_multi_vel(self,filenames,names,ds='v_24h',reindex=True,year=False):
        """ Load specified datasets into one pandas data frame.
        
        Parameters
        ----------
        filenames : list of path and name of MAT file from which to load data.
        name : list of names of Series.
        ds : name of dataset, i.e. v_24h or v_6h.
        reindex : if True, reindex to hierarchical year,doy index.
        year: if not false, return just data for year specified. Ignored if 
        reindex=False.
        
        Returns
        -------
        Velocity data as pandas DataFrame object in d.
        
        """
        store = {}
        
        for f,n in zip(filenames,names):
            print 'Loading ' + n
            self.load_vel(f,n,ds=ds,reindex=reindex,year=year)
            store[n] = deepcopy(self.d)
            print '...done.'
            print ' '
        self.d = pd.DataFrame(store)
               
    
    
    def reindex_year_doy(self,doy=None,start_year=2009,save=False):
        """ Convert incremental multi-year DOY to hierarchical year, DOY. 
        
        Parameters
        ----------
        doy : list of multi-year DOYs to convert. None to use self.d
        start_year : when the mult-year DOYs begin
        save : if True, apply to self.d, if False then return tuple
        Returns
        -------
        tuple : ( list(year), list(fractional doy) ) if save==True
        boolean : True to indicate completion if save==False
        
        Examples
        --------
        >>> reindex_year_doy([367.1,367.12])
        ( ())
        
        """
        if doy == None:
            doy = np.array(self.d.index)
        elif len(doy) < 2:
            doy = np.array([doy])
        else:
            doy = np.array(doy)
            
        #if not isiterable(doy):
        #    'Input data are not iterable (e.g. list or array).'
        #    raise ValueError
            
        years = np.zeros((doy.size))
        end_doy = doy[-1]
        end_year = 0
        current_year = start_year
        add_on = 1
        while 1:
            
            if not cal.isleap(current_year):
                add_on += 365
                if end_doy <= add_on:
                    end_year = current_year
                    break
            else:
                add_on += 366
                if end_doy <= add_on:
                    end_year = current_year
                    break
                
            current_year += 1
        
        #print 'End year' +  str(end_year)
            
        #end_year = start_year + int(math.floor(doy[-1] / 365))
        #leap_test = start_year + int(math.floor(doy[-1] / 366))
    
        #if leap_test < end_year:
        #    print 'leap_test true'
        #    end_year = leap_test
       # print end_year
        
        doy_count = 1
        reidx_doy = deepcopy(doy)
        for year in range(start_year,end_year+1):
            #print year
            if cal.isleap(year):
                add_on = 366
            else:
                add_on = 365
                
            if year == end_year:
                end = len(doy)
            else:
                end = np.where(doy >= doy_count+add_on)[0][0]
                
            start = np.where(doy >= doy_count)[0][0] 
            #print start
            #print end
            #print doy_count
            if doy_count + add_on < doy[start]:
                # skip year with no data - e.g. lev7/8, which are still indexed from 2009
                doy_count += add_on
                continue
            else:
                #print 'Found ' + str(year) + ' data, begins ' + str(start) + ', ends ' + str(end) + '.'
                if len(reidx_doy) > 1:
                    years[start:end] = int(year)
                    reidx_doy[start:end] -= doy_count-1
                else:
                    years = np.array([int(year)])
                    reidx_doy -= doy_count - 1
                doy_count += add_on
                
        if save == False: 
            if len(reidx_doy) > 1:               
                return ([years.tolist(),reidx_doy.tolist()]) 
            else:
                print 'here'
                return ([years[0],reidx_doy[0]])
        elif save == True:
            self.d.index = pd.Index([years.tolist(),reidx_doy.tolist()])
            return True
        else:
            print 'Unknown save option.'
            return False
        
        
    def reindex_datetime(self,doy=False,year=False):
        """ Reindex self.d with python datetimes.
        
        By default the function calculates new index based on existing year-doy
        hierarchical index. If this index is not present a ValueError will be
        returned.
        
        Parameters
        ----------
        If doy and year not specified, uses values of d.index.
        doy : Iterable of fractional dates of year (1-366 only)
        year : Iterable of years corresponding to each position in doy
        
        Returns
        -------
        True or False
        Reindexes self.d
        
        Examples
        --------
        >>> reindex_datetime()
        True
        
        >>> reindex_datetime(doy=[123.5,124.5],year=[2012,2012])
        True
        
        """
        if self.subset_year != None:
            # Data are read into object but have been subsetted to specific year
            year = [self.subset_year] * len(self.d)
            doy = self.d.index.tolist()
        elif doy == False and len(self.d.index[0]) == 2:
            # unzip hierarchical index
            # http://paddy3118.blogspot.co.uk/2007/02/unzip-un-needed-in-python.html
            aslist = self.d.index.tolist()
            year,doy = zip(*aslist)
            year = [int(round(y,0)) for y in year]
        else:
            print 'self.d.index is not yet converted to year-doy index.'
            raise ValueError
        
        if (doy != False and year == False) or (doy == False and year != False):
            print 'Both doy and year must be specified.'
            raise ValueError

                               
        index_as_dt = doy2datetime(year,doy)
        self.d.index = pd.Index(index_as_dt)
 

class Geod:
    
    def __init__(self):
        self.keys = ['YY','DOY','Seconds','Latitude','Longitude','Height','SigN',
                    'SigE','SigH','RMS','Hash','Atm','AtmPlusMinus','Fract DOY','Epoch',
                    'HashBF','NotF','North','East','Up','Set']
    
              
    def load_geod(self,filename,start_year,columns=False):
        """ Load in GPS GEOD file.
        
        Provide filename of file and start_year of data (to allow indexing of
        timestamp).
        
        Optionally, provide a list of column labels to read in - the column
        labels can be found in Geod.keys.
        
        """
        if columns == False:
            data = pd.read_csv(filename,sep=' ',header=None,names=self.keys)
        else:
            # Ensure we get Fract DOY so that we can complete indexing.
            if 'Fract DOY' not in columns: columns.append('Fract DOY')
            data = pd.read_csv(filename,sep=' ',header=None,names=self.keys,
                           usecols=columns)                   
        
        # Make time series index  
        api = gpsm()                 
        reindexed_doy = api.reindex_year_doy(doy=data['Fract DOY'].tolist(),start_year=start_year,save=False)
        as_datetime = doy2datetime(reindexed_doy[0],reindexed_doy[1])
        data.index = as_datetime
        # Drop the fract doy, no longer needed
        data = data.drop('Fract DOY',1)
        columns.remove('Fract DOY')
        
        return data
        
        