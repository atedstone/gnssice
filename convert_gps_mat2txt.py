# -*- coding: utf-8 -*-
"""
Convert GPS matlab data files to text files.

Created on Thu Feb 21 11:42:24 2013

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import numpy as np
import scipy.io
import gc

vel_f = ['lev0vel_2009_autumn2012_UTC-2_20121116_opt.mat',
         'lev1vel_2009_autumn2012_UTC-2_20121122_opt.mat',
         'lev2vel_2009_autumn2012_UTC-2_20121116_opt.mat',
         'lev3vel_2009_autumn2012_UTC-2_20121116_opt.mat',
         'lev4vel_2009_autumn2012_UTC-2_20121116_opt.mat',
         'lev5vel_2009_autumn2012_UTC-2_20121116_opt.mat',
         'lev6vel_2009_autumn2012_UTC-2_20121116_opt.mat']
         
for filename in vel_f:
    
    data = scipy.io.loadmat(filename)
    
    # Extract dataset and cull the random extra column of zeros
    v_6h = data['v_6h'][:,0:2]
    
    filename_noext = filename[0:-4]
    save_fname = filename_noext + '_6h.txt'
    np.savetxt(save_fname,v_6h,fmt=('%.6f','%.3f'))
    
    data = []
    gc.collect()
             
         
         
