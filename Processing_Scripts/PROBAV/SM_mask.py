# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 09:22:53 2018

@author: tih
"""

import h5py
import numpy as np
L2_file = ‘C:\PROBAV_L2A_20150506_085613_3_1KM_V001.HDF5’
h5 = h5py.File (L2_file,'r+')
SM = h5['/LEVEL2A/QUALITY/SM'].value
#Evaluate the three least significant bits for ‘clear’, ‘shadow’,
#‘undefined’, ‘cloud’, and ‘snow/ice’ and assign the outcome to variables
clr = np.where((SM&1 == 0) & (SM&2 == 0) & (SM&4 == 0))
shw = np.where((SM&1 != 0) & (SM&2 == 0) & (SM&4 == 0))
und = np.where((SM&1 == 0) & (SM&2 != 0) & (SM&4 == 0))
cld = np.where((SM&1 != 0) & (SM&2 != 0) & (SM&4 == 0))
ice = np.where((SM&1 == 0) & (SM&2 == 0) & (SM&4 != 0))
