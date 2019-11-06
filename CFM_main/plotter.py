#!/usr/bin/env python

# Copyright (C) 2019 C. Max Stevens <maxstev@uw.edu>
#
# Distributed under terms of the MIT license.

import matplotlib.pyplot as plt 
import h5py as h5
import os
import sys
import numpy as np
import scipy as sp
import pickle
import seaborn as sns 
sns.set()

'''
Example of plotting results from the CFM.
'''

rfolder = '/PATH/TO/RESULTS/FOLDER' # alter this to point to the results folder.
rfile = 'CFMresults.hdf5'
fn = os.path.join(rfolder,rfile)
f = h5.File(fn,'r')

timesteps = f['depth'][1:,0]
stps = len(timesteps)
depth = f['depth'][1:,1:]
density = f['density'][1:,1:]
temperature = f['temperature'][1:,1:]
dip_all = f['DIP'][:,:]
f.close()

f1,a1 = plt.subplots()
a1.plot(density[-1,:],depth[-1,:])
a1.invert_yaxis()
a1.grid(True)
f1.savefig('Example_DepthDensity.png')

