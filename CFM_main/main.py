#!/usr/bin/env python
"""
main.py
=======
The core file to run the CFM from the command line. This file creates classes 
for the spin up and main runs, and it runs the model. This can be bypassed if
you write your own script to call firn_density_spin and firn_density_nospin.py

Copyright © 2021 C. Max Stevens

Distributed under terms of the MIT license.
"""

import sys
import os
# from firn_density_spin import FirnDensitySpin
from firn_density_nospin import FirnDensityNoSpin
import time
import json
import shutil
import RCMpkl_to_spin as RCM

__author__ = "C. Max Stevens, Vincent Verjans, Brita Horlings, Annika Horlings, Jessica Lundin"
__license__ = "MIT"
__version__ = "1.1.0"
__maintainer__ = "Max Stevens"
__email__ = "maxstev@umd.edu"
__status__ = "Production"


if __name__ == '__main__':

    if len(sys.argv) >= 2:
        configName = os.path.join(os.path.dirname(__file__), sys.argv[1])
        print(configName)
    else:
        print('No .json configuration file specified. Exiting.')
        sys.exit()
        # configName = os.path.join(os.path.dirname(__file__), 'generic.json')

    with open(configName, "r") as f:
        jsonString = f.read()
        c = json.loads(jsonString)

    tic=time.time()

    print("")
    print("-----------------------------------------------------------------------")
    print("-----------------------------------------------------------------------")
    print("<<<<<<<< Running the Community Firn Model (CFM), Version %s >>>>>>>>" %__version__)
    print("<<<<<<<< Please cite your use:                                 >>>>>>>>")
    print("<<<<<<<< https://doi.org/10.5194/gmd-13-4355-2020              >>>>>>>>")
    print("<<<<<<<< Distributed under terms of the MIT license.           >>>>>>>>")
    print("<<<<<<<< Please consider telling us that you are using the CFM >>>>>>>>")
    print("<<<<<<<< (it helps to keep the project going!)                 >>>>>>>>")
    print("<<<<<<<< Questions/comments to Max Stevens: maxstev@umd.edu    >>>>>>>>")
    print("-----------------------------------------------------------------------")
    print("-----------------------------------------------------------------------")
    print("")
    
    if '-n' in sys.argv:
        NewSpin = True
    elif 'NewSpin' in c:
        NewSpin = c['NewSpin']
    else:
        NewSpin = False

    if 'input_type' not in c:
        c['input_type'] = "csv"

    if c['input_type'] == 'dataframe':
        pkl_name = os.path.join(c['InputFileFolder'],c['DFfile'])
        timeres = c['DFresample']
        climateTS, stepsperyear, depth_S1, depth_S2, desired_depth = RCM.makeSpinFiles(pkl_name,timeres = timeres, melt = c['MELT'])
    else:
        climateTS = None

    firn = FirnDensityNoSpin(configName, climateTS = climateTS, NewSpin = NewSpin)
    firn.time_evolve()

    shutil.copy(configName,c['resultsFolder'])
    
    print('run time =' , time.time()-tic , 'seconds')
