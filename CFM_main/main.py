#!/usr/bin/env python
'''This file creates classes for the spin up and main runs and runs the model''' 

import sys
import os
from firn_density_spin import FirnDensitySpin
from firn_density_nospin import FirnDensityNoSpin
import time
import json
import shutil

__author__ = "C. Max Stevens, Vincent Verjans, Brita Horlings, Annikah Horlings, Jessica Lundin"
__license__ = "MIT"
__version__ = "1.0.5"
__maintainer__ = "Max Stevens"
__email__ = "maxstev@uw.edu"
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
    print("<<<<<<<< Developed at the University of Washington             >>>>>>>>")
    print("<<<<<<<< Please cite your use!                                 >>>>>>>>")
    print("<<<<<<<< Distributed under terms of the MIT license.           >>>>>>>>")
    print("-----------------------------------------------------------------------")
    print("-----------------------------------------------------------------------")
    print("")
    
    if os.path.isfile(c['resultsFolder']+'/'+c['spinFileName']) and '-n' not in sys.argv:
        print('Skipping Spin-Up run;', c['resultsFolder']+'/'+c['spinFileName'], 'exists already')
        try:
            os.remove(c['resultsFolder']+'/'+c['resultsFileName'])
            print('deleted', c['resultsFolder']+'/'+c['resultsFileName'])
        except:
            pass
        firn = FirnDensityNoSpin(configName)
        firn.time_evolve()

    else:

        firnS = FirnDensitySpin(configName)
        firnS.time_evolve()

        firn = FirnDensityNoSpin(configName)
        firn.time_evolve()

    shutil.copy(configName,c['resultsFolder'])
    
    print('run time =' , time.time()-tic , 'seconds')
