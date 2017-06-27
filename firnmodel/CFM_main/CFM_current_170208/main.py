import sys
import os
from string import join
from firn_density_spin import FirnDensitySpin
from firn_density_nospin import FirnDensityNoSpin
import time
import json

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        configName = os.path.join(os.path.dirname(__file__), sys.argv[1])
    else:
        configName = os.path.join(os.path.dirname(__file__), 'generic.json')

    with open(configName, "r") as f:
        jsonString = f.read()
        c = json.loads(jsonString)

    # if '-s' in sys.argv:
    #     spin = 'on'
    #     print "spin on"
    # else:
    #     spin = 'off'
    #     print "spin off"

    # configSpin = {
    #     'on'  : FirnDensitySpin,
    #     'off' : FirnDensityNoSpin,
    # }
    
    tic=time.time()
    # try:
    # firn = configSpin[spin](configName)
    # except KeyError:
        # sys.exit("Error")
    print ""
    print "<<<<<<<< Running the Community Firn Model (CFM), Version 1.0 >>>>>>>>"
    print "<<<<<<<< Developed at the University of Washington           >>>>>>>>"
    print "<<<<<<<< Please cite your use!                               >>>>>>>>"
    print ""
    
    if os.path.isfile(c['resultsFolder']+'/'+c['spinFileName']) and '-n' not in sys.argv:
        print 'Skipping Spin-Up run;', c['resultsFolder']+'/'+c['spinFileName'], 'exists already'
        try:
            os.remove(c['resultsFolder']+'/'+c['resultsFileName'])
            print 'deleted', c['resultsFolder']+'/'+c['resultsFileName']
        except:
            pass
        firn = FirnDensityNoSpin(configName)
        firn.time_evolve()

    else:

        firnS = FirnDensitySpin(configName)
        firnS.time_evolve()

        firn = FirnDensityNoSpin(configName)
        firn.time_evolve()
    
    print 'run time =' , time.time()-tic , 'seconds'
