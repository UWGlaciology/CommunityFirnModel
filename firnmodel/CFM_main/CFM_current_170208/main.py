import sys
import os
from string import join
from firn_density_spin import FirnDensitySpin
from firn_density_nospin import FirnDensityNoSpin
import time

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        configName = os.path.join(os.path.dirname(__file__), sys.argv[1])
    else:
        configName = os.path.join(os.path.dirname(__file__), 'generic.json')

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

    firnS = FirnDensitySpin(configName)
    firnS.time_evolve()

    firn = FirnDensityNoSpin(configName)
    firn.time_evolve()
    
    print 'run time =' , time.time()-tic , 'seconds'
