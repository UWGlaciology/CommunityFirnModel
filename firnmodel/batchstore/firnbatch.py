'''
Created on Oct 23, 2013

@author: Jessica
'''

import firnmodel as frn
import os
import logging

for f in os.listdir('configs'):
    
    if f[0] == ".":
        continue
    
    #logging.basicConfig(filename='RUNDETAILS.log',level=logging.DEBUG,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    #console = logging.StreamHandler()
    #console.setLevel(logging.INFO)
    ## set a format which is simpler for console use
    #formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    ## tell the handler to use this format
    #console.setFormatter(formatter)
    ## add the handler to the root logger
    #logging.getLogger('').addHandler(console)
    
    print f

    configName = os.path.join('configs', f)
    frn.runModel(configName,True)
    frn.runModel(configName,False) 
