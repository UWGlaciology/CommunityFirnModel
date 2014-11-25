'''
Created on Oct 23, 2013

@author: Jessica
'''

import firnmodel_v03_firnmice as frn
import os
import logging
import json

#thenames=["Barnola1991","Arthern2010"]
thenames=["Arthern2010T"]

for f in thenames:
    
  #   if f[0] == ".":
#         continue
    
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
    
#     jsonFile = open("config_test_input2.json", "r")
#     data = json.load(jsonFile)
#     jsonFile.close()
# 
#     tmp = data["physRho"]
#     data["physRho"] = f
#     data["resultsFolder"] = f
# 
#     jsonFile = open("config_test_input2.json", "w+")
#     jsonFile.write(json.dumps(data))
#     jsonFile.close()
    
    configName = "config_test_input3.json"
    frn.runModel(configName,True)
    frn.runModel(configName,False) 
