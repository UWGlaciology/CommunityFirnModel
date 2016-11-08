
import firnmodel_v03 as frn
import os
import logging
import json

thenames=["Barnola1991","Arthern2010T","Arthern2010S","HLdynamic","HLSigfus","Li2004","Li2011","Helsen2008","Morris2014"]
# thenames=["Arthern2010T"]

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
    
    jsonFile = open("config_test_input_AGU.json", "r")
    data = json.load(jsonFile)
    jsonFile.close()
    
    re="AGUresults/" + f
    
    tmp = data["physRho"]
    data["physRho"] = f
    data["resultsFolder"] = re

    jsonFile = open("config_test_input_AGU.json", "w+")
    jsonFile.write(json.dumps(data))
    jsonFile.close()
    
    configName = "config_test_input_AGU.json"
    frn.runModel(configName,True)
    frn.runModel(configName,False) 
