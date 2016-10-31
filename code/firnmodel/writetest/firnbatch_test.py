import firnmodel_v05 as frn
import os
import logging
import json
import numpy as np
import sys

# thenames=["Barnola1991","Arthern2010T","Arthern2010S","HLdynamic","HLSigfus","Li2004","Li2011","Helsen2008","Morris2014"]
# thenames=["HLdynamic","HLSigfus","Helsen2008","Arthern2010S","Arthern2010T","Simonsen2013","Ligtenberg2011","Barnola1991","KuipersMunneke2015","Li2011"]
# thenames=["Helsen2008","Arthern2010S","Arthern2010T","Simonsen2013","Ligtenberg2011","Barnola1991","KuipersMunneke2015","Li2011"]
# thenames=["Arthern2010S"]

# 
# for mm in thenames:
#     print "mm=",mm
# 
#     jsonFile = open("genericDC.json", "r")
#     data = json.load(jsonFile)
#     jsonFile.close()
# 
#     re='domeC_' + mm
# 
# 
# #     tmp = data["physRho"]
#     data["physRho"] = mm
#     data["resultsFolder"] = re
# 
# 
# 
#     jsonFile = open("genericDC.json", "w+")
#     jsonFile.write(json.dumps(data,sort_keys=True, indent=4, separators=(',', ': ')))
# #     ,sort_keys=True, indent=4, separators=(',', ': '))
#     jsonFile.close()
configName = sys.argv[1]
print configName
frn.runModel(configName,True)
frn.runModel(configName,False)
    
        
        
