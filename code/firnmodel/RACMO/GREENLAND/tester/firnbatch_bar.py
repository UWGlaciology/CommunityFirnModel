
import firnmodel_v04_racmo as frn
import os
import logging
import json


for mm in range(11):
    print "mm=",mm
    for nn in range(11):
    
        jsonFile = open("racmo_config_bar.json", "r")
        data = json.load(jsonFile)
        jsonFile.close()
    
        resultsfolder='bar_%s_%s' % (mm,nn)
        in_temp='t2m_%s_%s.csv' % (mm,nn)
        in_smb='smb_%s_%s.csv' % (mm,nn)
    
        #re="AGUresults/" + f
    
    #     tmp = data["physRho"]
        data["InputFileNameTemp"] = in_temp
        data["InputFileNamebdot"] = in_smb
        data["resultsFolder"] = resultsfolder

        jsonFile = open("racmo_config_bar.json", "w+")
        jsonFile.write(json.dumps(data))
        jsonFile.close()
    
        configName = "racmo_config_bar.json"
        frn.runModel(configName,True)
        frn.runModel(configName,False) 
