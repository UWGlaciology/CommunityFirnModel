import os
import numpy as np
import sys
from shutil import copyfile
import sys
# from string import join
from firn_density_spin import FirnDensitySpin
from firn_density_nospin import FirnDensityNoSpin
import time
import json

'''
to run: parallel –a varden.txt –a modellist_ens.txt python firnbatch_RACMO_summit_ens.py

'''

# thenames=["Barnola1991","Arthern2010T","Arthern2010S","HLdynamic","HLSigfus","Li2004","Li2011","Helsen2008","Morris2014"]
# thenames=["HLdynamic","HLSigfus","Helsen2008","Arthern2010S","Arthern2010T","Simonsen2013","Ligtenberg2011","Barnola1991","KuipersMunneke2015","Li2011","Goujon2003"]
# thenames=["HLSigfus","Helsen2008","Arthern2010S","Arthern2010T","Simonsen2013","Ligtenberg2011","Barnola1991","KuipersMunneke2015","Li2011"]
# thenames=["Morris2014"]

# thefiles=["constant600","loop600","randspin600_1","randspin600_2","randspin600_3","randspin600_4","randspin600_5","randspin600_6","randspin600_7","randspin600_8","randspin600_8","randspin600_10"]
# thefiles=["randspin600_1","randspin600_2","randspin600_3","randspin600_4","randspin600_5","randspin600_6","randspin600_7","randspin600_8","randspin600_9","randspin600_10"]
# thefiles=["randspin600_11","randspin600_12","randspin600_13","randspin600_14","randspin600_15","randspin600_16","randspin600_17","randspin600_18","randspin600_19","randspin600_20"]

nn=sys.argv[1] # variable density choice
mm=sys.argv[2] # physics/model choice
# print nn
connm='RACMOjson/RACMO_Summit_config_'+nn+'_'+mm+'_ens.json'
copyfile('RACMO_Summit_config_master_varden.json', connm)

# for mm in thenames:
    
#     print "mm=",mm
# print connm
jsonFile = open(connm, "r")
data = json.load(jsonFile)
jsonFile.close()

re="RACMOresults_ens_all/r" + nn + "/" + mm #results folder
# tein="inputdata/RACMOinput/Summit/Summit_tskin_RACMO_" + nn + ".csv"
# smbin="inputdata/RACMOinput/Summit/Summit_smb_RACMO_" + nn + ".csv"

#     tmp = data["physRho"]
data["physRho"] = mm
data["resultsFolder"] = re
# data["InputFileNameTemp"] = tein
# data["InputFileNamebdot"] = smbin
data["srho_type"] = nn

if mm=="Arthern2010T":
    data["physGrain"] = True
else:
    data["physGrain"] = False

jsonFile = open(connm, "w+")
jsonFile.write(json.dumps(data,sort_keys=True, indent=4, separators=(',', ': ')))
#     ,sort_keys=True, indent=4, separators=(',', ': '))
jsonFile.close()
#     print re
#     os.mkdir(re)
configName = connm

firnS = FirnDensitySpin(configName)
firnS.time_evolve()

firn = FirnDensityNoSpin(configName)
firn.time_evolve()
    
        
        
