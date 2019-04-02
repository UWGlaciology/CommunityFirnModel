import os
import numpy as np
import sys
from shutil import copyfile, move
import sys
# from string import join
from firn_density_spin import FirnDensitySpin
from firn_density_nospin import FirnDensityNoSpin
import time
import json


'''
to run: parallel –a inputlist_ens_short.txt –a modellist_ens.txt -a sitelist_short.txt python firnbatch_RACMO_summit_ens.py
'''

dtype = 'RACMO'
# site = 'DYE2'

nn=sys.argv[1] # run ID
# nn='0'
mm=sys.argv[2] # physics
site = sys.argv[3]
# mm='HLdynamic'
# print nn
jdir = 'jsonstore/' + dtype + '/' + site
try:
	os.makedirs(jdir)
except:
	pass

connm=os.path.join(jdir,dtype+'_'+site+'_config_'+nn+'_'+mm+'.json')
# copyfile(dtype+'_'+site+'_config_master.json', connm)
# if nn=='daily':
copyfile('RACMO_firncover_config_master_daily.json', connm)
# else:
# copyfile('RACMO_firncover_config_master.json', connm)

jsonFile = open(connm, "r")
data = json.load(jsonFile)
jsonFile.close()

re=dtype+"results_ens_all190402/"+site+"/r" + nn + "/" + mm #results folder
tein=site+"_tskin_"+dtype+"_" + nn + "_s.csv"
smbin=site+"_acc_"+dtype+"_" + nn + "_s.csv"
meltin=site+"_melt_"+dtype+"_" + nn + "_s.csv"

#     tmp = data["physRho"]
data["InputFileFolder"] = os.path.join("inputdata", dtype+'input',site)
data["physRho"] = mm
data["resultsFolder"] = re
data["InputFileNameTemp"] = tein
data["InputFileNamebdot"] = smbin
# try:
data["InputFileNamemelt"] = meltin


if nn=='daily':
	data["stpsPerYearSpin"]=365.0
	data["stpsPerYear"]=365.0
	data["HbaseSpin"]=2970.0
else:
	data["stpsPerYearSpin"]=12.0
	data["stpsPerYear"]=12.0
	data["HbaseSpin"]=2800.0
# except:
# 	print("cannot find melt; model runs without melt")
# 	data["MELT"] = False

if mm=="Arthern2010T":
    data["physGrain"] = True
else:
    data["physGrain"] = False

if mm=='Helsen2008':
	data['MELT']=False
elif mm=='Li2011':
	data['MELT']=False
elif mm=='Goujon2003':
	data['MELT']=False

if site=='Summit':
	data["rhos0"] = 300.0
else:
	data["rhos0"] = 330.0

# if site=='CRAWFORD':
# 	data['MELT']=False
# else:
# 	data["MELT"] = True

jsonFile = open(connm, "w+")
jsonFile.write(json.dumps(data,sort_keys=True, indent=4, separators=(',', ': ')))
#     ,sort_keys=True, indent=4, separators=(',', ': '))
jsonFile.close()
#     print re
#     os.mkdir(re)
configName = connm
# try:
firnS = FirnDensitySpin(configName)
firnS.time_evolve()

firn = FirnDensityNoSpin(configName)
firn.time_evolve()

move(connm,os.path.join(re,dtype+'_'+site+'_config_'+nn+'_'+mm+'_ens.json'))
# except:
# 	print('!!!!! error', site, nn, mm)



    
        
        
