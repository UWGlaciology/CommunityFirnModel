import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.sparse.linalg as splin
from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr
from scipy.integrate import cumtrapz
import math
import ModelParameters.Gasses as MPG
import ModelParameters.Sites as MPS
import ModelParameters.Plotting as plots
import ModelParameters.Diffusivity as MPD
import ModelParameters.density as MPRHO
import csv
import json

# Set path to find all files to import and set output location
spot = os.path.dirname(sys.argv[0]) #Add Folder
print spot 
# os.chdir(spot) #change pwd to location of firnmodel.py
sys.path.insert(1,spot)
ResultsPlace=os.path.join(spot,'Results')
DataPathUser=os.path.join(spot,'DataImport/user_input')

def ClimateHistoryGenerator(Climate_config):

    GasPath = os.path.join(DataPathUser, 'GasHistory.csv') #Save location.
    
    with open(Climate_config, "r") as f:
        json_data=f.read()
        cc = json.loads(json_data)
        
    f_clim=np.loadtxt(os.path.join(DataPathUser,'Clim.csv'),delimiter=',',skiprows=0)
    
    time=f_clim[:,0]
    GAS1=np.ones(len(time))
    GAS2=np.ones(len(time))
    
    with open(GasPath, "w") as f:
        writer = csv.writer(f)
        writer.writerow(time)
        writer.writerow(GAS1)
        writer.writerow(GAS2)
    
    
if __name__ == "__main__":
    
    Climate_config = os.path.join(os.path.dirname(__file__),'ClimateConfig.json')
    
    ClimateHistoryGenerator(Climate_config)