'''
UW Community Firn Model: Firn Evolution Module
@authors: Jessica Lundin, Max Stevens, Paul Harris, Will Leahy, Michael Yoon, Huong Vo, Ed Waddington 
'''

import csv
import json
import sys
import math
import numpy as np
import scipy.io
import logging
from scipy import interpolate
from scipy.sparse import spdiags
import scipy.sparse.linalg as splin
#from plot import plotData
from shutil import rmtree
import matplotlib.pyplot as plt
import os
from string import join
import shutil
import time
#import data_interp as IntpData
#import netCDF4 as nc

spot = os.path.dirname(sys.argv[0]) #Add Folder


def startlogger():
    logging.basicConfig(filename='RUNDETAILS.log',level=logging.DEBUG,filemode='w',format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    
def endlogger():
    console=logging.getLogger()
    console.handlers[0].stream.close()
    console.removeHandler(console.handlers[0])

def DefineGlobals():
    global sPerYear,rhoi,rhoiMgm,rho1,rho2,R,g,KtoC,atmosP
    sPerYear = 31557600.0
    rhoi = 917.0
    rhoiMgm = 0.917
    rho1 = 550.0
    rho2 = 815.0
    R = 8.314
    g = 9.8
    KtoC = 273.15
    atmosP = 101325

def HerronLangwayAnalytic(c,h,THL,AHL): #uses m w.e. a^-1
    """Model steady-state firn density and age profiles and bubble close-off
    config -- firnmod.config.Config
    return -- age, rho (density) for steady-state dynamics
    """
    #print 'AHL=%s' % AHL 
    hSize = np.size(h)      
    rhoc = 0.550
    rhos = c['rhos0']/1000.0
    # Eventually make this into 1 'if else' conditional statement
    #if (c['userInput']):
    #    A = input_bdot[0] * rhoiMgm
    #else:
        #A = c['bdot0']*rhoiMgm#* 0.917 'bdot0' is in ice equivalent (I think), but H&L runs using W.E. so we multiply by 0.917
    
    A=AHL*rhoiMgm
    ### Right now only takes the first index in input_temp
    ### Going to have to add more info to this... AKA make sure the step years actually get the right year in the array
    #if (c['userInput']):
    #    k0 = 11.0 * np.exp(-10160/(R * c['Ts0'])) #Max Changed this. Look for an easy way to get the first value from input file
    #    k1 = 575.0 * np.exp(-21400/(R * c['Ts0']))
    ## Herron and Langway k0 = (6a), k1 = (6b)
    #else:
    #    k0 = 11.0 * np.exp(-10160/(R*c['Ts0'] ))
    #    k1 = 575.0 * np.exp(-21400/(R*c['Ts0'] ))
    k0 = 11.0 * np.exp(-10160/(R*THL ))
    k1 = 575.0 * np.exp(-21400/(R*THL ))    
    
# depth of critical density, eqn 8 from Herron and Langway
    h0_55 = 1/(rhoiMgm*k0) * (np.log(rhoc/(rhoiMgm-rhoc))-np.log(rhos/(rhoiMgm-rhos)))
    Z0 = np.exp(rhoiMgm*k0*h + np.log(rhos/(rhoiMgm-rhos)))
    # The boundary from zone 1 to zone 2 = t0_55
    t0_55 = 1/(k0*A) * np.log((rhoiMgm-rhos)/(rhoiMgm-rhoc ))
    rho_h0 = (rhoiMgm* Z0)/(1+Z0)
    if np.max(rho_h0) >= rhoiMgm:
        t0 = np.zeros(hSize)
        for jj in xrange(hSize):
            if rho_h0[jj]<=rhoiMgm-0.001:
                t0[jj] = (1/(k0*A)*np.log((rhoiMgm-rhos)/(rhoiMgm-rho_h0[jj])))
                jj_max = jj
            else:
                t0[jj] = (t0[jj_max])
        
    else:
        t0 = 1/(k0*A)*np.log((rhoiMgm-rhos)/(rhoiMgm-rho_h0))
    
    Z1 = np.exp(rhoiMgm*k1*(h-h0_55)/np.sqrt(A) + np.log(rhoc/(rhoiMgm-rhoc)))
    Z = np.concatenate((Z0[h<h0_55], Z1[h>h0_55]))
    rho_h = (rhoiMgm * Z)/(1+Z)
    tp = np.ones(hSize)
    for j in xrange(hSize):
        if rho_h[j]<rhoiMgm-0.01:
            tp[j] = 1/(k1*np.sqrt(A)) * np.log((rhoiMgm-rhoc)/(rhoiMgm-rho_h[j]))+ t0_55
            jMax = j
        else:
            tp[j] = tp[jMax]
    # Zone 1 and Zone 2 repsectively
    age = np.concatenate((t0[h<h0_55], tp[h>h0_55]))*sPerYear
    rho = rho_h*1000
    return age, rho

def solver(a_U,a_D,a_P,b): #function for solving matrix problem
    nz=np.size(b)
    
    Diags = (np.append([a_U,-a_P],[a_D],axis=0))
    cols=np.array([1, 0, -1])
       
    big_A=spdiags(Diags,cols,nz,nz,format='csc')
    big_A=big_A.T
            
    rhs=-b
    phi_t=splin.spsolve(big_A,rhs)
    return phi_t

def transient_solve_TR(z_edges_vec,z_P_vec,nt,dt,Gamma_P,phi_0,nz_P,nz_fv,phi_s):
    '''transient 1-d diffusion finite volume method'''
    #this is solving the heat evolution in the firn
    phi_t =phi_0
    
    for i_time in xrange(nt): #why is this looping? 1/7/16 
        Z_P = z_P_vec

        dZ = np.concatenate(([1],np.diff(z_edges_vec),[1]))
    
        dZ_u = np.diff(Z_P)
        dZ_u = np.append(dZ_u[0], dZ_u)

        dZ_d = np.diff(Z_P)
        dZ_d = np.append(dZ_d,dZ_d[-1])

        f_u = np.append(0, (1 -(z_P_vec[1:] - z_edges_vec)/dZ_u[1:]))
        f_d = np.append(1 - (z_edges_vec - z_P_vec[0:-1])/dZ_d[0:-1], 0)


        Gamma_U = np.append(Gamma_P[0], Gamma_P[0:-1] )
        Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])

    
        Gamma_u =  1/ ( (1 - f_u)/Gamma_P + f_u/Gamma_U )
        Gamma_d =  1/ ( (1 - f_d)/Gamma_P + f_d/Gamma_D )
    
        S_C=0
        S_C=S_C*np.ones(nz_P)
    
        b_0 = S_C*dZ
    
        D_u = (Gamma_u / dZ_u)
        D_d = (Gamma_d / dZ_d)
        a_U = D_u 
        a_D = D_d 

        a_P_0 = dZ/dt
            
        a_P =  a_U + a_D + a_P_0
       
        bc_u_0 = phi_s
        bc_type = 1.
        bc_u   = np.concatenate(([ bc_u_0], [bc_type]))
        bc_d_0 = 0 
        bc_type = 2
        bc_d   = np.concatenate(([ bc_d_0 ], [ bc_type ]))
        b = b_0 + a_P_0*phi_t       
    
        #Upper boundary
        a_P[0] = 1 
        a_U[0] = 0
        a_D[0] = 0
        b[0]=bc_u[0]
    
        #Down boundary
        a_P[-1] = 1 
        #a_U[-1] = 0
        a_D[-1] = 0
        a_U[-1] = 1
        b[-1]=dZ_u[-1]*bc_d[0]
    
        phi_t = solver(a_U,a_D,a_P,b)
    #         logging.debug("i: %d, phi_t: %s", i_time, str(phi_t))
    #         phi.np.append(phi_t)
        a_P=a_U+a_D+a_P_0
    return phi_t

def runModel(configName,spin):
    "Runs firnmodel with an input json file."
    tic=time.time()
    print 'spin=',spin
    
#     startlogger()
    ### LOGGING ###
    logging.basicConfig(filename='RUNDETAILS.log',level=logging.DEBUG,filemode='w',format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    if spin:
       logging.info('Spin Up initiated')        
    elif not spin:
       logging.info('Model run initiated')            
    logging.info("Model configName = %s, spin = %r" % (configName, spin))
    #####
    
    DefineGlobals()
    
    with open(configName, "r") as f:
        
        jsonString = f.read()
        c = json.loads(jsonString)
    
    logging.info("The physics are %s" % c['physRho'])
    
    ResultsPlace=os.path.join(spot,'results',c['resultsFolder'])
    DataPath = os.path.join(spot,c['InputFileFolder'])
    
    ##### Import data (ice-core or synthetic) - must be .csv with first row year and second row temp/bdot
    # if c['InputDataType']=='synthetic':
    
    FIDtemp=os.path.join(DataPath,c['InputFileNameTemp'])
    data_temp=np.genfromtxt(FIDtemp, delimiter=',')
    
#         FIDtempspin=os.path.join(DataPath,'temp_spin.csv')
#         data_temp_spin=np.genfromtxt(FIDtempspin, delimiter=',')
    
    input_year_temp=data_temp[0,:]
    input_temp=data_temp[1,:]
    if input_temp[0]<0.0:
        input_temp=input_temp+273.15
#         input_temp_mean=np.mean(input_temp)
    
#         input_year_temp_spin=data_temp_spin[0,:]
#         input_temp_spin=data_temp_spin[1,:]
#         if input_temp_spin[0]<0.0:
#             input_temp_spin=input_temp_spin+273.15
#         input_temp_mean_spin=np.mean(input_temp_spin)
    
    input_temp_mean=np.mean(input_temp[0:10])
#         input_temp_mean=input_temp[0]

    
    FIDbdot=os.path.join(DataPath,c['InputFileNamebdot'])
    data_bdot=np.genfromtxt(FIDbdot, delimiter=',')
    
#         FIDbdotspin=os.path.join(DataPath,'smb_spin.csv')
#         data_bdot_spin=np.genfromtxt(FIDbdotspin, delimiter=',')
    
    input_year_bdot=data_bdot[0,:]
    input_bdot=data_bdot[1,:]
    input_bdot[input_bdot<0]=0.001 #change any negative smb to a small positive.
#         input_bdot_mean=np.mean(input_bdot)
    input_bdot_mean=np.mean(input_bdot[0:10])
#         input_bdot_mean=input_bdot[0]

#         input_year_bdot_spin=data_bdot_spin[0,:]
#         input_bdot_spin=data_bdot_spin[1,:]
#         input_bdot_spin[input_bdot_spin<0]=0.001 #change any negative smb to a small positive.
#         input_bdot_mean_spin=np.mean(input_bdot_spin)


    print "temp_mean=",input_temp_mean
    print "bdot_mean=",input_bdot_mean
   
    ##### import racmo data    
    # elif c['InputDataType']=='racmo':
    #     logging.info('racmo input')
    #     FIDtemp=os.path.join(spot,c['InputFileNameTemp'])
    #     FIDbdot=os.path.join(spot,c['InputFileNamebdot'])
    #     nc_s = scipy.io.netcdf.netcdf_file(FIDbdot, 'r')
    #     nc_t = scipy.io.netcdf.netcdf_file(FIDtemp, 'r')
    #     # for the SMB file:
    #     smb = nc_s.variables['smb'][:].copy()    
    #     lat_smb = nc_s.variables['lat'][:].copy()
    #     lon_smb = nc_s.variables['lon'][:].copy()
    #     time_smb = nc_s.variables['time'][:].copy()
    #     input_year_bdot=time_smb
    #     input_bdot1=smb[:,5,5]
    #     input_bdot=input_bdot1*365./0.917/1000 #convert to m WE per year.
    #     input_bdot_mean=np.mean(input_bdot)

    #     # temperature files:
    #     t2m = nc_t.variables['t2m'][:].copy()
    #     lat_t2m = nc_t.variables['lat'][:].copy()
    #     lon_t2m = nc_t.variables['lon'][:].copy()
    #     time_t2m = nc_t.variables['time'][:].copy()
    #     input_year_temp=time_t2m
    #     input_temp=t2m[:,5,5]
    #     if input_temp[0]<0.0:
    #         input_temp=input_temp+273.15
            
    #     input_temp_mean=np.mean(input_temp)
    #     nc_s.close()
    #     nc_t.close()
        
    #########

    
    if spin:    
        
        if os.path.exists(c['resultsFolder']):
            rmtree(c['resultsFolder'])
        os.makedirs(c['resultsFolder'])

        gridLen = int((c['H']-c['HbaseSpin'])/(input_bdot_mean/c['stpsPerYearSpin']))
        gridheight = np.linspace(c['H'],c['HbaseSpin'],gridLen)
        z = c['H']-gridheight
        dz = np.diff(z)
        dz = np.append(dz,dz[-1])
        dx = np.ones((gridLen)) #assume box width of 1 unit. To include stress must update this for horizontal longitudinal stress
        print "gridlen=",gridLen
#         if c['InputDataType']=='racmo':
        THL=input_temp_mean
        AHL=input_bdot_mean           
#         else:
#             THL=input_temp[0]
#             AHL=input_bdot[0]

        age, rho = HerronLangwayAnalytic(c,z,THL,AHL) # initial density profile in spin up
                       
        if c['physGrain']:
            #r2 = np.linspace(c['r2s0'], (6*c['r2s0']), gridLen)
            if c['calcGrainSize']:
                r02 = -2.42e-9*(c['Ts0'])+9.46e-7 #m^2, Gow 1967
                r2 = r02*np.ones(gridLen)
            else:
                r2 = np.linspace(c['r2s0'], (6*c['r2s0']), gridLen) 
        
        # initial temperature profile
        Tz = input_temp_mean*np.ones(gridLen)
        agez0 = age
        rhoz0 = rho
        z0 = z
        Tz0 = Tz
        if c['physGrain']:
            r20 = r2
#         if c["physRho"]=='Morris2014':
#             Hx=np.zeros(gridLen) #Hx is temperature history
#             Hx=Hx+np.exp(-QMorris/(R*Tz))*dt #is this dt correct?            
        ##### use specified spin up length
#         years = c['yearSpin'] #should not need user input years here - just specify in the json
#         years = len(input_year_bdot_spin)
        ##### use auto-spin up length
        zz=np.min(z[rho>850.0]) #minimum depth > 850
        years = int(zz/AHL)+100 #this spins up so that a new layer reaches the minimum depth where density is >850.
#         years = 800
#         years = 10
        stp = int(years *c['stpsPerYearSpin'])
#         print 'stp=',stp
        modeltime=np.linspace(0,years,stp+1)
        
        
#         yr_start=max(input_year_temp_spin[0],input_year_bdot_spin[0])
#         yr_end=min(input_year_temp_spin[-1],input_year_bdot_spin[-1])
#         #print yr_start
#         #print yr_end
#         years = (yr_end-yr_start)*1.0 #add *1.0 to make float
#         #years = input_year[-1] - input_year[0] + 1
#         #stp = int(years * (years / len(input_year)))
#         stp = int(years *c['stpsPerYearSpin']) #Make sure that stpsPerYear is set properly in config.
#         modeltime=np.linspace(yr_start,yr_end,stp+1)
        
        
        
        
        logging.info("Spin up time = %s" % years)
        #Dcon = c['D_surf']*np.ones(gridLen)

    else: #not spin
        
        densityPath = os.path.join(c['resultsFolder'], 'densitySpin.csv')
        tempPath =    os.path.join(c['resultsFolder'], 'tempSpin.csv')
        agePath =     os.path.join(c['resultsFolder'], 'ageSpin.csv')
        depthPath =   os.path.join(c['resultsFolder'], 'depthSpin.csv')
        if c['bdot_type'] == 'mean':
            bdotPath= os.path.join(c['resultsFolder'], 'bdot_meanSpin.csv')
        if c["physGrain"]:
            r2Path = os.path.join(c['resultsFolder'], 'r2Spin.csv')
        if c["physRho"]=='Morris2014':
            HxPath = os.path.join(c['resultsFolder'], 'HxSpin.csv')
        
        initDepth =   np.genfromtxt(depthPath, delimiter = ',')
        initAge =     np.genfromtxt(agePath, delimiter = ',')
        initDensity = np.genfromtxt(densityPath, delimiter = ',')
        initTemp =    np.genfromtxt(tempPath, delimiter = ',')
        if c["physGrain"]:
            initr2 =    np.genfromtxt(r2Path, delimiter = ',')
        if c["physRho"]=='Morris2014':
            initHx = np.genfromtxt(HxPath, delimiter = ',')
        if c['bdot_type'] == 'mean':
            initbdot_mean=np.genfromtxt(bdotPath,delimiter=',') 
        z   = initDepth[1:]
        gridLen = np.size(z)
        age = initAge[1:]
        rho = initDensity[1:]
        Tz  = initTemp[1:]
        if c["physGrain"]:
            r2=initr2[1:]
        if c["physRho"]=='Morris2014':
            Hx=initHx[1:]      
        dz = np.diff(z)
        dz = np.append(dz,dz[-1])
        print "dz=",dz
        dx = np.ones(gridLen)

        yr_start=max(input_year_temp[0],input_year_bdot[0])
        yr_end=min(input_year_temp[-1],input_year_bdot[-1])
        years = (yr_end-yr_start)*1.0 #add *1.0 to make float
        #years = input_year[-1] - input_year[0] + 1
        #stp = int(years * (years / len(input_year)))
        stp = int(years *c['stpsPerYear']) #Make sure that stpsPerYear is set properly in config.
        modeltime=np.linspace(yr_start,yr_end,stp+1)
        logging.info("Model run time is %f" % years)
        #else:
        #    years = c['years']
        #    stp = int(years *c['stpsPerYear']) #how many steps there are in the model run
        #    modeltime=np.linspace(0,years,stp)
        
        ##### TWRITE: need to clean this up or make it easier to specify
#         TWrite = np.concatenate((xrange(0,110,10),xrange(101,150,1),xrange(155,250,5),xrange(260,2010,10))) #set times at which to write data. This line is firnmice.
#         #TWrite = (np.arange(0,2005,5)) #set times at which to write data
#         TWrite=TWrite[1:]
        
        inte=1 # how often the data should be written
        TWrite = modeltime[inte::inte] # vector of times data will be written. Can customize here.
#         TWst=min(min(np.where(modeltime>1958.0)))
# #         TWrite = modeltime[TWst::inte] # vector of times data will be written. Can customize here. 
#         TWrite = np.concatenate((modeltime[0:TWst:12], modeltime[TWst::inte])) # vector of times data will be written. Can customize here.
        ######
        
        z0 = z
        agez0 = age
        rhoz0 = rho
        Tz0 = Tz
        if c["physGrain"]:
            r20 = r2
        Dcon = c['D_surf']*np.ones(gridLen) #initialize a grid of Diffusivity constant
        compaction_rate = np.zeros(gridLen-1)
        
    # mass and stress
    mass = rho*dz
    sigma = mass * dx * g
    sigma = sigma.cumsum(axis = 0)


    # define the time domain
    timTotal = years*sPerYear #total time of model run, seconds
    dt = timTotal/stp #time step size, seconds
    print "dt=",dt
    
    fT10m = interpolate.interp1d(z,Tz) #temp at 10m depth
    T10m = fT10m(10)

    if spin:
        if c['stpsPerYearSpin'] or (stp / years) >= 1.:
            Ts = input_temp[0]*np.ones(stp) #vector of temperatures through time
        
        else: #if T<1, introduce a seasonal temperature cycle
            TPeriod = c['yearSpin']
            #   Temperature calculation from Anais Orsi, supplementry material, AGU 2012
            #   http://www.agu.org/journals/gl/gl1209/2012GL051260/supplement.shtml
#             Ts = input_temp[0] + c['TAmp']*(np.cos(2*np.pi*np.linspace(0,TPeriod,stp))+0.3*np.cos(4*np.pi*np.linspace(0,TPeriod,stp)))
            Ts = input_temp_mean
        t = 1.0 / c['stpsPerYearSpin']
        T_mean = T10m
        

#         bdotSec0 = input_bdot[0]/sPerYear/c['stpsPerYearSpin']
        bdotSec0 = input_bdot_mean/sPerYear/c['stpsPerYearSpin']
        print "Ts=",Ts[0]
        print "bdotSec0=",bdotSec0
#         bdotSec0 = input_bdot[0]/sPerYear/(stp / years)
        #print input_bdot[0],stp,years
        bdotSec = bdotSec0*np.ones(stp)
#         bdotSec = input_bdot_spin/sPerYear/c['stpsPerYearSpin']
        
        mass_sum=mass.cumsum(axis = 0)
        bdot_mean = np.concatenate(([mass_sum[0]/(rhoi*sPerYear)],mass_sum[1:]/(age[1:]*rhoi/t)))
        #print 'bdotSec=',bdotSec[-1]
        rhos0=c['rhos0']*np.ones(stp)
        #D_surf=c['D_surf']*np.ones(stp)
        #### bdot mean method 2 #####
#         bdot_mean0 = bdotSec0 * np.ones(gridLen)
#         bdot_mean = bdot_mean0
         
    else: #not spin
        #t = 1.0 / (stp / years)
        t = 1.0/c['stpsPerYear']
        #print 't = %s' % t
        TPeriod = years

        #Ts = input_temp[:]
        #bdotSec = input_bdot[:]/sPerYear/(stp / years)
        
        Ts=np.interp(modeltime,input_year_temp,input_temp) #interpolate to model time
#         T_mean = np.mean(Ts)*np.ones_like(Ts) # 4/13/15 this might be wrong - but works as long at the Ts is actually the mean annual temp.
        T_mean = T10m
#         if t < 1.0: #add seasonal signal
#             Ts = Ts + c['TAmp']*(np.cos(2*np.pi*np.linspace(0,TPeriod,stp+1))+0.3*np.cos(4*np.pi*np.linspace(0,TPeriod,stp+1)))
        
        
        bdot=np.interp(modeltime,input_year_bdot,input_bdot) #bdot is m ice equivalent/year. multiply by 0.917 for W.E. or 917.0 for kg/year

        
        bdotSec = bdot/sPerYear/(stp / years) #this is the accumulation rate vector at each time step (i.e. if time step is 0.5 year, it is bdot rate in m I.E./(0.5years))
        print "Ts=",Ts[0],Ts[-1]
        print "bdotSec0=",bdotSec[0],bdotSec[-1]
        
        #3/23/15: bdotSec is ice equivalent. Any model that is water equivalent needs to be multiplied by 0.917
        
        #### b_dot mean method 2 ####
#         bdot_mean0 = bdotSec[0] * np.ones(gridLen) #first value in bdotSec is the spin up value (steady state)
#         bdot_mean = bdot_mean0
#         bdot_mean_sum = bdot_mean * age
#         #tmean_vec = np.ones(gridLen)
        mass_sum=mass.cumsum(axis = 0)
        if c['bdot_type'] == 'mean':
            bdot_mean=initbdot_mean[1:]
#             print 'spin bdot mean=',bdot_mean[1:6]
            
            
#             bdot_mean1 = np.concatenate(([mass_sum[0]/(rhoi*sPerYear)],mass_sum[1:]/(age[1:]*rhoi/t)))
#             print 'calc bdot mean=',bdot_mean1[1:6] 

    # Eventually want to get these also under 'user_input' as lines 4 and 5 of csv file
        rhos0=c['rhos0']
        rhos0 = rhos0*np.ones(stp)
        #rhos0 = rhos0*np.ones(stp) + np.random.randint(-100,150,stp)
            
        D_surf=c['D_surf']*np.ones(stp) #this is a diffusivity tracker: D_surf can be smaller to make layers with lower/higher diffusivity (some fraction of the number calculated using parameterization
    
    

    
    #print '422mass=', mass[0:5]
    #print '422age=',age[0:5]
    #print 'masssum=', mass_sum[0:2]
    #print 'bdotmean SPIN', bdot_mean[0:5]
    
    if not spin:
        spacewriteint=c['spacewriteint'] #how many of the spatial nodes to write to file.
        rho_time = np.append(modeltime[0],rho[0::spacewriteint])
        Tz_time = np.append(modeltime[0],Tz[0::spacewriteint])
        age_time = np.append(modeltime[0],age[0::spacewriteint]/sPerYear)
        z_time = np.append(modeltime[0],z[0::spacewriteint])
        comprate_time = np.append(modeltime[0],compaction_rate[0::spacewriteint])
        D_time = np.append(modeltime[0],Dcon[0::spacewriteint])
        stress_time = np.append(modeltime[0],sigma[0::spacewriteint])
        Clim_time = np.append(modeltime[0],[bdot[0],Ts[0]]) #not sure if bdot or bdotSec
        bdot_time = np.append(modeltime[0],bdot_mean[0::spacewriteint])
        if c["physGrain"]:
            r2_time = np.append(modeltime[0],r2[0::spacewriteint])
        if c["physRho"]=='Morris2014':
            Hx_time = np.append(modeltime[0],Hx[0::spacewriteint])
         
        '''initialize files to write in time loop'''
        densityPath = os.path.join(c['resultsFolder'], 'density.csv')
        tempPath = os.path.join(c['resultsFolder'], 'temp.csv')
        agePath = os.path.join(c['resultsFolder'], 'age.csv')
        depthPath = os.path.join(c['resultsFolder'], 'depth.csv')
        compratePath = os.path.join(c['resultsFolder'], 'comprate.csv')
        DconPath = os.path.join(c['resultsFolder'], 'Dcon.csv')
        stressPath = os.path.join(c['resultsFolder'], 'stress.csv')
        ClimPath = os.path.join(c['resultsFolder'], 'Clim.csv')
        bdotPath = os.path.join(c['resultsFolder'], 'bdot_mean.csv')
        if c["physGrain"]:
            r2Path = os.path.join(c['resultsFolder'], 'r2.csv')
        if c["physRho"]=='Morris2014':
            HxPath = os.path.join(c['resultsFolder'], 'Hx.csv')
            
        with open(densityPath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(rho_time)
        with open(tempPath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(Tz_time)
        with open(agePath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(age_time)
        with open(depthPath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(z_time)
        with open(compratePath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(comprate_time)        
        with open(DconPath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(D_time)
        with open(stressPath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(stress_time)            
        with open(ClimPath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(Clim_time)
        with open(bdotPath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(bdot_time)        
        if c["physGrain"]:
            with open(r2Path, "w") as f:
                writer = csv.writer(f)
                writer.writerow(r2_time) 
        if c["physRho"]=='Morris2014':
            with open(HxPath, "w") as f:
                writer = csv.writer(f)
                writer.writerow(Hx_time) 
                   

#        # initialize grain growth
#        if c['physGrain']:
##             r2 = np.linspace(c['r2s0'], (6*c['r2s0']), gridLen)
#            if c['calcGrainSize'] == 'True':
#                r02 = -2.42e-9*(c['Ts0'])+9.46e-7 #m^2, Gow 1967
#                r2 = r02*np.ones(gridLen)
#            else:
#                r2 = np.linspace(c['r2s0'], (6*c['r2s0']), gridLen) 
#            r2_time = np.concatenate(([0], r2))
#            r2Path = os.path.join(c['resultsFolder'], 'r2.csv')
#            with open(r2Path, "a") as f:
#                writer = csv.writer(f)
#                writer.writerow(r2_time)
                
        
        bcoAgeMartAll = []
        bcoDepMartAll = []
        bcoAge815All = []
        bcoDep815All = []
        LIZAgeAll = []
        LIZDepAll = []
        intPhiAll = []
        #tortAll = []
        
        ### Initial BCO,LIZ, and DIP ###
        #Vc = (6.95e-4)*T10m-0.043 #Martinerie et al., 1994, Eq. 2: critical pore volume at close off
        #bcoMartRho = ((Vc+1/(rhoi*(1e-3)))**-1)*1000 # Martinerie density at close off
        bcoMartRho = 1/( 1/(917.0) + T10m*6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        bcoAgeMart = min(age[rho>=bcoMartRho])/sPerYear # close-off age from Martinerie
        bcoDepMart = min(z[rho>=(bcoMartRho)])
        bcoAgeMartAll.append(bcoAgeMart) #age at the 815 density horizon
        bcoDepMartAll.append(bcoDepMart) #this is the 815 close off depth
        
        # bubble close-off age and depth assuming rho_crit = 815kg/m^3
        bcoAge815 = min(age[rho>=(rho2)])/sPerYear #close-off age where rho=815 kg m^-3
        bcoDep815 = min(z[rho>=(rho2)]) #depth of 815 horizon
        bcoAge815All.append(bcoAge815) #age at the 815 density horizon
        bcoDep815All.append(bcoDep815) #this is the 815 close off depth
        
        ### Lock-in depth and age
        LIZMartRho = bcoMartRho - 14.0 #LIZ depth (Blunier and Schwander, 2000)      
        LIZAgeMart = min(age[rho>LIZMartRho])/sPerYear #lock-in age
        LIZDepMart = min(z[rho>=(LIZMartRho)]) #lock in depth
        #phiLockIn=0.15
        #lockIn = min(age[phiClosed[phiClosed<phi]/phi[phiClosed<phi] >= phiLockIn])/sPerYear #lock-in age Old, from Jessica's original. Too shallow, likely.       
        LIZAgeAll.append(LIZAgeMart)
        LIZDepAll.append(LIZDepMart)
        
        ### Porosity, including DIP
        phi = 1-rho/rhoi # total porosity
        phi[phi<=0]=1e-16 
        phiC = 1-bcoMartRho/rhoi; #porosity at close off
        #phiClosed = np.ones(gridLen)
        phiClosed = 0.37*phi*(phi/phiC)**-7.6 #Closed porosity, from Goujon. See Buizert thesis (eq. 2.3) as well
        

                
        phiOpen = phi - phiClosed #open porosity
        phiOpen[phiOpen<=0] = 1.e-10 #don't want negative porosity.
                   
        intPhi = np.sum(phi * dz) #depth-integrated porosity
        intPhiAll.append(intPhi)
            
        bcoPath = os.path.join(c['resultsFolder'], 'BCO.csv')
        lidPath = os.path.join(c['resultsFolder'], 'LID.csv')
        intPhiPath = os.path.join(c['resultsFolder'], 'porosity.csv')
        


    del_s=-30.*np.ones(stp)
    del_z=del_s[0]*np.ones(gridLen)
    tic2=time.time()
    
    for ii in xrange(stp): #start main time-stepping loop
#         if spin:
# #             stp = int(years *c['stpsPerYearSpin'])
#             mtime=ii
        
#         if not spin:
        mtime=modeltime[ii] #placeholder for writing data.
#         if mtime in range(0,400,20):
#             print 'time=',mtime
        if spin:
            rho_old=rho
        #print 'ii=%s' % ii
                                              
        if c['physRho']=='HLdynamic': #uses m w.e. a^-1
            Q1=10160.
            Q2=21400.
            k1=11.0
            k2=575.0
            aHL=1.0
            bHL=0.5
#             print mtime
            A = bdotSec[ii]*(1/t)*sPerYear*rhoiMgm #A from the input json file is m ice equivalent.
#             print A
            drho_dt = np.zeros(gridLen)
            dr_dt=drho_dt
            dr_dt[rho<rho1] = k1*np.exp(-Q1/(R*Tz[rho<rho1]))*(rhoiMgm-rho[rho<rho1]/1000)*np.power(A,aHL)*1000/sPerYear
            drho_dt[rho<rho1] = dr_dt[rho<rho1]

            dr_dt[rho>=rho1] = k2*np.exp(-Q2/(R*Tz[rho>=rho1]))*(rhoiMgm-rho[rho>=rho1]/1000)*np.power(A,bHL)*1000/sPerYear
            drho_dt[rho>=rho1] = dr_dt[rho>=rho1]                    
        
        elif c['physRho']=='HLSigfus': #uses m w.e. a^-1
            Q1=10160.
            Q2=21400.
            k1=11.0
            k2=575.0
            aHL=1.0
            
            A = bdotSec[ii]*(1/t)*sPerYear*rhoiMgm
            drho_dt = np.zeros(gridLen)
            #if max(rho)>rho1:
            f550 = interpolate.interp1d(rho,sigma)
            sigma550 = f550(rho1)
            rhoDiff = (rhoiMgm-rho/1000)
            #rhoDiff[rhoDiff<=0] = 1.0e-16

            drho_dt[rho<rho1] = k1 * np.exp(-Q1/(R*Tz[rho<rho1]))*(rhoiMgm-rho[rho<rho1]/1000)*np.power(A,aHL)*1000/sPerYear
            k = np.power(k2*np.exp(-Q2/(R*Tz[rho>=rho1])),2)/sPerYear
            sigmaDiff = (sigma[rho>=rho1]-sigma550)
            #if rhoDiff[rho>=rho1]>0: #Max is not sure why this was coded in.
            drho_dt[(rho>=rho1) & (rho<rhoi)] = k * (sigmaDiff*rhoDiff[rho>=rho1]) / (g * np.log((rhoiMgm-rho1/1000)/(rhoDiff[(rho>=rho1) & (rho<rhoi)])))
            #else:
            drho_dt[(rho>=rho1) & (rho>=rhoi)] = 0
                        
        elif c['physRho'] =='Li2004': #Li and Zwally (2011), Equation from Arthern, 2010 (eq. 2): not sure where Rob got that? Uses m W.E./year for bdot
              #### Needs to have the vapor flux coded in if we want to use these physics properly.
            drho_dt = np.zeros(gridLen)
            #for j in xrange(gridLen):
            dr_dt = (rhoi-rho)*bdotSec[ii]*(1/t)*sPerYear*rhoiMgm*(139.21-0.542*T10m)*8.36*(KtoC-Tz)**-2.061
            drho_dt = dr_dt/sPerYear
        
        elif c['physRho'] =='Li2011': #uses m water equivalent per year. Note that some T in the paper uses K, some uses C.
            dr_dt = np.zeros(gridLen)
#             TmC=T_mean[ii]-273.15
            TmC=T10m-273.15
#             A = bdotSec[ii]*(1/t)*sPerYear*rhoiMgm
            A = np.mean(bdotSec)*(1/t)*sPerYear*rhoiMgm
            beta1 = -9.788 + 8.996 * A - 0.6165 * TmC # uses A, which is instant bdot. Needs to be switched to a long-term accumulation rate, especially to run small time steps. 
            beta2 = beta1/(-2.0178 + 8.4043 * A - 0.0932 * TmC) # Pay attention to units in paper: T is in deg C.
            dr_dt = np.zeros(gridLen)
            if c['bdot_type'] == 'instant':
                dr_dt[rho<=rho1] = (rhoi-rho[rho<=rho1])*A*beta1*8.36*(KtoC-Tz[rho<=rho1])**-2.061
                dr_dt[rho>rho1] = (rhoi-rho[rho>rho1])*A*beta2*8.36*(KtoC-Tz[rho>rho1])**-2.061
            elif c['bdot_type'] == 'mean':
                dr_dt[rho<=rho1] = (rhoi-rho[rho<=rho1])*  bdot_mean[rho<=rho1]*(1/t)*sPerYear*0.917*beta1*8.36*(KtoC-Tz[rho<=rho1])**-2.061
                dr_dt[rho>rho1] = (rhoi-rho[rho>rho1])*  bdot_mean[rho>rho1]*(1/t)*sPerYear*0.917*beta2*8.36*(KtoC-Tz[rho>rho1])**-2.061                            
            drho_dt = dr_dt/sPerYear        
                                                                          
        elif c['physRho'] =='Helsen2008': #Equation taken from Arthern et al., 2010. b_dot is in m WE per year. 
            #dr_dt = np.zeros(gridLen)
            if c['bdot_type'] == 'instant':
                dr_dt = (rhoi-rho) * bdotSec[ii]*0.917*(1/t)*sPerYear * (76.138-0.28965*T10m)*8.36*(KtoC-Tz)**-2.061
            elif c['bdot_type'] == 'mean':
                dr_dt = (rhoi-rho) * bdot_mean*0.917*(1/t)*sPerYear * (76.138-0.28965*T10m)*8.36*(KtoC-Tz)**-2.061
            drho_dt = dr_dt/sPerYear
        
        elif c['physRho'] =='Arthern2010S': #this is the steady-state solution. b_dot is in kg/year.
            ar1=0.07
            ar2=0.03
            Ec=60.0e3
            Eg=42.4e3
            dr_dt = np.zeros(gridLen)
            if c['bdot_type'] == 'instant':
                dr_dt[rho<rho1] = (rhoi-rho[rho<rho1])*ar1*(bdotSec[ii]*(1/t)*sPerYear)*917.0*g*np.exp(-Ec/(R*Tz[rho<rho1])+Eg/(R*T10m))
                dr_dt[rho>=rho1] = (rhoi-rho[rho>=rho1])*ar2*(bdotSec[ii]*(1/t)*sPerYear)*917.0*g*np.exp(-Ec/(R*Tz[rho>=rho1])+Eg/(R*T10m))            
            elif c['bdot_type'] == 'mean':
                dr_dt[rho<rho1] = (rhoi-rho[rho<rho1])*ar1*bdot_mean[rho<rho1]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho<rho1])+Eg/(R*T10m))
                dr_dt[rho>=rho1] = (rhoi-rho[rho>=rho1])*ar2*bdot_mean[rho>=rho1]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho>=rho1])+Eg/(R*T10m))
            
            drho_dt = dr_dt/sPerYear
            
        elif c['physRho']=='Arthern2010T': #full transient physics: uses stress
            kc1=9.2e-9
            kc2=3.7e-9
            Ec=60.0e3
            if not c['physGrain']:
                print "Grain growth should be on for Arthern Transient"           
            dr_dt = np.zeros(gridLen)
            dr_dt[rho<rho1] = kc1*(rhoi-rho[rho<rho1])*np.exp(-Ec/(R*Tz[rho<rho1]))*sigma[rho<rho1]/(r2[rho<rho1])
            dr_dt[rho>=rho1] = kc2*(rhoi-rho[rho>=rho1])*np.exp(-Ec/(R*Tz[rho>=rho1]))*sigma[rho>=rho1]/(r2[rho>=rho1])
            drho_dt = dr_dt#/sPerYear
            
        elif c['physRho'] =='Simonsen2013': # b_dot is in kg/year.
            ar1=0.07
            ar2=0.03
            Ec=60.0e3
            Eg=42.4e3
            F0=0.68 # firnmice value?
            F1=1.03 # firnmice value?
#             F0=0.8 # Simonsen's recommended (email correspondence)
#             F1=1.25 # Simonsen's recommended (email correspondence)
            
            dr_dt = np.zeros(gridLen)
            if c['bdot_type'] == 'instant':
                gamma=61.7/((bdotSec[ii]*(1/t)*sPerYear*917.0)**(0.5))*np.exp(-3800./(R*T10m))
                dr_dt[rho<rho1] = F0*(rhoi-rho[rho<rho1])*ar1*(bdotSec[ii]*(1/t)*sPerYear)*917.0*g*np.exp(-Ec/(R*Tz[rho<rho1])+Eg/(R*T10m))
                dr_dt[rho>=rho1] = F1*gamma*(rhoi-rho[rho>=rho1])*ar2*(bdotSec[ii]*(1/t)*sPerYear)*917.0*g*np.exp(-Ec/(R*Tz[rho>=rho1])+Eg/(R*T10m))
            elif c['bdot_type'] == 'mean':
                gamma=61.7/((bdot_mean[rho>=rho1]*(1/t)*sPerYear*917.0)**(0.5))*np.exp(-3800./(R*T10m))
                dr_dt[rho<rho1] = F0*(rhoi-rho[rho<rho1])*ar1*bdot_mean[rho<rho1]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho<rho1])+Eg/(R*T10m))
                dr_dt[rho>=rho1] = F1*gamma*(rhoi-rho[rho>=rho1])*ar2*bdot_mean[rho>=rho1]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho>=rho1])+Eg/(R*T10m))  
            drho_dt = dr_dt/sPerYear

        elif c['physRho'] =='Ligtenberg2011': #b_dot is in mm W.E. per year.
            ar1=0.07
            ar2=0.03
            Ec=60.0e3
            Eg=42.4e3
            dr_dt = np.zeros(gridLen)
            if c['bdot_type'] == 'instant':
                M_0=1.435-0.151*np.log(bdotSec[ii]*sPerYear*1e3*0.917)
                M_1=2.366-0.293*np.log(bdotSec[ii]*sPerYear*1e3*0.917)
                M_0=np.max(M_0,0.25)
                M_1=np.max(M_1,0.25)
                dr_dt[rho<rho1] = (rhoi-rho[rho<rho1])*M_0*ar1*bdotSec[ii]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho<rho1])+Eg/(R*T10m))
                dr_dt[rho>=rho1] = (rhoi-rho[rho>=rho1])*M_1*ar2*bdotSec[ii]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho>=rho1])+Eg/(R*T10m))                           
            elif c['bdot_type'] == 'mean':
                M_0=1.435-0.151*np.log(bdot_mean[rho<rho1]*sPerYear*917.0)
                M_1=2.366-0.293*np.log(bdot_mean[rho>=rho1]*sPerYear*917.0)
                M_0=np.max(M_0,0.25)
                M_1=np.max(M_1,0.25)
                dr_dt[rho<rho1] = (rhoi-rho[rho<rho1])*M_0*ar1*bdot_mean[rho<rho1]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho<rho1])+Eg/(R*T10m))
                dr_dt[rho>=rho1] = (rhoi-rho[rho>=rho1])*M_1*ar2*bdot_mean[rho>=rho1]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho>=rho1])+Eg/(R*T10m))
            drho_dt = dr_dt/sPerYear
            
        elif c['physRho'] =='KuipersMunneke2015': #b_dot is in mm W.E. per year.
            ar1=0.07
            ar2=0.03
            Ec=60.0e3
            Eg=42.4e3
            dr_dt = np.zeros(gridLen)
            if c['bdot_type'] == 'instant':
                M_0=1.042-0.0916*np.log(bdotSec[ii]*sPerYear*1e3*0.917)
                M_1=1.734-0.2039*np.log(bdotSec[ii]*sPerYear*1e3*0.917)            
                M_0=np.max(M_0,0.25)
                M_1=np.max(M_1,0.25)
                print 'M_0=',M_0
                print 'M_1=',M_1
                dr_dt[rho<rho1] = (rhoi-rho[rho<rho1])*M_0*ar1*bdotSec[ii]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho<rho1])+Eg/(R*T10m))
                dr_dt[rho>=rho1] = (rhoi-rho[rho>=rho1])*M_1*ar2*bdotSec[ii]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho>=rho1])+Eg/(R*T10m))                           
            elif c['bdot_type'] == 'mean':
                M_0=1.042-0.0916*np.log(bdot_mean[rho<rho1]*sPerYear*917.0)
                M_1=1.734-0.2039*np.log(bdot_mean[rho>=rho1]*sPerYear*917.0)
                M_0=np.max(M_0,0.25)
                M_1=np.max(M_1,0.25)
                dr_dt[rho<rho1] = (rhoi-rho[rho<rho1])*M_0*ar1*bdot_mean[rho<rho1]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho<rho1])+Eg/(R*T10m))
                dr_dt[rho>=rho1] = (rhoi-rho[rho>=rho1])*M_1*ar2*bdot_mean[rho>=rho1]*(1/t)*sPerYear*917.0*g*np.exp(-Ec/(R*Tz[rho>=rho1])+Eg/(R*T10m))
            drho_dt = dr_dt/sPerYear             
        
        elif c['physRho'] =='Spencer2001':   #does not work   # Uncommented out lines 464 - 475
            pass
            #C11 = 3.38e9
            #C12 = 46.8e3
            #C13 = 0.000121
            #C14 = -0.689
            #C15 = 0.149
            #C21 = 9.06e8
            #C22 = 41.0e3
            #C23 = 0.0856
            #C24 = -1.05
            #C25 = -0.0202
            #C31 = 1.38e7
            #C32 = 30.1e3
            #C33 = 0.284
            #C34 = -0.0734
            #C35 = 0.00322

#             drho_dt = np.zeros(gridLen)
#             for j in xrange(gridLen):
#                 if rho[j]<rho1:
#                     dr_dt = rho[j]*C11/sPerYear*np.exp(-C12/(R*Tz[j]))*(1-rho[j]/rhoi)*((1-(1-rho[j]/rhoi)**C13)**C14)*sigma[j]**C15
#                     drho_dt[j] = dr_dt/sPerYear
#                 elif rho[j]<rho2:
#                     dr_dt = rho[j]*C21/sPerYear*np.exp(-C22/(R*Tz[j]))*(1-rho[j]/rhoi)*((1-(1-rho[j]/rhoi)**C23)**C24)*sigma[j]**C25
#                     drho_dt[j] = dr_dt/sPerYear
#                 else:
#                     dr_dt = rho[j]*C31/sPerYear*np.exp(-C32/(R*Tz[j]))*(1-rho[j]/rhoi)*((1-(1-rho[j]/rhoi)**C33)**C34)*sigma[j]**C35
#                     drho_dt[j] = dr_dt/sPerYear

        elif c['physRho']=='Goujon2003_1':
            top2m=np.nonzero(z<=1.)
            gg=9.82
            rho[top2m]=c['rhos0']
            
            sigma_MPa=sigma/(1.0e6)
            sigma_bar=sigma/(1.0e5)
            Qgj=60.0e3
            n=3.0 
            #drho_dt = np.zeros(gridLen)
            dDdt = np.zeros(gridLen)
            # rhoi2=917.0
            rhoi2cgs=.9165*(1.-1.53e-4*(T10m-273.15)) #from Anais
            rhoi2=rhoi2cgs*1000.0
            if ii==0:
                print 'rhoi2=', rhoi2
            D = rho/rhoi2
            
            Dm23=0.9 #transition from zone 2 to 3
            
            D0 = 0.00226 * T10m + 0.03 #from Arnaud, not Goujon (Goujon uses C, not K)
            if D0>0.59: #Transition has to be less than D=0.6
                D0=0.59
                print 'D0 at max of 0.59'
            Dms = D0 + 0.009 #D0 + epsilon factor, maximum is 0.599
            # D0=Dms
            Dmsrho=Dms*rhoi2 #density of change
            if ii==stp-1:
                print 'Dmsrho=', Dmsrho
            ind1=np.argmax(D>=Dms) #first index where the density is greater than or equal to the transition
            if ii<=2:
                print 'ind1=',ind1
            # ind1=ind1
            # if ii<=2:
                # print 'ind1_2=',ind1
            Dm = D[ind1] #actual transition density. Use this so that the transition falls at a node and avoids the sinularity
            if ii<=2:
                print 'Ds=',D[ind1-2:ind1+3]
            # print 'Dms=',Dms
            
            A = 7.89e3*np.exp(-Qgj/(R*Tz))*1.0e-3 #A given in MPa^-3 s^-1, Goujon uses bar as pressure unit
            if ii<=2:
                print ii
                print 'BIGA=',A
                print 'D0=',D0
                print 'Dm=',Dm
            ccc = 15.5 # no units given
            Z0g=110.2*D0**3.-148.594*D0**2.+87.6166*D0-17. #from Anais' code
            
            lp=(D/D0)**(1.0/3.0)
            Zg=Z0g+ccc*(lp-1.0)
            if ii<=2:
                print 'Zg=',Zg
            lpp_n=(4.0*Z0g * (lp-1.0)**2.0 * (2.0*lp+1.0) + ccc * (lp-1.0)**3.0 * (3.0*lp + 1.0))
            lpp_d=(12.0*lp * (4.0*lp - 2.0*Z0g*(lp-1.0) - ccc*(lp-1.0)**2.0))
            lpp=lp + (lpp_n/lpp_d)
            a = (np.pi/(3.0*Zg * lp**2.0)) * (3.0 * (lpp**2.0 - 1.0)*Z0g + lpp**2.0*ccc*(2.0*lpp-3.0)+ccc)

            # lpp = lp + ((4.0*Z0g*(lp-1.0)**(2.0) * (2.0*lp+1.0) + ccc*(lp-1.0)**3.0 * (3.0*lp+1.0)) / (12.0*lp * (4.0 * lp - 2.0 * Z0g * (lp-1.0) - ccc*(lp-1.0)**2.0)))
            # a = (np.pi / (3.0*Zg*lp**2.0)) * (3.0*(lpp**2.0 - 1.0)*Z0g + lpp**2.0 *ccc*(2.0*lpp - 3.0) + ccc)
            sigmastar = (4.0*np.pi*sigma_bar)/(a*Zg*D);
            # gamma=(5.3*A[ind1] * (Dm**2*D0)**(1.0/3.0) * (a[ind1]/np.pi)**(1.0/2.0) * (sigmastar[ind1]/3.0)**n) / ((sigma_bar[ind1]/(Dm**2))*(1-(5.0/3.0*Dm)));
            gamma_an=(5.3*A[ind1] * (Dms**2*Dms)**(1.0/3.0) * (a[ind1]/np.pi)**(1.0/2.0) * (sigmastar[ind1]/3.0)**n) / ((sigma_bar[ind1]/(Dms**2))*(1-(5.0/3.0*Dms))); 
            # print gamma_an
            gamma=0.455/sPerYear
            # gamma=0.5/sPerYear
            if spin and ii<=100:
                gamma=gamma_an/1.3
            else:
                gamma=gamma_an
            if gamma<=0:
                print 'GAMMA less than ZERO:',ii
                
                
            # if ii==stp-1:
            if ii<=3:
                print ii
                print 'lpp_n=',lpp_n
                print 'lpp_d=',lpp_d
                print 'lp=',lp
                print 'D0=',D0
                print 'A=',A
                print 'AAA=',a
                print 'gamma=', gamma
                
            if ii==stp-1:
                print 'GAMMAfinal=',gamma
            

            # dDdt[D<=Dm]=gamma*(sigma_bar[D<=Dm]/(D[D<=Dm])**2.0)*(1.0-(5.0/3.0)*D[D<=Dm])

            dDdt[D<=Dm]=gamma*(sigma_bar[D<=Dm])*(1.0-(5.0/3.0)*D[D<=Dm])/((D[D<=Dm])**2.0)
            dDdt[D>Dm]=5.3*A[D>Dm]* (((D[D>Dm]**2.0)*D0)**(1/3.)) * (a[D>Dm]/np.pi)**(1.0/2.0) * (sigmastar[D>Dm]/3.0)**n
            
            print 'dDdt_ind=',dDdt[ind1-3:ind1+3]
            # dDdt[D<=Dms]=gamma*(sigma_bar[D<=Dms])*(1.0-(5.0/3.0)*D[D<=Dms])/((D[D<=Dms])**2.0)
            # dDdt[D>Dms]=5.3*A[D>Dms]* (((D[D>Dms]**2.0)*D0)**(1/3.)) * (a[D>Dms]/np.pi)**(1.0/2.0) * (sigmastar[D>Dms]/3.0)**n
            
            # print dDdt
            
            rhoC = rho2 #should be Martinerie density
            frho2 = interpolate.interp1d(rho,sigma_bar) #jessica did this... not sure why/if it is right. Max, 12/4/15. It is a function...
            sigmarho2 = frho2(rhoC) #pressure at close off
            # print 'sigrho2=', sigmarho2

            ind2=np.argmax(D>=Dm23)
            
            # sigma_b = sigmarho2 * (D*(1-rhoC/rhoi)) / (rhoC/rhoi*(1-D))
            # sigma_b = (sigma_MPa[ind2] * (D*(1-Dm23)) / (Dm23*(1-D)))/10. #works for Ex2
            sigma_b = ((sigma_bar[ind2]+ atmosP/1.0e5) * (D*(1-Dm23)) / (Dm23*(1-D)))
            # print sigma_b
            sigmaEff = (sigma_bar + atmosP/1.0e5 - sigma_b)#*1000
            sigmaEff[sigmaEff<=0]=1.0e-9
            
            ind2=np.argmax(D>=Dm23)
            ind3=ind2+10
            # if ii==10:
            #     print 'sigMpa=',sigma_MPa[ind2:ind3]
            #     print 'sigb=',sigma_b[ind2:]
            #     print 'atmosP=',atmosP/1.0e5            
            #     print 'sigEff=',sigmaEff[ind2:ind3]
                    
            # dDdt[D>Dm23] = 2.*A[D>Dm23] * ( (D[D>Dm23]*(1-D[D>Dm23])) / (1-(1-D[D>Dm23])**(1/n))**n ) * (2*sigmaEff[D>Dm23]/n)**3.0
            
            Ad = 1.2e3*np.exp(-Qgj/(R*Tz))#*1.0e-3
            T34=0.98
            # dDdt[D>T34] = 9/4*Ad[D>T34]*(1-D[D>T34])*sigmaEff[D>T34]
            
            if ii<=3:
                print 'minSEF=',np.min(sigmaEff)
                print 'maxD=',np.max(A)
                print 'dDdt=',ii, dDdt
                
            if ii==stp-1:
                print 'minSEF=',np.min(sigmaEff)
                print 'maxD=',np.max(A)
                print 'dDdt=',ii, dDdt
            drho_dt=dDdt*rhoi2#/1000.
            drho_dt[top2m]=0.0


        elif c['physRho']=='Goujon2003':
            top2m=np.nonzero(z<=1.1)
            gg=9.82
            rho[top2m]=c['rhos0']
            
            sigma_MPa=sigma/(1.0e6)
            sigma_bar=sigma/(1.0e5)
            Qgj=60.0e3
            n=3.0 
            dDdt = np.zeros(gridLen)
            # rhoi2=917.0
            rhoi2cgs=.9165*(1.-1.53e-4*(T10m-273.15)) #from Anais
            rhoi2=rhoi2cgs*1000.0
            
            D = rho/rhoi2
            
            Dm23=0.9 #transition from zone 2 to 3
            rho23=Dm23*rhoi2

            
            D0 = 0.00226 * T10m + 0.03 #from Arnaud, not Goujon (Goujon uses C, not K)
            if D0>0.59: #Transition has to be less than D=0.6
                D0=0.59
            Dms = D0 + 0.009 #D0 + epsilon factor, maximum is 0.599
            Dmsrho=Dms*rhoi2 #density of change
            print 'Dmsrho', Dmsrho

            ind1=np.argmax(D>=Dms) #first index where the density is greater than or equal to the transition
            print 'ind1', ind1
            Dm = D[ind1] #actual transition density. Use this so that the transition falls at a node and avoids the sinularity
            Dmrho=Dm*rhoi2
            print 'Dm', Dm
            print 'Dmrho', Dmrho
            print 'DD', D[ind1-1:ind1+2]          
            print 'DD2', D[ind1-1:ind1+2]*rhoi2
         
            A = 7.89e3*np.exp(-Qgj/(R*Tz))*1.0e-3 #A given in MPa^-3 s^-1, Goujon uses bar as pressure unit
            ccc = 15.5 # no units given
            Z0g=110.2*D0**3.-148.594*D0**2.+87.6166*D0-17. #from Anais' code           
            lp=(D/D0)**(1.0/3.0)
            Zg=Z0g+ccc*(lp-1.0)
            lpp_n=(4.0*Z0g * (lp-1.0)**2.0 * (2.0*lp+1.0) + ccc * (lp-1.0)**3.0 * (3.0*lp + 1.0))
            lpp_d=(12.0*lp * (4.0*lp - 2.0*Z0g*(lp-1.0) - ccc*(lp-1.0)**2.0))
            lpp=lp + (lpp_n/lpp_d)
            a = (np.pi/(3.0*Zg * lp**2.0)) * (3.0 * (lpp**2.0 - 1.0)*Z0g + lpp**2.0*ccc*(2.0*lpp-3.0)+ccc)
            sigmastar = (4.0*np.pi*sigma_bar)/(a*Zg*D);
#             gamma=(5.3*A[ind1] * (Dms**2*Dms)**(1.0/3.0) * (a[ind1]/np.pi)**(1.0/2.0) * (sigmastar[ind1]/3.0)**n) / ((sigma_bar[ind1]/(Dms**2))*(1-(5.0/3.0*Dms))); 

            if ii==0:
                gamma=0.5/sPerYear
                gamma_old=gamma
            else:
                gamma=gamma_old
            print 'gamma before', gamma
                
            dDdt[D<=Dm]=gamma*(sigma_bar[D<=Dm])*(1.0-(5.0/3.0)*D[D<=Dm])/((D[D<=Dm])**2.0)
            dDdt[D>Dm]=5.3*A[D>Dm]* (((D[D>Dm]**2.0)*D0)**(1/3.)) * (a[D>Dm]/np.pi)**(1.0/2.0) * (sigmastar[D>Dm]/3.0)**n 
            
            gam_div=1.005
            
            # iterate to increase gamma first if not in steady state    
            if ii!=0 and dDdt[ind1]<=dDdt[ind1+1] and gamma!=gamma_old2: #and dDdt_old[ind1]!=dDdt_old[ind1]:
                while dDdt[ind1]<dDdt[ind1+1]:
                    gamma=gamma*(gam_div)
                    print "g1", gamma
                    dDdt[D<=Dm]=gamma*(sigma_bar[D<=Dm])*(1.0-(5.0/3.0)*D[D<=Dm])/((D[D<=Dm])**2.0)
                    dDdt[D>Dm]=5.3*A[D>Dm]* (((D[D>Dm]**2.0)*D0)**(1/3.)) * (a[D>Dm]/np.pi)**(1.0/2.0) * (sigmastar[D>Dm]/3.0)**n  
            # then iterate to find the maximum value of gamma that will make a continuous drho/dt
            while dDdt[ind1]>=dDdt[ind1+1]:
                gamma=gamma/gam_div
                print "g2", gamma
                dDdt[D<=Dm]=gamma*(sigma_bar[D<=Dm])*(1.0-(5.0/3.0)*D[D<=Dm])/((D[D<=Dm])**2.0)
                dDdt[D>Dm]=5.3*A[D>Dm]* (((D[D>Dm]**2.0)*D0)**(1/3.)) * (a[D>Dm]/np.pi)**(1.0/2.0) * (sigmastar[D>Dm]/3.0)**n                      
            raw_input('press enter')
            print 'gamma after ', gamma
            gamma_old2=gamma_old
            gamma_old=gamma
            
            rhoC = rho2 #should be Martinerie density
            frho2 = interpolate.interp1d(rho,sigma_bar) #jessica did this... not sure why/if it is right. Max, 12/4/15. It is a function...
            sigmarho2 = frho2(rhoC) #pressure at close off
#             print 'sigrho2=', sigmarho2

            ind2=np.argmax(D>=Dm23)
            
#             sigma_b = sigmarho2 * (D*(1-rhoC/rhoi)) / (rhoC/rhoi*(1-D))
#             sigma_b = (sigma_MPa[ind2] * (D*(1-Dm23)) / (Dm23*(1-D)))/10. #works for Ex2
#             sigma_b = ((sigma_bar + atmosP/1.0e5) * (D*(1-Dm23)) / (Dm23*(1-D)))
            sigma_b = ((atmosP/1.0e5) * (D*(1-Dm23)) / (Dm23*(1-D)))
#             print sigma_b
            sigmaEff = (sigma_bar + atmosP/1.0e5 - sigma_b)#*1000
            sigmaEff[sigmaEff<=0]=1.0e-9
            
            ind2=np.argmax(D>=Dm23)
            ind3=ind2+10
             

                    
            dDdt[D>Dm23] = 2.*A[D>Dm23] * ( (D[D>Dm23]*(1-D[D>Dm23])) / (1-(1-D[D>Dm23])**(1/n))**n ) * (2*sigmaEff[D>Dm23]/n)**3.0
            
            Ad = 1.2e3*np.exp(-Qgj/(R*Tz))*1.0e-1
            T34=0.98
            if ii==1:
                print 'rho23=',rho23
                print 'rho34=',T34*rhoi2 
            dDdt[D>T34] = 9/4*Ad[D>T34]*(1-D[D>T34])*sigmaEff[D>T34]

            dDdt_old=dDdt
            drho_dt=dDdt*rhoi2#/1000.
            drho_dt[top2m]=0.0
            print ii




        elif c['physRho']=='Morris2014': #uses stress. Need to choose physics for zone 2. 
            QMorris=112.0e3
#             kMorris=11.0
            kMorris=1.2 # m w.e., morris page 361
            rhoW=1000.0
            if spin and ii==0:
                Hx=np.zeros(gridLen) #Hx is temperature history
                Hx=Hx+np.exp(-QMorris/(R*Tz))*dt #is this dt correct?
          
            drho_dt = np.zeros(gridLen)
            #dr_dt = drho_dt
#             drho_dt[rho<=rho1] = (kMorris/(rhoW*g)) * ((rhoi-rho[rho<=rho1])/rho[rho<=rho1]) * 1/Hx[rho<=rho1] * np.exp(-QMorris/(R*Tz[rho<=rho1]))*sigma[rho<=rho1]
            drho_dt[rho<=rho1] = (kMorris/(rhoW*g)) * ((rhoi-rho[rho<=rho1])) * 1/Hx[rho<=rho1] * np.exp(-QMorris/(R*Tz[rho<=rho1]))*sigma[rho<=rho1]
            #print 'max=', np.max(drho_dt)
            
            ########### Choose which zone 2 physics you want
            
            #HL Dynamic
            Q2=21400.
            k2=575.0
            bHL=0.5
            A = bdotSec*(1/t)*sPerYear*rhoiMgm #A from the input json file is m ice equivalent.
            drho_dt[rho>=rho1] = k2*np.exp(-Q2/(R*Tz[rho>=rho1]))*(rhoiMgm-rho[rho>=rho1]/1000)*np.power(A[ii],bHL)*1000/sPerYear
            
            #Arthern
            #ar2=0.03
            #Ec=60.0e3
            #Eg=42.4e3
            #drho_dt[rho>rho1] = ((rhoi-rho[rho>=rho1])*ar2*bdotSec[ii]*(1/t)*sPerYear*rho[rho>rho1]*g*np.exp(-Ec/(R*Tz[rho>rho1])+Eg/(R*T10m)))/sPerYear
            
            #HL Sigfus         
            #f550 = interpolate.interp1d(rho,sigma)
            #sigma550 = f550(rho1)
            #rhoDiff = (rhoiMgm-rho/1000)                                     
            #k = np.power(k2*np.exp(-Q2/(R*Tz[rho>rho1])),2)/sPerYear
            #sigmaDiff = (sigma[rho>rho1]-sigma550)
            #sigmaDiff[sigmaDiff<0]=0 #Stops model from running away - spin up artifact.
            #dd=(k*(sigmaDiff*rhoDiff[rho>rho1]) / (g*np.log((rhoiMgm-rho1/1000)/(rhoDiff[rho>rho1]))))
            ##print 'sig1 = ', dd[0]             
            #drho_dt[rho>rho1] = (k*(sigmaDiff*rhoDiff[rho>rho1]) / (g*np.log((rhoiMgm-rho1/1000)/(rhoDiff[rho>rho1])))) #use H&L Sigfus physics for zone 2.
            
            #Li 2011
            #TmC=T10m-273.15
            #A = bdotSec[ii]*(1/t)*sPerYear*rhoiMgm
            #beta1 = -9.788 + 8.996*A - 0.6165*TmC 
            #beta2 = beta1/(-2.0178 + 8.4043*A - 0.0932*TmC) # this is the one from the paper. Does not work with scaled beta1. Pay attention to units in paper.
            #drho_dt[rho>rho1] = ((rhoi-rho[rho>rho1])*A*beta2*8.36*(KtoC-Tz[rho>rho1])**-2.061)/sPerYear

            ##############
            
            Hx_new=0
            Hx = np.concatenate(([Hx_new],Hx[0:-1]))            
            Hx=Hx+np.exp(-QMorris/(R*Tz))*dt #just an initial order of magnitude value for Hx              
                    
        elif c['physRho']=='Barnola1991': #uses m W.E. (zone 1) and stress (zone 2)
            Q1=10160.
            k1=11.0
            aHL=1.0
            alphaBarnola=-37.455
            betaBarnola=99.743
            deltaBarnola=-95.027
            gammaBarnola=30.673
            A0b =  2.54e4
            n=3.0
            QBarnola=60.0e3
            # print ii
            rho[rho>rhoi]=rhoi #The Barnola model will go a fraction over the ice density (order 10^-3), so this stops that.
            drho_dt = np.zeros(gridLen)
            D = rho/rhoi
            nBa = n*np.ones(gridLen)
            A0 = A0b*np.ones(gridLen)/1.e18 #this is for the n=3 region.
            
            #Vc = (6.95e-4)*T10m-0.043
            #rhoCRhoi = ((Vc+1/(rhoi*(1e-3)))**-1)*1000/rhoi
            #        
            #sigma_b = atmosP*(D*(1-rhoCRhoi))/(rhoCRhoi*(1-D)) # This should be bubble pressure?
            #
            #sigmadiff = sigma - sigma_b
            #print sigmadiff
            #sigmaEff[sigmaEff<0]=1.e-15
            sigmaEff=sigma
            #nBa[sigmaEff<0.1e6]=1.
            #A0[sigmaEff<0.1e6]=A0b/1.e9#for the n=1 region
            
            #this is Max's vectorized version
            # zone 2   
            fe = 10.0**(alphaBarnola*(rho[rho<=800.]/1000)**3.+betaBarnola*(rho[rho<=800.]/1000)**2.+deltaBarnola*rho[rho<=800.]/1000+gammaBarnola)
            drho_dt[rho<=800.] = rho[rho<=800.]*A0[rho<=800.]*np.exp(-QBarnola/(R*Tz[rho<=800.]))*fe*(sigmaEff[rho<=800.]**nBa[rho<=800.])
            
            # zone 1
            A = bdotSec*(1/t)*sPerYear*rhoiMgm
            drho_dt[rho<rho1] = dr_dt = k1 * np.exp(-Q1/(R*Tz[rho<rho1]))*(rhoiMgm-rho[rho<rho1]/1000)*np.power(A[ii],aHL)*1000/sPerYear
            
            # zone 3
            fs = (3./16.)*(1-rho[rho>800.]/rhoi)/(1-(1-rho[rho>800.]/rhoi)**(1./3.))**3.
            drho_dt[rho>800.] = rho[rho>800.]*A0[rho>800.]*np.exp(-QBarnola/(R*Tz[rho>800.]))*fs*(sigmaEff[rho>800.]**nBa[rho>800.])                    
        
        else:
            print 'Error: you need to choose model physics in json (check spelling)'
            sys.exit()
        
        #update the density using explicit method
        rho = rho + dt*drho_dt
        
#         dthickness=mass/(dt*drho_dt)        
        
        #update the age (seconds)
#         age = np.insert(age,0,0) #insert new surface value
#         age= age[:gridLen]+dt #remove bottom element
        ageold=age
        age = np.concatenate(([0],age[:-1]))+dt

        #Heat Diffusion # is this lagrangian?
        if c['heatDiff'] == 'on' :
            nz_P=len(z)
            nz_fv = nz_P-2
            nt = 1
            
            z_edges_vec = z[1:-2]+dz[2:-1]/2
            z_edges_vec = np.concatenate(([z[0]],z_edges_vec, [z[-1]]))
            z_P_vec = z
            rhoP = rho
            phi_s = Ts[ii]
            phi_0 = Tz
            #print Tz

            K_ice = 9.828*np.exp(-0.0057*phi_0)
            K_firn = K_ice * (rhoP / 1000) ** (2 - 0.5 * (rhoP / 1000)) #
            c_firn = 152.5 + 7.122 * phi_0 # specific heat of ice
#             c_firn = c_ice * rho #schwander - specific heat of firn
            Gamma_P = K_firn/(c_firn * rho)

            Tz = transient_solve_TR(z_edges_vec,z_P_vec,nt,dt,Gamma_P, phi_0,nz_P,nz_fv,phi_s)
            Tz = np.concatenate(([Ts[ii]],Tz[:-1]))
            
            fT10m = interpolate.interp1d(z,Tz) #temp at 10m depth
            T10m = fT10m(10)
            
        #update the length of the boxes
        z_old=z
        #dzNew = bdotSec[ii]*rhoi/rhos0[ii]*dt #these lines worked on 3/23
        dzNew = bdotSec[ii]*rhoi/rhos0[ii]*sPerYear
        dz = mass/rho*dx
        dz = np.concatenate(([dzNew],dz[:-1]))
#         dz =dz.transpose()
        z = dz.cumsum(axis = 0)
        
        # update z and rho vectors
        z = np.concatenate(([0],z[:-1]))         
        rho = np.concatenate(([rhos0[ii]],rho[:-1]))
        
        zdiffnew=(z[1:]-z[1])
        zdiffold=(z_old[0:-1]-z_old[0])
        
        compaction_rate=(zdiffold-zdiffnew)/dt*sPerYear #this is cumulative compaction rate in m/yr from 0 to the node specified in depth
        
        if not spin:
            Dcon = np.concatenate(([D_surf[ii]],Dcon[:-1])) #Dcon is a diffusivity tracker. Just creating a vector at each time step that is some fraction which is multiplied by diffusivity in firnair.

        #update mass
        massNew =bdotSec[ii]*sPerYear*rhoi
        mass = np.concatenate(([massNew],mass[:-1]))
        mass_sum=mass.cumsum(axis=0)
        #if ii==2:
        #    print "mass=" , mass[0:4]

        sigma = mass * g * dx # kg*m*s^-2 --> Pa
        sigma = sigma.cumsum(axis = 0)
        
        ### bdot_mean method 1: divide mass by time ####
        bdot_mean = np.concatenate(([mass_sum[0]/(rhoi*sPerYear)],mass_sum[1:]*t/(age[1:]*rhoi)))
#         if ii==5:
#             print 'bdotmean=',bdot_mean[0:5]
        # bdot_mean is in units of m ice equivalent per time step in units of seconds. To get m water, multiply by 0.917. To get kg, multiply by 917.0.
        
#         if ii==2:
#             print "mass_sum2=", mass_sum[0:2]
#             print "bdotmean=" , bdot_mean[0:4]        
        #bdot_mean=mass_sum/age
        ################
        
        ### bdot mean method 2 ####
#         bdot_meanNew=((bdot_mean*ageold/dt+bdotSec[ii])*dt)
#         bdot_meanNew=bdot_meanNew#/(age)
#         bdot_mean = (np.concatenate(([bdotSec[ii]],bdot_meanNew[:-1]/age[1:])))
#         if ii==1:
#             print 'age=',age[0:2]
#             print 'ageold=',ageold[0:2]
#             print 'dt=',dt
#             print 'bdotmean=',bdot_meanNew[0:2]
#         #print bdot_mean[0:5]
        ##############  
          
                
        #Update grain growth
        if c['physGrain']: #this is grain growth
            kgr=1.3e-7 #grain growth rate from Arthern (2010)
            Eg=42.4e3
            # Temp should be Temp at 10m depth.
            #dr2Dt = kgr * np.exp(-Eg/(R*T10m))
            dr2Dt = kgr * np.exp(-Eg/(R*Tz)) #units are m^2/s
            r2 = r2 + dr2Dt * dt #* c['stpsPerYear'] #does this need the c['stps']? 11/19/14
            #Add time elapsed to beginning of row for output file
            r2_time = np.concatenate(([t*ii + 1], r2))
            if c['calcGrainSize']: #This uses surface temperature to get an initial grain size.
                r2 = np.concatenate(([-2.42e-9*Ts[ii]+9.46e-7], r2[:-1]))
            else:
#                 r2 = np.concatenate(([-2.42e-9*(c['Ts0'])+9.46e-7], r2[:-1]))
                r2 = np.concatenate(([(0.1e-3)**2], r2[:-1]))

        if spin and ii == (stp-1):
            rho_time = np.concatenate(([t*ii + 1],rho))
            Tz_time = np.concatenate(([t*ii + 1],Tz))
            age_time = np.concatenate(([t*ii + 1],age))
            z_time = np.concatenate(([t*ii + 1],z))
            if c['bdot_type'] == 'mean':
                bdot_time = np.concatenate(([t*ii + 1],bdot_mean))
            if c['physGrain']:
                r2_time = np.concatenate(([t*ii + 1],r2))
            if c["physRho"]=='Morris2014':
                Hx_time = np.concatenate(([t*ii + 1],Hx))
                
            '''initialize files to write in time loop'''
            densityPath = os.path.join(c['resultsFolder'], 'densitySpin.csv')
            tempPath = os.path.join(c['resultsFolder'], 'tempSpin.csv')
            agePath = os.path.join(c['resultsFolder'], 'ageSpin.csv')
            depthPath = os.path.join(c['resultsFolder'], 'depthSpin.csv')
            if c['bdot_type'] == 'mean':
                bdotPath = os.path.join(c['resultsFolder'], 'bdot_meanSpin.csv')
            if c['physGrain']:
                r2Path = os.path.join(c['resultsFolder'], 'r2Spin.csv')
            if c["physRho"]=='Morris2014':
                HxPath = os.path.join(c['resultsFolder'], 'HxSpin.csv')
                
            with open(densityPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(rho_time)
            with open(tempPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(Tz_time)
            with open(agePath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(age_time)
            with open(depthPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(z_time)                
            if c['bdot_type'] == 'mean':
                with open(bdotPath, "a") as f:
                    writer = csv.writer(f)
                    writer.writerow(bdot_time)              
            if c['physGrain']:
                print "PhysGrain is on"
                with open(r2Path, "a") as f:
                    writer = csv.writer(f)
                    writer.writerow(r2_time)  
            if c["physRho"]=='Morris2014':
                with open(HxPath, "a") as f:
                    writer = csv.writer(f)
                    writer.writerow(Hx_time)   
                                             
#         elif not spin and (t*ii)%1 == 0:

        #elif not spin and [True for jj in TWrite if jj == t*ii+1] == [True]:
        elif not spin and [True for jj in TWrite if jj == mtime] == [True]: #write model output to file
            rho_time = np.append(mtime,rho[0::spacewriteint])
            Tz_time = np.append(mtime,Tz[0::spacewriteint])
            age_time = np.append(mtime,age[0::spacewriteint]/sPerYear)
            z_time = np.append(mtime,z[0::spacewriteint])
            comprate_time =np.append(mtime,compaction_rate[0::spacewriteint])
            Dcon_time = np.append(mtime,Dcon[0::spacewriteint])
            stress_time = np.append(mtime,sigma[0::spacewriteint])
            if c['physGrain']:
                r2_time = np.append(mtime,r2[0::spacewriteint])
            if c["physRho"]=='Morris2014':
                Hx_time = np.append(mtime,Hx[0::spacewriteint])
            Clim_time = np.append(mtime,[bdot[ii],Ts[ii]])
            bdot_time = np.append(mtime,bdot_mean[0::spacewriteint])
                
            with open(densityPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(rho_time)
            with open(tempPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(Tz_time)
            with open(agePath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(age_time)
            with open(depthPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(z_time)
            with open(compratePath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(comprate_time)
            with open(DconPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(Dcon_time)
            with open(stressPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(stress_time)
            with open(ClimPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(Clim_time)
            with open(bdotPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(bdot_time)
            if c['physGrain']:
                with open(r2Path, "a") as f:
                    writer = csv.writer(f)
                    writer.writerow(r2_time)
            if c["physRho"]=='Morris2014':                
                with open(HxPath, "a") as f:
                    writer = csv.writer(f)
                    writer.writerow(Hx_time)         
                
            ### BCO,LIZ, and DIP ###
#             Vc = (6.95e-4)*T10m-0.043 #Martinerie et al., 1994, Eq. 2: critical pore volume at close off
#             bcoMartRho = ((Vc+1/(rhoi*(1e-3)))**-1)*1000 # Martinerie density at close off
            bcoMartRho = 1/( 1/(917.0) + T10m*6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
            bcoAgeMart = min(age[rho>=bcoMartRho])/sPerYear # close-off age from Martinerie
            bcoDepMart = min(z[rho>=(bcoMartRho)])
            bcoAgeMartAll.append(bcoAgeMart) #age at the 815 density horizon
            bcoDepMartAll.append(bcoDepMart) #this is the 815 close off depth
            
            ### bubble close-off age and depth assuming rho_crit = 815kg/m^3
            bcoAge815 = min(age[rho>=(rho2)])/sPerYear #close-off age where rho=815 kg m^-3
            bcoDep815 = min(z[rho>=(rho2)]) #depth of 815 horizon
            bcoAge815All.append(bcoAge815) #age at the 815 density horizon
            bcoDep815All.append(bcoDep815) #this is the 815 close off depth
            
            ### Lock-in depth and age
            LIZMartRho = bcoMartRho - 14.0 #LIZ depth (Blunier and Schwander, 2000)      
            LIZAgeMart = min(age[rho>LIZMartRho])/sPerYear #lock-in age
            LIZDepMart = min(z[rho>=(LIZMartRho)]) #lock in depth
            phiLockIn=0.15 #relic?
            lockIn = min(age[phiClosed[phiClosed<phi]/phi[phiClosed<phi] >= phiLockIn])/sPerYear #lock-in age Old, from Jessica's original. Too shallow, likely.       
            LIZAgeAll.append(LIZAgeMart)
            LIZDepAll.append(LIZDepMart)
            
            ### Porosity, including DIP
            phi = 1-rho/rhoi # total porosity
            phi[phi<=0]=1e-16 
            phiC = 1-bcoMartRho/rhoi; #porosity at close off
            phiClosed = np.ones(gridLen)
            phiClosed = 0.37*phi*(phi/phiC)**-7.6 #Closed porosity, from Goujon. See Buizert thesis (eq. 2.3) as well
                    
            phiOpen = phi - phiClosed #open porosity
            phiOpen[phiOpen<=0] = 1.e-10 #don't want negative porosity.
                    
            intPhi = np.sum(phi * dz) #depth-integrated porosity
            intPhiAll.append(intPhi)      
            
        if not spin and ii==10:
            loop10time=time.time()-tic2
            elapsed=time.time()-tic
            elapsed=elapsed*(stp/10)
            elapsed_min=elapsed/60.
            mins=np.floor(elapsed_min)
            secs=(elapsed_min-mins)*60 
            print "estimated model run time is %s minutes %s seconds" % (mins, secs)
            ####### end time stepping loop
                    
    #print 'dt = ', dt
    elapsed=time.time()-tic
    elapsed_min=elapsed/60.
    mins=np.floor(elapsed_min)
    secs=(elapsed_min-mins)*60
       

    if spin:
        logging.info("spin up took %s minutes %s seconds" % (mins, secs))
   
    if not spin: #BCO, LID, DIP writer
        with open(bcoPath, "w") as f: #write BCO.csv file. rows are: time, BCO age (mart), BCO depth (mart),BCO age (815), BCO depth (815) 
            writer = csv.writer(f)
            writer.writerow(np.append(modeltime[0],TWrite[:len(bcoAge815All)])) 
            writer.writerow(bcoAgeMartAll)
            writer.writerow(bcoDepMartAll)
            writer.writerow(bcoAge815All)
            writer.writerow(bcoDep815All)
        with open(lidPath, "w") as f: #write LIZ information. Rows: time, age, dep
            writer = csv.writer(f)
            writer.writerow(np.append(modeltime[0],TWrite[:len(LIZDepAll)]))
            writer.writerow(LIZAgeAll)
            writer.writerow(LIZDepAll)
        with open(intPhiPath, "w") as f: #depth-integrated porosity
            writer = csv.writer(f)
            writer.writerow(np.append(modeltime[0],TWrite[:len(intPhiAll)]))
            writer.writerow(intPhiAll)
                
        logging.debug("spin up years = %s" % c['yearSpin'])  
        logging.debug("spin up steps per year = %s" % c["stpsPerYearSpin"])
        logging.debug("model run years = %s" % c["years"])
        logging.debug("model run steps per year = %s" % c["stpsPerYear"])
        logging.info("model run took %s minutes %s seconds" % (mins, secs))
        logpath = os.path.join(spot,c['resultsFolder'])    
        shutil.copy(os.path.join(spot,'RUNDETAILS.log'),logpath)
#         endlogger()
#         os.remove('RUNDETAILS.log')

    if c['plotting'] == 'on' and spin == 0:
        if os.path.exists(c['plotsFolder']):
            rmtree(c['plotsFolder'])
        os.makedirs(c['plotsFolder'])   
        plotData(c['physGrain'], c['resultsFolder'], c['plotsFolder'])

        
if __name__ == '__main__':
    
    startlogger()
    
    if len(sys.argv) >= 2 and '-s' not in sys.argv:
        configName = os.path.join(os.path.dirname(__file__), sys.argv[1])
    elif len(sys.argv) >= 2 and '-s' in sys.argv:
        configName = os.path.join(os.path.dirname(__file__), sys.argv[1])
    else:
#         configName = os.path.join(os.path.dirname(__file__), 'configs', 'configLi1.json')
        configName = os.path.join(os.path.dirname(__file__), 'generic.json')

    if '-s' in sys.argv:
        spin = True
    else:
        spin = False
        
#     spin = True

# This right now is for testing purposes.
    with open(configName, "r") as f:
        jsonString = f.read()
        c = json.loads(jsonString)

        
        #input_year, input_temp, input_bdot = input_data.retrieve_data(c["InputFileName"])


    runModel(configName,spin)

    
    