'''
Created on Jun 1, 2013
@author: Jessica Lundin
'''

import csv
import json
import sys
import math
import numpy as np
import logging
from scipy import interpolate
from scipy.sparse import spdiags
import scipy.sparse.linalg as splin
from plot import plotData
from shutil import rmtree
import matplotlib.pyplot as plt
import os
from string import join
import shutil
import time
# import input_data

#global input_temp
#global input_year
#global input_AR


def HerronLangwayAnalytic(c,h,THL,AHL):
    """Model steady-state firn density and age profiles and bubble close-off
    config -- firnmod.config.Config
    return -- age, rho (density) for steady-state dynamics
    """
    print 'AHL=%s' % AHL 
    hSize = np.size(h)      
    rhoc = 0.550
    rhos = c['rhos0']/1000.0
    # Eventually make this into 1 'if else' conditional statement
    #if (c['userInput']):
    #    A = input_AR[0] * c['rhoiMgm']
    #else:
        #A = c['bdot0']*c['rhoiMgm']#* 0.917 'bdot0' is in ice equivalent (I think), but H&L runs using W.E. so we multiply by 0.917
    A=AHL*c['rhoiMgm']
    ### Right now only takes the first index in input_temp
    ### Going to have to add more info to this... AKA make sure the step years actually get the right year in the array
    #if (c['userInput']):
    #    k0 = 11.0 * np.exp(-10160/(c['R'] * c['Ts0'])) #Max Changed this. Look for an easy way to get the first value from input file
    #    k1 = 575.0 * np.exp(-21400/(c['R'] * c['Ts0']))
    ## Herron and Langway k0 = (6a), k1 = (6b)
    #else:
    #    k0 = 11.0 * np.exp(-10160/(c['R']*c['Ts0'] ))
    #    k1 = 575.0 * np.exp(-21400/(c['R']*c['Ts0'] ))
    k0 = 11.0 * np.exp(-10160/(c['R']*THL ))
    k1 = 575.0 * np.exp(-21400/(c['R']*THL ))    
    
# depth of critical density, eqn 8 from Herron and Langway
    h0_55 = 1/(c['rhoiMgm']*k0) * (np.log(rhoc/(c['rhoiMgm']-rhoc))-np.log(rhos/(c['rhoiMgm']-rhos)))
    Z0 = np.exp(c['rhoiMgm']*k0*h + np.log(rhos/(c['rhoiMgm']-rhos)))
    # The boundary from zone 1 to zone 2 = t0_55
    t0_55 = 1/(k0*A) * np.log((c['rhoiMgm']-rhos)/(c['rhoiMgm']-rhoc ))
    rho_h0 = (c['rhoiMgm']* Z0)/(1+Z0)
    if np.max(rho_h0) >= c['rhoiMgm']:
        t0 = np.zeros(hSize)
        for jj in xrange(hSize):
            if rho_h0[jj]<=c['rhoiMgm']-0.001:
                t0[jj] = (1/(k0*A)*np.log((c['rhoiMgm']-rhos)/(c['rhoiMgm']-rho_h0[jj])))
                jj_max = jj
            else:
                t0[jj] = (t0[jj_max])
        
    else:
        t0 = 1/(k0*A)*np.log((c['rhoiMgm']-rhos)/(c['rhoiMgm']-rho_h0))
    
    Z1 = np.exp(c['rhoiMgm']*k1*(h-h0_55)/np.sqrt(A) + np.log(rhoc/(c['rhoiMgm']-rhoc)))
    Z = np.concatenate((Z0[h<h0_55], Z1[h>h0_55]))
    rho_h = (c['rhoiMgm'] * Z)/(1+Z)
    tp = np.ones(hSize)
    for j in xrange(hSize):
        if rho_h[j]<c['rhoiMgm']-0.01:
            tp[j] = 1/(k1*np.sqrt(A)) * np.log((c['rhoiMgm']-rhoc)/(c['rhoiMgm']-rho_h[j]))+ t0_55
            jMax = j
        else:
            tp[j] = tp[jMax]
    # Zone 1 and Zone 2 repsectively
    age = np.concatenate((t0[h<h0_55], tp[h>h0_55]))*c['sPerYear']
    rho = rho_h*1000
    return age, rho

def solver(a_U,a_D,a_P,b):
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
    
    for i_time in xrange(nt):
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

    tic=time.time()
    "Runs firnmodel with an input json file."
    count = 0
    logging.getLogger()
    logging.basicConfig(filename='RUNDETAILS.log',level=logging.DEBUG,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    #logging.getLogger('').addHandler(console)
    
    if spin:
        logging.info('Spin Up initiated')
        
    elif not spin:
        logging.info('Model run initiated')  
          
    logging.info("Model configName = %s, spin = %r" % (configName, spin))
    
    with open(configName, "r") as f:
        
        jsonString = f.read()
        c = json.loads(jsonString)
    
    logging.info("The physics are %s" % c['physRho'])
    
    if c["userInput"]:
    # Right now, hard coded for only 3 variables: year, temperature, Accumulation Rate
        FIDtemp=c["InputFileNameTemp"]
        userdata_temp=np.genfromtxt(FIDtemp, delimiter=',')
        input_year_temp=userdata_temp[0,:]
        input_temp=userdata_temp[1,:]
        
        FIDaccu=c["InputFileNameAccu"]
        userdata_accu=np.genfromtxt(FIDaccu, delimiter=',')
        input_year_accu=userdata_accu[0,:]
        input_AR=userdata_accu[1,:]            
    
    if spin:    
        
        if os.path.exists(c['resultsFolder']):
            rmtree(c['resultsFolder'])
        os.makedirs(c['resultsFolder'])

        if (c['userInput']):
            gridLen = int((c['H']-c['HbaseSpin'])/(input_AR[0]/c['stpsPerYearSpin']))
        else:
            gridLen = int((c['H']-c['HbaseSpin'])/(c['bdot0']/c['stpsPerYearSpin']))
        gridheight = np.linspace(c['H'],c['HbaseSpin'],gridLen)
        z = c['H']-gridheight
        dz = np.diff(z)
        dz = np.append(dz,dz[-1])
        dx = np.ones((gridLen)) #assume box width of 1 unit. To include stress must update this for horizontal longitudinal stress
        
        if (c['userInput']):
            THL=input_temp[0]
            AHL=input_AR[0]
        
        else:
            THL=c['Ts0']
            AHL=c['bdot0']
        
        age, rho = HerronLangwayAnalytic(c,z,THL,AHL)
        
        if c['physGrain']:
#             r2 = np.linspace(c['r2s0'], (6*c['r2s0']), gridLen)
            if c['calcGrainSize']:
                r02 = -2.42e-9*(c['Ts0'])+9.46e-7 #m^2, Gow 1967
                r2 = r02*np.ones(gridLen)
            else:
                r2 = np.linspace(c['r2s0'], (6*c['r2s0']), gridLen) 
        
        
#         age = np.zeros(gridLen)
#         rho = c['rhos0']*np.ones(gridLen)
        # initial temperature profile
        if (c['userInput']):
            Tz = input_temp[0]*np.ones(gridLen)
        else:
            Tz = c['Ts0']*np.ones(gridLen) #init Temp profile
        agez0 = age
        rhoz0 = rho
        z0 = z
        Tz0 = Tz
        if c['physGrain']:
            r20 = r2
        #if (c['userInput']):
        #    years = c['yearSpin']
        #    year = input_year[-1] - input_year[0] + 1
        #    stp = int(years * (year / len(input_year))) #this makes the time step consistent with the input data
        #else:
        #    years = c['yearSpin']
        #    stp = int(years *c['stpsPerYearSpin'])
        
        years = c['yearSpin'] #should not need user input years here - just specify in the json
        stp = int(years *c['stpsPerYearSpin'])+1
        
        #Dcon = c['D_surf']*np.ones(gridLen)

    else: #not spin
        
        densityPath = os.path.join(c['resultsFolder'], 'densitySpin.csv')
        tempPath =    os.path.join(c['resultsFolder'], 'tempSpin.csv')
        agePath =     os.path.join(c['resultsFolder'], 'ageSpin.csv')
        depthPath =   os.path.join(c['resultsFolder'], 'depthSpin.csv')
        if c["physGrain"]:
            r2Path = os.path.join(c['resultsFolder'], 'r2Spin.csv')
        
        initDepth =   np.genfromtxt(depthPath, delimiter = ',')
        initAge =     np.genfromtxt(agePath, delimiter = ',')
        initDensity = np.genfromtxt(densityPath, delimiter = ',')
        initTemp =    np.genfromtxt(tempPath, delimiter = ',')
        if c["physGrain"]:
            initr2 =    np.genfromtxt(r2Path, delimiter = ',')
          
        z   = initDepth[1:]
        gridLen = np.size(z)
        age = initAge[1:]
        rho = initDensity[1:]
        Tz  = initTemp[1:]
        if c["physGrain"]:
            r2=initr2  
        dz = np.diff(z)
        dz = np.append(dz,dz[-1])
        dx = np.ones(gridLen)
        if (c['userInput']):
            yr_start=max(input_year_temp[0],input_year_accu[0])
            yr_end=min(input_year_temp[-1],input_year_accu[-1])
            print yr_start
            print yr_end
            years = yr_end-yr_start
            
            #years = input_year[-1] - input_year[0] + 1
            #stp = int(years * (years / len(input_year)))
            stp = int(years *c['stpsPerYear']) #Make sure that stpsPerYear is set properly in config.
            print stp
            modeltime=np.linspace(yr_start,yr_end,stp+1)
            #modeltime=modeltime-min(modeltime)
            print 'modeltime= %s' % modeltime
        else:
            years = c['years']
            stp = int(years *c['stpsPerYear']) #how many steps there are in the model run
            modeltime=np.linspace(0,years,stp+1)
        TWrite = np.concatenate((xrange(0,110,10),xrange(101,150,1),xrange(150,250,5),xrange(250,2010,10))) #set times at which to write data. This line is firnmice.
        #TWrite = (np.arange(0,2005,5)) #set times at which to write data
        
        inte=4 # how often the data should be written
        #TWrite = modeltime[inte::inte] # vector of times data will be written. Can customize here. 
        #print 'Twrite = %s' %TWrite
        z0 = z
        agez0 = age
        rhoz0 = rho
        Tz0 = Tz
        if c["physGrain"]:
            r20 = r2
        Dcon = c['D_surf']*np.ones(gridLen) #initialize a grid of Diffusivity constant

        
    # mass and stress
    mass = rho*dz
    sigma = mass * dx * c['g']
    sigma = sigma.cumsum(axis = 0)

    # define the time domain
    timTotal = years*c['sPerYear']
    dt = timTotal/stp
    
    fT10m = interpolate.interp1d(z,Tz) #temp at 10m depth
    T10m = fT10m(10)



    if spin:
        if c['stpsPerYearSpin'] or (stp / years) == 1.:
            if (c['userInput']):
                Ts = input_temp[0]*np.ones(stp)
            else:
                Ts = c['Ts0']*np.ones(stp)
        else:
            TPeriod = c['yearSpin']
            #   Temperature calculation from Anais Orsi, supplementry material, AGU 2012
            #   http://www.agu.org/journals/gl/gl1209/2012GL051260/supplement.shtml
            if (c['userInput']):
                Ts = input_temp[0] + c['TAmp']*(np.cos(2*np.pi*np.linspace(0,TPeriod,stp))+0.3*np.cos(4*np.pi*np.linspace(0,TPeriod,stp)))
            else:
                Ts = c['Ts0'] + c['TAmp']*(np.cos(2*np.pi*np.linspace(0,TPeriod,stp))+0.3*np.cos(4*np.pi*np.linspace(0,TPeriod,stp)))
        
        if (c['userInput']):
            print 'stp=%s,years=%s' %(stp,years)
            #t = (1.0 / (stp/years))
            t = 1.0 / c['stpsPerYearSpin']
        else:
            t = 1.0 / c['stpsPerYearSpin']
        if (c['userInput']):
            bdotSec0 = input_AR[0]/c['sPerYear']/c['stpsPerYearSpin']
        else:
            bdotSec0 = c['bdot0']/c['sPerYear']/c['stpsPerYearSpin']
        bdotSec = bdotSec0*np.ones(stp)
        rhos0=c['rhos0']*np.ones(stp)
        #D_surf=c['D_surf']*np.ones(stp)
        
    else: #not spin
        if (c['userInput']):
            #t = 1.0 / (stp / years)
            t = 1.0/c['stpsPerYear']
            print 't = %s' % t
            TPeriod = years
        else:
            t = 1.0/c['stpsPerYear']
            TPeriod = c['years']
#         if c['BCtemp'] == 'stepChange':
#             Ts = c['Ts0']*np.ones(stp)
#             Ts[99*c['stpsPerYear']:] = Ts[99*c['stpsPerYear']:]+c['Tpert']
        # Instead of storing the 'Ts0', it takes the input data and makes a new vector of the changing temperatures
        if (c['userInput']):
            #Ts = input_temp[:]
            #bdotSec = input_AR[:]/c['sPerYear']/(stp / years)
            
            Ts=np.interp(modeltime,input_year_temp,input_temp) #interpolate to model time
            bdot=np.interp(modeltime,input_year_accu,input_AR)

            bdotSec = bdot/c['sPerYear']/(stp / years)
#           rhos0 =
#           D_surf =
        else:
            Ts = c['Ts0']*np.ones(stp)
        # if there is a StepChange, add the perturbation

            if c['BCtemp'] == 'stepChange':
                #   Temperature calculation from Anais Orsi, supplementry material, AGU 2012
                #   http://www.agu.org/journals/gl/gl1209/2012GL051260/supplement.shtml
                Ts = c['Ts0']*np.ones(stp)
                Ts[99*c['stpsPerYear']:] = Ts[99*c['stpsPerYear']:]+c['Tpert'] #indexing starts at zero, so 99 gives a change at year 100.
                #Ts[99*c['stpsPerYear']:] = c['Tpert'] #use this line if you want the value specified in the json file to be the actual value
            if t < 1.0: #add seasonal signal
                    Ts = Ts + c['TAmp']*(np.cos(2*np.pi*np.linspace(0,TPeriod,stp))+0.3*np.cos(4*np.pi*np.linspace(0,TPeriod,stp)))

            # Accumulation Rate being added in, vector bdot
            bdotSec0 = c['bdot0']/c['sPerYear']/c['stpsPerYear']
            bdotSec = bdotSec0*np.ones(stp)
            if c['BCbdot'] == 'stepChange':
                bdotPert = c['bdotPert']/c['sPerYear']/c['stpsPerYear']
                bdotSec[99*c['stpsPerYear']:] = bdotSec[99*c['stpsPerYear']:]+bdotPert
                #bdotSec[99*c['stpsPerYear']:] = bdotPert #use this line if you want the value specified in the json file to be the actual value, rather than the size of the pert.

    # Eventually want to get these also under 'user_input' as lines 4 and 5 of csv file
        rhos0=c['rhos0']
        rhos0 = rhos0*np.ones(stp)
        if c['BCrhos0'] == 'stepChange':
            rhos0Pert = c['rhos0Pert']
            rhos0[99*c['stpsPerYear']:] = rhos0[99*c['stpsPerYear']:]+rhos0Pert
            
        D_surf=c['D_surf']*np.ones(stp) #this is a diffusivity tracker: D_surf can be smaller to make layers with lower/higher diffusivity (some fraction of the number calculated using parameterization
        if c['BC_D_surf'] == 'stepChange':
            D_surf_Pert = c['D_surf_Pert']
            D_surf[99*c['stpsPerYear']:] = D_surf[99*c['stpsPerYear']:]+D_surf_Pert
    
    if not spin:
        rho_time = np.append(modeltime[0],rho)
        Tz_time = np.append(modeltime[0],Tz)
        age_time = np.append(modeltime[0],age)
        z_time = np.append(modeltime[0],z)
        D_time = np.append(modeltime[0],Dcon)
        if (c['userInput']):
            Clim_time = np.append(modeltime[0],[bdot[0],Ts[0]]) #not sure if bdot or bdotSec
        if c["physGrain"]:
            r2_time = np.append(modeltime[0],r2)
         
        '''initialize files to write in time loop'''
        densityPath = os.path.join(c['resultsFolder'], 'density.csv')
        tempPath = os.path.join(c['resultsFolder'], 'temp.csv')
        agePath = os.path.join(c['resultsFolder'], 'age.csv')
        depthPath = os.path.join(c['resultsFolder'], 'depth.csv')
        DconPath = os.path.join(c['resultsFolder'], 'Dcon.csv')
        if (c['userInput']):        
            ClimPath = os.path.join(c['resultsFolder'], 'Clim.csv')
        if c["physGrain"]:
            r2Path = os.path.join(c['resultsFolder'], 'r2.csv')
            
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
        with open(DconPath, "w") as f:
            writer = csv.writer(f)
            writer.writerow(D_time)
        if (c['userInput']):            
            with open(ClimPath, "w") as f:
                writer = csv.writer(f)
                writer.writerow(Clim_time)
        if c["physGrain"]:
            with open(r2Path, "w") as f:
                writer = csv.writer(f)
                writer.writerow(r2_time)            

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
        #bcoMartRho = ((Vc+1/(c['rhoi']*(1e-3)))**-1)*1000 # Martinerie density at close off
        bcoMartRho = 1/( 1/(917.0) + T10m*6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
        bcoAgeMart = min(age[rho>=bcoMartRho])/c['sPerYear'] # close-off age from Martinerie
        bcoDepMart = min(z[rho>=(bcoMartRho)])
        bcoAgeMartAll.append(bcoAgeMart) #age at the 815 density horizon
        bcoDepMartAll.append(bcoDepMart) #this is the 815 close off depth
        
        # bubble close-off age and depth assuming rho_crit = 815kg/m^3
        bcoAge815 = min(age[rho>=(c['rho2'])])/c['sPerYear'] #close-off age where rho=815 kg m^-3
        bcoDep815 = min(z[rho>=(c['rho2'])]) #depth of 815 horizon
        bcoAge815All.append(bcoAge815) #age at the 815 density horizon
        bcoDep815All.append(bcoDep815) #this is the 815 close off depth
        
        ### Lock-in depth and age
        LIZMartRho = bcoMartRho - 14.0 #LIZ depth (Blunier and Schwander, 2000)      
        LIZAgeMart = min(age[rho>LIZMartRho])/c['sPerYear'] #lock-in age
        LIZDepMart = min(z[rho>=(LIZMartRho)]) #lock in depth
        #lockIn = min(age[phiClosed[phiClosed<phi]/phi[phiClosed<phi] >= c['phiLockIn']])/c['sPerYear'] #lock-in age Old, from Jessica's original. Too shallow, likely.       
        LIZAgeAll.append(LIZAgeMart)
        LIZDepAll.append(LIZDepMart)
        
        ### Porosity, including DIP
        phi = 1-rho/c['rhoi'] # total porosity
        phi[phi<=0]=1e-16 
        print "min phi= ", np.min(phi)
        phiC = 1-bcoMartRho/c['rhoi']; #porosity at close off
        print "phiC=", phiC
        #phiClosed = np.ones(gridLen)
        phiClosed = 0.37*phi*(phi/phiC)**-7.6 #Closed porosity, from Goujon. See Buizert thesis (eq. 2.3) as well
        
        #for jj in range(gridLen): #porosity, from Goujon et al. 2003 (good explanation in Buizert thesis, eq. 2.3). Not sure why this is looped.
        #    if phi[jj]>phiC:
        #        phiClosed[jj] = c['gamma'] * phi[jj]*(phi[jj]/phiC)**-7.6
        #        jjMax = jj
        #    else:
        #        phiClosed[jj] = phiClosed[jjMax] #phiClosed is closed porosity
                
        phiOpen = phi - phiClosed #open porosity
        phiOpen[phiOpen<=0] = 1.e-10 #don't want negative porosity.
                   
        intPhi = np.sum(phi * dz) #depth-integrated porosity
        intPhiAll.append(intPhi)
            
        bcoPath = os.path.join(c['resultsFolder'], 'BCO.csv')
        lidPath = os.path.join(c['resultsFolder'], 'LID.csv')
        intPhiPath = os.path.join(c['resultsFolder'], 'porosity.csv')
     
    for i in xrange(stp): #start main time-stepping loop
        if not spin:
            mtime=modeltime[i] #placeholder for writing data.
#         print 'i=%s' % i
        
#         if c['physRho']=='HLdynamic':
#             A = bdotSec*(1/t)*c['sPerYear']*c['rhoiMgm'] #A from the input json file is m ice equivalent.
#             drho_dt = np.zeros(gridLen)
#             for j in xrange(gridLen):
#                 #print 'A=%s' % A
#                 if rho[j]<c['rho1']:
#                     dr_dt = c['k1']*np.exp(-c['Q1']/(c['R']*Tz[j]))*(c['rhoiMgm']-rho[j]/1000)*np.power(A[i],c['aHL'])*1000/c['sPerYear']
#                     drho_dt[j] = dr_dt
#                 else:
#                     dr_dt = c['k2']*np.exp(-c['Q2']/(c['R']*Tz[j]))*(c['rhoiMgm']-rho[j]/1000)*np.power(A[i],c['bHL'])*1000/c['sPerYear']
#                     drho_dt[j] = dr_dt
                    
                    
        if c['physRho']=='HLdynamic':
            A = bdotSec*(1/t)*c['sPerYear']*c['rhoiMgm'] #A from the input json file is m ice equivalent.
            drho_dt = np.zeros(gridLen)
            dr_dt=drho_dt
            dr_dt[rho<c['rho1']] = c['k1']*np.exp(-c['Q1']/(c['R']*Tz[rho<c['rho1']]))*(c['rhoiMgm']-rho[rho<c['rho1']]/1000)*np.power(A[i],c['aHL'])*1000/c['sPerYear']
            drho_dt[rho<c['rho1']] = dr_dt[rho<c['rho1']]

            dr_dt[rho>=c['rho1']] = c['k2']*np.exp(-c['Q2']/(c['R']*Tz[rho>=c['rho1']]))*(c['rhoiMgm']-rho[rho>=c['rho1']]/1000)*np.power(A[i],c['bHL'])*1000/c['sPerYear']
            drho_dt[rho>=c['rho1']] = dr_dt[rho>=c['rho1']]                    
        
        elif c['physRho']=='HLSigfus':
            A = bdotSec*(1/t)*c['sPerYear']*c['rhoiMgm']
            drho_dt = np.zeros(gridLen)
            #if max(rho)>c['rho1']:
            f550 = interpolate.interp1d(rho,sigma)
            sigma550 = f550(c['rho1'])
            rhoDiff = (c['rhoiMgm']-rho/1000)

            drho_dt[rho<c['rho1']] = c['k1']*np.exp(-c['Q1']/(c['R']*Tz[rho<c['rho1']]))*(c['rhoiMgm']-rho[rho<c['rho1']]/1000)*np.power(A[i],c['aHL'])*1000/c['sPerYear']
            k = np.power(c['k2']*np.exp(-c['Q2']/(c['R']*Tz[rho>=c['rho1']])),2)/c['sPerYear']
            sigmaDiff = (sigma[rho>=c['rho1']]-sigma550)
            #if rhoDiff[rho>=c['rho1']]>0: #Max is not sure why this was coded in.
            drho_dt[rho>=c['rho1']] = k*(sigmaDiff*rhoDiff[rho>=c['rho1']])/(c['g']*np.log((c['rhoiMgm']-c['rho1']/1000)/(rhoDiff[rho>=c['rho1']])))
            #else:
                #drho_dt[rho>=c['rho1']] = 0
                        
        elif c['physRho'] =='Li2004': #Equation from Li Zwally, 2004 ...10/21/14 note: Max not sure if this is working properly
            drho_dt = np.zeros(gridLen)
            #for j in xrange(gridLen):
            dr_dt = (c['rhoi']-rho)*bdotSec[i]*(1/t)*c['sPerYear']*c['rhoiMgm']*(139.21-0.542*T10m)*8.36*(c['KtoC']-Tz)**-2.061
            drho_dt = dr_dt/c['sPerYear']
        
        elif c['physRho'] =='Li2011':
            dr_dt = np.zeros(gridLen)
            TmC=T10m-273.15
            A = bdotSec[i]*(1/t)*c['sPerYear']*c['rhoiMgm']
            beta1 = -9.788 + 8.996*A - 0.6165*TmC 
                #beta1 = -9.788 + 8.996*A*100 - 0.6165*T10m #scaled: A*100 so accumulation is cm/year, instead of m/year. Seems to work?
            beta2 = beta1/(-2.0178 + 8.4043*A - 0.0932*TmC) # this is the one from the paper. Does not work with scaled beta1. Pay attention to units in paper.
                #beta2 = (-9.788+8.966*A-0.6165*T10m)/(-2.0178 + 8.4043*A - 0.0932*T10m) #the one that seems to work: the published beta1, not the beta1 you need to use in model.
            dr_dt[rho<=c['rho1']] = (c['rhoi']-rho[rho<=c['rho1']])*A*beta1*8.36*(c['KtoC']-Tz[rho<=c['rho1']])**-2.061
            dr_dt[rho>c['rho1']] = (c['rhoi']-rho[rho>c['rho1']])*A*beta2*8.36*(c['KtoC']-Tz[rho>c['rho1']])**-2.061
            drho_dt = dr_dt/c['sPerYear']        
                                                                          
        elif c['physRho'] =='Helsen2008': #Equation from Arthern et al., 2010 ...10/21/14 Max not sure if this is working properly.
            #dr_dt = np.zeros(gridLen)
            dr_dt = (c['rhoi']-rho)*bdotSec[i]*(1/t)*c['sPerYear']*(76.138-0.28965*Ts[i])*8.36*(c['KtoC']-Tz)**-2.061
            drho_dt = dr_dt/c['sPerYear']
        
        elif c['physRho'] =='Arthern2010S': #this is the steady-state solution
            dr_dt = np.zeros(gridLen)
            dr_dt[rho<c['rho1']] = (c['rhoi']-rho[rho<c['rho1']])*c['ar1']*bdotSec[i]*(1/t)*c['sPerYear']*rho[rho<c['rho1']]*c['g']*np.exp(-c['Ec']/(c['R']*Tz[rho<c['rho1']])+c['Eg']/(c['R']*T10m))
            dr_dt[rho>=c['rho1']] = (c['rhoi']-rho[rho>=c['rho1']])*c['ar2']*bdotSec[i]*(1/t)*c['sPerYear']*rho[rho>=c['rho1']]*c['g']*np.exp(-c['Ec']/(c['R']*Tz[rho>=c['rho1']])+c['Eg']/(c['R']*T10m))
            drho_dt = dr_dt/c['sPerYear']
            
        elif c['physRho']=='Arthern2010T':
            if not c['physGrain']:
                print "Grain growth should be on for Arthern Transient"           
            dr_dt = np.zeros(gridLen)
            kc1=9.2e-9
            kc2=3.7e-9
            thresh=c['rho2']
            dr_dt[rho<thresh] = kc1*(c['rhoi']-rho[rho<thresh])*np.exp(-c['Ec']/(c['R']*Tz[rho<thresh]))*sigma[rho<thresh]/(r2[rho<thresh])
            dr_dt[rho>=thresh] = kc2*(c['rhoi']-rho[rho>=thresh])*np.exp(-c['Ec']/(c['R']*Tz[rho>=thresh]))*sigma[rho>=thresh]/(r2[rho>=thresh])
            drho_dt = dr_dt#/c['sPerYear']
        
        elif c['physRho'] =='Spencer2001':      # Uncommented out lines 464 - 475
            pass

#             drho_dt = np.zeros(gridLen)
#             for j in xrange(gridLen):
#                 if rho[j]<c['rho1']:
#                     dr_dt = rho[j]*c['C11']/c['sPerYear']*np.exp(-c['C12']/(c['R']*Tz[j]))*(1-rho[j]/c['rhoi'])*((1-(1-rho[j]/c['rhoi'])**c['C13'])**c['C14'])*sigma[j]**c['C15']
#                     drho_dt[j] = dr_dt/c['sPerYear']
#                 elif rho[j]<c['rho2']:
#                     dr_dt = rho[j]*c['C21']/c['sPerYear']*np.exp(-c['C22']/(c['R']*Tz[j]))*(1-rho[j]/c['rhoi'])*((1-(1-rho[j]/c['rhoi'])**c['C23'])**c['C24'])*sigma[j]**c['C25']
#                     drho_dt[j] = dr_dt/c['sPerYear']
#                 else:
#                     dr_dt = rho[j]*c['C31']/c['sPerYear']*np.exp(-c['C32']/(c['R']*Tz[j]))*(1-rho[j]/c['rhoi'])*((1-(1-rho[j]/c['rhoi'])**c['C33'])**c['C34'])*sigma[j]**c['C35']
#                     drho_dt[j] = dr_dt/c['sPerYear']
# 
#         elif c['physRho']=='Goujon2003':
#             drho_dt = np.zeros(gridLen)
#             D = rho/c['rhoi']
#             rhoC = c['rho2'] #should be Martinerie density
#             if max(rho)>c['rho2']:
#                 frho2 = interpolate.interp1d(rho,sigma)
#                 sigmarho2 = frho2(rhoC)
#             sigma_b = sigmarho2*(D*(1-rhoC/c['rhoi']))/(rhoC/c['rhoi']*(1-D))
#             sigmaEff = sigma + c['atmosP'] - sigma_b
#             for j in xrange(gridLen):
#                 if D[j]<0.6:
#                     gamma = 1e-14
#                     drho_dt[j]=(gamma*(sigma[j]/D[j]**2)*(1-(5/3)*D[j]))*c['rhoi']
#                 elif D[j]<0.9:
#                     A = 7.89e-15*np.exp(-c['Qgj']/(c['R']*Tz))
#                     Z0 = 2.0
#                     c = 15.5
#                     D0 = 0.00226*(Ts[i]-c['KtoC'])+0.647
#                     lp = (D[j]/D0)**(1.0/3.0)
#                     Z = Z0+c*(lp-1.0)
# #                     lpp = lp + (4.0*Z0*(lp-1.)**2.*(2.*lp+1.)+c(lp-1.)**3.*(3.*lp+1.))/(12.0*lp*(4.0*lp-2.*Z0*(lp-1.)-c*(lp-1.)**2.0))
# #                     a = (np.pi/(2*Z))
#                     a = 1.0
#                     sigmaStar = (4.0*np.pi*sigma[j])/(a*Z*D[j])
#                     drho_dt[j] = c['rhoi']*5.3*A[j]*((D[j]**2)*D0)** (1.0/3.0) * (a/np.pi)**0.5 *(sigmaStar/3.0)**c['n']
#                 elif D[j]<0.95:
#                     A = 7.89e-15*np.exp(-c['Qgj']/(c['R']*Tz))
#                     drho_dt[j] = (2*A[j]*((D[j]*(1-D[j]))/(1-(1-D[j])**(1/c['n']))**c['n']) * (2*sigmaEff[j]/c['n'])**c['n'])*c['rhoi']
#                 else:
#                     A = 1.2e-3*np.exp(-c['Qgj']/(c['R']*Tz))
#                     drho_dt[j] = 9/4*A[j]*(1-D[j])*sigmaEff[j]*c['rhoi']

        elif c['physRho']=='Morris2014':
            if spin and i==0:
                Hx=np.zeros(gridLen) #Hx is temperature history
                Hx=Hx+np.exp(-c['QMorris']/(c['R']*Tz))*dt #is this dt correct?
          
            
            #Hx_new=0
            #Hx = np.concatenate(([Hx_new],Hx[0:-1]))            
            #Hx=Hx+np.exp(-c['QMorris']/(c['R']*Tz))*dt #just an initial order of magnitude value for Hx  
            #
            drho_dt = np.zeros(gridLen)
            #dr_dt = drho_dt
            drho_dt[rho<=c['rho1']] = (c['kMorris']/(c['rhoW']*c['g'])) * ((c['rhoi']-rho[rho<=c['rho1']])/rho[rho<=c['rho1']]) * 1/Hx[rho<=c['rho1']] * np.exp(-c['QMorris']/(c['R']*Tz[rho<=c['rho1']]))*sigma[rho<=c['rho1']]
            #print 'max=', np.max(drho_dt)
            
            ### Choose which zone 2 physics you want
            
            #HL Dynamic
            #A = bdotSec*(1/t)*c['sPerYear']*c['rhoiMgm'] #A from the input json file is m ice equivalent.
            #drho_dt[rho>=c['rho1']] = c['k2']*np.exp(-c['Q2']/(c['R']*Tz[rho>=c['rho1']]))*(c['rhoiMgm']-rho[rho>=c['rho1']]/1000)*np.power(A[i],c['bHL'])*1000/c['sPerYear']
            
            #Arthern
            #drho_dt[rho>c['rho1']] = ((c['rhoi']-rho[rho>=c['rho1']])*c['ar2']*bdotSec[i]*(1/t)*c['sPerYear']*rho[rho>c['rho1']]*c['g']*np.exp(-c['Ec']/(c['R']*Tz[rho>c['rho1']])+c['Eg']/(c['R']*T10m)))/c['sPerYear']
            
            #HL Sigfus         
            f550 = interpolate.interp1d(rho,sigma)
            sigma550 = f550(c['rho1'])
            rhoDiff = (c['rhoiMgm']-rho/1000)                                     
            k = np.power(c['k2']*np.exp(-c['Q2']/(c['R']*Tz[rho>c['rho1']])),2)/c['sPerYear']
            sigmaDiff = (sigma[rho>c['rho1']]-sigma550)
            sigmaDiff[sigmaDiff<0]=0 #Stops model from running away - spin up artifact.
            dd=(k*(sigmaDiff*rhoDiff[rho>c['rho1']]) / (c['g']*np.log((c['rhoiMgm']-c['rho1']/1000)/(rhoDiff[rho>c['rho1']]))))
            #print 'sig1 = ', dd[0]             
            drho_dt[rho>c['rho1']] = (k*(sigmaDiff*rhoDiff[rho>c['rho1']]) / (c['g']*np.log((c['rhoiMgm']-c['rho1']/1000)/(rhoDiff[rho>c['rho1']])))) #use H&L Sigfus physics for zone 2.
            
            ##Li 2011
            #TmC=T10m-273.15
            #A = bdotSec[i]*(1/t)*c['sPerYear']*c['rhoiMgm']
            #beta1 = -9.788 + 8.996*A - 0.6165*TmC 
            #beta2 = beta1/(-2.0178 + 8.4043*A - 0.0932*TmC) # this is the one from the paper. Does not work with scaled beta1. Pay attention to units in paper.
            #drho_dt[rho>c['rho1']] = ((c['rhoi']-rho[rho>c['rho1']])*A*beta2*8.36*(c['KtoC']-Tz[rho>c['rho1']])**-2.061)/c['sPerYear']

#            rho[rho<c['rho1']]
#                    #equation 5 from Morris and Wingham
#                    dr_dt = c['kMorris']/(c['rhoW']*c['g'])*(c['rhoi']-rho[j])/rho[j]*1/H[j]*np.exp(-c['QMorris']/(R*Tz[j]))*sigma[j]
##                     dr_dt = (c['rhoi']-rho[j])*c['ar1']*bdotSec[i]*(1/t)*c['sPerYear']*rho[j]*c['g']*np.exp(-c['Ec']/(c['R']*Tz[j])+c['Eg']/(c['R']*T10m))
#                    drho_dt[j] = dr_dt
#                else: #what is this bit?
#                    pass
##                     dr_dt = (c['rhoi']-rho[j])*c['ar2']*bdotSec[i]*(1/t)*c['sPerYear']*rho[j]*c['g']*np.exp(-c['Ec']/(c['R']*Tz[j])+c['Eg']/(c['R']*T10m))
##                     drho_dt[j] = dr_dt/c['sPerYear']

            #Hx = Hx + np.exp(-c['QMorris']/(c['R']*Tz))*dt #equation 7 from Morris and Wingham
            #Hx_new=np.exp(-c['QMorris']/(c['R']*Tz[0]))
            #Hx = np.concatenate(([Hx_new],Hx[0:-1]))
            Hx_new=0
            Hx = np.concatenate(([Hx_new],Hx[0:-1]))            
            Hx=Hx+np.exp(-c['QMorris']/(c['R']*Tz))*dt #just an initial order of magnitude value for Hx  
            


                    
        elif c['physRho']=='Barnola1991':
            # print i
            rho[rho>c['rhoi']]=c['rhoi'] #The Barnola model will go a fraction over the ice density (order 10^-3), so this stops that.
            drho_dt = np.zeros(gridLen)
            D = rho/c['rhoi']
            nBa = c['n']*np.ones(gridLen)
            A0b =  c['A0Barnola']
            A0 = A0b*np.ones(gridLen)/1.e18 #this is for the n=3 region.
            #Vc = (6.95e-4)*T10m-0.043
            #rhoCRhoi = ((Vc+1/(c['rhoi']*(1e-3)))**-1)*1000/c['rhoi']
            #        
            #sigma_b = c['atmosP']*(D*(1-rhoCRhoi))/(rhoCRhoi*(1-D)) # This should be bubble pressure?
            #
            #sigmadiff = sigma - sigma_b
            #print sigmadiff
            #sigmaEff[sigmaEff<0]=1.e-15
            sigmaEff=sigma
            #nBa[sigmaEff<0.1e6]=1.
            #A0[sigmaEff<0.1e6]=A0b/1.e9#for the n=1 region
            
            #this is Max's vectorized version
            # zone 2   
            fe = 10.0**(c['alphaBarnola']*(rho[rho<=800.]/1000)**3.+c['betaBarnola']*(rho[rho<=800.]/1000)**2.+c['deltaBarnola']*rho[rho<=800.]/1000+c['gammaBarnola'])
            drho_dt[rho<=800.] = rho[rho<=800.]*A0[rho<=800.]*np.exp(-c['QBarnola']/(c['R']*Tz[rho<=800.]))*fe*(sigmaEff[rho<=800.]**nBa[rho<=800.])
            
            # zone 1
            A = bdotSec*(1/t)*c['sPerYear']*c['rhoiMgm']
            drho_dt[rho<c['rho1']] = dr_dt = c['k1']*np.exp(-c['Q1']/(c['R']*Tz[rho<c['rho1']]))*(c['rhoiMgm']-rho[rho<c['rho1']]/1000)*np.power(A[i],c['aHL'])*1000/c['sPerYear']
            
            # zone 3
            fs = (3./16.)*(1-rho[rho>800.]/c['rhoi'])/(1-(1-rho[rho>800.]/c['rhoi'])**(1./3.))**3.
            drho_dt[rho>800.] = rho[rho>800.]*A0[rho>800.]*np.exp(-c['QBarnola']/(c['R']*Tz[rho>800.]))*fs*(sigmaEff[rho>800.]**nBa[rho>800.])                    

        
        elif c['physRho']=='BarnolaJ':
            drho_dt = np.zeros(gridLen)
            D = rho/c['rhoi']
            nBa = c['n']*np.ones(gridLen)
            A0 =  c['A0Barnola']*np.ones(gridLen)
            Vc = (6.95e-4)*T10m-0.043
            rhoCRhoi = ((Vc+1/(c['rhoi']*(1e-3)))**-1)*1000/c['rhoi']
            sigma_b = np.zeros(gridLen)
            for j in xrange(gridLen):
                if D[j]<1:
                    sigma_b[j] = c['atmosP']*(D[j]*(1-rhoCRhoi))/(rhoCRhoi*(1-D[j]))
                    jMax = j
                else:
                    sigma_b[j] = sigma_b[jMax]
            sigmaEff = sigma - sigma_b
            for j in xrange(gridLen):
                if sigmaEff[j]<0:
                    sigmaEff[j]=0
                if sigmaEff[j]<0.1e6:
                    nBa[j] = 1.
                    A0[j] = A0[j]/1e8
                else:
                    A0[j] = A0[j]/1e18
                    
                if rho[j]<c['rho1']:
                    A = bdotSec*(1/t)*c['sPerYear']*c['rhoiMgm']
                    drho_dt[j] = dr_dt = c['k1']*np.exp(-c['Q1']/(c['R']*Tz[j]))*(c['rhoiMgm']-rho[j]/1000)*np.power(A[i],c['aHL'])*1000/c['sPerYear']
                elif rho[j]<800:
                    fe = 10.0**(c['alphaBarnola']*(rho[j]/1000)**3.+c['betaBarnola']*(rho[j]/1000)**2.+c['deltaBarnola']*rho[j]/1000+c['gammaBarnola'])
                    drho_dt[j] = rho[j]*A0[j]*np.exp(-c['QBarnola']/(c['R']*Tz[j]))*fe*(sigmaEff[j]**nBa[j])
                elif rho[j]<c['rhoi']:
#                     print (1.-(1.-rho/c['rhoi']))
                    denom = (1.-(1.-rho/c['rhoi'])**(1./3.))
                    #print denom
                    for jj in xrange(gridLen):
                        if math.isnan(denom[jj]):
                            denom[jj] = 1
                    fs = 3./16.*(1.- rho/c['rhoi'])/denom**3
                    drho_dt[j] = rho[j]*A0[j]*np.exp(-c['QBarnola']/(c['R']*Tz[j]))*fs[j]*(sigmaEff[j]**nBa[j])
        
        
        
        else:
            print 'Error: you need to choose model physics in json (check spelling)'
            sys.exit()
        
        #update the density using explicit method
        rho = rho + dt*drho_dt
        rho[rho>=c['rhoi']] = c['rhoi']-1.0e-10
#         plt.plot(rho,drho_dt)
#         plt.ylim(0,8e-7)
#         plt.show()

        #update the age (seconds)
#         age = np.insert(age,0,0) #insert new surface value
#         age= age[:gridLen]+dt #remove bottom element
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
            phi_s = Ts[i]
            phi_0 = Tz
            #print Tz

            K_ice = 9.828*np.exp(-0.0057*phi_0)
            K_firn = K_ice * (rhoP / 1000) ** (2 - 0.5 * (rhoP / 1000)) #
            c_firn = 152.5 + 7.122 * phi_0 # specific heat of ice
#             c_firn = c_ice * rho #schwander - specific heat of firn
            Gamma_P = K_firn/(c_firn * rho)


            Tz = transient_solve_TR(z_edges_vec,z_P_vec,nt,dt,Gamma_P, phi_0,nz_P,nz_fv,phi_s)
            Tz = np.concatenate(([Ts[i]],Tz[:-1]))
            
        #update the length of the boxes
        dzNew = bdotSec[i]*c['rhoi']/rhos0[i]*dt
        dz = mass/rho*dx
        dz = np.concatenate(([dzNew],dz[:-1]))
#         dz =dz.transpose()
        z = dz.cumsum(axis = 0)
        
        #max added:
        z = np.concatenate(([0],z[:-1])) 
        
        #if c['physRho']=='Morris2014':
        #    print 'i=', i
        #    Hx = Hx + np.exp(-c['QMorris']/(c['R']*Tz))*dt #equation 7 from Morris and Wingham
        #    Hx_new=np.exp(-c['QMorris']/(c['R']*Tz[0]))
        #    Hx = np.concatenate(([Hx_new],Hx[:,-1]))
        ####

        rho = np.concatenate(([rhos0[i]],rho[:-1]))
        
        if not spin:
            Dcon = np.concatenate(([D_surf[i]],Dcon[:-1])) #Dcon is a diffusivity tracker. Just creating a vector at each time step that is some fraction which is multiplied by diffusivity in firnair.

        #update mass
        massNew =bdotSec[i]*c['sPerYear']*c['rhoi']
        mass = np.concatenate(([massNew],mass[:-1]))

        sigma = mass * c['g'] * dx
        sigma = sigma.cumsum(axis = 0)
        

        
        #Update grain growth
        if c['physGrain']: #this is grain growth
            # Temp should be Temp at 10m depth.
            #dr2Dt = c['kg'] * np.exp(-c['Eg']/(c['R']*T10m))
            dr2Dt = c['kg'] * np.exp(-c['Eg']/(c['R']*Tz))
            r2 = r2 + dr2Dt * dt #* c['stpsPerYear'] #does this need the c['stps']? 11/19/14
            #Add time elapsed to beginning of row for output file
            r2_time = np.concatenate(([t*i + 1], r2))
            if c['calcGrainSize']: #This uses surface temperature to get an initial grain size. 
                r2 = np.concatenate(([-2.42e-9*Ts[i]+9.46e-7], r2[:-1]))
                if not spin:
                    print "grain size on", r2[0]
                    
            else:
                r2 = np.concatenate(([-2.42e-9*(c['Ts0'])+9.46e-7], r2[:-1]))

        if spin and i == (stp-1):
#             print t*i + 1
            rho_time = np.concatenate(([t*i + 1],rho))
            Tz_time = np.concatenate(([t*i + 1],Tz))
            age_time = np.concatenate(([t*i + 1],age))
            z_time = np.concatenate(([t*i + 1],z))
            r2_time = np.concatenate(([t*i + 1],r2))
            '''initialize files to write in time loop'''
            densityPath = os.path.join(c['resultsFolder'], 'densitySpin.csv')
            tempPath = os.path.join(c['resultsFolder'], 'tempSpin.csv')
            agePath = os.path.join(c['resultsFolder'], 'ageSpin.csv')
            depthPath = os.path.join(c['resultsFolder'], 'depthSpin.csv')
            if c['physGrain']:
                r2Path = os.path.join(c['resultsFolder'], 'r2Spin.csv')
                
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
            if c['physGrain']:
                with open(r2Path, "a") as f:
                    writer = csv.writer(f)
                    writer.writerow(r2)                
#         elif not spin and (t*i)%1 == 0:

        elif not spin and [True for jj in TWrite if jj == t*i+1] == [True]:
        #elif not spin and [True for jj in TWrite if jj == mtime] == [True]:
            rho_time = np.append(mtime,rho)
            Tz_time = np.append(mtime,Tz)
            age_time = np.append(mtime,age)
            z_time = np.append(mtime,z)
            Dcon_time = np.append(mtime,Dcon)
            if c['physGrain']:
                r2_time = np.append(mtime,r2)
            if (c['userInput']):
                Clim_time = np.append(mtime,[bdot[i],Ts[i]])
                
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
            with open(DconPath, "a") as f:
                writer = csv.writer(f)
                writer.writerow(Dcon_time)
            if (c['userInput']):
                with open(ClimPath, "a") as f:
                    writer = csv.writer(f)
                    writer.writerow(Clim_time)
            if c['physGrain']:
                with open(r2Path, "a") as f:
                    writer = csv.writer(f)
                    writer.writerow(r2_time)                
                          
                
            ### BCO,LIZ, and DIP ###
            #Vc = (6.95e-4)*T10m-0.043 #Martinerie et al., 1994, Eq. 2: critical pore volume at close off
            #bcoMartRho = ((Vc+1/(c['rhoi']*(1e-3)))**-1)*1000 # Martinerie density at close off
            bcoMartRho = 1/( 1/(917.0) + T10m*6.95E-7 - 4.3e-5) # Martinerie density at close off; see Buizert thesis (2011), Blunier & Schwander (2000), Goujon (2003)
            bcoAgeMart = min(age[rho>=bcoMartRho])/c['sPerYear'] # close-off age from Martinerie
            bcoDepMart = min(z[rho>=(bcoMartRho)])
            bcoAgeMartAll.append(bcoAgeMart) #age at the 815 density horizon
            bcoDepMartAll.append(bcoDepMart) #this is the 815 close off depth
            
            # bubble close-off age and depth assuming rho_crit = 815kg/m^3
            bcoAge815 = min(age[rho>=(c['rho2'])])/c['sPerYear'] #close-off age where rho=815 kg m^-3
            bcoDep815 = min(z[rho>=(c['rho2'])]) #depth of 815 horizon
            bcoAge815All.append(bcoAge815) #age at the 815 density horizon
            bcoDep815All.append(bcoDep815) #this is the 815 close off depth
            
            ### Lock-in depth and age
            LIZMartRho = bcoMartRho - 14.0 #LIZ depth (Blunier and Schwander, 2000)      
            LIZAgeMart = min(age[rho>LIZMartRho])/c['sPerYear'] #lock-in age
            LIZDepMart = min(z[rho>=(LIZMartRho)]) #lock in depth
            #lockIn = min(age[phiClosed[phiClosed<phi]/phi[phiClosed<phi] >= c['phiLockIn']])/c['sPerYear'] #lock-in age Old, from Jessica's original. Too shallow, likely.       
            LIZAgeAll.append(LIZAgeMart)
            LIZDepAll.append(LIZDepMart)
            
            ### Porosity, including DIP
            phi = 1-rho/c['rhoi'] # total porosity
            #phi[phi<=0]=1e-16 
            phiC = 1-bcoMartRho/c['rhoi']; #porosity at close off
            #phiClosed = np.ones(gridLen)
            phiClosed = 0.37*phi*(phi/phiC)**-7.6 #Closed porosity, from Goujon. See Buizert thesis (eq. 2.3) as well
            
            #for jj in range(gridLen): #porosity, from Goujon et al. 2003 (good explanation in Buizert thesis, eq. 2.3). Not sure why this is looped.
            #    if phi[jj]>phiC:
            #        phiClosed[jj] = c['gamma'] * phi[jj]*(phi[jj]/phiC)**-7.6
            #        jjMax = jj
            #    else:
            #        phiClosed[jj] = phiClosed[jjMax] #phiClosed is closed porosity
                    
            phiOpen = phi - phiClosed #open porosity
            phiOpen[phiOpen<=0] = 1.e-10 #don't want negative porosity.
                    
            intPhi = np.sum(phi * dz) #depth-integrated porosity
            intPhiAll.append(intPhi)      
                        
                         

            
#    plt.figure(1)
#    p1,= plt.plot(age/c['sPerYear'],rho)
#    p2,= plt.plot(agez0/c['sPerYear'],rhoz0,'--r')
#    plt.ylabel('Density (kg m$^{-3}$)',fontsize=14)
#    plt.xlabel('Age (year)',fontsize=14)
#    plt.savefig('rho_dep.png')
  
#    plt.figure(2)
#    p1,= plt.plot(age/c['sPerYear'],z)
#    p2,= plt.plot(agez0/c['sPerYear'],z0,'--r')
#    plt.ylabel('Depth (m)',fontsize=14)
#    plt.xlabel('Age (year)',fontsize=14)
#    plt.savefig('rho_dep.png')
  
    #plt.figure(3)
    #p1,= plt.plot(rho,z)
    #p2,= plt.plot(rhoz0,z0,'--r')
    #plt.xlabel('Density (kg m$^{-3}$)',fontsize=14)
    #plt.ylabel('Depth (m)',fontsize=14)
    #ax=plt.gca()
    #ax.set_ylim(ax.get_ylim()[::-1])
#    plt.savefig('rho_dep.png')
     
#    plt.figure(4)
#    p1,= plt.plot(Tz-c['KtoC'],z)
#    p2,= plt.plot(Tz0-c['KtoC'],z0,'--r')
#    plt.xlabel('Temperature (K)',fontsize=14)
#    plt.ylabel('Depth (m)',fontsize=14)
#    ax=plt.gca()
#    ax.set_ylim(ax.get_ylim()[::-1])
#    plt.xlim(c['Ts0']-c['KtoC']-20,c['Ts0']-c['KtoC']+20)
     
#    plt.show()
    
#     if not spin:
#         plt.figure(5)
#         p1,= plt.plot(TWrite[:len(bco815All)], bco815All)
#     #     p2,= plt.plot(Tz0+c['KtoC'],z0,'--r')
#         plt.xlabel('BCO (Age)',fontsize=14)
#         plt.ylabel('Time (years)',fontsize=14)
#         ax=plt.gca()
#         ax.set_ylim(ax.get_ylim()[::-1])
#     
#         plt.figure(6)
#         p1,= plt.plot(TWrite[:len(bco815All)], bcoDepAll)
#     #     p2,= plt.plot(Tz0+c['KtoC'],z0,'--r')
#         plt.xlabel('BCO Depth (m)',fontsize=14)
#         plt.ylabel('Time (years)',fontsize=14)
#         ax=plt.gca()
#         ax.set_ylim(ax.get_ylim()[::-1])
#         
# #         plt.figure(7)
#         
# #         p1,= plt.plot(TWrite[:len(intPhiAll)], intPhiAll)
# #     #     p2,= plt.plot(Tz0+c['KtoC'],z0,'--r')
# #         plt.xlabel('BCO (Age)',fontsize=14)
# #         plt.ylabel('Time (years)',fontsize=14)
# #         ax=plt.gca()
# #         ax.set_ylim(ax.get_ylim()[::-1])
#     
#         plt.figure(8)
#         p1,= plt.plot(TWrite[:len(bco815All)], lockInAll)
#     #     p2,= plt.plot(Tz0+c['KtoC'],z0,'--r')
#         plt.xlabel('Lock In Age (m)',fontsize=14)
#         plt.ylabel('Time (years)',fontsize=14)
#         ax=plt.gca()
#         ax.set_ylim(ax.get_ylim()[::-1])
#     
    #plt.show()
    
    elapsed=time.time()-tic
    elapsed_min=elapsed/60.
    mins=np.floor(elapsed_min)
    secs=(elapsed_min-mins)*60    

    if spin:
        logging.info("spin up took %s minutes %s seconds" % (mins, secs))
   
    if not spin:
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
        logpath = os.path.join(os.getcwd(),c['resultsFolder'])    
        shutil.copy('RUNDETAILS.log',logpath)
        #os.remove('RUNDETAILS.log')

    if c['plotting'] == 'on' and spin == 0:
        if os.path.exists(c['plotsFolder']):
            rmtree(c['plotsFolder'])
        os.makedirs(c['plotsFolder'])   
        plotData(c['physGrain'], c['resultsFolder'], c['plotsFolder'])

        
if __name__ == '__main__':
    logging.basicConfig(filename='RUNDETAILS.log',level=logging.DEBUG,filemode='w',format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    
    if len(sys.argv) >= 2 and '-s' not in sys.argv:
        configName = os.path.join(os.path.dirname(__file__), sys.argv[1])
    elif len(sys.argv) >= 2 and '-s' in sys.argv:
        configName = os.path.join(os.path.dirname(__file__), sys.argv[1])
    else:
#         configName = os.path.join(os.path.dirname(__file__), 'configs', 'configLi1.json')
        configName = os.path.join(os.path.dirname(__file__), 'config_test_input.json')

    if '-s' in sys.argv:
        spin = True
    else:
        spin = False

# This right now is for testing purposes.
    with open(configName, "r") as f:
        jsonString = f.read()
        c = json.loads(jsonString)

        
        #input_year, input_temp, input_AR = input_data.retrieve_data(c["InputFileName"])


    runModel(configName,spin)

    
    