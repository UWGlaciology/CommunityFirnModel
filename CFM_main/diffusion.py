#!/usr/bin/env python
'''
code to handle diffusion (heat, enthalpy)
calls solver
Draws from Numerical Heat Transfer and Heat Flow (Patankar, 1980)
Enthalpy from Voller, Swaminathan, and Thomas (1990)
Isotope diffusion now has its own class.
'''

from solver import transient_solve_TR
# from solver import transient_solve_EN_old
# from solver import transient_solve_EN_new
from solver import transient_solve_EN, Marshall, apparent_heat
from constants import *
import numpy as np
from scipy import interpolate

def firnConductivity(self,iii,K_ice):
    '''
    Function to set the firn's thermal conductivity
    based on one of a number of parameterizations.

    Choose your favorite!
    Default is Calonne et al. 2019.
    References are provided at the end of this script.
    '''
    if self.c['conductivity']=='Calonne2019':  #Calonne et al. 2019
        rho_transition = 450.0 #[kg/m^3]
        a = 0.02 #[m^3/kg]
        theta       = 1 / (1 + np.exp(-2*a*(self.rho - rho_transition)))
        kref_firn   = 2.107 + 0.003618 * (self.rho - RHO_I)
        kref_snow   = 0.024 - 1.23e-4 * self.rho + 2.5e-6 * self.rho**2
        kref_i      = 2.107 # [W/m/K]
        kref_a      = 0.024 # [W/m/K] 
        K_air       = kref_a # use this for now; at some point find equation for T-dependence of air
        K_firn      = (1-theta) * K_ice*K_air/(kref_i*kref_a) * kref_snow + theta * K_ice/kref_i * kref_firn # equation 5
    elif self.c['conductivity']=='Schwander':
        K_firn  = K_ice * (self.rho/RHO_I) ** (2 - 0.5 * (self.rho/RHO_I))    # Schwander 1997, eq. A11
    elif self.c['conductivity']=='Yen_fixed':
        K_firn  = 2.22362 * (self.rho / 1000)**1.885                          # Yen 1981, eq 34 w/ fixed K_ice (original)
    elif self.c['conductivity']=='Yen_var':
        K_firn  = K_ice * (self.rho / 1000)**1.885                            # Yen 1981, modified for variable K_ice
    elif self.c['conductivity']=='Anderson':
        K_firn  = 0.021 + 2.5 * (self.rho/1000.)**2                           # Anderson (1976)
    elif self.c['conductivity']=='Yen_b':
        K_firn  = 0.0688 * np.exp(0.0088*(self.Tz-273.15) + 4.6682*self.rho/1000) # Yen 1981, eq. 35.
    elif self.c['conductivity']=='Sturm':
        K_firn  = 0.138 - 1.01*(self.rho/1000) + 3.233*(self.rho/1000)**2     # Sturm, 1997.; rho < 0.6
    elif self.c['conductivity']=='VanDusen':
        K_firn  = 2.1e-2 + 4.2e-4 * self.rho + 2.2e-9 * (self.rho)**3         # Van Dusen 1929 (via C&P)
    elif self.c['conductivity']=='Schwerdtfeger':
        K_firn  = (2 * K_ice * self.rho) / (3*RHO_I - self.rho)               # Schwerdtfeger (via C&P)
    elif self.c['conductivity']=='Riche':
        K_firn  = 3.e-6 * self.rho**2 - 1.06e-5 * self.rho + 0.024            # Riche and Schneebeli 2013 eq. 10
    elif self.c['conductivity']=='Jiawen':
        K_firn  = 0.0784 + 2.697 * (self.rho/1000.)**2                        # Jiawen 1991 eq. 3 
    elif self.c['conductivity']=='Calonne2011':
        K_firn  = 0.024 - 1.23e-4 * self.rho + 2.5e-6 * self.rho**2           # Calonne et al. 2011
    elif self.c['conductivity'] =='mix':
        if iii==0:
            print('Mixed conductivity (dig into code for details')
        K_firn = np.zeros_like(self.rho)
        Kdict = {}
        Kdict['Sturm']      = 0.138 - 1.01*(self.rho/1000) + 3.233*(self.rho/1000)**2
        Kdict['Anderson']   = 0.021 + 2.5 * (self.rho/1000.)**2
        K_firn[self.z<0.2]  = Kdict['Sturm'][self.z<0.2]
        K_firn[self.z>=0.3] = Kdict['Anderson'][self.z>=0.3]
        Kcond = ((self.z>=0.2) & (self.z<0.3))
        K_firn[Kcond]       = (Kdict['Sturm'][Kcond] + Kdict['Anderson'][Kcond])/2
    else:
        if iii==0:
            print('Conductivity is not set to one of the values; using Calonne (2019)')
        # K_firn = 0.021 + 2.5 * (self.rho/1000.)**2                           # Anderson (1976)
        rho_transition = 450.0 #[kg/m^3]
        a = 0.02 #[m^3/kg]
        theta       = 1 / (1 + np.exp(-2*a*(self.rho - rho_transition)))
        kref_firn   = 2.107 + 0.003618 * (self.rho - RHO_I)
        kref_snow   = 0.024 - 1.23e-4 * self.rho + 2.5e-6 * self.rho**2
        kref_i      = 2.107 # [W/m/K]
        kref_a      = 0.024 # [W/m/K] 
        K_air       = kref_a # use this for now; at some point find equation for T-dependence of air
        K_firn      = (1-theta) * K_ice*K_air/(kref_i*kref_a) * kref_snow + theta * K_ice/kref_i * kref_firn # equation 5

    return K_firn
##########################

def heatDiff(self,iii):
    '''
    Heat diffusion function

    :param z:
    :param dz:
    :param Ts:
    :param rho:

    :returns self.Tz:
    :returns self.T10m:
    
    thermal diffusivity: alpha = K_firn / (rho*c_firn)
    '''

    nz_P            = len(self.z)
    nz_fv           = nz_P - 2
    nt              = 1

    z_edges_vec1 = self.z[0:-1] + np.diff(self.z) / 2
    z_edges_vec = np.concatenate(([self.z[0]], z_edges_vec1, [self.z[-1]]))
    z_P_vec     = self.z
    
    phi_s           = self.Tz[0]
    phi_0           = self.Tz

    K_ice           = 9.828 * np.exp(-0.0057 * phi_0) # thermal conductivity, Cuffey and Paterson, eq. 9.2 (Yen 1981)
    c_firn          = 152.5 + 7.122 * phi_0 # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
    # c_firn        = CP_I # If you prefer a constant specific heat.

    K_firn = firnConductivity(self,iii,K_ice) # thermal conductivity

    if self.c['MELT']:
        try:
            if self.c['LWCheat']=='lowK':
                K_firn[self.LWC>0]=K_firn[self.LWC>0]/1.e4
        except:
            pass

    Gamma_P         = K_firn

    tot_rho         = self.rho
    c_vol           = self.rho * c_firn

    self.Tz         = transient_solve_TR(z_edges_vec, z_P_vec, nt, self.dt[iii], Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol)

    self.T10m       = self.Tz[np.where(self.z>=10.0)[0][0]]

    if self.c['MELT']:
        if self.c['LWCheat']=='effectiveT':
            pass

        elif np.any(self.Tz>273.1500001):
            print(f'WARNING: TEMPERATURE EXCEEDS MELTING TEMPERATURE at {iii}')
            print('WARM TEMPERATURES HAVE BEEN SET TO 273.15; MODEL RUN IS CONTINUING')

        self.Tz[self.Tz>=273.15]=273.15

    return self.Tz, self.T10m

##########################
### end heat diffusion ###
##########################

def enthalpyDiff(self,iii):
    '''
    enthalpy diffusion function, new
    1/30/19 - method from Voller and Swaminathan
    LWC is in volume (m^3)
    thermal diffusivity: alpha = K_firn / (rho*c_firn)
    '''
    Tstart          = self.Tz.copy()
    nz_P            = len(self.z) - 1
    nz_fv           = nz_P - 2

    if np.any(self.LWC>0): # this behavior is depricated; keeping code for now. (6/16/21)
        nt = 10 # number of iterations for the solver
    else:
        nt = 1

    z_edges_vec1    = self.z[0:-1] + np.diff(self.z) / 2
    z_edges_vec     = np.concatenate(([self.z[0]], z_edges_vec1, [self.z[-1]]))
    z_P_vec         = self.z

    # z_P_vec = (self.z[1:]+self.z[:-1])/2
    # z_edges_vec = self.z

    phi_s           = self.Tz[0] - T_MELT # work in [C] so that reference Temperature is 0 for enthalpy
    phi_0           = self.Tz - T_MELT
  
    vol_ice     = self.mass / RHO_I     # volume of the ice portion of each volume
    vol_tot     = vol_ice + self.LWC    # total volume of ice and liquid in each volume
    mass_liq    = self.LWC * RHO_W_KGM  # mass of liquid water
    rho_liq_eff = mass_liq / self.dz      # effective density of the liquid portion
    tot_rho     = (self.mass + mass_liq) / self.dz # 'total' density of volume (solid plus liquid)
    g_liq_1     = self.LWC / vol_tot     # liquid volume fraction (of the material portion, porosity ignored)
    g_ice_1     = vol_ice / vol_tot     # solid/ice volume fraction 

    K_water = 0.55575                         # thermal conductivity, water (W/m/K)
    K_ice   = 9.828 * np.exp(-0.0057 * self.Tz) # thermal conductivity, ice (W/m/K), Cuffey and Paterson, eq. 9.2 (Yen 1981)
    # K_mix = g_liq_1*K_liq + g_ice_1*K_ice

    ### Specific Heats
    c_firn          = 152.5 + 7.122 * self.Tz # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
    # c_firn  = CP_I # If you prefer a constant specific heat
    c_ice = c_firn
    c_liq = 4219.9 # J/kg/K, taken from engineeringtoolbox.com. Ha!
    # c_vol = g_ice_1 * RHO_I * c_ice + g_liq_1 * RHO_W_KGM * c_liq #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)
    c_vol = (g_ice_1 * c_ice + g_liq_1 * c_liq) * tot_rho #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)

    ### Conductivity
    K_firn = firnConductivity(self,iii,K_ice) # thermal conductivity [W/m/K]

    K_liq = K_water * (rho_liq_eff/1000)**1.885 # I am assuming that conductivity of water in porous material follows a similar relationship to ice.
    K_eff = g_liq_1*K_liq + g_ice_1*K_firn # effective conductivity
    
    ICT = 0 #Iteration Count Threshold (deprecated)

    ### Total enthalpy before solver (for testing conservation)
    tot_heat_pre = np.sum(CP_I_kJ*self.mass*self.Tz + T_MELT*CP_W/1000*self.LWC*RHO_W_KGM + LF_I_kJ*self.LWC*RHO_W_KGM)

    lwc_old = self.LWC.copy()
    
    phi_ret, g_liq, count, iterdiff,g_sol   = transient_solve_EN(z_edges_vec, z_P_vec, nt, self.dt[iii], K_eff, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol, self.LWC, self.mass, self.dz,ICT,self.rho,iii)
    # phi_ret, g_liq, count, iterdiff,g_sol   = Marshall(z_edges_vec, z_P_vec, nt, self.dt[iii], K_eff, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol, self.LWC, self.mass, self.dz,ICT,self.rho,iii)
    
    ### Below for testing code where firn layers are the finite volumes.
    # phi_ret, g_liq, count, iterdiff,g_sol   = transient_solve_EN(z_edges_vec, z_P_vec, nt, self.dt[iii], K_eff[0:-1], phi_0[0:-1], nz_P, nz_fv, phi_s, tot_rho[0:-1], c_vol[0:-1], self.LWC[0:-1], self.mass[0:-1], self.dz[0:-1],ICT,self.rho[0:-1],iii)
    # phi_ret = np.append(phi_ret,phi_ret[-1])
    # g_liq = np.append(g_liq,g_liq[-1])
    # g_sol= np.append(g_sol,g_sol[-1])

    LWC_ret = g_liq * self.dz
    # self.LWC        = g_liq * vol_tot


    delta_mass_liq  = mass_liq - (LWC_ret * RHO_W_KGM)
    dml_sum = 0.0 

    self.LWC = LWC_ret.copy()
    self.Tz = phi_ret + 273.15
    self.T10m       = self.Tz[np.where(self.z>=10.0)[0][0]]

    ### Total enthalpy after solver (for testing conservation)
    tot_heat_post = np.sum(CP_I_kJ*self.mass*self.Tz + T_MELT*CP_W/1000*self.LWC*RHO_W_KGM + LF_I_kJ*self.LWC*RHO_W_KGM)

    if (np.abs(tot_heat_post-tot_heat_pre)/tot_heat_pre)>1e-2:
        print(f'change in enthalpy at iteration {iii}!')
        print('pre:', tot_heat_pre)
        print('post:', tot_heat_post)
        ediff = (tot_heat_post-tot_heat_pre)                
        print('difference (kJ):', (tot_heat_post-tot_heat_pre))
        print('difference %:', ediff/tot_heat_pre)

    if np.any(self.Tz>273.1500001):
        print('WARNING: TEMPERATURE EXCEEDS MELTING TEMPERATURE')
        print('Maximal temperature was:',np.max(self.Tz),' at layers:',np.where(self.Tz == np.max(self.Tz)))
        print('iii, modeltime', iii, self.modeltime[iii])
        print('WARM TEMPERATURES HAVE BEEN SET TO 273.15; MODEL RUN IS CONTINUING')
    self.Tz[self.Tz>=273.15]=273.15

    if np.any(delta_mass_liq<0):
        if np.any(np.abs(delta_mass_liq[delta_mass_liq<0])>1e-7):
            print('------')

            print('If you are seeing this message there was a liquid mass gain in diffusion.') 
            print('Please email maxstev@umd.edu so I can fix it.')

        dml_sum = np.sum(delta_mass_liq[delta_mass_liq<0])
    
    delta_mass_liq  = np.maximum(delta_mass_liq,0) # fix for numerical instabilities with small time steps.
    self.mass       = self.mass + delta_mass_liq
    self.rho        = self.mass/self.dz

    tot_mass_new = self.mass + self.LWC*1000

    return self.Tz, self.T10m, self.rho, self.mass, self.LWC, dml_sum

##############################
### end enthalpy diffusion ###
##############################

def heatDiff_highC(self,iii):

    '''
    IN DEVELOPMENT

    One way of dealing with liquid water in the firn
    is to just set the heat capacity to be very high. 
    '''
    if iii==0:
        print('WARNING: heatDiff_highC IS IN DEVELOPMENT')

    nz_P            = len(self.z)
    nz_fv           = nz_P - 2
    nt              = 1

    z_edges_vec1 = self.z[0:-1] + np.diff(self.z) / 2
    z_edges_vec = np.concatenate(([self.z[0]], z_edges_vec1, [self.z[-1]]))
    z_P_vec     = self.z

    dt_sub = self.dt[iii]

    phi_s           = self.Tz[0]
    phi_0           = self.Tz

    g_liq = self.LWC/self.dz
    H_L_liq = RHO_W_KGM*LF_I #volumetric latent enthalpy [J/m3]
    H_lat = g_liq*H_L_liq

    K_ice           = 9.828 * np.exp(-0.0057 * phi_0) # thermal conductivity, Cuffey and Paterson, eq. 9.2 (Yen 1981)
    c_firn          = 152.5 + 7.122 * phi_0 # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
    # c_firn        = CP_I # If you prefer a constant specific heat.
    
    K_firn = firnConductivity(self,iii,K_ice) # thermal conductivity

    Gamma_P         = K_firn
    tot_rho         = self.rho
    c_vol_0         = self.rho * c_firn 
    
    C_lat = H_lat/T_MELT

    c_vol = c_vol_0 + C_lat


    self.Tz        = transient_solve_TR(z_edges_vec, z_P_vec, nt, dt_sub, Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol)
    self.T10m       = self.Tz[np.where(self.z>=10.0)[0][0]]

    self.Tz, self.LWC, self.rho, self.mass, refrozen_mass = LWC_correct(self)

    if np.any(self.Tz>273.15001):
        print('Tz higher than 273.15001')
        iHT = np.where(self.Tz>273.15001)[0]
        print(f'iHT: {iHT}')
        print(f'HT: {self.Tz[iHT]}')

    self.Tz[self.Tz>=273.15]=273.15

    dml_sum = 0

    return self.Tz, self.T10m, self.rho, self.mass, self.LWC, dml_sum

##############################
### end highC diffusion ###
##############################

def heatDiff_Teff(self,iii):
    '''    
    IN DEVELOPMENT
    artificially set the temperature of volumes with liquid water 
    to be higher than T_melt 
    '''

    if iii==0:
        print('WARNING: heatDiff_Teff IS IN DEVELOPMENT')

    nz_P            = len(self.z)
    nz_fv           = nz_P - 2
    nt              = 1

    z_edges_vec1 = self.z[0:-1] + np.diff(self.z) / 2
    z_edges_vec = np.concatenate(([self.z[0]], z_edges_vec1, [self.z[-1]]))
    z_P_vec     = self.z

    Q = LF_I * self.LWC * RHO_W_KGM
    deltaT = Q / (self.mass*CP_I)

    # Tc = self.Tz - 273.15
    # c0 = deltaT>-Tc
    # deltaT[c0] = -Tc[c0]

    Tz_eff = self.Tz + deltaT
    
    # phi_s           = self.Tz[0]
    # phi_0           = self.Tz

    phi_s           = Tz_eff[0]
    phi_0           = Tz_eff

    K_ice           = 9.828 * np.exp(-0.0057 * self.Tz) # thermal conductivity, Cuffey and Paterson, eq. 9.2 (Yen 1981)
    c_firn          = 152.5 + 7.122 * phi_0 # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
    # c_firn        = CP_I # If you prefer a constant specific heat.
    lwc_layers  = np.where(self.LWC>0)[0]
    
    K_firn = firnConductivity(self,iii,K_ice) # thermal conductivity

    Gamma_P         = K_firn
    tot_rho         = self.rho
    c_vol           = self.rho * c_firn

    T_eff_new         = transient_solve_TR(z_edges_vec, z_P_vec, nt, self.dt[iii], Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol)

    excessT = np.maximum(0.0,(T_eff_new - T_MELT))
    LWC_new = (excessT * self.mass * CP_I)/ LF_I / 1000 # divide by 1000 to put in volume (m3)
    self.LWC = LWC_new
    T_eff_new[self.LWC>0] = T_MELT
    self.Tz = T_eff_new

    self.T10m       = self.Tz[np.where(self.z>=10.0)[0][0]]

    # self.Tz, self.LWC, self.rho, self.mass, refrozen_mass = LWC_correct(self)

    if np.any(self.Tz>273.15001):
        print('Tz higher than 273.15001')
        iHT = np.where(self.Tz>273.15001)[0]
        print(f'iHT: {iHT}')
        print(f'HT: {self.Tz[iHT]}')

    self.Tz[self.Tz>=273.15]=273.15

    dml_sum = 0

    return self.Tz, self.T10m, self.rho, self.mass, self.LWC, dml_sum

##############################
### end T-eff diffusion ###
##############################

def heatDiff_LWCcorr(self,iii, iters,correct_therm_prop):
    '''

    IN DEVELOPMENT  

    just run the heat diffusion as normal and then balance energy.  
    '''

    if iii==0:
        print('WARNING: heatDiff_LWCcorr IS IN DEVELOPMENT')

    nz_P            = len(self.z)
    nz_fv           = nz_P - 2
    nt              = 1

    z_edges_vec1 = self.z[0:-1] + np.diff(self.z) / 2
    z_edges_vec = np.concatenate(([self.z[0]], z_edges_vec1, [self.z[-1]]))
    z_P_vec     = self.z
    
    # iters = 24
    dt_sub = self.dt[iii]/iters

    for jj in range(iters):

        phi_s           = self.Tz[0]
        phi_0           = self.Tz

        ###
        if correct_therm_prop:
            vol_ice     = self.mass / RHO_I     # volume of the ice portion of each volume
            vol_tot     = vol_ice + self.LWC    # total volume of ice and liquid in each volume
            mass_liq    = self.LWC * RHO_W_KGM  # mass of liquid water
            rho_liq_eff = mass_liq / self.dz      # effective density of the liquid portion
            tot_rho     = (self.mass + mass_liq) / self.dz # 'total' density of volume (solid plus liquid)
            g_liq_1     = self.LWC / vol_tot     # liquid volume fraction (of the material portion, porosity ignored)
            g_ice_1     = vol_ice / vol_tot     # solid/ice volume fraction 

            K_water = 0.55575                         # thermal conductivity, water (W/m/K)
            K_ice   = 9.828 * np.exp(-0.0057 * self.Tz) # thermal conductivity, ice (W/m/K), Cuffey and Paterson, eq. 9.2 (Yen 1981)
            # K_mix = g_liq_1*K_liq + g_ice_1*K_ice

            c_firn          = 152.5 + 7.122 * self.Tz # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
            # c_firn  = CP_I # If you prefer a constant specific heat
            c_ice = c_firn
            c_liq = 4219.9 # J/kg/K, taken from engineeringtoolbox.com. Ha!
            # c_vol = g_ice_1 * RHO_I * c_ice + g_liq_1 * RHO_W_KGM * c_liq #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)
            c_vol = (g_ice_1 * c_ice + g_liq_1 * c_liq) * tot_rho #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)

            K_firn = firnConductivity(self,iii,K_ice) # thermal conductivity

            K_liq = K_water * (rho_liq_eff/1000)**1.885 # I am assuming that conductivity of water in porous material follows a similar relationship to ice.
            K_eff = g_liq_1*K_liq + g_ice_1*K_firn # effective conductivity
            Gamma_P = K_eff
        ###

        ###
        else:
            K_ice           = 9.828 * np.exp(-0.0057 * phi_0) # thermal conductivity, Cuffey and Paterson, eq. 9.2 (Yen 1981)
            c_firn          = 152.5 + 7.122 * phi_0 # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
            # c_firn        = CP_I # If you prefer a constant specific heat.
            
            K_firn = firnConductivity(self,iii,K_ice) # thermal conductivity

            Gamma_P         = K_firn
            tot_rho         = self.rho
            c_vol           = self.rho * c_firn


        self.Tz        = transient_solve_TR(z_edges_vec, z_P_vec, nt, dt_sub, Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol)
        self.T10m       = self.Tz[np.where(self.z>=10.0)[0][0]]

        self.Tz, self.LWC, self.rho, self.mass, refrozen_mass = LWC_correct(self)

        if np.any(self.Tz>273.15001):
            print('Tz higher than 273.15001')
            iHT = np.where(self.Tz>273.15001)[0]
            print(f'iHT: {iHT}')
            print(f'HT: {self.Tz[iHT]}')

        self.Tz[self.Tz>=273.15]=273.15

    dml_sum = 0

    return self.Tz, self.T10m, self.rho, self.mass, self.LWC, dml_sum

##############################
### end T-eff diffusion ###
##############################

def LWC_correct(self):
    '''
    *** TEST FUNCTION ***
    If there is LWC in a layer after temperature diffusion and the temperature
    is less than zero, one option is to just balance the energy to increase the
    temperature and lower the LWC. It isn't the best way to solve the problem 
    but it is one way. 

    This should be vectorized but that is not a priority.
    '''

    ind_wetcold = np.where((self.Tz<T_MELT) & (self.LWC>0))[0]
    refrozen_mass = np.zeros_like(self.rho)
    if ind_wetcold.size!=0:
        cold_content = CP_I * self.mass * (T_MELT - self.Tz) # [J]
        ### LWC is volume (m^3)
        heattofreeze = self.LWC * RHO_W_KGM * LF_I # [J]
        
        for kk in ind_wetcold:
            if cold_content[kk] < heattofreeze[kk]:
                # not enough cold content
                # temperature raised to T_MELT
                # some water refreeze to bring T to T_MELT

                self.Tz[kk] = T_MELT
                self.LWC[kk] = self.LWC[kk] - (cold_content[kk]/1000/LF_I)
                refrozen_mass[kk] = cold_content[kk]/LF_I
                self.mass[kk] = self.mass[kk] + refrozen_mass[kk]
                self.rho[kk] = self.mass[kk]/self.dz[kk]
                # self.LWC[kk] = self.LWC[kk] - (cold_content[kk]/1000/LF_I)
            else: #enough cold content, all LWC refreezes
                # Temperature is raised from refreezing
                refrozen_mass[kk] = self.LWC[kk] * 1000
                self.LWC[kk] = 0
                self.Tz[kk] = self.Tz[kk] + heattofreeze[kk]/CP_I/self.mass[kk]
                self.mass[kk] = self.mass[kk] + refrozen_mass[kk]
                self.rho[kk] = self.mass[kk]/self.dz[kk]
                if self.Tz[kk]>273.15001:
                    print(f'kk is {kk}, Tz is {self.Tz[kk]}')
                

                # self.Tz[kk] = self.Tz[kk] + (heattofreeze[kk]/1000/LF_I)
        # print(f'max Tz is {np.max(self.Tz)}')
        if np.any(self.LWC<0):
            print("negative LWC from correction (LWC_correct function)")
            self.LWC[self.LWC<0] = 0
        if np.any(self.Tz>273.15001):
            print("temps above T_MELT from correction (LWC_correct function)")
            print()
            # print(self.Tz[self.Tz>T_MELT])
        self.Tz[self.Tz>T_MELT] = T_MELT
    return self.Tz, self.LWC, self.rho, self.mass, refrozen_mass


#### Radiation penetration (in development)
# def rad_pen(self,E_rp):

#     def exco(rho):
#     return -0.0338*rho +33.54

#     k_ex = exco(self.rho)
#     c_firn    = 0.021 + 2.5 * (self.rho/1000.)**2
#     deltaE_layers = E_rp * np.exp(-k_ex*self.z)
#     deltaT = deltaE_layers/(self.mass*CP_I)
#     self.Tz[1:] = self.Tz[1:]+deltaT


'''
### References for conductivity parameterizations ###
# Also can refer to Physics of Glaciers, chapter 9.2

Anderson EA (1976) A point energy and mass balance model of a snow cover. (doi:10.1016/S0074-6142(99)80039-4)
Brandt RE and Warren SG (1997) Temperature measurements and heat transfer in near-surface snow at the South Pole. J. Glaciol. 43(144), 339–351
Calonne, N., Flin, F., Morin, S., Lesaffre, B., du Roscoat, S. R., & Geindreau, C. (2011). Numerical and experimental investigations of the effective thermal conductivity of snow. Geophysical Research Letters, 38, L23501. https://doi.org/10.1029/2011GL049234
Calonne, N., Milliancourt, L., Burr, A., Philip, A., Martin, C. L., Flin, F., & Geindreau, C. (2019). Thermal conductivity of snow, firn, and porous ice from 3-D image-based computations. Geophysical Research Letters, 46, 13,079–13,089. https://doi. org/10.1029/2019GL085228
Jiawen R, Dahe Q and Maohuan H (1991) THERMAL PROPERTIES AND TEMPERATURE DISTRIBUTION OF SNOW/FIRN ON THE LAW DOME ICE CAP, ANTARCTICA. Antarct. Res. 2(2), 38–46
Lüthi MP and Funk M (2001) Modelling heat flow in a cold, high-altitude glacier: Interpretation of measurements from Colle Gnifetti, Swiss Alps. J. Glaciol. 47(157), 314–324 (doi:10.3189/172756501781832223)
Riche F and Schneebeli M (2013) Thermal conductivity of snow measured by three independent methods and anisotropy considerations. Cryosphere 7(1), 217–227 (doi:10.5194/tc-7-217-2013)
Schwander J, Sowers T, Barnola J-M, Blunier T, Fuchs A and Malaizé B (1997) Age scale of the air in the summit ice: Implication for glacial-interglacial temperature change. J. Geophys. Res. Atmos. 102(D16), 19483–19493 (doi:10.1029/97JD01309)
Schwerdtfeger P (1963) Theoretical derivation of the thermal conductivity and diffusivity of snow. IAHS Publ 61, 75–81 http://iahs.info/uploads/dms/061007.pdf
Sturm M, Holmgren J, König M and Morris K (1997) The thermal conductivity of seasonal snow. J. Glaciol. 43(143), 26–41 (doi:10.1017/S0022143000002781)
Van Dusen MS (1929) Thermal conductivity of non-metallic solids. International critical tables of numerical data, physics, chemistry and technology. McGraw-Hill New York, 216–217
Yen Y-C (1981) Review of Thermal Properties of Snow, Ice, and Sea Ice. CRREL Rep. 81-10, 1–27 http://acwc.sdp.sirsi.net/client/search/asset/1005644

'''