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
from solver import transient_solve_EN
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
        K_firn  = 0.0688 * np.exp(0.0088*(phi_0-273.15) + 4.6682*self.rho/1000) # Yen 1981, eq. 35.
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

        elif np.any(self.Tz>273.15):
            print('WARNING: TEMPERATURE EXCEEDS MELTING TEMPERATURE')
            print('Maximum temperature was:',np.max(self.Tz),' at layers:',np.where(self.Tz == np.max(self.Tz)))
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
    nz_P            = len(self.z)
    nz_fv           = nz_P - 2

    if np.any(self.LWC>0):
        nt              = 3 # number of iterations for the solver
    else:
        nt = 1

    z_edges_vec1    = self.z[0:-1] + np.diff(self.z) / 2
    z_edges_vec     = np.concatenate(([self.z[0]], z_edges_vec1, [self.z[-1]]))
    z_P_vec         = self.z
    
    phi_s           = self.Tz[0]
    phi_0           = self.Tz

    H_ice = RHO_I * CP_I * (self.Tz - T_MELT)
    H_liq = RHO_W_KGM * LF_I
  
    vol_ice     = self.mass / RHO_I     # volume of the ice portion of each volume
    vol_tot     = vol_ice + self.LWC    # total volume of ice and liquid in each volume
    mass_liq    = self.LWC * RHO_W_KGM  # mass of liquid water
    rho_liq_eff = mass_liq/self.dz      # effective density of the liquid portion
    tot_rho     = (self.mass + mass_liq) / self.dz # 'total' density of volume (solid plus liquid)
    g_liq_1     = self.LWC / vol_tot     # liquid volume fraction (of the material portion, porosity ignored)
    g_ice_1     = vol_ice / vol_tot     # solid/ice volume fraction 
    g_liq       = self.LWC / self.dz    # liquid volume fraction, total volume
    g_ice       = vol_ice / self.dz     # solid/ice volume fraction, total volume

    # H_tot = H_ice * g_ice + H_liq * g_liq #total enthalpy (don't need?)

    K_water = 0.55575                         # thermal conductivity, water (W/m/K)
    K_ice   = 9.828 * np.exp(-0.0057 * phi_0) # thermal conductivity, ice (W/m/K), Cuffey and Paterson, eq. 9.2 (Yen 1981)
    # K_mix = g_liq_1*K_liq + g_ice_1*K_ice

    c_firn          = 152.5 + 7.122 * phi_0 # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
    # c_firn  = CP_I # If you prefer a constant specific heat
    c_ice = c_firn
    c_liq = 4219.9 # J/kg/K, taken from engineeringtoolbox.com. Ha!
    # c_vol = g_ice_1 * RHO_I * c_ice + g_liq_1 * RHO_W_KGM * c_liq #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)
    c_vol = (g_ice_1 * c_ice + g_liq_1 * c_liq) * tot_rho #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)

    K_firn = firnConductivity(self,iii,K_ice) # thermal conductivity
    K_liq = K_water * (rho_liq_eff/1000)**1.885 # I am assuming that conductivity of water in porous material follows a similar relationship to ice.
    K_eff = g_liq_1*K_liq + g_ice_1*K_firn # effective conductivity
    
    Gamma_P = K_eff

    deltaH = RHO_W_KGM * LF_I #Voller 1990, eq 11. (letting T_ref = T_melt) units: J/m^3

    self.Tz, g_liq   = transient_solve_EN(z_edges_vec, z_P_vec, nt, self.dt[iii], Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol, g_liq, deltaH)

    self.T10m       = self.Tz[np.where(self.z>=10.0)[0][0]]

    if np.any(self.Tz>273.16):
        print('WARNING: TEMPERATURE EXCEEDS MELTING TEMPERATURE')
        print('Maximal temperature was:',np.max(self.Tz),' at layers:',np.where(self.Tz == np.max(self.Tz)))
        print('iii, modeltime', iii, self.modeltime[iii])
        print('WARM TEMPERATURES HAVE BEEN SET TO 273.15; MODEL RUN IS CONTINUING')
        self.Tz[self.Tz>=273.15]=273.15

    ### sort out how much liquid has been frozen and update density and LWC.
    # if np.any(self.rho>917.0):
    #     print('high rho before reallocate',iii)
    lwc_old = self.LWC.copy()
    mass_old = self.mass.copy()
    rho_old = self.rho.copy()
    tot_mass_old = mass_old + lwc_old*1000

    self.LWC        = g_liq * vol_tot
    delta_mass_liq  = mass_liq - (self.LWC * RHO_W_KGM)
    if np.any(delta_mass_liq<0):
        print(self.modeltime[iii],'Fixing negative values of delta_mass_liq, min value:',min(delta_mass_liq))
        delta_mass_liq  = np.maximum(delta_mass_liq,0) # fix for numerical instabilities with small time steps.
    self.mass       = self.mass + delta_mass_liq
    self.rho        = self.mass/self.dz

    tot_mass_new = self.mass + self.LWC*1000

    return self.Tz, self.T10m, self.rho, self.mass, self.LWC

##############################
### end enthalpy diffusion ###
##############################

#### Radiation penetration (work in progress)
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