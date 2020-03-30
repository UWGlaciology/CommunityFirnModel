#!/usr/bin/env python
'''
code to handle diffusion (heat, enthalpy, isotope)
calls solver
Draws from Numerical Heat Transfer and Heat Flow (Patankar, 1980)

Isotope diffusion now has its own class.
'''

from solver import transient_solve_TR
# from solver import transient_solve_EN_old
# from solver import transient_solve_EN_new
from solver import transient_solve_EN
from constants import *
import numpy as np
from scipy import interpolate

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

    # z_edges_vec   = self.z[1:-2] + self.dz[2:-1] / 2
    # z_edges_vec   = np.concatenate(([self.z[0]], z_edges_vec, [self.z[-1]]))
    # z_P_vec       = self.z

    z_edges_vec1 = self.z[0:-1] + np.diff(self.z) / 2
    z_edges_vec = np.concatenate(([self.z[0]], z_edges_vec1, [self.z[-1]]))
    z_P_vec     = self.z
    
    phi_s           = self.Tz[0]
    phi_0           = self.Tz

    K_ice           = 9.828 * np.exp(-0.0057 * phi_0) # thermal conductivity, Cuffey and Paterson, eq. 9.2 (Yen 1981)
    c_firn          = 152.5 + 7.122 * phi_0 # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
    # c_firn        = CP_I # If you prefer a constant specific heat.

    ### Conductivity. Choose your favorite! ###
    # References are provided at the end of this script.

    # alpha_DAVIDCS_m2s = 18.3 / S_PER_YEAR # special for work Max was doing for David Clemens-Sewell

    # K_firn = alpha_DAVIDCS_m2s * self.rho * c_firn
    # K_firn[self.z>=20.0]    = K_ice[self.z>=20.0] * (self.rho[self.z>=20.0]/RHO_I) ** (2 - 0.5 * (self.rho[self.z>=20.0]/RHO_I))    # Schwander 1997, eq. A11

    # K_firn    = K_ice * (self.rho/RHO_I) ** (2 - 0.5 * (self.rho/RHO_I))    # Schwander 1997, eq. A11
    K_firn    = 2.22362 * (self.rho / 1000)**1.885                          # Yen 1981, eq 34 w/ fixed K_ice (original)
    # K_firn    = K_ice * (self.rho / 1000)**1.885                            # Yen 1981, modified for variable K_ice
    # K_firn    = 0.021 + 2.5 * (self.rho/1000.)**2                           # Anderson (1976)
    # K_firn    = 0.0688 * np.exp(0.0088*phi_0 + 4.6682*self.rho)             # Yen 1981, eq. 35.
    # K_firn    = 0.138 - 1.01*(self.rho/1000) + 3.233*(self.rho/1000)**2     # Sturm, 1997.; rho < 0.6
    # K_firn    = 2.1e-2 + 4.2e-4 * self.rho + 2.2e-9 * (self.rho)**3         # Van Dusen 1929 (via C&P)
    # K_firn    = (2 * K_ice * self.rho) / (3*RHO_I - self.rho)               # Schwerdtfeger (via C&P)
    # K_firn    = 3.e-6 * self.rho**2 - 1.06e-5 * self.rho + 0.024            # Riche and Schneebeli 2013 eq. 10
    # k_firn    = 0.0784 + 2.697 * (self.rho/1000.)**2                        # Jiawen 1991 eq. 3
    if self.c['MELT']:
        if self.c['LWCheat']=='lowK':
            K_firn[self.LWC>0]=K_firn[self.LWC>0]/1.e4

    Gamma_P         = K_firn # (new)
 
    # Gamma_P         = K_firn / (c_firn) # old)
    tot_rho         = self.rho
    c_vol = self.rho * c_firn
    # print('c_vol',c_vol)

    self.Tz         = transient_solve_TR(z_edges_vec, z_P_vec, nt, self.dt[iii], Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol)

    # fT10m           = interpolate.interp1d(self.z, self.Tz)                                 # temp at 10m depth
    # self.T10m       = fT10m(10)
    self.T10m       = self.Tz[np.where(self.z>=10.0)[0][0]] #Should be slightly faster than the interpolate.
    if self.c['MELT']:
        if self.c['LWCheat']=='effectiveT':
            pass
        elif np.any(self.Tz>273.15):
            print('WARNING: TEMPERATURE EXCEEDS MELTING TEMPERATURE')
            print('Maximal temperature was:',np.max(self.Tz),' at layers:',np.where(self.Tz == np.max(self.Tz)))
            print('WARM TEMPERATURES HAVE BEEN SET TO 273.15; MODEL RUN IS CONTINUING')
            self.Tz[self.Tz>=273.15]=273.15

    return self.Tz, self.T10m
### end heat diffusion
######################

def enthalpyDiff(self,iii):
    '''
    enthalpy diffusion function, new
    1/30/19 - method from Voller and Swaminathan
    LWC is in volume (m^3)
    thermal diffusivity: alpha = K_firn / (rho*c_firn)
    '''
    Tstart = self.Tz.copy()
    nz_P            = len(self.z)
    nz_fv           = nz_P - 2
    if np.any(self.LWC>0):
        nt              = 3
    else:
        nt = 1

    z_edges_vec1    = self.z[0:-1] + np.diff(self.z) / 2
    z_edges_vec     = np.concatenate(([self.z[0]], z_edges_vec1, [self.z[-1]]))
    z_P_vec         = self.z
    
    phi_s           = self.Tz[0]
    phi_0           = self.Tz

    H_ice = RHO_I * CP_I * (self.Tz - T_MELT)
    H_liq = RHO_W_KGM * LF_I
  
    vol_ice     = self.mass / RHO_I # volume of the ice portion of each volume
    vol_tot     = vol_ice + self.LWC # total volume of ice and liquid in each volume
    mass_liq    = self.LWC * RHO_W_KGM # mass of liquid water
    rho_liq_eff = mass_liq/self.dz # effective density of the liquid portion
    tot_rho     = (self.mass + mass_liq) / self.dz # 'total' density of volume (solid plus liquid)
    g_liq_1     = self.LWC / vol_tot # liquid volume fraction (of the material portion, porosity ignored)
    g_ice_1     = vol_ice / vol_tot # solid/ice volume fraction 
    g_liq       = self.LWC / self.dz # liquid volume fraction, total volume
    g_ice       = vol_ice / self.dz # solid/ice volume fraction, total volume

    # H_tot = H_ice * g_ice + H_liq * g_liq #total enthalpy (don't need?)

    K_water = 0.55575 # thermal conductivity of water (W/m/K)
    K_ice   = 9.828 * np.exp(-0.0057 * phi_0) # thermal conductivity (W/m/K), Cuffey and Paterson, eq. 9.2 (Yen 1981)
    # K_mix = g_liq_1*K_liq + g_ice_1*K_ice

    # c_firn          = 152.5 + 7.122 * phi_0 # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
    c_firn  = CP_I # If you prefer a constant specific heat
    c_ice = c_firn
    c_liq = 4219.9 # J/kg/K, taken from engineeringtoolbox.com. Ha!
    # c_vol = g_ice_1 * RHO_I * c_ice + g_liq_1 * RHO_W_KGM * c_liq #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)
    c_vol = (g_ice_1 * c_ice + g_liq_1 * c_liq) * tot_rho

    ### Conductivity. Choose your favorite! ###
    ### References are provided at the end of this script.

    # K_firn[self.z>=20.0]    = K_ice[self.z>=20.0] * (self.rho[self.z>=20.0]/RHO_I) ** (2 - 0.5 * (self.rho[self.z>=20.0]/RHO_I))    # Schwander 1997, eq. A11
    # K_firn    = K_ice * (self.rho/RHO_I) ** (2 - 0.5 * (self.rho/RHO_I))    # Schwander 1997, eq. A11
    # K_firn    = 2.22362 * (self.rho / 1000)**1.885                          # Yen 1981, eq 34 w/ fixed K_ice (original)
    # K_firn    = K_ice * (self.rho / 1000)**1.885                            # Yen 1981, modified for variable K_ice
    # K_firn    = 0.021 + 2.5 * (self.rho/1000.)**2                           # Anderson (1976)
    # K_firn    = 0.0688 * np.exp(0.0088*phi_0 + 4.6682*self.rho)             # Yen 1981, eq. 35.
    # K_firn    = 0.138 - 1.01*(self.rho/1000) + 3.233*(self.rho/1000)**2     # Sturm, 1997.; rho < 0.6
    # K_firn    = 2.1e-2 + 4.2e-4 * self.rho + 2.2e-9 * (self.rho)**3         # Van Dusen 1929 (via C&P)
    K_firn    = (2 * K_ice * self.rho) / (3*RHO_I - self.rho)               # Schwerdtfeger (via C&P)
    # K_firn    = 3.e-6 * self.rho**2 - 1.06e-5 * self.rho + 0.024            # Riche and Schneebeli 2013 eq. 10
    # k_firn    = 0.0784 + 2.697 * (self.rho/1000.)**2                        # Jiawen 1991 eq. 3
    
    K_liq = K_water * (rho_liq_eff/1000)**1.885 # I am assuming that conductivity of water in porous material follows a similar relationship to ice.

    K_eff = g_liq_1*K_liq + g_ice_1*K_firn # effective conductivity
    
    Gamma_P = K_eff

    deltaH = RHO_W_KGM * LF_I #Voller 1990, eq 11. (letting T_ref = T_melt) units: J/m^3

    self.Tz, g_liq   = transient_solve_EN(z_edges_vec, z_P_vec, nt, self.dt[iii], Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, c_vol, g_liq, deltaH)

    fT10m           = interpolate.interp1d(self.z, self.Tz)                                 # temp at 10m depth
    self.T10m       = fT10m(10)

    if np.any(self.Tz>273.16):
        print('WARNING: TEMPERATURE EXCEEDS MELTING TEMPERATURE')
        print('Maximal temperature was:',np.max(self.Tz),' at layers:',np.where(self.Tz == np.max(self.Tz)))
        print('iii, modeltime', iii, self.modeltime[iii])
        print('WARM TEMPERATURES HAVE BEEN SET TO 273.15; MODEL RUN IS CONTINUING')
        self.Tz[self.Tz>=273.15]=273.15

    ### sort out how much liquid has been frozen and update density and LWC.
    if np.any(self.rho>917.0):
        print('high rho before reallocate',iii)
    lwc_old = self.LWC.copy()
    mass_old = self.mass.copy()
    rho_old = self.rho.copy()

    self.LWC        = g_liq * vol_tot
    delta_mass_liq  = mass_liq - (self.LWC * RHO_W_KGM)
    if np.any(delta_mass_liq<0):
        print(self.modeltime[iii],'Fixing negative values of delta_mass_liq, min value:',min(delta_mass_liq))
        delta_mass_liq  = np.maximum(delta_mass_liq,0) # fix for numerical instabilities with small time steps.
    self.mass       = self.mass + delta_mass_liq
    self.rho        = self.mass/self.dz

    if np.any(self.rho>917.0):
    # if iii>12101:
        print('high rho after reallocate',iii)
        # ind9 = np.where(self.rho>917.0)[0]
        # ind9 = np.arange(0,6)
        # print('ind9',ind9)
        # print('Tstart',Tstart[0:5])
        # print('Tz',self.Tz[0:5])
        # print('dml',delta_mass_liq[ind9])
        # print('lwc',self.LWC[ind9])
        # print('dz',self.dz[ind9])
        # print('mass',self.mass[ind9])
        # print('rho',self.rho[ind9])
        # print('mass_liq',mass_liq[ind9])
        # print('lwc_old',lwc_old[ind9])
        # print('mass_old',mass_old[ind9])
        # print('rho_old', rho_old[ind9])
        # print('max dml', max(delta_mass_liq))
        # print('min dml', min(delta_mass_liq))
        # input('enter')



    return self.Tz, self.T10m, self.rho, self.mass, self.LWC
### end enthalpy diffusion


'''
### References for conductivity parameterizations ###
# Also can refer to Physics of Glaciers, chapter 9.2

Anderson EA (1976) A point energy and mass balance model of a snow cover. (doi:10.1016/S0074-6142(99)80039-4)
Brandt RE and Warren SG (1997) Temperature measurements and heat transfer in near-surface snow at the South Pole. J. Glaciol. 43(144), 339–351
Jiawen R, Dahe Q and Maohuan H (1991) THERMAL PROPERTIES AND TEMPERATURE DISTRIBUTION OF SNOW/FIRN ON THE LAW DOME ICE CAP, ANTARCTICA. Antarct. Res. 2(2), 38–46
Lüthi MP and Funk M (2001) Modelling heat flow in a cold, high-altitude glacier: Interpretation of measurements from Colle Gnifetti, Swiss Alps. J. Glaciol. 47(157), 314–324 (doi:10.3189/172756501781832223)
Riche F and Schneebeli M (2013) Thermal conductivity of snow measured by three independent methods and anisotropy considerations. Cryosphere 7(1), 217–227 (doi:10.5194/tc-7-217-2013)
Schwander J, Sowers T, Barnola J-M, Blunier T, Fuchs A and Malaizé B (1997) Age scale of the air in the summit ice: Implication for glacial-interglacial temperature change. J. Geophys. Res. Atmos. 102(D16), 19483–19493 (doi:10.1029/97JD01309)
Schwerdtfeger P (1963) Theoretical derivation of the thermal conductivity and diffusivity of snow. IAHS Publ 61, 75–81 http://iahs.info/uploads/dms/061007.pdf
Sturm M, Holmgren J, König M and Morris K (1997) The thermal conductivity of seasonal snow. J. Glaciol. 43(143), 26–41 (doi:10.1017/S0022143000002781)
Van Dusen MS (1929) Thermal conductivity of non-metallic solids. International critical tables of numerical data, physics, chemistry and technology. McGraw-Hill New York, 216–217
Yen Y-C (1981) Review of Thermal Properties of Snow, Ice, and Sea Ice. CRREL Rep. 81-10, 1–27 http://acwc.sdp.sirsi.net/client/search/asset/1005644

'''