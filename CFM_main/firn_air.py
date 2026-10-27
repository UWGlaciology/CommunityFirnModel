#!/usr/bin/env python

'''
firn_air.py
===========
Firn-air (gas) transport for the CFM: diffusion, gravitational and thermal
fractionation, and advection of trace gases and isotopes through the open pore
space of the firn column.

Copyright © 2021 C. Max Stevens

Distributed under terms of the MIT license.

Overview
--------
One ``FirnAir`` instance is created per gas (see ``gaschoice`` in the firn-air
config file). ``FirnDensityNoSpin`` builds the instances at initialization and
calls ``FirnAir.firn_air_diffusion`` once per time step, after the densification,
temperature, and (if used) meltwater updates, so the gas sees the current
density/temperature/grid. Each call:

    1. recomputes porosity (total, closed, open), the close-off depth, and the
       lock-in depth from the current density profile;
    2. estimates the air pressure perturbation caused by the change in open
       pore volume since the previous step (used by the 'Darcy' advection);
    3. computes the diffusivity profile from the parameterization chosen with
       ``Diffu_param``, plus an eddy diffusivity for the convective zone (and,
       optionally, the lock-in zone);
    4. solves the 1-D transient advection-diffusion equation for the gas
       concentration with ``solver.transient_solve_TR`` (finite volume,
       Patankar 1980), which includes gravitational and thermal fractionation
       source terms;
    5. advances the (simple) gas-age tracker.

Solved equation (schematically, per unit open pore volume)::

    d(phi_op * C)/dt = d/dz[ (D + D_eddy) * phi_op * dC/dz ] - d/dz[ w * phi_op * C ] + S

where C is the gas concentration (a ratio to the surface value; the model is
linear in C), phi_op is the open porosity, D the molecular diffusivity, w the
air velocity relative to the firn, and S the gravitational + thermal source
term. The upper boundary is a fixed concentration (the atmospheric history
``Gs``); the lower boundary is zero-gradient.

Configuration
-------------
Firn air is switched on with ``"FirnAir": true`` in the main CFM .json file; the
file named by ``"AirConfigName"`` (default ``AirConfig.json``) holds the keys
read here (stored as ``self.cg``)::

    gaschoice            list of gas names; each must be implemented in ``gasses()``
    gravity              'on' / 'off': gravitational fractionation
    thermal              'on' / 'off': thermal fractionation (requires gravity 'on')
    runtype              'transient' or 'steady' (steady uses ``steady_T`` for
                         the close-off density and diffusivity temperature, and
                         ``steady_bdot``; the main CFM then forces Herron-Langway
                         physics and instantaneous accumulation)
    steady_T, steady_bdot   steady-state temperature [K] and accumulation [m ice eq. / yr]
    Diffu_param          diffusivity parameterization: 'Severinghaus', 'Schwander',
                         'Freitag', 'Adolph', 'Witrant', 'Battle', or 'Christo'
                         (see ``FirnAir.diffusivity``)
    lockin               if true, include eddy diffusivity in the lock-in zone
                         and shut off molecular diffusion below the lock-in depth
    ConvectiveZoneDepth  e-folding depth of the near-surface eddy diffusivity [m]
    advection_type       'zero', 'Darcy', or 'Christo' (see ``solver.w``)
    outputs              which firn-air fields the main model writes

Units
-----
Depths in m, density in kg m^-3, temperature in K, pressure in Pa, diffusivities
in m^2 s^-1, time steps in s (the model time vector is in years). Molar masses
are in kg mol^-1.

Notes
-----
Several constants are site-specific or empirical placeholders and are flagged
where they occur: the reference air pressure ``p_a`` (7.0e4 Pa), the fixed
gas-age value assigned at the lock-in depth (15 yr), and the mean close-off
density relation (Martinerie et al.).

References
----------

* Adolph, A. C., Albert, M. R. (2014). Gas diffusivity and permeability through the
  firn column at Summit, Greenland. The Cryosphere, 8, 319-328.
* Battle et al. (``Battle`` diffusivity; the code does not record the full citation).
* Blunier, T., Schwander, J. (2000). Gas enclosure in ice: age difference and
  fractionation (lock-in density, ``LIDRho = rho_co - 14``).
* Buizert, C. PhD thesis, University of Copenhagen (the ``Christo`` options;
  ``gasses`` relative diffusivities and the CO2 free-air diffusivity, App. A3).
* Cuffey, K. M., Paterson, W. S. B. (2010). The Physics of Glaciers (``Severinghaus``).
* Freitag, J. et al. (2002) (``Freitag`` diffusivity, Eq. 15, and the Darcy permeability).
* Goujon, C., Barnola, J.-M., Ritz, C. (2003). J. Geophys. Res., 108(B7), 2366
  (closed-porosity relation, Eqs. 9-10).
* Grachev, A. M., Severinghaus, J. P. (2003). Thermal diffusion factors (``omega``
  for d15N2 and d40Ar).
* Kawamura, K. et al. (2006). Surface eddy diffusivity of the convective zone.
* Martinerie, P. et al. Mean close-off density as a function of temperature.
* Patankar, S. V. (1980). Numerical Heat Transfer and Fluid Flow.
* Schwander, J. (1988). Gas diffusion in firn (``Schwander`` diffusivity).
* Witrant, E. et al. (2012). Atmos. Chem. Phys., 12, 11465-11483 (``Witrant``).

'''


import numpy as np 
from solver import solver
from solver import transient_solve_TR
# from Gasses import gasses 
# from Sites import sites 
from reader import read_input
import json
import scipy.interpolate as interpolate
from constants import *
import os

class FirnAir:
    '''
    Class for modeling firn-air transport of a single gas or isotope ratio.

    One instance is created per gas in ``gaschoice``; ``FirnDensityNoSpin`` keeps
    them in the dictionary ``self.FA`` keyed by gas name and calls
    ``firn_air_diffusion`` each time step. Gas concentrations are tracked as a
    ratio to the surface value, so the surface boundary condition is the
    (normalized) atmospheric history ``Gs`` -- ones for the isotope ratios
    'd15N2' and 'd40Ar', which evolve only through fractionation.

    Attributes
    ----------
    cg : dict
        firn-air configuration (contents of the file named by ``AirConfigName``)
    Gs : array of floats
        surface gas concentration at each model time step (ratio to reference)
    Gz : array of floats
        gas concentration profile (same ratio convention as ``Gs``)
    Tz : array of floats
        firn temperature profile (unit: K)
    z : array of floats
        depth of each model node (unit: m)
    rho : array of floats
        firn density profile (unit: kg m^-3)
    bdot : array of floats
        accumulation rate time series (unit: m ice eq. per year)
    gaschoice : str
        name of the gas (must be implemented in ``gasses``)
    gam_x : float
        free-air diffusivity of this gas relative to CO2 (unitless)
    M : float
        molar mass of the gas (unit: kg mol^-1)
    deltaM : float
        molar mass of the gas minus that of air (unit: kg mol^-1); drives
        gravitational fractionation
    d_0 : float
        free-air diffusivity of CO2 at the surface temperature and ``p_a``
        (unit: m^2 s^-1)
    omega : float
        thermal diffusion sensitivity of the gas (unit: K^-1); drives thermal
        fractionation, zero for gases without a thermal term
    czd : float
        e-folding depth of the convective-zone eddy diffusivity (unit: m)
    p_a : float
        reference (surface) air pressure (unit: Pa); fixed at 7.0e4 Pa
    air_pressure : array of floats
        air-pressure perturbation relative to hydrostatic, from the change in
        open pore volume over the last step (unit: Pa)
    air_pressure_base : array of floats
        hydrostatic air pressure at initialization (unit: Pa)
    rho_co, por_co : float
        density (kg m^-3) and porosity at close-off
    bcoRho : float
        mean close-off density, from the Martinerie relation (unit: kg m^-3)
    LIDRho : float
        density at the lock-in depth (unit: kg m^-3)
    por_tot, por_cl, por_op : arrays of floats
        total, closed, and open porosity profiles (unitless)
    air_volume, air_volume_old : arrays of floats
        open pore volume of each layer at the current and previous step (unit: m)
    z_co, LID : float
        close-off depth and lock-in depth (unit: m); set each step
    diffu, d_eddy : arrays of floats
        molecular and eddy diffusivity profiles (unit: m^2 s^-1)
    gas_age : array of floats
        gas-age tracker (unit: years). Simplified: ages accumulate by one time
        step per step and are fixed at 15 yr at and below the lock-in depth
        (a placeholder, not a calculated gas age).

    Returns
    -------
    Gz : array of floats
        (from ``firn_air_diffusion``) updated gas concentration profile
    '''

    def __init__(self,air_config,Gs,z,modeltime,Tz, rho, dz, gaschoice, bdot):
        '''
        Sets up the porosity, pressure, and concentration fields for one gas.
        The initial concentration is the first value of the surface history
        everywhere in the column.
        :param air_config: firn-air configuration dictionary (see module docstring for keys)
        :param Gs: time series of surface gas concentration, one value per model time step
        :param z: depth of each model node (unit: m)
        :param modeltime: model time vector (unit: years); not used by the class
        :param Tz: initial firn temperature profile (unit: K)
        :param rho: initial firn density profile (unit: kg m^-3)
        :param dz: thickness of each model layer (unit: m)
        :param gaschoice: name of the gas; must be implemented in ``gasses``
        :param bdot: accumulation rate time series (unit: m ice eq. per year)
        '''

        self.cg = air_config
        self.Gs = Gs

        self.Gz             = self.Gs[0]*np.ones_like(z)
        self.Tz             = Tz
        self.z              = z
        self.rho            = rho
        self.bdot           = bdot
        self.gaschoice      = gaschoice

        # self.gam_x, self.M, self.deltaM, self.d_0, self.omega = gasses(self.cg['gaschoice'], Tz[0], P_0, M_AIR)
        self.gam_x, self.M, self.deltaM, self.d_0, self.omega = gasses(gaschoice, Tz[0], P_0, M_AIR)

        self.czd            = self.cg['ConvectiveZoneDepth']
        self.p_a            = 7.0e4
        # self.air_pressure     = self.p_a * np.ones_like(z)
        self.air_pressure   = self.p_a * np.exp(M_AIR*GRAVITY*z/(R*self.Tz))
        self.air_pressure_base = np.copy(self.air_pressure)
        # self.air_pressure_0 = np.copy(self.air_pressure)

        self.rho_co, self.por_co, self.por_tot, self.por_cl, self.por_op, self.bcoRho, self.LIDRho = self.porosity()
        
        self.air_volume     = self.por_op * dz

        self.gas_age = np.zeros_like(self.z)
    ### end __init__ ###
    ####################

    def diffusivity(self):
        '''
        Calculates the molecular and eddy diffusivity profiles for the current step.

        Molecular diffusivity of the gas in firn comes from the parameterization
        chosen with ``Diffu_param``. Each is a function of the free-air
        diffusivity ``d_0`` and ``gam_x`` (gas relative to CO2), temperature,
        and porosity (and, for 'Witrant', the accumulation rate)::

            'Severinghaus'  Cuffey and Paterson relation, depends on density, T, p_a
            'Schwander'     Schwander (1988), Eq. 2, depends on total porosity
            'Freitag'       Freitag (2002), Eq. 15, depends on open porosity
            'Adolph'        Adolph and Albert (2014), via a permeability-porosity fit
            'Witrant'       Witrant et al. (2012); transitions smoothly to low
                            diffusivity past a threshold that depends on bdot
            'Battle'        Battle et al., depends on open porosity
            'Christo'       profile read from a data file (c_diffu_NEEM.txt)
                            at a hardcoded absolute path -- not usable
                            elsewhere without editing this file

        An eddy diffusivity is then added: exponentially decaying from the
        surface with e-folding depth ``ConvectiveZoneDepth`` (surface value
        2.4e-5 m^2/s, Kawamura et al. 2006). If ``lockin`` is true, the eddy
        diffusivity between the lock-in depth and the close-off depth is set
        to the molecular value there, and the molecular diffusivity below the
        lock-in depth is set to ~zero (diffusion then proceeds through the
        eddy term only).

        Finally, diffusivity is set to 1e-40 (effectively zero) wherever
        density exceeds the mean close-off density, and the first three layers
        at that floor are forced to it: a single layer of low diffusivity is
        numerically unstable.

        Uses ``self.z_co``, ``self.LID``, and (for 'Witrant') ``self.bdot_t``,
        which ``firn_air_diffusion`` sets immediately beforehand.

        :return diffu: molecular diffusivity profile (unit: m^2/s, type: array)
        :return d_eddy: eddy diffusivity profile (unit: m^2/s, type: array)
        '''

        ## Constants
        d_eddy_sc   = self.d_0 #Eddy diffusivity in the convective zone
        d_ind       = np.min(np.where(self.z>self.z_co)) #indices of nodes past close-off depth
        if self.cg['runtype']=='steady':
            Tz_d = self.cg['steady_T'] * np.ones_like(self.Tz)
        else:
            Tz_d = self.Tz
        
        ##### Parameterizations for diffusivity #####
        if self.cg['Diffu_param'] == "Severinghaus": # Use Severinghaus relationship from Cuffey and Paterson
            diffu_full      = self.gam_x * self.d_0 * ((P_0/self.p_a) * (Tz_d/273.15)**1.85 * (2.00 * (1-(self.rho/RHO_I))-0.167))  

        elif self.cg['Diffu_param'] == "Schwander": # Use Schwander 1988, Eq. 2 Diffusivity (does not work very well) use 4e2 for d_0
            k_sch           = P_0 / self.p_a * (Tz_d / 253.16)**1.85 # Constant given in Schwander
            diffu_full      = self.gam_x * k_sch * (23.7 * self.por_tot - 2.84) / (1000**2) # 1/1000**2 is unit conversion.                 
        
        elif self.cg['Diffu_param'] == "Freitag": # Use Freitag, 2002, Eq 15 Diffusivity use 9e2 for d_0
            diffu_full      = 1.0 * self.gam_x * self.d_0 * self.por_op ** 2.1  

        elif self.cg['Diffu_param'] == "Adolph": # Use Adolph and Albert, 2014
            perm            = 10**-7.29 * self.por_op**3.71 
            diffu_full      = self.d_0 * self.gam_x * (perm/(10**-6.973))**(1/2.79)         
        
        elif self.cg['Diffu_param']=="Witrant":    ### Use Witrant, 2012
            diffu_full      = self.d_0 * (2.5 * self.por_op - 0.31) * (Tz_d / 273.15)**(1.8) * P_0 / self.p_a
            if self.bdot_t<0.1:
                # d_thresh = 3.16e-8
                d_thresh = 1.
            else:
                # d_thresh = 3.16e-6
                d_thresh = 100.

            thresh_ind = np.where(diffu_full*S_PER_YEAR<=d_thresh)[0]
            z_thresh = self.z[thresh_ind[0]]

            diffu_t = (d_thresh - (1e-2)) / (1 + np.exp(50*(self.z - (self.z_co + z_thresh)/2)/self.z_co)) + (1e-2)

            diffu_full[thresh_ind] = diffu_t[thresh_ind]/S_PER_YEAR

            diffu_full = self.gam_x * diffu_full
        
        elif self.cg['Diffu_param']=="Battle":
            diffu_full      = (self.gam_x * 1.16 * (23.7 * self.por_op - 2.84)) / (24*3600)
        
        elif self.cg['Diffu_param']=="Christo":
            pp              = "/Users/maxstev/Documents/Grad_School/Research/FIRN/CFM/CommunityFirnModel/gasmodel/DataImport"
            diffu_data      = np.loadtxt(os.path.join(pp,'c_diffu_NEEM.txt'))
            h               = diffu_data[:,0]
            diffu_full_data = self.gam_x*self.d_0*diffu_data[:,1]

            diffu_full      = np.interp(self.z,h,diffu_full_data)
        
        
                        
        ###### Add in high diffusivity in convective zone and low diffusivity below LID     
        ### Convective zone###
        d_eddy          = np.zeros(np.size(diffu_full))
        ind             = np.nonzero(self.z<self.czd)
        d_eddy_surf     = 2.426405E-5 #Kawamura, 2006
        # H_scale       = czd
        d_eddy_up       = d_eddy_surf * np.exp(-1*self.z/self.czd)
        
        #### Lock-in zone physics ###
        if self.cg['lockin']:
            LI_ind = np.where(diffu_full<=1e-9)[0][0]

            ind             = np.flatnonzero(self.z>self.LID)
            ind2            = np.flatnonzero(self.z<self.z_co)
            ind3            = np.intersect1d(ind,ind2)
            ind4 = np.where(self.z<=self.LID)[0]
            # ind5 = np.where(self.z>=self.LID)[0][0]

            # diffu_full[ind4] = diffu_full[ind4]-diffu_full[ind4[-1]] #re-scale diffusivity so it becomes zero at LID

            d_eddy[ind3]    = diffu_full[ind3]      # set eddy diffusivity in LID equal to diffusivity at LID
            diffu_full[ind] = 1e-20                 # set molecular diffusivity equal to zero for "main" diffusivity after LID - eddy diffusivity term drives diffusion below
            d_eddy          = d_eddy + d_eddy_up #make eddy diffusivity vector have convective and lock-in zone values
            
        else:
            # diffu_full[np.flatnonzero(self.z>self.z_co)] = 1.0e-40
            d_eddy          = d_eddy_up #eddy diffusivity includes only convective zone
         
        # ij = np.where(self.z>self.LID)[0][0]
        # diffu_full 

        # diffu_full[self.z>self.z_co]=1e-40
        diffu_full[self.rho>self.bcoRho]=1e-40

        diffu_full[diffu_full<=0] = 1.e-40

        i4=np.where(diffu_full==1.e-40)[0][0]
        diffu_full[i4:i4+3]=1e-40 # a single layer of low diffusivity is numerically unstable, so ensure that low-diffusivity is at least 3 layers thick.

        diffu=diffu_full
        
        return diffu , d_eddy #, dd
    ### end diffusivity ###
    #######################

    def porosity(self): #,rho,T
        '''
        Calculates the total, closed, and open porosity, and the close-off and
        lock-in densities, from the current density profile.

        The mean close-off density ``bcoRho`` follows the Martinerie
        relation, 1/rho_co = 1/RHO_I + 6.95e-7*T - 4.3e-5, evaluated at the
        temperature of the first node deeper than 20 m (or at ``steady_T`` for
        ``runtype == 'steady'``). Closed porosity follows Goujon et al. (2003),
        Eqs. 9-10 (alpha = 0.37), and is capped at the total porosity. Open
        porosity is total minus closed, and floored at 1e-25 to avoid
        division by zero in the solver. The local close-off density ``rho_co``
        is the density of the first node with (essentially) zero open porosity;
        the lock-in density is ``rho_co - 14`` kg/m^3 (Blunier and Schwander,
        2000).

        Raises an IndexError if the model domain does not extend below 20 m or
        past close-off (the lookups assume both exist).
        :return rho_co: density at close-off (unit: kg/m^3, type: float)
        :return por_co: mean close-off porosity (type: float)
        :return por_tot: total porosity profile (type: array)
        :return por_cl: closed porosity profile (type: array)
        :return por_op: open porosity profile (type: array)
        :return bcoRho: mean close-off density (unit: kg/m^3, type: float)
        :return LIDRho: density at the lock-in depth (unit: kg/m^3, type: float)
        '''
        
        indT=np.where(self.z>20)[0][0]
        if self.cg['runtype']=='steady':
            self.bcoRho         = 1/( 1/(RHO_I) + self.cg['steady_T'] * 6.95E-7 - 4.3e-5) # Martinerie mean close-off density
            # self.bcoRho = 825.3 
        else:
            self.bcoRho         = 1/( 1/(RHO_I) + self.Tz[indT] * 6.95E-7 - 4.3e-5) # Martinerie mean close-off density 
            # self.bcoRho = 825.3
        
        ### Porosity, from Goujon et al., 2003, equations 9 and 10
        self.por_tot        = 1-self.rho/RHO_I # Total porosity
        self.por_co         = 1 - self.bcoRho/RHO_I # mean close-off porosity
        alpha               = 0.37 # constant determined in Goujon
        self.por_cl         = np.zeros_like(self.por_tot)

        self.por_cl[self.por_tot>0]     = alpha*self.por_tot[self.por_tot>0]*(self.por_tot[self.por_tot>0]/self.por_co)**(-7.6)

        ind                 = self.por_cl>self.por_tot # where (if anywhere) closed porosity exceeds total porosity
        self.por_cl[ind]    = self.por_tot[ind] # fix por_cl in that case so it cannot exceed
        self.por_op         = self.por_tot - self.por_cl # Open Porosity

        co_ind              = np.where(self.por_op<=1e-10)[0][0]
        self.rho_co         = self.rho[co_ind]
        self.LIDRho         = self.rho_co - 14 #LID depth (Blunier and Schwander, 2000)

        op_i = np.where(self.por_op<=0)[0][0]
        co_i = np.where(self.rho>=self.rho_co)[0][0]
        ct_i = np.where(self.por_cl>=self.por_tot)[0][0]
        li_i = np.where(self.rho>=self.LIDRho)[0][0]

        # print('first full close', self.z[ct_i])
        # print('open',self.z[op_i])
        # print('open',self.rho[op_i])
        # print('bco',self.z[co_i])
        # print('bco',self.rho[co_i])
        # print('porosity',self.por_tot[co_i])
        # print('LID z',self.z[li_i])
        # print('LID rho',self.rho[li_i])
        # print('liz width',self.z[co_i]-self.z[li_i])

        self.por_op[self.por_op<=0] = 1.0e-25
        
        return self.rho_co, self.por_co, self.por_tot, self.por_cl, self.por_op, self.bcoRho, self.LIDRho
    ### end porosity ###
    ####################

    def firn_air_diffusion(self,AirParams,iii):

        '''
        Advances the gas concentration profile by one model time step.

        In order: copies the per-step model state into the instance, updates
        porosity and the close-off (``z_co``) and lock-in (``LID``) depths,
        estimates the air-pressure perturbation from the change in open pore
        volume (``air_volume_old / air_volume``), computes diffusivity, and
        solves the advection-diffusion equation with
        ``solver.transient_solve_TR`` (gravitational and thermal fractionation
        per ``gravity`` and ``thermal``; advection per ``advection_type``).
        The upper boundary is the surface concentration ``phi_s = Gz[0]`` and
        the lower boundary is zero-gradient. After the solve the profile is
        shifted down one node and the new surface value ``Gs[iii]`` is
        inserted at the top (matching the model's layer-advection step).

        The gas-age tracker ages every non-zero entry by one time step, shifts
        the profile with the firn, and sets the lock-in node and any
        unassigned node below the lock-in depth to 15 yr (placeholder).

        :param AirParams: per-step model state; each key is set as an instance
            attribute. Keys: 'Tz' (K), 'rho' (kg/m^3), 'dt' (s), 'z' (m),
            'rhos0' (kg/m^3; unused here), 'dz_old' and 'dz' (m; thickness
            before/after this step's densification), 'rho_old' (kg/m^3), and
            'w_firn' (firn downward velocity, m/s)
        :param iii: current time step index (indexes ``Gs`` and ``bdot``)
        :return Gz: updated gas concentration profile (type: array)
        :return diffu: molecular diffusivity profile (unit: m^2/s, type: array)
        :return w_p: air velocity from advection, at node centers (unit: m/s, type: array)
        :return gas_age: simplified gas-age tracker (unit: years, type: array)
        '''


        for k,v in list(AirParams.items()):
            setattr(self,k,v)

        nz_P        = len(self.z)

        z_edges1 = self.z[0:-1] + np.diff(self.z) / 2
        z_edges = np.concatenate(([self.z[0]], z_edges1, [self.z[-1]]))
        z_P_vec     = self.z
        # if (iii>100 and iii<110):
        #     print('z_P',z_P_vec[0:5])
        #     print('len p', len(z_P_vec))
        #     print('edges',z_edges[0:5])
        #     print('len edges',len(z_edges))
        # input('enter')
        # phi_s     = self.Ts[iii]
        phi_s       = self.Gz[0]
        phi_0       = self.Gz

        # K_ice     = 9.828 * np.exp(-0.0057 * phi_0)
        # K_firn    = K_ice * (self.rho / 1000) ** (2 - 0.5 * (self.rho / 1000))

        self.rho_co, self.por_co, self.por_tot, self.por_cl, self.por_op, self.bcoRho, self.LIDRho = self.porosity()
        
        # self.air_pressure_old     = np.copy(self.air_pressure)
        porosity_old            = (RHO_I-self.rho_old)/RHO_I

        por_co1              = 1 - self.rho_co/RHO_I # Porosity at close-off
        alpha               = 0.37 # constant determined in Goujon
        por_cl1              = np.zeros_like(porosity_old)
        por_cl1[porosity_old>0] = alpha*porosity_old[porosity_old>0]*(porosity_old[porosity_old>0]/por_co1)**(-7.6)
        ind                 = por_cl1>porosity_old
        por_cl1[ind]         = porosity_old[ind]
        por_op_old          = porosity_old - por_cl1 # Open Porosity

        self.air_volume_old     = por_op_old * self.dz_old
        # self.air_volume_old   = np.copy(self.air_volume)

        self.air_volume         = self.por_op * self.dz
        volfrac = self.air_volume_old / self.air_volume
        # volfrac = np.concatenate(([volfrac[0]],volfrac))
        self.air_pressure       = (self.p_a*np.exp(M_AIR*GRAVITY*self.z/(R*self.Tz))) * volfrac - (self.p_a*np.exp(M_AIR*GRAVITY*self.z/(R*self.Tz))) # assume air pressure is atmos in entire column


        self.pressure_grad      = np.gradient(self.air_pressure,self.z) 
        self.z_co               = min(self.z[self.rho>=(self.bcoRho)]) #close-off depth; bcoRho is mean close off density
        # self.z_co               = min(self.z[self.rho>=(self.rho_co)]) #close-off depth; bcoRho is close off density
        # self.z_co               = self.z[np.where(self.rho>=self.rho_co)[0][0]] #close-off depth; bcoRho is close off density
        self.LID                = min(self.z[self.rho>=(self.LIDRho)]) #lock in depth; LIDRho is lock-in density
        # print(self.LID)
        self.bdot_t = self.bdot[iii]

        self.diffu, self.d_eddy = self.diffusivity()
        # ij = np.where(self.diffu<=1e-20)[0][0]
        # print('Diffu zero',self.z[ij])

        airdict = {
            'd_eddy':           self.d_eddy,
            'por_op':           self.por_op,
            'Tz':               self.Tz,
            'deltaM':           self.deltaM,
            'omega':            self.omega,
            'dz':               self.dz,
            'rho':              self.rho,
            'gravity':          self.cg['gravity'],
            'thermal':          self.cg['thermal'],
            'air_pressure':     self.air_pressure,
            'pressure_grad':    self.pressure_grad,
            'z':                self.z, 
            'dt':               self.dt,
            'z_co':             self.z_co,
            'p_a':              self.p_a,
            'por_tot':          self.por_tot,
            'por_cl':           self.por_cl,
            'w_firn':           self.w_firn,
            'advection_type':   self.cg['advection_type'],
            'iii': iii,
            'gaschoice': self.gaschoice
            }

        msk = np.where(self.z>self.z_co)[0][0]
        c_vol = np.ones_like(self.rho) # Just a filler here to make the diffusion work.
        
        self.Gz, w_p = transient_solve_TR(z_edges, z_P_vec, self.dt, self.diffu, phi_0, nz_P, phi_s, c_vol, airdict)
        self.Gz = np.concatenate(([self.Gs[iii]], self.Gz[:-1]))

        ind_LID = np.where(self.z>=self.LID)[0]
        self.gas_age[self.gas_age>0] = self.gas_age[self.gas_age>0] + self.dt/S_PER_YEAR
        self.gas_age = np.concatenate(([0], self.gas_age[:-1]))
        self.gas_age[ind_LID[0]] = 15.0
        ii2 = np.where((self.z>self.LID) & (self.gas_age==0))[0]
        self.gas_age[ii2] = 15

        # print('diffu',self.diffu[0:10])

        return self.Gz, self.diffu, w_p, self.gas_age
    ### end firn_air_diffusion ###
    ##############################

GASES_AVAILABLE = ('CO2', 'CH4', 'd15N2', 'd40Ar', 'SF6', 'C14', 'C13', 'CFC11', 'CFC12',
                   'C13_CFC12', 'CC14', 'CFC113', 'CFC115', 'R134a', 'CH3CCl3',
                   'HCFC22', 'C13_CH4', 'FOG') # keep in sync with the branches in gasses()

def gasses(gaschoice, T, p_a, M_air):
    
    '''
    Returns the gas-specific constants used by the firn-air model.

    Implemented gases: CO2, CH4, d15N2, d40Ar, SF6, C14, C13, CFC11, CFC12,
    C13_CFC12, CC14, CFC113, CFC115, R134a, CH3CCl3, HCFC22, C13_CH4, FOG.
    ('FOG' is a synthetic test tracer: CO2-like with a 100 yr decay constant.)
    Free-air diffusivities relative to CO2 (``gam_x``) are the literature
    values compiled in the Buizert thesis; ``omega`` (thermal sensitivity) is
    nonzero only for d15N2 and d40Ar (Grachev and Severinghaus, 2003; the
    d40Ar value is scaled by its 4 amu mass difference).

    The free-air CO2 diffusivity returned as ``d_0`` is computed from
    temperature and pressure (Buizert thesis, Appendix A3):
    D = 5.75e-10 * T^1.81 * (101325 / p_a) [m^2/s].

    An unrecognized ``gaschoice`` raises a ValueError listing the available
    gases (``GASES_AVAILABLE``). Known limitation: the radioactive decay
    constant defined for C14 and FOG is not returned or used.
    :param gaschoice: name of the gas species
    :param T: temperature (unit: K)
    :param p_a: atmospheric pressure (unit: Pa)
    :param M_air: mean molar mass of air (unit: kg/mol)
    :return gam_x: free-air diffusivity relative to CO2 (type: float)
    :return M: molar mass of the gas (unit: kg/mol, type: float)
    :return deltaM: molar mass of the gas minus that of air (unit: kg/mol, type: float)
    :return d_0: free-air diffusivity of CO2 at T and p_a (unit: m^2/s, type: float)
    :return omega: thermal diffusion sensitivity (unit: K^-1, type: float)
    '''

    #d_0 = 5.e2 # Free air diffusivity, CO2, m**2/yr Schwander, 1988 reports 7.24 mm**2/s =379 m**2/yr
    d_0 = 1.6e-5 # m^2/s :wikipedia value. changed 9/27/13  Schwander, 1988 reports 7.24 mm**2/s = 7.24e-6 m**2/yr

    D_ref_CO2 = 5.75E-10*T**1.81*(101325/p_a) #Christo Thesis, appendix A3
    print('gas choice is ', gaschoice)
    
    if gaschoice == 'CO2':
        gam_x = 1. #free-air diffusivity relative to CO2. Unitless (gamma in Buizert thesis, page 13).
        M = 44.01e-3 # molecular mass, kg/mol
        decay = 0.
        omega = 0.0
        
    elif gaschoice == 'CH4':
        gam_x = 1.367
        M = 16.04e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'd15N2':
        
        gam_x = 1.275*0.9912227 # not sure of the origin here... Christo's model?
        #gam_x = 
        M = 1.E-3 + M_air
        decay = 0.
        omega = 0.0147/1000 # Grachev and Severinghaus, 2003 

    elif gaschoice == 'SF6':
        gam_x = 0.554
        M = 146.06e-3
        decay = 0.
        omega = 0.
        
    elif gaschoice == 'C14':
        gam_x = 0.991368
        M = 46.01e-3
        decay = 1./8267.
        omega = 0.
        
    elif gaschoice == 'C13':
        gam_x = 0.9955648
        M = 45.01e-3
        decay = 0.
        omega = 0.
        
    elif gaschoice == 'CFC11':
        gam_x = 0.525
        M = 137.37e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'CFC12':
        gam_x = 0.596
        M = 120.91e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'C13_CFC12':
        gam_x = 0.59552
        M = 121.91e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'CC14':
        gam_x = 0.470
        M = 153.82e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'CFC113':
        gam_x = 0.453
        M = 187.38e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'CFC115':
        gam_x = 0.532
        M = 154.47e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'R134a':
        gam_x = 0.630
        M = 102.03e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'CH3CCl3':
        gam_x = 0.485
        M = 133.40e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'HCFC22':
        gam_x = 0.710
        M = 86.47e-3
        decay = 0.
        omega = 0.

    elif gaschoice == 'C13_CH4':
        gam_x = 1.340806
        M = 17.04e-3
        decay = 0.
        omega = 0.
        
    elif gaschoice == 'd40Ar':
        gam_x = 1.21
        M = 4.e-3 + M_air
        decay = 0.
        omega = 0.00985/1000. * 4 # Grachev and Severinghaus, 2003

    elif gaschoice == 'FOG':
        gam_x = 1.0
        M = 44e-3
        decay = 1./100.
        omega = 0.

    else:
        raise ValueError(f"Unknown gas '{gaschoice}' in the firn-air config. "
                         f"Available gases: {', '.join(GASES_AVAILABLE)}")
            
    deltaM = (M - M_air) #delta molecular mass from CO2.
    #gam_x = D_gas #* D_ref_CO2
    d_0=D_ref_CO2
                   
    return gam_x, M, deltaM, d_0, omega
    ### end gasses ###
    ##################
