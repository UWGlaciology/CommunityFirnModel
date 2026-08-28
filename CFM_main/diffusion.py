#!/usr/bin/env python
'''
diffusion.py

Functions for handling firn temperature evolution via diffusion, including
standard heat diffusion (no liquid water) and refreezing/enthalpy diffusion
(when liquid water is present in the firn column). Calls solver functions
defined in solver.py.

Draws from:
    - Numerical Heat Transfer and Fluid Flow (Patankar, 1980) for the
      finite volume discretization.
    - Voller, Swaminathan, and Thomas (1990) for the enthalpy formulation
      of phase change.

Grid convention:
    - self.z: layer edges [m] as tracked by the CFM elsewhere in the model.
    - z_P: layer (finite volume) centers, computed here from self.z.
    - z_edges: layer edges, padded with one dummy point past self.z[-1] so
      that the last real layer has a well-defined volume.

Main entry points (called each time step from the CFM's main loop):
    - heatDiff(self, iii): standard heat diffusion, no liquid water present.
      Calls transient_solve_TR in solver.py.
    - refreezeDiff(self, iii, _solver=...): heat diffusion with refreezing,
      used when liquid water (LWC) is present in the firn. Dispatches to one
      of several solver functions in solver.py (see solver_map inside the
      function) for comparison/testing purposes.

Other functions:
    - firnConductivity(self, iii, K_ice): computes firn thermal conductivity
      using one of several parameterizations, set via self.c['conductivity'].
      References for each parameterization are noted inline.
    - total_enthalpy(), total_h2o_mass(): diagnostic helpers for checking
      energy/mass conservation.

Isotope diffusion has its own dedicated class/module.

Note: refreezeDiff's default solver argument must match a key in solver_map;
see solver.py's module header for a description of the available solvers.
'''

from solver import (
    transient_solve_TR,
    transient_solve_enthalpy,
    transient_solve_ahc,
    transient_solve_decp,
    transient_solve_ncz
)
from constants import *
import numpy as np

def firnConductivity(self, iii, K_ice):
    '''
    Compute firn thermal conductivity using a selectable parameterization.

    The parameterization is chosen via self.c['conductivity']; see the
    elif chain in this function for the full list of supported string
    values. If self.c['conductivity'] does not match any known option,
    falls back to the Calonne et al. (2019) parameterization and prints
    a warning (once, at iii==0).

    :param iii: current model time step index (used only to control
        one-time print statements on the first step)
    :param K_ice: thermal conductivity of ice [W/m/K] at each layer,
        typically temperature-dependent (see caller, e.g. Cuffey and
        Paterson eq. 9.2 / Yen 1981)

    :return: K_firn, firn thermal conductivity [W/m/K] at each layer,
        same shape as self.rho

    Supported self.c['conductivity'] options (see inline comments for
    full references)::

        'Calonne2019' (default/fallback), 'Schwander', 'Yen_fixed',
        'Yen_var', 'Anderson', 'Yen_b', 'Sturm', 'VanDusen',
        'Schwerdtfeger', 'Riche', 'Jiawen', 'Calonne2011', 'mix'
        ('mix' blends Sturm and Anderson by depth; see inline code for
        the exact depth thresholds used)

    NOTE: the 'Calonne2019' branch and the fallback (else) branch contain
        duplicated code (identical formula, copy-pasted). Consider
        refactoring the fallback to simply call/reuse the 'Calonne2019'
        branch's logic to avoid the two implementations silently
        diverging if one is edited without the other.
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

def heatDiff(self, iii):
    '''
    Standard heat diffusion step (no liquid water / no phase change).

    Sets up the finite volume grid (layer centers and edges, with a dummy
    point appended past the last layer so it has a well-defined volume),
    computes temperature-dependent ice conductivity and firn conductivity
    (via firnConductivity), specific heat, and volumetric heat capacity,
    then calls transient_solve_TR in solver.py to advance temperature by
    one time step. Also enforces an upper bound of 273.15 K (melting point)
    on the result, printing a warning if this bound was exceeded before
    clamping.

    Use refreezeDiff instead of this function when liquid water (LWC) is
    present in the firn column, since that requires latent heat handling
    that this function does not provide.

    An older/alternate version of this logic may exist as heatDiffOLD()
    (see comment in original source); not covered by this docstring.

    :param iii: current model time step index (used to index self.dt,
        and to control one-time print statements via firnConductivity)

    :return: self.Tz (updated temperature profile [K]), self.T10m
        (temperature at 10 m depth [K], or NaN if the firn column is
        shallower than 10 m)

    Reference: Patankar (1980); Cuffey and Paterson (2010), eqs. 9.1-9.2
        for ice conductivity/specific heat.

    NOTE: as of this writing, this function contains a couple of no-op
        self-assignments (K_firn = K_firn; phi_0 = phi_0) left over from
        editing; safe to remove.
    '''

    nz_P            = len(self.z) #- 1
    nz_fv           = nz_P - 2 # this does not get used

    ### Add a dummy point at the end so that the fields at z[-1] have a volume associated with them
    z_dummy = np.zeros(len(self.z)+1)
    z_dummy[:-1] = self.z
    z_dummy[-1] = self.z[-1]+np.diff(self.z)[-1]
 
    ### z_P is layer centers
    z_P = (z_dummy[1:] + z_dummy[:-1])/2 # this assumes that self.z are the edges of the firn layers; this gets the centers of the layers for finite volume solver
    ### z_edges is layer edges
    z_edges = z_dummy

    phi_s           = self.Tz[0]
    phi_0           = self.Tz

    K_ice           = 9.828 * np.exp(-0.0057 * phi_0) # thermal conductivity, Cuffey and Paterson, eq. 9.2 (Yen 1981)
    K_firn = firnConductivity(self,iii,K_ice)

    c_firn          = 152.5 + 7.122 * phi_0 # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
    # c_firn        = CP_I # If you prefer a constant specific heat.

    if self.c['MELT']:
        try:
            if self.c['LWCheat']=='lowK':
                K_firn[self.LWC>0]=K_firn[self.LWC>0]/1.e4
        except:
            pass

    Gamma_P         = K_firn
    c_vol           = self.rho * c_firn

    self.Tz         = transient_solve_TR(z_edges, z_P, self.dt[iii], Gamma_P, phi_0, nz_P, phi_s, c_vol)

    try:
        self.T10m       = self.Tz[np.where(self.z>=10.0)[0][0]]
    except:
        self.T10m = np.nan

    if self.c['MELT']:
        if self.c['LWCheat']=='effectiveT':
            pass

        elif np.any(self.Tz>273.1500001):
            print(f'WARNING: TEMPERATURE EXCEEDS MELTING TEMPERATURE at {iii}')
            print('WARM TEMPERATURES HAVE BEEN SET TO 273.15; MODEL RUN IS CONTINUING')

        self.Tz[self.Tz>=273.15]=273.15

    return self.Tz, self.T10m

##########################
### end heatDiff ###
##########################

##########################
### Meltwater refreezing methods
##########################

# def total_enthalpy(T_C, th_solid, th_liquid, z_edges):
#     dZ = np.diff(z_edges)
#     return np.sum(enthalpy_of(T_C, th_solid, th_liquid) * dZ)   # [J/m2]

# def total_h2o_mass(th_solid, th_liquid, z_edges):
#     dZ = np.diff(z_edges)
#     return np.sum((th_solid + th_liquid) * dZ)                  # [kg/m2]

def refreezeDiff(self, iii, solver_name='transient_solve_enthalpy'):
    '''
    Heat diffusion step with refreezing (liquid water present in firn).

    Sets up the finite volume grid (same convention as heatDiff), converts
    temperature to deg C (fusion at T=0), computes effective thermal
    conductivity as a liquid/solid volume-fraction-weighted mix, then
    dispatches to one of several refreezing solver functions in solver.py
    (see solver_map below) to advance temperature, solid mass, and liquid
    mass by one time step. Which solver is used is controlled by the
    _solver argument, primarily to support side-by-side comparison of
    solver methods during development (see solver.py module header).

    After the solver call, converts results back to Kelvin/density/LWC,
    clamps any temperature exceeding melting point back to 273.15 K, and
    runs several diagnostic checks: liquid-mass-gain detection (flags
    layers where diffusion appears to have created liquid water), an
    overshoot-clamp report (only for solvers that return claw_mushy/
    claw_dry -- currently only transient_solve_enthalpy), and an overall
    mass conservation check (pre- vs. post-solve total H2O mass).

    :param iii: current model time step index (used to index self.dt,
        for print diagnostics, and passed through to firnConductivity)
    :param _solver: string key selecting which solver function to call;
        must match a key in solver_map (currently 'transient_solve_enthalpy',
        'transient_solve_ahc', or 'transient_solve_decp').
    :return: self.Tz (updated temperature [K]), self.T10m (temperature at
        10 m depth [K] or None), self.rho (updated density [kg/m3]),
        self.mass (updated solid mass [kg/m2] per layer), self.LWC
        (updated liquid water volume [m3] per layer), dml_sum (sum of any
        negative liquid-mass changes flagged as unexpected losses [kg/m2];
        0.0 if none detected)

    Reference: Voller, Swaminathan, and Thomas (1990) for enthalpy method;
        see individual solver docstrings in solver.py for method-specific
        references.

    NOTE: several intermediate quantities computed in this function
        (vol_total, mass_total, rho_total, and the commented-out c_vol
        block) are not currently used downstream; candidates for removal
        or for completing an intended (but currently incomplete) enthalpy
        conservation check (tot_heat_pre is computed but no tot_heat_post
        comparison exists yet).
    '''

    solver_map = {
        'transient_solve_enthalpy': transient_solve_enthalpy,
        'transient_solve_ahc': transient_solve_ahc,
        'transient_solve_decp': transient_solve_decp,
        'transient_solve_ncz': transient_solve_ncz,
    }
    
    solver = solver_map[solver_name]


    ### Grid:
    z_dummy = np.zeros(len(self.z)+1) # Add a dummy point at the end so that the fields at z[-1] have a volume associated with them 
    z_dummy[:-1] = self.z
    z_dummy[-1] = self.z[-1]+np.diff(self.z)[-1] 
    z_P = (z_dummy[1:] + z_dummy[:-1])/2 # this assumes that self.z are the edges of the firn layers; this gets the centers of the layers for finite volume solver
    z_edges = z_dummy
    ###

    T_s           = self.Tz[0] - T_MELT # Surface work in [C] so that reference Temperature is 0 for enthalpy
    TzC           = self.Tz - T_MELT
  
    vol_solid     = self.mass / RHO_I     # volume of the ice portion of each volume
    vol_liquid     = self.LWC
    # vol_total     = vol_solid + vol_liquid    # total volume of ice and liquid in each layer (porosity ignored)

    mass_liquid    = vol_liquid * RHO_W_KGM  # mass of liquid water
    mass_liquid_start = mass_liquid.copy()
    mass_solid    = self.mass
    # mass_total    = mass_solid + mass_liquid    

    th_liquid      = mass_liquid / self.dz      # th_wat = φ_wat · ρ_wat — liquid water mass per total volume
    th_liquid_old = th_liquid.copy()
    th_solid      = self.rho                   # th_solid = φ_ice · ρ_ice — ice mass per total volume

    # rho_total     = (self.mass + mass_liquid) / self.dz # 'total' density of volume (solid plus liquid)

    ### claude code calls this phi. It is volume fraction
    g_liquid    = th_liquid / RHO_W_KGM # liquid volume fraction (of the material portion, porosity ignored)
    g_solid     = th_solid / RHO_I     # solid/ice volume fraction 

    ### Specific Heats
    ### not sure c_vol gets used
    # c_firn          = 152.5 + 7.122 * self.Tz # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
    # c_firn  = CP_I # If you prefer a constant specific heat
    # c_solid = c_firn
    
    # c_liquid = 4219.9 # J/kg/K, taken from engineeringtoolbox.com. Ha!
    
    # # c_vol = g_solid * RHO_I * c_solid + g_liquid * RHO_W_KGM * c_liquid #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)
    # c_vol = (g_solid * c_solid + g_liquid * c_liquid) * rho_total #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)
    #######

    ### Conductivity
    K_solid   = 9.828 * np.exp(-0.0057 * self.Tz) # thermal conductivity of ice (W/m/K), Cuffey and Paterson, eq. 9.2 (Yen 1981)
    K_firn = firnConductivity(self,iii,K_solid) # thermal conductivity [W/m/K]
    
    K_water = 0.55575                         # thermal conductivity, water (W/m/K)
    # K_liquid = K_water * (th_wat/1000)**1.885 # assume that conductivity of water in porous material follows a similar relationship to ice.
    K_liquid = K_water
    K_eff = g_liquid * K_liquid + g_solid * K_firn # effective conductivity
    ###
    
    ### Total enthalpy/mass before solver (for testing conservation)
    # tot_heat_pre = np.sum(CP_I_kJ * self.mass * self.Tz + T_MELT * CP_W/1000 * vol_liquid * RHO_W_KGM + LF_I_kJ * vol_liquid * RHO_W_KGM)
    tot_mass_pre = np.sum(self.mass + vol_liquid*1000)
    lwc_old = vol_liquid.copy()

    ### call the solver
    solver_out = solver(z_edges, z_P, self.dt[iii], K_eff, TzC, th_liquid, th_solid, iii)
    ###
    self.total_count += solver_out['count']

    self.Tz         = solver_out['TzC_return'] + 273.15
    self.rho        = solver_out['th_solid']
    _mass_liquid    = solver_out['th_liquid'] * self.dz
    self.LWC        = _mass_liquid / RHO_W_KGM 
    self.mass       = self.rho * self.dz
    try:
        self.T10m       = self.Tz[np.where(self.z>=10.0)[0][0]]
    except:
        self.T10m = None

    gain = (self.LWC * RHO_W_KGM) - mass_liquid_start      # >0 = liquid created
    gain_layers = np.where(gain > 1e-6)[0]
    if gain_layers.size:
        print(f"[iii={iii}] liquid gained in {gain_layers.size} layers, "
            f"max {gain.max():.3e} kg/m2")
        had_liq = mass_liquid_start[gain_layers] > 0
        print(f"init mass: {mass_liquid_start[gain_layers].max():.3e}")
        print(f"   of those, initially mushy: {had_liq.sum()}, "
            f"initially dry: {(~had_liq).sum()}")
        print(f"   T at those layers (C): {solver_out['TzC_return'][gain_layers][:5]}")

    if 'claw_mushy' in solver_out.keys():
        CLAW_TOL = 5e-2   # kg/m2 per step; tune to your noise floor
        claw_total = solver_out['claw_mushy'] + solver_out['claw_dry']

        if claw_total > CLAW_TOL:
            total_lwc_mass = np.sum(lwc_old * RHO_W_KGM)   # kg/m2, pre-solve liquid mass
            rel_claw = claw_total / total_lwc_mass if total_lwc_mass > 0 else np.inf

            print(f"[iii={iii}] clamp clawed back liquid: "
                f"mushy={solver_out['claw_mushy']:.3e}, dry={solver_out['claw_dry']:.3e} kg/m2 "
                f"(TOTAL {claw_total:.3e})")
            print(f"    relative to total column LWC ({total_lwc_mass:.3e} kg/m2): "
                f"{rel_claw*100:.3f}%")
    
    delta_mass_liquid  = mass_liquid_start - (self.LWC * RHO_W_KGM)
    dml_sum = 0.0 

    if np.any(self.Tz>273.1500001):
        print('WARNING: TEMPERATURE EXCEEDS MELTING TEMPERATURE')
        print('Maximal temperature was:',np.max(self.Tz),' at layers:',np.where(self.Tz == np.max(self.Tz)))
        print('iii, modeltime', iii, self.modeltime[iii])
        print('WARM TEMPERATURES HAVE BEEN SET TO 273.15; MODEL RUN IS CONTINUING')
    self.Tz[self.Tz>=273.15]=273.15

    if np.any(delta_mass_liquid<0):
        if np.any(np.abs(delta_mass_liquid[delta_mass_liquid<0])>1e-5):
            print('------')
            print('If you are seeing this message there was a liquid mass gain in diffusion.') 
            print('Please email maxstev@umd.edu so I can fix it.')
        dml_sum = np.sum(delta_mass_liquid[delta_mass_liquid<0])

    tot_mass_post = np.sum(self.mass + _mass_liquid)

    if np.abs((tot_mass_post-tot_mass_pre)/tot_mass_pre)>1e-3: # flag if there is larger than 0.1% difference
        print(f'change in mass (enthalpy solver) at iteration {iii}!')
        print('pre:', tot_mass_pre)
        print('post:', tot_mass_post)

    return self.Tz, self.T10m, self.rho, self.mass, self.LWC, dml_sum


########################
### end refreezeDiff ###
########################

def enthalpyDiff_old(self, iii):
    '''
    Legacy enthalpy diffusion function (used until mid-July 2026).

    Superseded by refreezeDiff (which dispatches to transient_solve_enthalpy,
    transient_solve_ahc, or transient_solve_decp). Retained for
    testing/comparison purposes only; calls transient_solve_EN_old in
    solver.py. Not part of the current default CFM time-stepping pipeline.

    Method: Voller and Swaminathan (1991)/Voller, Swaminathan, and Thomas
    (1990) enthalpy formulation. LWC is tracked in volume [m^3].
    Thermal diffusivity: alpha = K_firn / (rho * c_firn).

    Grid convention: layers = volumes (finite volume centers/edges), same
    as heatDiff/refreezeDiff.

    :param iii: current model time step index (used to index self.dt, for
        firnConductivity's one-time print, and for diagnostic messages)

    :return: self.Tz (updated temperature [K]), self.T10m (temperature at
        10 m depth [K] or None), self.rho (updated density [kg/m3]),
        self.mass (updated solid mass [kg] per layer), self.LWC (updated
        liquid water volume [m3] per layer), dml_sum (sum of any negative
        liquid-mass changes flagged as unexpected losses [kg]; 0.0 if none
        detected)

    NOTE: sets nt=10 if any LWC>0 else nt=1, intending to control solver
        iteration count -- however, transient_solve_EN_old no longer uses
        this argument (iteration count is controlled by its own max_iter
        parameter instead). This nt logic is currently inert; see
        transient_solve_EN_old's docstring for details.

    NOTE: contains several commented-out alternate calculations (e.g., for
        c_vol, K_liq) preserved from earlier development; left as-is since
        this function is for legacy comparison testing rather than active
        development.
    '''

    Tstart          = self.Tz.copy()

    # T_old = self.Tz.copy() # initial temperature profile

    ### Add a dummy point at the end so that the fields at z[-1] have a volume associated with them 
    z_dummy = np.zeros(len(self.z)+1)
    z_dummy[:-1] = self.z
    z_dummy[-1] = self.z[-1]+np.diff(self.z)[-1]
 
    z_P = (z_dummy[1:] + z_dummy[:-1])/2 # this assumes that self.z are the edges of the firn layers; this gets the centers of the layers for finite volume solver
    z_edges = z_dummy

    phi_s           = self.Tz[0] - T_MELT # work in [C] so that reference Temperature is 0 for enthalpy
    phi_0           = self.Tz - T_MELT
  
    vol_ice     = self.mass / RHO_I     # volume of the ice portion of each volume
    vol_tot     = vol_ice + self.LWC    # total volume of ice and liquid in each volume
    mass_liq    = self.LWC * RHO_W_KGM  # mass of liquid water
    rho_liq_eff = mass_liq / self.dz      # effective density of the liquid portion
    # tot_rho     = (self.mass + mass_liq) / self.dz # 'total' density of volume (solid plus liquid)
    g_liq_1     = self.LWC / vol_tot     # liquid volume fraction (of the material portion, porosity ignored)
    g_ice_1     = vol_ice / vol_tot     # solid/ice volume fraction 

    K_water = 0.55575                         # thermal conductivity, water (W/m/K)
    K_ice   = 9.828 * np.exp(-0.0057 * self.Tz) # thermal conductivity, ice (W/m/K), Cuffey and Paterson, eq. 9.2 (Yen 1981)
    # K_mix = g_liq_1*K_liq + g_ice_1*K_ice

    ### Specific Heats
    # c_firn          = 152.5 + 7.122 * self.Tz # specific heat, Cuffey and Paterson, eq. 9.1 (page 400)
    c_firn  = CP_I # If you prefer a constant specific heat
    c_ice = c_firn
    c_liq = 4219.9 # J/kg/K, taken from engineeringtoolbox.com. Ha!
    # c_vol = g_ice_1 * RHO_I * c_ice + g_liq_1 * RHO_W_KGM * c_liq #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)
    c_vol = (g_ice_1 * c_ice + g_liq_1 * c_liq) * tot_rho #Voller eq. 10., the 'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)

    # tot_rho = (mass_solid + mass_liquid) / dz # 'total' density of volume (solid plus liquid)
    # vol_tot = vol_ice + vol_liquid
    # g_liq   = vol_liquid / vol_tot  # liquid volume fraction
    # g_ice   = vol_ice / vol_tot     # ice volume fraction 
    # c_liq = 4219.9 # J/kg/K
    # c_ice = 2097.0 # J/kg/K
    # c_vol   = (g_ice * c_ice + g_liq * c_liq) * tot_rho #'volume-averaged specific heat of mixture', or rho * cp. (so really heat capacity)

    ### Conductivity
    K_firn = firnConductivity(self,iii,K_ice) # thermal conductivity [W/m/K]

    K_liq = K_water * (rho_liq_eff/1000)**1.885 # I am assuming that conductivity of water in porous material follows a similar relationship to ice.
    K_eff = g_liq_1*K_liq + g_ice_1*K_firn # effective conductivity

    ### Total enthalpy/mass before solver (for testing conservation)
    tot_heat_pre = np.sum(CP_I_kJ*self.mass*self.Tz + T_MELT*CP_W/1000*self.LWC*RHO_W_KGM + LF_I_kJ*self.LWC*RHO_W_KGM)
    tot_mass_pre = np.sum(self.mass + self.LWC*1000)

    lwc_old = self.LWC.copy()

    phi_ret, g_liq, count, iterdiff,g_sol   = transient_solve_EN_old(z_edges, z_P, self.dt[iii], K_eff, phi_0, phi_s, self.LWC, self.mass, self.dz, iii)

    LWC_ret = g_liq * self.dz
    # self.LWC        = g_liq * vol_tot

    delta_mass_liq  = mass_liq - (LWC_ret * RHO_W_KGM)
    dml_sum = 0.0 

    self.LWC = LWC_ret.copy()
    self.Tz = phi_ret + 273.15
    try:
        self.T10m       = self.Tz[np.where(self.z>=10.0)[0][0]]
    except:
        self.T10m = None

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

    tot_mass_post = np.sum(self.mass + self.LWC*1000)

    if np.abs((tot_mass_post-tot_mass_pre)/tot_mass_pre)>1e-3: # flag if there is larger than 0.1% difference
        print(f'change in mass (enthalpy solver) at iteration {iii}!')
        print('pre:', tot_mass_pre)
        print('post:', tot_mass_post)

    return self.Tz, self.T10m, self.rho, self.mass, self.LWC, dml_sum

##############################
### end enthalpy diffusion ###
##############################

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