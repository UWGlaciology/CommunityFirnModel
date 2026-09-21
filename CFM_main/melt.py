#!/usr/bin/env python
from constants import *
import numpy as np
import time
import sys

from diffusion import heatDiff

from darcy_funcs import hydrconducsat_Calonne
from darcy_funcs import vG_Yama_params
from darcy_funcs import phead_vG
from darcy_funcs import krel_vG
from darcy_funcs import thetae_update
from darcy_funcs import thetaeff_equaliser
from darcy_funcs import dfdg_derivative
from darcy_funcs import runoffZuoOerlemans
from darcy_funcs import runoffDarcy
from darcy_funcs import flux_bisection
from darcy_funcs import flux_newtonraphson

'''
Functions to handle meltwater percolation through the firn column.

Two percolation schemes are provided:

- bucket(): a discrete "bucket" scheme, treating each firn layer as a
  reservoir with a finite storage capacity (refreezing + irreducible
  retention). Computationally efficient; melt is applied as a single
  instantaneous pulse per model time step. See markdown reference doc
  for full details.

- darcyscheme(): a more physically detailed scheme solving for actual
  water flux between layers via an unsaturated Darcy-flow formulation
  (Hirashima et al. 2010), with van Genuchten hydraulic properties and
  its own adaptive internal sub-stepping. Requires the helper functions
  defined in darcy_funcs.py (imported below). See markdown reference doc
  for full details.
'''

#############
def _debug_report(debug, label, value):
    """
    Print a labeled diagnostic value if debug output is enabled.

    Used throughout `bucket()` to optionally trace intermediate LWC/runoff
    bookkeeping values (e.g. blocked LWC at various stages) without
    cluttering the main routine with conditional print statements.

    Parameters
    ----------
    debug : bool
        If False, this function does nothing.
    label : str
        Short description of the value being reported.
    value : float
        The diagnostic value to print.
    """
    if debug:
        print(f'  [debug] {label}: {value}')

#########################
### end _debug_report ###
#########################

def _available_pore_space(rho, dz, rho_i=RHO_I, cap_density=RHO_I):
    """
    Compute the pore space available for liquid water storage, capped so
    filling it with water cannot push a layer's density above `cap_density`.

    Implements the "saturated water content" of Eq. 9, Wever et al. (2014)
    (see also discussion in Yamaguchi et al. (2010)).

    Parameters
    ----------
    rho : ndarray
        layer density [kg m-3]. May exceed rho_i (e.g. a hypothetical
        "potential" density after refreezing); such layers get zero porosity.
    dz : ndarray
        layer thickness [m].
    rho_i : float, optional
        Ice density used in phi = (rho_i - rho) / rho_i. Default RHO_I.
        A value fractionally above 917 avoids porosity==0 exactly at
        rho==917 due to floating-point round-off.
    cap_density : float, optional
        Density [kg m-3] beyond which pore space is capped (default RHO_I).

    Returns
    -------
    phivol_av : ndarray
        Available pore space per layer [m].
    """
    phi = np.zeros_like(rho)
    below = rho < rho_i
    phi[below] = (rho_i - rho[below]) / rho_i
    phivol_av = phi * dz * (RHO_I / RHO_W_KGM)

    ilim = np.where(rho + phivol_av * RHO_W_KGM / dz > cap_density)[0]
    if ilim.size > 0:
        phivol_av[ilim] = np.maximum(
            dz[ilim] * (cap_density - rho[ilim]) / RHO_W_KGM, 0.
        )
    return phivol_av
#################################
### end _available_pore_space ###
#################################

def _irreducible_lwc(rho, phivol_av, RhoImp, rho_i=RHO_I,
                      coleou_lesaffre=True, irr_val=0.02):
    """
    Compute the irreducible liquid water content (LWC) each layer can hold
    against gravity drainage.

    Two formulations are supported:

    1. Coléou & Lesaffre (1998), via Langen et al. (2017) Eqs. 3-4:
       irreducible water mass fraction depends on density, then converted
       to a volume fraction of available pore space.
    2. Fixed fraction `irr_val` of available pore space (simpler, user-set).

    layers at or above `RhoImp` (the impermeability density threshold) are
    assigned zero irreducible LWC, since they are treated as solid ice
    lenses incapable of holding retained water.

    Parameters
    ----------
    rho : ndarray
        layer density [kg m-3].
    phivol_av : ndarray
        Available pore space per layer [m] (see `_available_pore_space`).
    RhoImp : float
        Density threshold [kg m-3] above which layers are impermeable.
    rho_i : float, optional
        Ice density [kg m-3]. Default RHO_I.
    coleou_lesaffre : bool, optional
        If True, use the Coléou & Lesaffre (1998) formulation. If False,
        use the fixed fraction `irr_val`. Default True.
    irr_val : float, optional
        Fixed irreducible water fraction of pore space, used only if
        `coleou_lesaffre` is False. Default 0.02.

    Returns
    -------
    LWCirr : ndarray
        Irreducible LWC per layer [m].
    """
    if not coleou_lesaffre:
        LWCirr = irr_val * phivol_av
        LWCirr[rho >= RhoImp] = 0.
        return LWCirr

    wmi = 0.057 * (rho_i - rho) / rho + 0.017
    wmi[rho >= RhoImp] = 0.

    swi       = np.zeros_like(wmi)
    imsk      = rho < rho_i
    swi[imsk] = (wmi[imsk] / (1 - wmi[imsk])) * rho_i * rho[imsk] \
                / (RHO_W_KGM * (rho_i - rho[imsk]))

    LWCirr = phivol_av * swi
    LWCirr[rho >= RhoImp] = 0.
    return LWCirr
############################
### end _irreducible_lwc ###
############################

def _find_impermeable_layers(rho, dz, RhoImp, ThickImp=0., DownToIce=False):
    """
    Identify impermeable layers (and top-of-lens indices) for percolation
    schemes.

    Three modes, controlled by `DownToIce` and `ThickImp`:

    - DownToIce=True: all layers are permeable until the last layer with
      rho < RhoImp; that layer and everything below it is impermeable
      (water bypasses shallower ice lenses, only trapped once it reaches
      the "ice sheet" proper). `imptop` is the first index of this block.
    - DownToIce=False, ThickImp>0: contiguous runs of layers with
      rho >= RhoImp are ice lenses. A lens is impermeable if its total
      thickness is >= ThickImp, or if it is the bottom-most lens in the
      domain (always impermeable regardless of thickness).
    - DownToIce=False, ThickImp==0: every layer with rho >= RhoImp is
      impermeable, treated as one lens per contiguous run.

    The last layer in the domain is always forced impermeable.

    Parameters
    ----------
    rho : ndarray
        layer density [kg m-3].
    dz : ndarray
        layer thickness [m].
    RhoImp : float
        Density threshold for impermeability [kg m-3].
    ThickImp : float, optional
        Minimum ice-lens thickness [m] for impermeability (ignored if
        DownToIce is True). Default 0.
    DownToIce : bool, optional
        See above. Default False.

    Returns
    -------
    imp : ndarray of int
        Sorted, unique indices of impermeable layers (always includes the
        last layer in the domain).
    imptop : ndarray of int
        Sorted index of the shallowest layer of each impermeable lens/block.
    """
    nnd = len(rho)

    if DownToIce:
        below_thresh = np.where(rho < RhoImp)[0]
        if below_thresh.size > 0:
            imp = np.arange(below_thresh[-1] + 1, nnd)
        else:
            imp = np.arange(0, nnd)

        imptop = imp[:1] if imp.size > 0 else np.array([], dtype=int)

    else:
        lens_top = np.array([ii for ii in range(1, nnd)
                              if rho[ii] >= RhoImp and rho[ii - 1] < RhoImp])
        lens_bot = np.array([ii for ii in range(0, nnd - 1)
                              if rho[ii] >= RhoImp and rho[ii + 1] < RhoImp])
        if rho[0] >= RhoImp:
            lens_top = np.append(0, lens_top)
        lens_bot = np.append(lens_bot, nnd - 1)

        lens_top = lens_top.astype(int)
        lens_bot = lens_bot.astype(int)

        imp    = np.array([], dtype=int)
        imptop = np.array([], dtype=int)
        for ii in range(len(lens_top)):
            lens_thickness = dz[lens_top[ii]:lens_bot[ii] + 1].sum()
            is_bottom_lens = (ii == len(lens_top) - 1)
            if lens_thickness >= ThickImp or is_bottom_lens:
                imp    = np.append(imp, np.arange(lens_top[ii], lens_bot[ii] + 1))
                imptop = np.append(imptop, lens_top[ii])

    imp     = imp.astype(int)
    imptop  = imptop.astype(int)

    if imp.size == 0:
        imp     = np.array([nnd - 1])
        imptop  = imp.copy()

    return np.unique(imp), np.sort(imptop)

###################################
### end _find_impermeable_layers ###
###################################

def bucket(self, iii):
    """
    Route meltwater and rain through the firn column via a bucket-style
    percolation scheme, updating density, temperature, and liquid water
    content (LWC) for one model time step.

    Physical process (in order of operations):
      1. Surface layers are melted according to `self.snowmeltSec[iii]`,
         and the grid is regridded (layers shifted upward, new layers
         appended at the base to retain the same number).
      2. Storage capacity of each layer is computed as the sum of:
           - refreezing capacity (from cold content, limited by pore space)
           - irreducible-retention capacity (Coleou & Lesaffre (1998) or
             fixed-fraction formulation)
      3. Ice lenses are identified as impermeable barriers (by density
         and/or thickness threshold, or via DownToIce mode).
      4. Surface liquid input (melt + rain) is distributed downward into
         available storage capacity until blocked by an impermeable layer
         or storage is exhausted.
      5. LWC in excess of irreducible retention is redistributed into any
         deeper storage capacity, or blocked above impermeable barriers.
      6. Refreezing is applied (in two passes: once after redistribution,
         once more after ponding/runoff, to catch cold content freed up
         by redistribution).
      7. Water blocked by impermeable barriers either ponds in pore space
         above the barrier (if Ponding=True) or runs off; a user-set
         fraction can also run off directly without attempting storage.
      8. Optionally, lateral runoff is applied via the Zuo & Oerlemans
         (1996) parameterization.

    Mass conservation is checked twice (after refreeze pass 1 and after
    refreeze pass 2) and a warning is printed if violated beyond tolerance.

    User-configurable options are read from `self.c`; if missing, defaults
    are used with a warning (see USER CHOICES block below).

    In Zuo and Oerlemans runoff parameterization, surface slope is inferred 
    to be rise/run (m/m, dimensionless). Zuo and Oerlemans (1996) Eq.(22) 
    constants (c1, c2, c3) require a small (~0.001-0.1) dimensionless slope 
    value for the exponential term to behave sensibly; the paper itself 
    doesn't give an explicit formula for computing S, only that "the surface 
    slope S is derived from a fit to the altitude profile along the GIMEX 
    transect" (inferred from differential GPS).
    NOTE: rise/run (m/m) is our best-supported inference from context, not an
    explicitly stated definition in the paper -- treat with appropriate caution.

    Originally coded by Vincent Verjans; edited by Max S.

    Parameters
    ----------
    iii : int
        Current model time step index.

    Returns
    -------
    tuple
        (rho, age, dz, Tz, r2, z, mass, dzn, LWC, meltgridtrack,
         refrozentot, runofftot, dh_melt) — all updated in place on
         `self` and also returned explicitly for the caller to reassign.
    """
    debug = self.c.get('debug_bucket', False)

    ### USER CHOICES ###
    try:
        ColeouLesaffre      = self.c['ColeouLesaffre']
        IrrVal              = self.c['IrrVal'] if not ColeouLesaffre else 0.
        RhoImp              = self.c['RhoImp']
        DownToIce           = self.c['DownToIce']
        ThickImp            = self.c['ThickImp'] if not DownToIce else 0.
        Ponding             = self.c['Ponding']
        DirectRunoff        = self.c['DirectRunoff']
        RunoffZuoOerlemans  = self.c['RunoffZuoOerlemans']
        Slope               = self.c['Slope']
    except KeyError:
        print('You should add the new melt variables to your .json. '
              'See melt.py and example.json. Using defaults.')
        ColeouLesaffre      = True
        IrrVal              = 0.
        RhoImp              = 830.
        DownToIce           = False
        ThickImp            = 0.1
        Ponding             = False
        DirectRunoff        = 0.0
        RunoffZuoOerlemans  = False
        Slope               = 0.1
    
    try:
        keep_firnthickness  = self.c['keep_firnthickness']
    except KeyError:
        print("keep_firnthickness not in .json; setting to False")
        keep_firnthickness  = False
    ### END USER CHOICES ###

    ### Determine mass of melted firn ###
    melt_volume_IE = self.snowmeltSec[iii] * S_PER_YEAR     # [m ie]
    melt_volume_WE = melt_volume_IE * RHO_I_MGM             # [m we]
    melt_mass      = melt_volume_WE * RHO_W_KGM             # [kg]

    # total_liquid_mass_start = np.sum(melt_mass) + np.sum(self.LWC * RHO_W_KGM)

    nnd       = len(self.z)          # number of layers
    rhoi      = RHO_I
    runofftot = 0.              # initialise total runoff [m we]
    mass_sum  = np.cumsum(self.mass)   # cumulative mass [kg]

    ### Melting of surface layers ###
    i_surf_new = np.where(mass_sum > melt_mass)[0][0]   # Index that becomes the new surface node after melting (partially-melted node)
    n_melted   = i_surf_new + 1                                 # number of layers melted

    ### Partially melted layer properties ###
    pm_mass = mass_sum[i_surf_new] - melt_mass      # remaining mass
    pm_dz   = pm_mass / self.rho[i_surf_new]               # remaining thickness
    pm_lwc  = self.LWC[i_surf_new] / self.dz[i_surf_new] * pm_dz  # LWC of the pm layer
    pm_Tz   = T_MELT

    dzo     = self.dz.copy()

    if i_surf_new > 0:
        dh_melt = -1 * (np.sum(dzo[0:i_surf_new]) + (dzo[i_surf_new] - pm_dz))
    else:
        dh_melt = -1 * (dzo[i_surf_new] - pm_dz)
    
    avg_dh_melted = -1 * dh_melt / n_melted

    ### Liquid water input at the surface ###
    liq_input_mass = max(
        melt_mass + (np.sum(self.LWC[0:i_surf_new + 1]) - pm_lwc) * RHO_W_KGM, 0
    )  # avoid negative input due to numerical round-off
    liq_input_vol = liq_input_mass / RHO_W_KGM

    try:
        liq_input_vol = liq_input_vol + self.rainSec[iii] * S_PER_YEAR * RHO_I_MGM  # [m]
    except (AttributeError, IndexError):
        pass  # no rain input provided

    liqmcinit = pm_lwc + np.sum(self.LWC[i_surf_new + 1:]) + liq_input_vol  # mass conservation check

    ### Regridding ###
    if melt_mass > 0:
        if i_surf_new > 0:
            self.rho        = np.concatenate((self.rho[i_surf_new:-1], self.rho[-1] * np.ones(n_melted)))
            self.bdot_mean  = np.concatenate((self.bdot_mean[i_surf_new:-1], self.bdot_mean[-1] * np.ones(n_melted)))
            self.age        = np.concatenate((self.age[i_surf_new:-1], self.age[-1] * np.ones(n_melted)))
            self.Dcon       = np.concatenate((self.Dcon[i_surf_new:-1], self.Dcon[-1] * np.ones(n_melted)))
            self.dzn        = np.concatenate((np.zeros(n_melted), self.dz[1:]))[0:self.compboxes]
            if self.r2 is not None:
                self.r2     = np.concatenate((self.r2[i_surf_new:-1], self.r2[-1] * np.ones(n_melted)))
        else:
            self.dzn        = self.dz[0:self.compboxes]  # avoids bug from undefined self.dzn

        self.LWC = np.concatenate(([pm_lwc], self.LWC[i_surf_new + 1:-1], self.LWC[-1] * np.ones(n_melted)))

        
        if keep_firnthickness:
            nb_th   = np.maximum(avg_dh_melted, self.dz[-1])
            self.dz = np.concatenate(([pm_dz], self.dz[i_surf_new + 1:-1], nb_th * np.ones(n_melted)))
        else:
            self.dz = np.concatenate(([pm_dz], self.dz[i_surf_new + 1:-1], self.dz[-1] * np.ones(n_melted)))

        self.Tz     = np.concatenate(([pm_Tz], self.Tz[i_surf_new + 1:-1], self.Tz[-1] * np.ones(n_melted)))
        self.z      = self.dz.cumsum(axis=0)
        self.z      = np.concatenate(([0], self.z[:-1]))
        self.mass   = self.rho * self.dz

        if self.doublegrid:
            meltgridtrack = np.concatenate((self.gridtrack[i_surf_new:-1], self.gridtrack[-1] * np.ones(n_melted)))
        else:
            meltgridtrack = np.zeros(nnd)
    else:
        meltgridtrack   = self.gridtrack
        self.dzn        = self.dz[0:self.compboxes]
    ### end regridding ###

    ### Calculate excessive LWC (above irreducible holding capacity) ###
    phivol_av   = _available_pore_space(self.rho, self.dz, rho_i=rhoi, cap_density=RHO_I)
    LWCirr      = _irreducible_lwc(self.rho, phivol_av, RhoImp, rho_i=rhoi,
                               coleou_lesaffre=ColeouLesaffre, irr_val=IrrVal)
    LWC_excess  = np.maximum(self.LWC - LWCirr, 0)  # LWC beyond irreducible content [m]

    _debug_report(debug, 'LWC_excess total [kg]', np.sum(LWC_excess) * RHO_W_KGM)
    _debug_report(debug, 'LWC_irr total [kg]', np.sum(LWCirr) * RHO_W_KGM)

    ### Refreezing capacity from cold content, and retention capacity after refreezing ###
    ### 'excess' LWC = beyond irreducible value; 'additional' LWC = beyond what
    ### is currently present, up to the irreducible value.
    ### (In theory, refr_cap should be 0 for layers with existing LWC, except the uppermost.)
    cold_content            = CP_I * self.mass * (T_MELT - self.Tz)       # cold content [J]
    refr_cap_from_cc        = (cold_content / LF_I) / RHO_W_KGM           # refreeze capacity from cold content [m we]
    refr_cap_from_cc_supp   = np.maximum(0, refr_cap_from_cc - self.LWC)  # capacity for LWC beyond what's present [m we]
    refr_cap                = np.minimum(refr_cap_from_cc, phivol_av)     # total refreeze capacity [m we]
    
    LWC_to_ice    = RHO_W_KGM / RHO_I * self.LWC                      # volume existing LWC would take if frozen
    refr_vol_supp = np.maximum(0, phivol_av - LWC_to_ice)             # volume available for additional LWC [m]
    refr_cap_supp = np.minimum(refr_cap_from_cc_supp, refr_vol_supp)  # refreeze capacity for additional LWC [m we]

    ### Potential density/porosity after maximal refreezing, to compute retention capacity ###
    rho_pot       = (self.mass + refr_cap * RHO_W_KGM) / self.dz        # potential density after refreezing [kg m-3]
    phivol_av_pot = _available_pore_space(rho_pot, self.dz, rho_i=rhoi, cap_density=RHO_I)
    LWCirr_pot    = _irreducible_lwc(rho_pot, phivol_av_pot, RhoImp, rho_i=rhoi,
                                   coleou_lesaffre=ColeouLesaffre, irr_val=IrrVal)

    LWC_unf     = np.maximum(0, self.LWC - refr_cap)   # unfrozen LWC remaining after refreeze [m]
    retcap_supp = np.maximum(0, LWCirr_pot - LWC_unf)  # retention capacity for additional LWC [m]
    stcap       = refr_cap_supp + retcap_supp          # total storage capacity for additional LWC [m]

    ### Ice lens algorithm: find impermeable layers ###
    imp, _ = _find_impermeable_layers(self.rho, self.dz, RhoImp,
                                      ThickImp=ThickImp, DownToIce=DownToIce)

    stcap[imp] = 0.                 # zero storage capacity for impermeable layers
    stcap_cum  = np.cumsum(stcap)   # cumulative storage capacity

    ### Store surface melt according to stcap of layers from surface to bottom ###
    LWCblocked = np.zeros(nnd)

    if liq_input_vol > 0:
        ### i_in_bot: Bottom-most layer index that receives distributed surface liquid input (either where 
        ###        cumulative storage capacity meets liq_in_vol, or capped by the nearest impermeable barrier)
        if stcap_cum[-1] >= liq_input_vol: # There is room in the column to accomodate all the liquid input
            i_in_bot = np.where(stcap_cum >= liq_input_vol)[0][0]  
        else:
            i_in_bot = nnd - 1  # not enough capacity anywhere; fall through to bottom

        if i_in_bot >= imp[0]: # Impermeable barrier (or lack of pore space) prevents full distribution  
            i_in_bot             = max(0, imp[0] - 1)
            storageinp           = np.concatenate((stcap[0:i_in_bot + 1], np.zeros(nnd - i_in_bot - 1)))
            LWCblocked[i_in_bot] = liq_input_vol - np.sum(storageinp)  # excess water blocked above barrier
        else: # No impermeable barrier reached; distribute by storage capacity
            if i_in_bot == 0:
                storageinp = np.concatenate(([liq_input_vol], np.zeros(nnd - 1)))
            else:
                storageinp = np.concatenate(
                    (stcap[0:i_in_bot], [liq_input_vol - stcap_cum[i_in_bot - 1]], np.zeros(nnd - i_in_bot - 1))
                )
    elif liq_input_vol == 0:
        storageinp = np.zeros(nnd)
    else:
        print('Negative liquid input! Check your inputs. Exiting.')
        sys.exit()

    _debug_report(debug, 'LWCblocked after initial storage (C1) [kg]',
                  np.sum(LWCblocked) * RHO_W_KGM)

    stcap_remaining = stcap - storageinp  # updated storage capacity

    ### Set LWC_excess in impermeable layers as blocked LWC ###
    indsblc             = np.intersect1d(np.where(LWC_excess > 0)[0], imp)
    LWCblocked[indsblc] += LWC_excess[indsblc]
    self.LWC[indsblc]   -= LWC_excess[indsblc]
    LWC_excess[indsblc] = 0.

    _debug_report(debug, 'LWCblocked after excess-in-impermeable transfer (C2) [kg]',
                  np.sum(LWCblocked) * RHO_W_KGM)

    ### Distribute LWC_excess into layers with storage capacity, or block above impermeable barriers ###
    ### j_scan: The moving "cursor" index — starts at the uppermost unprocessed node with excess LWC, 
    ###         and advances downward through the column as each segment is resolved
    ### j_target: The node the current segment's excess LWC is being routed to — either the next node 
    ###           with storage capacity, or the node just above an impermeable barrier if one is hit first
    LWC_redist      = np.copy(self.LWC)  # LWC will be modified by LWC_excess transfers
    storage_excess  = np.zeros(nnd)      # LWC stored in each layer via this redistribution

    if np.any(LWC_excess) > 0:
        tostore     = 0                           # running total of excess LWC awaiting storage/blocking in the current segment
        inds_ex     = np.where(LWC_excess > 0)[0] # layers with excess LWC
        ind_ex_bot  = inds_ex[-1]                 # bottom-most layer with excess LWC
        j_scan      = inds_ex[0]                  # start from uppermost layer with excess LWC

        if np.any(stcap_remaining > 0):
            ind_st_bot = np.where(stcap_remaining > 0)[0][-1]  # bottom-most layer that can still store LWC_excess
        else:
            ind_st_bot = 0

        if np.any(stcap_remaining[1:] > 0) and (ind_st_bot > j_scan): 
            ### There is storage capacity in the column, deeper than j_scan
            ### Loop terminates once all excess layers are processed AND nothing remains to be stored/blocked
            while (j_scan <= ind_ex_bot) or (tostore > 0):
                if (np.where(stcap_remaining[j_scan:] > 0)[0]).size > 0:
                    j_target = j_scan + np.where(stcap_remaining[j_scan:] > 0)[0][0]  # next layer with storage capacity
                else:
                    ### No remaining storage capacity below j_scan: block excess above each
                    ### underlying impermeable barrier (mirrors the "no storage" branch below)
                    inds_ex_deep = inds_ex[inds_ex >= j_scan]
                    
                    for jj2 in inds_ex_deep:
                        matches = np.where(imp >= jj2)[0]
                        if matches.size > 0:
                            j_target = imp[matches[0]] - 1
                        else:
                            _debug_report(True, 'WARNING: no impermeable layer found at/after index', jj2)
                            j_target = nnd - 1
                        
                        LWCblocked[j_target] += LWC_excess[jj2]
                        LWC_redist[jj2]       = LWCirr[jj2]
                    break  # exit the while loop

                ### Compare depth of nearest impermeable barrier at/below j_scan vs. j_target:
                ### if the barrier is deeper than j_target, water can reach j_target unobstructed.
                if imp[np.where(imp >= j_scan)[0][0]] > j_target:
                    ### j_scan and j_target are not separated by an impermeable barrier
                    tostore                           += np.sum(LWC_excess[j_scan:j_target + 1]) # accumulate excess along the way
                    LWC_redist[j_scan:j_target + 1]    = np.minimum(LWC_redist[j_scan:j_target + 1], LWCirr[j_scan:j_target + 1]) # trim excess down to irreducible value
                    storage_excess[j_target]           = min(stcap_remaining[j_target], tostore) # j_target stores as much of tostore as it can hold
                    tostore                           -= storage_excess[j_target] # remainder (if any) carries forward
                    j_scan                             = j_target + 1 # advance scan index past the layer just filled
                    if j_scan >= ind_st_bot:
                        ### No further storage capacity remains in the column: block whatever is left over
                        j_target             = imp[np.where(imp >= j_scan)[0][0]] - 1
                        LWCblocked[j_target] += tostore
                        tostore              = 0.
                else:
                    ### Impermeable barrier exists between j_scan and j_target
                    j_target                         = imp[np.where(imp >= j_scan)[0][0]] - 1 # barrier is redefined as the node just above it
                    tostore                         += np.sum(LWC_excess[j_scan:j_target + 1])
                    LWC_redist[j_scan:j_target + 1]  = np.minimum(LWC_redist[j_scan:j_target + 1], LWCirr[j_scan:j_target + 1]) # all accumulated excess is blocked here, since it can't pass the barrier
                    LWCblocked[j_target]            += tostore
                    tostore                          = 0.
                    
                    if j_target < ind_ex_bot:
                        j_scan = inds_ex[np.where(inds_ex > j_target)[0][0]] # jump to next layer with excess LWC, past the barrier
                    else:
                        j_scan = ind_ex_bot + 1  # all excess layers handled; terminate while loop

            _debug_report(debug, 'LWCblocked after redistribution loop (C3) [kg]',
                          np.sum(LWCblocked) * RHO_W_KGM)
        else:
            # No usable storage capacity in the column: block all excess above nearest barrier
            for j_scan in inds_ex:
                matches = np.where(imp >= j_scan)[0]
                if matches.size > 0:
                    j_target = imp[matches[0]] - 1
                else:
                    _debug_report(True, 'WARNING: no impermeable layer found at/after index', j_scan)
                    j_target = nnd - 1
                
                LWCblocked[j_target] += LWC_excess[j_scan]
                LWC_redist[j_scan]    = LWCirr[j_scan]

            _debug_report(debug, 'LWCblocked, no storage capacity available (C4) [kg]',
                          np.sum(LWCblocked) * RHO_W_KGM)

    ### Combine storage from surface-input distribution (storageinp, Section 3) with
    ### storage from this excess-LWC redistribution (storage_excess) into final LWC profile
    storagetot = storageinp + storage_excess
    LWC_redist = LWC_redist + storagetot

    ### Refreezing (pass 1) ###
    ### Refreeze either all available liquid (if cold content allows) or as much
    ### as capacity permits. refr_cap already reflects both the energy budget
    ### (cold content) and the physical pore-space limit (Section 2).
    freeze      = np.minimum(LWC_redist, refr_cap)  # refreezing per layer [m we]
    self.mass   = self.mass + RHO_W_KGM * freeze    # update mass [kg]
    self.LWC    = LWC_redist - freeze        # update LWC [m], remaining unfrozen water
    self.rho    = self.mass / self.dz               # update density [kg m-3]
    ### NOTE: dz is unchanged here -- refreezing increases density in place,
    ### consistent with the fixed-volume assumption used throughout this scheme
    ### (see _available_pore_space docstring).

    latheat        = freeze * RHO_W_KGM * LF_I     # latent heat released [J]
    cold_content  -= latheat                       # remaining cold content [J]
    ### Clean up tiny negative cold_content from floating-point round-off
    ### (should be exactly >=0 in theory, since freeze was capped by refr_cap,
    ### which was itself derived from cold_content).
    cold_content[(cold_content < 0) & (cold_content > -1e-9)] = 0  # clean numerical noise

    refrozentot         = np.sum(freeze)                    # total refrozen water [m we]
    deltaT              = cold_content / (CP_I * self.mass) # remaining temperature deficit below T_MELT
    ### Only layers that actually refroze something get a temperature update --
    ### layers with freeze==0 had no phase change, so Tz is left untouched.
    self.Tz[freeze > 0] = T_MELT - deltaT[freeze > 0]

    ### Dry cold-firn sanity check ###
    ### Physically, a layer colder than T_MELT should hold no liquid water --
    ### any LWC there should already have been refrozen above. This check
    ### catches (a) negligible numerical residue, which is zeroed out silently,
    ### and (b) any genuine leftover LWC in a cold layer, which indicates a bug
    ### upstream and is reported for debugging.
    coldlayers = np.where(self.Tz < T_MELT)[0]
    
    if np.all(self.LWC[coldlayers] < 1e-9):
        # Residual LWC is negligible (floating-point noise) -- safe to zero out
        self.LWC[coldlayers] = 0.
    if np.any(self.LWC[coldlayers] > 0.):
        # Genuine problem: some cold layer still holds non-negligible LWC
        # after refreezing. This should not happen -- flag loudly regardless
        # of the debug flag, since it indicates a mass/energy conservation bug.
        print('#############')
        print('Problem: water content in a cold layer (bucket, post-refreeze #1)')
        print(f'iii: {iii}')
        xx = np.where((self.LWC > 0) & (self.Tz < T_MELT))[0]
        print(f'Layer depths: {self.z[xx]}')
        print(f'Layer LWC: {self.LWC[xx]}')
        print(f'Layer T: {self.Tz[xx]}')
        print(f'Layer rho: {self.rho[xx]}')
        print('#############')

    ### Direct runoff of a fraction of blocked LWC (user choice) ###
    ### DirectRunoff represents water assumed to escape laterally immediately,
    ### bypassing any attempt at ponding above the barrier.
    runofftot   = runofftot + DirectRunoff * np.sum(LWCblocked)
    LWCblocked  = (1 - DirectRunoff) * LWCblocked

    _debug_report(debug, 'runofftot after direct-runoff fraction (R1) [kg]',
                  runofftot * RHO_W_KGM)

    if np.any(LWCblocked > 0):
        if Ponding:
            ### Ponding: water blocked above an impermeable ice lens fills
            ### available pore space upward through the firn, forming a
            ### perched saturated zone (a "perched water table")
            phiempty        = self.dz * (rhoi - self.rho) / RHO_W_KGM - self.LWC
            phiempty[imp]   = 0.               # no ponding directly in impermeable layers
            phiempty[self.rho > RhoImp] = 0.   # handles layers that became impermeable via refreezing

            ### Process blocked layers from deepest to shallowest: ponding above
            ### a barrier fills contiguous space upward, so resolving the
            ### deepest barrier first avoids conflicts if multiple blocked
            ### layers are stacked in the same column.
            for kk in np.flip(np.where(LWCblocked > 0)[0]):
                ### Cumulative empty pore space from kk upward, reversed so
                ### index 0 = kk itself, increasing as we go shallower
                phiempty_cumf = np.cumsum(np.flip(phiempty[0:kk + 1]))

                if phiempty_cumf[-1] >= LWCblocked[kk]:
                    ifill = np.where(phiempty_cumf > LWCblocked[kk])[0][0]  # layer[kk-ifill] holds the remainder
                else:
                    ### Not enough pore space even up to the surface: excess runs off
                    ifill           = kk
                    runofftot       = runofftot + LWCblocked[kk] - phiempty_cumf[kk]
                    LWCblocked[kk]  = phiempty_cumf[kk]

                if ifill == 0:
                    ### The blocked layer itself has enough empty space; no need to pond upward
                    self.LWC[kk] += LWCblocked[kk]
                    phiempty[kk] -= LWCblocked[kk]
                else:
                    ### Fill layers [kk-ifill+1, kk] completely (using all their empty space),
                    ### then deposit the remainder into layer kk-ifill
                    self.LWC[kk - ifill + 1:kk + 1] += phiempty[kk - ifill + 1:kk + 1]
                    LWCblocked[kk]                  -= np.sum(phiempty[kk - ifill + 1:kk + 1])
                    phiempty[kk - ifill + 1:kk + 1]  = 0.
                    self.LWC[kk - ifill]            += LWCblocked[kk]
                    phiempty[kk - ifill]            -= LWCblocked[kk]

                    if np.any(self.LWC < 0):
                        self.LWC[self.LWC < 0] = 0.0 # guard against floating-point overshoot

                LWCblocked[kk] = 0. # LWCblocked[kk] fully accommodated (stored or run off)

            _debug_report(debug, 'runofftot after ponding (R2, Ponding=True) [kg]',
                          runofftot * RHO_W_KGM)
            _debug_report(debug, 'LWCblocked after ponding (should be ~0) [kg]',
                          np.sum(LWCblocked) * RHO_W_KGM)
        else:
            ### No ponding: any water blocked above an impermeable barrier
            ### runs off outright, with no attempt to store it as a perched
            ### saturated zone.
            runofftot   = runofftot + np.sum(LWCblocked)
            LWCblocked  = 0 * LWCblocked # LWCblocked fully accounted for (all became runoff)

            _debug_report(debug, 'runofftot after blocked-LWC runoff (R2, Ponding=False) [kg]',
                          runofftot * RHO_W_KGM)

    ### Zuo and Oerlemans (1996) lateral runoff routine ###
    ### Optional additional runoff pathway, applied after refreezing and after
    ### the direct-runoff/ponding logic above. Operates on whatever LWC remains
    ### in self.LWC at this point (already post-refreeze, post-ponding).
    if RunoffZuoOerlemans:
        ### Recompute pore space and irreducible LWC using post-refreezing
        ### self.rho (density has changed since Section 2/4's calculations,
        ### so these are deliberately recomputed here rather than reused).
        phivol_av   = _available_pore_space(self.rho, self.dz, rho_i=rhoi, cap_density=RHO_I)
        LWCirr      = _irreducible_lwc(self.rho, phivol_av, RhoImp, rho_i=rhoi,
                                   coleou_lesaffre=ColeouLesaffre, irr_val=IrrVal)
        LWC_rfZO    = np.maximum(0, self.LWC - LWCirr)  # mobile LWC subject to lateral runoff [m]

        if np.any(LWC_rfZO > 0):
            indsrfZO    = np.where(LWC_rfZO > 0)[0]
            c1zuo       = 1.5 * 24 * 3600    # Zuo & Oerlemans (1996) constant [s]
            c2zuo       = 25. * 24 * 3600    # Zuo & Oerlemans (1996) constant [s]
            c3zuo       = 140.                # Zuo & Oerlemans (1996) constant [/]
            
            ### Characteristic drainage timescale, Eq.(22): steeper Slope -> smaller
            ### tstar -> faster drainage.
            tstar = c1zuo + c2zuo * np.exp(-1 * c3zuo * Slope)

            ### Forward-Euler discretization of the loss term (-W/t*) from the
            ### governing ODE, Eq.(21): dW/dt = Pw - W/t*. This is not a direct
            ### quote of Eq.(21) itself (which describes production + loss for
            ### the accumulated meltwater W); here we only apply the loss term
            ### to the currently mobile LWC over one model timestep.
            rfZO           = np.zeros(nnd)
            rfZO[indsrfZO] = self.dt[iii] * LWC_rfZO[indsrfZO] / tstar  # [m]

            self.LWC  = self.LWC - rfZO
            runofftot = runofftot + np.sum(rfZO)

            _debug_report(debug, 'runofftot after Zuo-Oerlemans runoff (R3) [kg]',
                          runofftot * RHO_W_KGM)

    ### Mass conservation check ###
    ### Total liquid water should be conserved across: what remains liquid,
    ### what has refrozen so far, and what has run off, relative to the
    ### liquid input computed near the start of the function (liqmcinit).
    liqmcfinal = np.sum(self.LWC) + refrozentot + runofftot
    if abs(liqmcfinal - liqmcinit) > 1e-3:
        print(f'Mass conservation error (bucket, check 1) at step {iii}\n'
              f'    Init: {liqmcinit} m\n    Final: {liqmcfinal} m')

    ### Refreezing (pass 2) ###
    ### Ponding, direct runoff, and lateral runoff (Sections 6-7) can move LWC
    ### into layers that were not accounted for in refreeze pass 1 (Section 5),
    ### and those layers may still have leftover cold content available. This
    ### second pass catches any additional refreezing made possible by that
    ### redistribution.
    ### NOTE: cap_density was 916.99 in earlier versions. 
    phivol_av = _available_pore_space(self.rho, self.dz, rho_i=rhoi, cap_density=RHO_I)

    ### Cold content is recomputed fresh here (not reused from Section 5's
    ### cold_content, which has already been partially consumed by pass 1).
    cold_content_new      = CP_I * self.mass * (T_MELT - self.Tz)       # cold content [J], post-redistribution
    refr_cap_from_cc_new  = (cold_content_new / LF_I) / RHO_W_KGM       # refreeze capacity from cold content [m we]
    refr_cap_new          = np.minimum(refr_cap_from_cc_new, phivol_av) # total refreeze capacity [m we]

    freeze      = np.minimum(self.LWC, refr_cap_new) # refreezing per layer, this pass [m we]
    self.mass   = self.mass + RHO_W_KGM * freeze
    self.LWC    = self.LWC - freeze
    self.rho    = self.mass / self.dz

    latheat           = freeze * RHO_W_KGM * LF_I
    cold_content_new -= latheat
    cold_content_new[(cold_content_new < 0) & (cold_content_new > -1e-9)] = 0 # clean numerical noise

    refrozentot         = refrozentot + np.sum(freeze) # running total across both refreeze passes
    
    deltaT              = cold_content_new / (CP_I * self.mass)
    self.Tz[freeze > 0] = T_MELT - deltaT[freeze > 0] # running total across both refreeze passes  

    ### Dry cold-firn sanity check ###
    ### Same rationale as post-refreeze #1: a layer colder than T_MELT should
    ### hold no liquid water after refreezing. Negligible residue is zeroed
    ### out silently; genuine leftover LWC indicates a conservation bug and is
    ### reported regardless of the debug flag.
    coldlayers = np.where(self.Tz < T_MELT)[0]

    if np.all(self.LWC[coldlayers] < 1e-9):
        self.LWC[coldlayers] = 0.

    if np.any(self.LWC[coldlayers] > 0.):
        print('Problem: water content in a cold layer (bucket, post-refreeze #2)')
        xx = np.where((self.LWC > 0) & (self.Tz < T_MELT))[0]
        print(f'Layer depths: {self.z[xx]}')
        print(f'Layer LWC: {self.LWC[xx]}')
        print(f'Layer T: {self.Tz[xx]}')
        print(f'Layer rho: {self.rho[xx]}')

    self.rho[self.rho > RHO_I] = RHO_I  # clip numerical overshoot above ice density

    ### Mass conservation check 2 ###
    ### Final check, using a tighter tolerance than check 1 since this is after
    ### all redistribution, ponding, and refreezing passes are complete.
    liqmcfinal = np.sum(self.LWC) + refrozentot + runofftot
    if abs(liqmcfinal - liqmcinit) > 1e-5:
        print(f'Mass conservation error (bucket, check 2) at step {iii}\n'
              f'    Init: {liqmcinit} m\n    Final: {liqmcfinal} m')

    return self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, \
        self.dzn, self.LWC, meltgridtrack, refrozentot, runofftot, dh_melt

##################
### end bucket ###
##################

##########################
def darcyscheme(self,iii):
    '''
    Darcy-flow percolation scheme for meltwater routing through firn.

    A more physically detailed alternative to the bucket scheme (`bucket()`),
    this routine solves for actual water flux between adjacent firn layers
    using an unsaturated Darcy-flow formulation, following Hirashima et al.
    (2010), with van Genuchten unsaturated-flow parameters (grain-size- and
    density-dependent, via Yamaguchi et al. and Calonne et al.
    parameterizations) and an upwind/downwind conductivity scheme at layer
    interfaces (Szymkiewicz, 2009).

    Because Darcy flow can require much finer time resolution than a single
    outer model time step to remain numerically stable, this scheme runs its
    own internal adaptive sub-stepping loop, advancing from timer=0 to
    timer=self.dt[iii] in increments of dtsub (bounded by dtmin/dtmax and
    adjusted based on how close the resolved flux came to its equilibrium
    limit each sub-step).

    At each sub-step, in order:
      1. Surface nodes are melted (mass supplied at a constant rate,
         meltflux_mass, over the sub-step), and the grid is regridded, similar
         to the surface-melting logic in `bucket()`.
      2. Ice lenses are (re-)identified as impermeable barriers, only when the
         density profile has changed since the last sub-step (a caching
         optimization, since this scan is relatively expensive to repeat
         every sub-step).
      3. Node density is artificially capped at `rholim_dcy` for the purposes
         of flow calculations only, giving ice lenses a small residual
         permeability rather than an abrupt, numerically unstable wall of
         zero porosity. Any leftover LWC in these density-modified nodes is
         forced to run off at the very end of the routine, since their true
         density does not actually support long-term water storage.
      4. Liquid input (melt + rain) for this sub-step is distributed
         downward into available accommodation space at the surface,
         analogous to `bucket()`'s surface-input distribution, but based
         purely on physical pore space rather than combined
         refreeze+retention storage capacity.
      5. Effective saturation, van Genuchten shape parameters, and hydraulic
         conductivity (saturated and relative) are computed for each node.
      6. Nodes that cannot have any downward outflow are identified: dry
         nodes, impermeable nodes, nodes directly above an impermeable
         barrier, and nodes sitting in a fully-saturated column directly
         above a barrier (a perched water table with no remaining pore
         space to push water into via Darcy flux).
      7. For each remaining candidate interface, the limiting flux ("qlim",
         Hirashima et al. 2010 Eq. 20) is solved for iteratively -- an
         analytic first guess (or a warm-started guess from the previous
         sub-step's flux, if flux has been stable), refined via bisection or
         Newton-Raphson depending on how far the initial guess is from
         equilibrium.
      8. The actual water flux transferred this sub-step is computed via an
         exponential relaxation toward the limiting flux (Hirashima et al.
         2010 Eq. 23), applied to update LWC. Lateral runoff is optionally
         computed for water with no downward outflow path.
      9. Refreezing is applied based on each node's cold content and
         available pore volume.
      10. The sub-step size is adapted based on how close the resolved flux
          came to the limiting flux, then clamped to [dtmin, dtmax] and to
          not overshoot the outer model time step.

    Configuration options (read from `self.c`; defaults used with a warning
    if missing -- see USER CHOICES block below):

    dtmin_darcy, dtmax_darcy, dtsub_darcy : float [s]
        Minimum, maximum, and initial internal sub-step size.
    RhoImp_darcy : float [kg m-3]
        Density threshold above which a node is considered part of an
        impermeable ice lens.
    ThickImp_darcy : float [m]
        Minimum ice-lens thickness required for a lens to be impermeable.
    lat_runoff_darcy : bool
        If True, compute lateral runoff for nodes with available water but
        no downward outflow path.
    slope_proxy_darcy : float [m/m], dimensionless
        Surface slope used in the lateral-runoff parameterization. As with
        `bucket()`'s analogous `Slope` option, this is best understood as a
        dimensionless rise/run value, inferred from context rather than an
        explicit definition in the cited source -- see inline comments near
        its use for details.
    rholim_dcy : float [kg m-3]
        Density cap used only for Darcy flow calculations (not the true
        density), allowing ice lenses some residual permeability rather
        than an abrupt wall of zero porosity. Must be greater than
        RhoImp_darcy.
    ColeouLesaffre, IrrVal : bool, float
        Same meaning and same config keys as in `bucket()` -- governs
        whether irreducible water content is computed via the Coleou &
        Lesaffre (1998) density-based formulation or a fixed pore-space
        fraction. Shared between both percolation schemes for consistency.

    References
    ----------
    Calonne, N. et al. -- saturated hydraulic conductivity parameterization
        (`hydrconducsat_Calonne`).
    Coleou, C. and Lesaffre, B. (1998). Irreducible water saturation in
        snow. Annals of Glaciology, 26, 64-68.
    Cuffey, K.M. and Paterson, W.S.B. (2010). The Physics of Glaciers, 4th
        ed. -- specific heat of ice parameterization, Eq. 9.1.
    Hirashima, H. et al. (2010). Numerical modeling of liquid water
        movement through layered snow based on new measurements of the
        water retention curve. Cold Regions Science and Technology, 64(2),
        94-103. -- limiting-flux (qlim) formulation, Eqs. 1, 5, 11, 20, 23.
    Szymkiewicz, R. (2009) -- upwind/downwind conductivity selection at
        layer interfaces (implemented as `bigkedg` in this function).
    Yamaguchi, S. et al. (2010, 2012) -- van Genuchten parameter
        parameterization (`vG_Yama_params`) and saturated water content
        discussion.

    Parameters
    ----------
    iii : int
        Current model time step index.

    Returns
    -------
    tuple
        (rho, age, dz, Tz, r2, z, mass, dzn, LWC, meltgridtrack, refr_tot,
         runofftot) -- all updated in place on `self` and also returned
         explicitly for the caller to reassign.

    Notes
    -----
    Originally coded by Vincent Verjans. As of this documentation pass, the
    configuration options above were newly wired to `self.c` (previously
    hardcoded); two operator-precedence bugs in the flux-solving loop and
    the post-loop sanity checks were identified and corrected (see version
    history / code review notes for details).
    '''

    ticdarcy = time.time()
    timetot  = self.dt[iii]      # total duration to be covered by the Darcy routine
    dtsub    = 60                # [s] duration of Darcy time steps, adjusted iteratively

    ### USER CHOICES ###
    try:
        dtmin           = self.c['dtmin_darcy']         # [s] minimal time step for the Darcy routine
        dtmax           = self.c['dtmax_darcy']         # [s] maximal time step for the Darcy routine
        dtsub           = self.c['dtsub_darcy']         # [s] starting time step for the Darcy routine
        RhoImp          = self.c['RhoImp_darcy']        # density from which a layer is impermeable [kg m-3]
        ThickImp        = self.c['ThickImp_darcy']      # minimum ice-lens thickness for impermeability [m]
        lat_runoff      = self.c['lat_runoff_darcy']    # compute lateral runoff above impermeable layers [True/False]
        slope_proxy     = self.c['slope_proxy_darcy']   # proxy for the slope [dimensionless]
        rholim_dcy      = self.c['rholim_dcy']          # max density used in Darcy calcs to allow some permeability [kg m-3]
        eps_cvg         = self.c['eps_cvg_darcy']       # convergence criterion for equilibrium head [m]
        runoff_method   = self.c['runoff_method_darcy'] # 'ZuoOerlemans' or 'Darcy'
        ColeouLesaffre  = self.c['ColeouLesaffre']       # reuses the same flag as bucket()
        IrrVal          = self.c['IrrVal'] if not ColeouLesaffre else 0.
    except KeyError:
        print('You should add the new Darcy-scheme variables to your .json. '
              'See melt.py and example.json. Using defaults.')
        dtmin           = 60
        dtmax           = 3600
        dtsub           = 60.
        RhoImp          = 873.
        ThickImp        = 0.5
        lat_runoff      = True
        slope_proxy     = 0.02
        rholim_dcy      = 910.
        eps_cvg         = 0.1e-3
        runoff_method   = 'Darcy'
        ColeouLesaffre  = True
        IrrVal          = 0.02
    ### END USER CHOICES ###

    ### Surface fluxes ###
    melt_mass_tot = self.snowmeltSec[iii] * S_PER_YEAR * RHO_I       # total melt over the Darcy routine [kg]
    meltflux_mass = melt_mass_tot / timetot                          # melt flux throughout the Darcy routine [kg s-1]
    try:
        rain_vol_tot = self.rainSec[iii] * S_PER_YEAR * RHO_I_MGM    # total rain [m we]
    except (AttributeError, IndexError):
        rain_vol_tot = 0.
    rainflux_vol = rain_vol_tot / timetot   # rain flux throughout the Darcy routine [m we s-1]

    eps_cvg  = 0.1e-3 #convergence criterion for equilibrium head when solving for qlim of Hirashima et al. (2010) [m]

    ### Firn variables ###
    ncv            = len(self.z)         # number of nodes (number of control volumes)
    initial_lwc    = np.copy(self.LWC)   # LWC before Darcy scheme
    phi            = (RHO_I-self.rho)/RHO_I  # update porosity
    phi            = np.maximum(0,phi)   # avoid numerical errors of very small negative phi
    cp_i           = 152.5+7.122*self.Tz # specific heat of ice [J kg-1 K-1] Cuffey and Paterson 2010 (9.1)
    rg             = np.sqrt(self.r2)    # grain radius [m]
    runofftot      = 0.                  # total runoff over the entire Darcy routine
    refr_tot       = 0.                  # total refreezing over the entire Darcy routine
    dltz           = np.append(self.dz[0:-1]/2+self.dz[1:]/2,self.dz[-1]/2) # distance between centres of nodes

    if self.doublegrid: #if we have doublegrid: need to adjust gridtrack
        meltgridtrack  = np.copy(self.gridtrack) #prepare gridtrack adjusted for melting
    elif self.doublegrid==False:
        meltgridtrack = np.zeros_like(self.dz) # just return a zero array

    timer      = 0                  # timer of the Darcy routine
    rho_lens0  = np.copy(self.rho)  # initialise rho used for lenses at the Darcy time step
    glwflux_d2 = np.zeros(ncv-1)    # glw flux computed at Darcy step -2
    glwflux_d1 = np.zeros(ncv-1)    # glw flux computed at Darcy step -1
    if RhoImp>rholim_dcy:
        print('RhoImp must be below rholim_dcy in Darcy scheme, exiting')
        sys.exit()

    while timer < timetot:  
        ### Melt upper nodes ###
        if meltflux_mass > 0:
            meltstep_mass = meltflux_mass * dtsub             # melt mass over dtsub time interval [kg]
            mass_sum      = np.cumsum(self.mass)              # depth-cumulated mass (local; not self.mass_sum)
            ind1          = np.where(mass_sum >= meltstep_mass)[0][0]  # bottom-most node affected by melt
            pm_mass       = mass_sum[ind1] - meltstep_mass     # mass of partially melted node
            pm_dz         = pm_mass / self.rho[ind1]           # thickness of partially melted node
            pm_lwc        = self.LWC[ind1] / self.dz[ind1] * pm_dz  # LWC of partially melted node
            lwc_p         = np.sum(self.LWC[0:ind1 + 1]) - pm_lwc  # LWC of melted firn contributing to percolation
            n_mlt         = ind1 + 1                           # number of nodes melted, including the partial node

            self.dz  = np.concatenate(([pm_dz], self.dz[ind1 + 1:-1], self.dz[-1] * np.ones(n_mlt)))
            self.dzn = np.concatenate((np.zeros(n_mlt), self.dz[1:]))
            self.dzn = self.dzn[0:self.compboxes]
            self.z   = np.cumsum(self.dz)
            self.z   = np.append(0, self.z[0:-1])
            dltz     = np.append(self.dz[0:-1] / 2 + self.dz[1:] / 2, self.dz[-1] / 2)  # distance between node centres
            self.rho = np.append(self.rho[ind1:-1], self.rho[-1] * np.ones(n_mlt))
            phi      = (RHO_I - self.rho) / RHO_I              # update porosity
            phi      = np.maximum(0, phi)                      # avoid numerical errors from very small negative phi
            self.mass = self.dz * self.rho
            self.age  = np.append(self.age[ind1:-1], self.age[-1] * np.ones(n_mlt))
            self.LWC  = np.concatenate(([pm_lwc], self.LWC[ind1 + 1:-1], self.LWC[-1] * np.ones(n_mlt)))
            self.Tz   = np.append(self.Tz[ind1:-1], self.Tz[-1] * np.ones(n_mlt))
            cp_i      = 152.5 + 7.122 * self.Tz                # specific heat of ice [J kg-1 K-1] (Cuffey & Paterson 2010, Eq. 9.1)
            self.r2   = np.append(self.r2[ind1:-1], self.r2[-1] * np.ones(n_mlt))
            rg        = np.sqrt(self.r2)                       # grain radius [m]

            if self.doublegrid:
                meltgridtrack = np.concatenate((meltgridtrack[ind1:-1], meltgridtrack[-1] * np.ones(n_mlt)))

            # Liquid water input at the surface node
            liq_input = meltstep_mass / RHO_W_KGM + lwc_p + rainflux_vol * dtsub  # [m we]

        else:
            # Only rain as liquid water input; no melting this sub-step
            n_mlt         = 0
            self.dzn      = self.dz[0:self.compboxes]
            liq_input     = rainflux_vol * dtsub  # [m we]
            meltstep_mass = 0.
         
        ### Spot ice lenses and evaluate their thickness ###
        if timer == 0 or (n_mlt > 1 or self.rho[0] >= RhoImp
                           or len(rho_lens0[rho_lens0 >= RhoImp]) != len(self.rho[self.rho >= RhoImp])):
            # Only re-run the (relatively expensive) lens-detection scan when the
            # density profile has changed since the last sub-step, or on the
            # first sub-step. Otherwise, reuse imp/imptop from the previous sub-step.
            imp, imptop = _find_impermeable_nodes(self.rho, self.dz, RhoImp, ThickImp=ThickImp, DownToIce=False)

        rho_lens0 = np.copy(self.rho)  # rho snapshot used to detect lens-distribution changes at the next sub-step

        ### Define all Darcy-scheme variables ###
        rhodcy = np.minimum(self.rho, rholim_dcy)   # cap density to allow residual percolation through ice lenses
        rhomod = np.where(self.rho != rhodcy)[0]    # nodes whose density was modified for flow purposes only
        phidcy = (RHO_I - rhodcy) / RHO_I           # porosity computed using the capped (rhodcy) density

        theta_w = self.LWC / self.dz                # volumetric water content

        # Available pore space and irreducible LWC, computed via the same
        # helpers used in bucket(), for consistency between percolation schemes.
        phivol_av_dcy = _available_pore_space(rhodcy, self.dz, rho_i=RHO_I, cap_density=RHO_I)
        LWC_i = _irreducible_lwc(rhodcy, phivol_av_dcy, RhoImp, rho_i=RHO_I,
                                  coleou_lesaffre=ColeouLesaffre, irr_val=IrrVal)

        theta_s = phidcy * RHO_I_MGM                # theta at saturation (Yamaguchi 2010: ~10% pore space filled with air)
        theta_s = theta_s - 1e-6                     # numerical adjustment to avoid rho slightly above RHO_I
        theta_s[theta_s < 1e-6] = 1e-9               # zero porosity in ice layers can cause numerical problems

        theta_i = LWC_i / self.dz                    # irreducible water content as a volumetric fraction
        theta_i = np.minimum(theta_i, theta_s - 1e-9)  # limited to available porosity (numerical safety)
        theta_i[self.rho >= RhoImp] = 0.             # zero irreducible water content for ice lenses

        LWC_i = theta_i * self.dz                    # recompute LWC_i after the theta_s-based clamp above

        LWCav = np.maximum(0, self.LWC - LWC_i)      # LWC available for water flow
        LWCacm = self.dz * theta_s - self.LWC        # extra LWC that can be accommodated in each layer
        LWCacm[imp] = 0.                             # no water storage in impermeable ice lenses

        ### Distribute liquid water input at the surface into available accommodation space ###
        i_input_bottom = np.where(np.cumsum(LWCacm) >= liq_input)[0][0]   # deepest node receiving liquid input

        if i_input_bottom >= imptop[0]:
            # Input would reach an impermeable lens: cap distribution to the node just above it
            i_input_bottom = imptop[0] - 1

        if i_input_bottom == -1:
            # Impermeable lens sits right at the surface: all input runs off immediately
            runoff0 = 1 * liq_input
        elif i_input_bottom == 0:
            # Only the surface node receives input, limited by its accommodation space
            self.LWC[0] = self.LWC[0] + min(liq_input, LWCacm[0])
            runoff0 = liq_input - min(liq_input, LWCacm[0])   # input not accommodated by the surface node runs off
        elif i_input_bottom > 0:
            # Nodes [0:i_input_bottom] are filled completely; node i_input_bottom takes the remainder
            liq_input_remaining = liq_input - np.cumsum(LWCacm)[i_input_bottom - 1]
            self.LWC[0:i_input_bottom] = self.LWC[0:i_input_bottom] + LWCacm[0:i_input_bottom]
            self.LWC[i_input_bottom] = self.LWC[i_input_bottom] + min(liq_input_remaining, LWCacm[i_input_bottom])
            runoff0 = liq_input_remaining - min(liq_input_remaining, LWCacm[i_input_bottom])

        ### Recompute volumetric water content and accommodation space to reflect the updated LWC
        theta_w = self.LWC / self.dz
        LWCav = np.maximum(0, self.LWC - LWC_i)      # LWC available for water flow
        LWCacm = self.dz * theta_s - self.LWC        # remaining accommodation space per node
        LWCacm[imp] = 0.                              # no water storage in impermeable ice lenses
        
        ### Effective water saturation, Hirashima (2010) Eq. (5) -- rescales
        ### volumetric water content onto a 0-1 scale between the irreducible
        ### content (theta_i, "0") and saturated content (theta_s, "1").
        theta_e = (theta_w - theta_i) / (theta_s - theta_i)
        stab_e  = 1e-9                                  # stabilization bound for theta_e
        theta_e = np.maximum(stab_e, theta_e)            # avoid non-positive effective saturation
        theta_e = np.minimum(1 - stab_e, theta_e)        # avoid effective saturation equal to 1

        avG, nvG, mvG = vG_Yama_params(rg, rhodcy)       # van Genuchten shape parameters (alpha, n, m)
        bigk_s = hydrconducsat_Calonne(rg, rhodcy)       # hydraulic conductivity at saturation [m s-1]
        bigk_r = krel_vG(mvG, theta_e)                   # relative hydraulic conductivity
        bigk   = bigk_r * bigk_s                          # hydraulic conductivity, Hirashima (2010) Eq. (11) [m s-1]
        bigk_d = np.append(bigk[1:], 0)                  # hydraulic conductivity staggered down (next node's value)

        hd    = phead_vG(avG, nvG, mvG, theta_e)          # pressure head [m]
        dlthd = np.append(np.diff(hd), 0)                 # head difference between node and underlying node
        dhdz  = dlthd / dltz                              # d(head)/d(z); positive => absolute head increases downward

        ### Upwind/downwind conductivity selection at each interface, following
        ### Szymkiewicz (2009) Eq. (3c): the sign of the *total* head gradient
        ### (pressure head gradient + gravity, i.e. dhdz+1) determines whether
        ### flow is driven downward (use the upper node's conductivity, bigk)
        ### or upward against gravity (use the lower node's conductivity, bigk_d).
        ### This avoids the numerical instability of a naive central-conductivity
        ### scheme. Note: Szymkiewicz (2009) defines potential head as negative;
        ### what matters here is only the sign of q in their Eq. (1).
        bigkedg = np.zeros(ncv)
        bigkedg[dhdz + 1 >= 0] = bigk[dhdz + 1 >= 0]      # downward total-head gradient
        bigkedg[dhdz + 1 < 0]  = bigk_d[dhdz + 1 < 0]     # upward total-head gradient

        ### Determine flow at interfaces of volumes, following Hirashima et al. (2010) ###
        ## Detect all nodes that cannot have outflow ##
        indsdry = np.where(theta_e <= 1e-3)[0]   # essentially dry nodes: no mobile water, no outflow possible

        if np.any(theta_e > 1e-3):
            i0dry = np.where(theta_e > 1e-3)[0][-1] + 1   # depth below which all nodes are dry
        else:
            i0dry = 1   # only surface influx in the surface node

        imp_d = imp - 1   # nodes with an impermeable volume directly below them (also blocked from outflow)

        ## Find continuous saturated stacks directly above an impermeable boundary ##
        # A node sitting in a fully-saturated (theta_e >= 0.95) column immediately
        # above an ice lens has no available pore space to push water further into
        # via Darcy flux -- physically, this is a perched water table that has
        # already filled up to that point. Any additional water reaching these
        # nodes must be handled via ponding/lateral runoff elsewhere, not through
        # the normal flux solver below.
        satnofl = []   # indices of nodes in a saturated, no-outflow stack above a lens

        for lens_top_idx in imptop:                 # check above each impermeable lens
            if lens_top_idx <= i0dry:                # only relevant where water is percolating
                i_scan_up = lens_top_idx - 1         # start just above the lens
                reached_surface = (i_scan_up < 0)    # already at the surface?

                # Walk upward while each node remains fully saturated
                while (not reached_surface) and theta_e[i_scan_up] >= 0.95:
                    satnofl.append(i_scan_up)
                    i_scan_up -= 1
                    if i_scan_up < 0:
                        reached_surface = True   # surface reached; terminate the scan

        satnofl = np.array(satnofl)     # convert to array
        satnofl_d = satnofl - 1         # nodes directly above a saturated-no-outflow node (also blocked)

        # i00: all node indices that must have zero outflow --
        # dry nodes, impermeable nodes, nodes just above an impermeable barrier,
        # nodes just above a saturated-no-outflow stack, and the domain's bottom boundary.
        i00 = np.unique(np.concatenate((indsdry, imp, imp_d, satnofl_d, [ncv - 1]))).astype(int)

        # i11: candidate nodes (within the wet region, 0 to i0dry) that CAN have
        # outflow > 0 -- everything not excluded by i00. These are passed to the
        # iterative flux solver below.
        i11 = [ii for ii in range(i0dry) if ii not in i00]

        ### Iterative guesses to determine qlim of Hirashima et al. (2010) ###
        glw      = np.zeros(ncv - 1)      # guess of downward-transported LWC at each interface (qlim in Hirashima 2010)
        glwc     = np.copy(self.LWC)      # working LWC profile, updated as glw is resolved interface-by-interface
        glwcacm  = np.copy(LWCacm)        # working accommodation space, updated alongside glwc
        gtheta_e = np.copy(theta_e)       # working effective saturation, updated alongside glwc
        ghd      = np.zeros(ncv - 1)      # working pressure head, updated alongside glwc

        # Resolve interfaces bottom-up: solving the deepest interface first, then
        # updating glwc/glwcacm, ensures each shallower interface's calculation
        # already reflects the flux decisions made below it.
        for j1 in np.flip(i11):
            # qlim must equalize hd[j1] and hd[j1+1], assuming inflow from j1-1 has not yet occurred
            inds11 = np.array([j1, j1 + 1])   # indices on each side of the interface of interest

            # --- Initial guess for glw[j1] ---
            # Defensive abs() on the right-hand side guards against a slightly
            # negative flux arising from numerical round-off (e.g. glwcacm
            # marginally negative); fluxes are otherwise non-negative by design.
            if timer == 0 or abs(glwflux_d1[j1] - glwflux_d2[j1]) > 0.1 * abs(glwflux_d1[j1]):
                # Flux has changed significantly since the last sub-step (not in
                # steady state): guess via the analytic saturation-equalizer.
                glw[j1] = thetaeff_equaliser(theta_i[inds11], theta_s[inds11], glwc[inds11], self.dz[inds11])
                glw[j1] = max(glw[j1], 0)                                    # avoid negative guess values
                glw[j1] = min(glw[j1], min(LWCav[j1], glwcacm[j1 + 1]))      # cannot exceed available/accommodatable LWC
            else:
                # Flux is stable: warm-start from the previous sub-step's flux rate
                glw[j1] = dtsub * glwflux_d1[j1]
                glw[j1] = min(glw[j1], min(LWCav[j1], glwcacm[j1 + 1]))

            # --- Refine the guess and evaluate the equilibrium residual ---
            gtheta_e[inds11] = thetae_update(glw[j1], theta_i[inds11], theta_s[inds11], glwc[inds11], self.dz[inds11])
            ghd[inds11] = phead_vG(avG[inds11], nvG[inds11], mvG[inds11], gtheta_e[inds11])   # pressure head [m]
            f_eq = ghd[j1] - ghd[j1 + 1] - dltz[j1]   # Hirashima (2010) Eq. (20), evaluated at this interface

            if abs(f_eq) > 1.:
                # Guess is far from equilibrium: use the more robust bisection method
                glw[j1] = flux_bisection(glw[j1], LWCav, glwcacm, theta_i[inds11], theta_s[inds11],
                                          glwc[inds11], self.dz[inds11], avG[inds11], nvG[inds11],
                                          mvG[inds11], eps_cvg)
            else:
                # Guess is close to equilibrium: use the faster-converging Newton-Raphson method
                glw[j1] = flux_newtonraphson(glw[j1], LWCav, glwcacm, theta_i[inds11], theta_s[inds11],
                                              glwc[inds11], self.dz[inds11], avG[inds11], nvG[inds11],
                                              mvG[inds11], eps_cvg)

            # Update the working LWC profile and accommodation space to reflect
            # this interface's resolved flux, before moving to the interface above.
            glwc = self.LWC + np.append(0, glw) - np.append(glw, 0)
            glwcacm = self.dz * theta_s - glwc
            glwcacm[imp] = 0.   # no water storage in impermeable ice lenses

        glwflux_d2 = np.copy(glwflux_d1)   # shift: previous "last step" flux becomes "two steps ago"
        glwflux_d1 = glw / dtsub           # save this sub-step's resolved flux, for use as next step's warm-start
        
        ### Compute the water fluxes ###
        qlim = np.append(glw, 0)          # limiting flux at each interface, Hirashima (2010) Eq. (20) [m] (0 outflow from last node)
        q0   = bigkedg * (dhdz + 1)       # instantaneous Darcy flux under current conditions, Hirashima (2010) Eq. (1) [m s-1]
        ifl  = np.logical_and(qlim > 1e-20, q0 > 1e-20)   # interfaces with meaningful (non-zero) flow

        qstep = np.zeros(ncv)             # total flux transferred over this sub-step, per interface
        # Exponential relaxation toward the equilibrium-limiting flux qlim, following
        # Hirashima (2010) Eq. (23). This avoids overshooting equilibrium that a naive
        # q0*dtsub estimate could produce, without requiring an infinitesimally small
        # sub-step.
        qstep[ifl] = qlim[ifl] * (1 - np.exp(-dtsub * q0[ifl] / qlim[ifl]))   # [m we]
        qstep[-1] = 0                      # lower boundary condition: no outflow from the last node
        qstep = np.minimum(qstep, LWCav)   # safety clamp: flow cannot exceed available water

        lwcin  = np.append(0., qstep[0:-1])   # inflow to each node (no upper boundary flux into the surface node)
        lwcout = np.copy(qstep)               # outflow from each node [m we]

        if lat_runoff:
            ### Lateral runoff for water that is available but has no downward outflow
            ### path (sitting directly above an impermeable barrier or a saturated,
            ### no-outflow stack).
            ii_rf = np.concatenate((imp_d, satnofl_d)).astype(int)
            ii_rf = np.intersect1d(ii_rf, np.where(LWCav > 0)[0])

            if runoff_method == 'ZuoOerlemans':
                rfout = runoffZuoOerlemans(dtsub, slope_proxy, LWCav, ii_rf)
            elif runoff_method == 'Darcy':
                rfout = runoffDarcy(dtsub, slope_proxy, bigk, ii_rf)
            else:
                print(f"Unrecognized runoff_method_darcy: '{runoff_method}'. "
                      "Must be 'ZuoOerlemans' or 'Darcy'. Exiting.")
                sys.exit()

            rfout = np.minimum(rfout, LWCav - lwcout) # runoff cannot push total outflow beyond LWCav
            lwcout += rfout
            runoff1 = np.sum(rfout) # total lateral (non-surface) runoff this sub-step [m we]
        else:
            runoff1 = 0.        

        self.LWC = self.LWC + lwcin - lwcout

        if timer + dtsub == timetot:
            # At the very end of the Darcy routine: any LWC remaining in nodes whose
            # density was artificially capped (rhomod, Section 4) must run off, since
            # their true density does not actually support holding that water -- the
            # earlier permeability relaxation was a numerical convenience for flow
            # calculations only, not a physical storage allowance.
            runoff1 = runoff1 + np.sum(self.LWC[rhomod])
            self.LWC[rhomod] = 0.

        runofftot = runofftot + runoff0 + runoff1   # accumulate surface-input runoff + this sub-step's outflow runoff
        
        ### Refreezing ###
        cold_content = cp_i * self.mass * (T_MELT - self.Tz)          # cold content of nodes [J]
        refr_pot_ht  = cold_content / LF_I                             # refreezing potential from cold content [kg]
        refr_pot_v   = RHO_I_MGM * phi * self.dz * RHO_W_KGM - 1e-6    # refreezing potential from available volume [kg] (numerical safety margin)
        refr_pot     = np.minimum(refr_pot_ht, refr_pot_v)             # refreezing potential per node [kg]

        refr = np.minimum(self.LWC * RHO_W_KGM, refr_pot)              # refreezing per node [kg]
        refr[refr < 0] = 0                                              # avoid negative refreezing

        self.LWC  = np.maximum(0, self.LWC - refr / RHO_W_KGM)         # liquid mass loss [m we] (avoids numerical rounding errors)
        self.mass = self.mass + refr                                    # solid mass gain [kg]
        self.rho  = self.mass / self.dz                                 # update density [kg m-3]
        phi       = (RHO_I - self.rho) / RHO_I                          # update porosity
        phi       = np.maximum(0, phi)                                  # avoid numerical errors from very small negative phi

        latheat      = refr * LF_I                                      # latent heat released by refreezing [J]
        cold_content = cold_content - latheat                           # remaining cold content [J]
        self.Tz      = T_MELT - cold_content / (cp_i * self.mass)       # updated node temperatures [K]
        cp_i         = 152.5 + 7.122 * self.Tz                          # specific heat of ice [J kg-1 K-1] (Cuffey & Paterson 2010, Eq. 9.1)

        refr_tot = refr_tot + np.sum(refr) / RHO_W_KGM                  # accumulate total refreezing [m we]

        timer = timer + dtsub   # advance the sub-step timer

        ### Adapt the sub-step size ###
        if np.any(ifl):
            qratio = max(qstep[ifl] / qlim[ifl])   # largest ratio of actual to limiting flux, among flowing interfaces
        else:
            qratio = 0.   # no interface had outflow > 0

        if qratio < 0.25:
            dtsub = 1.25 * dtsub   # flow well below the limit: grow the sub-step
        if qratio > 0.75:
            dtsub = 0.75 * dtsub   # flow near/at the limit: shrink the sub-step to avoid overshoot

        dtsub = max(dtsub, dtmin)             # keep sub-step within configured bounds
        dtsub = min(dtsub, dtmax)
        dtsub = min(dtsub, timetot - timer)   # ensure the final sub-step lands exactly on timetot
        
    ### Sanity checks ###
    ### Water balance check: total liquid remaining + total refrozen + total
    ### runoff should equal the total input (melt + initial LWC + rain).
    water_balance_residual = (
        np.sum(self.LWC) + refr_tot + runofftot
        - (melt_mass_tot / RHO_W_KGM + np.sum(initial_lwc) + rain_vol_tot)
    )
    if abs(water_balance_residual) > 1e-12:
        print('Liquid water loss/gain, amount:', water_balance_residual)

    if np.any(self.Tz > T_MELT):
        print('Max Tz:', np.max(self.Tz))

    ### Check cold layers are dry ###
    coldlayers = np.where(self.Tz < T_MELT)[0]
    if np.all(self.LWC[coldlayers] < 1e-9):
        self.LWC[coldlayers] = 0.
    if np.any(self.LWC[coldlayers] > 0.):
        print('Problem: water content in a cold layer')

    if time.time() - ticdarcy >= 10:
        print(f'{iii} CFM time: {self.modeltime[iii]}')
        print(f'Darcy scheme run time: {np.around(time.time() - ticdarcy, 2)}')

    return self.rho, self.age, self.dz, self.Tz, self.r2, self.z, self.mass, \
        self.dzn, self.LWC, meltgridtrack, refr_tot, runofftot

#######################
### end darcyscheme ###
#######################