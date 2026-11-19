#!/usr/bin/env python

'''
vapor.py

Vapor-phase mass transport within the firn pore space, driven by
temperature-gradient-induced saturation vapor density gradients
(Arthern, Winebrenner, and Waddington 2000, Eqs. 1, 5, 6; flux form
adapted from Whillans and Grootes, 1985).

This module computes the net mass flux into each layer from vapor diffusion
[kg m-2 s-1]. Vapor transport changes the MASS of a layer at fixed volume (it
fills or empties existing pore space), so firn_density_nospin.py adds it to
self.mass and recomputes rho at fixed dz. It is NOT a densification term and
must not be added to drho_dt (which is a density change at fixed mass, i.e. a
volume change).

Because rho_v = rho_v(Tz) depends only on the (already-updated) temperature
field, this term requires no implicit solve of its own -- it is evaluated
directly, reusing the finite-volume face-conductance machinery in solver.py
(_fv_geometry), but applied to a KNOWN field rather than an unknown one being
solved for. The flux is computed once for each face and used by both layers
that share it, so the column total changes only by the fluxes at the top and
bottom (both zero for now): mass is conserved to rounding.

The diffusion of water vapor in firn is also what the isotope diffusion
(isotopeDiffusion.py, Johnsen et al. 2000) is, and the two share their
ingredients, which are defined here: the saturation vapor pressure
(saturation_vapor_pressure), the diffusivity of vapor in free air
(vapor_diffusivity) and the tortuosity of firn (firn_tortuosity_inverse; the
effective pore space is the porosity times 1/tortuosity). The isotope diffusion
does not include the vapor transport (the mass that the vapor moves carries no
isotopic composition yet); the model prints a warning if both are on.
'''

import numpy as np
from solver import _fv_geometry
from constants import *


M_WATER = 0.018   # molar mass of water [kg mol-1]


def saturation_vapor_pressure(Tz, atmosphere='air'):
    '''
    Saturation vapor pressure over ice [Pa]. Shared by the vapor transport and by the isotope diffusion
    (isotopeDiffusion.py), so that both use the same one.

    'air' (Earth): 3.454e12 exp(-6133 / T), as in Johnsen et al. (2000).
    'CO2' (Mars): Arthern, Winebrenner, and Waddington (2000), Eq. 6, 611 Pa at the triple point (273.16 K)
    and a constant latent heat of sublimation (2.83e6 J kg-1). The two agree to 0.35% for 223 to 273 K.

    :param Tz: temperature [K]
    :param atmosphere: 'air' or 'CO2'
    :return: p_sat [Pa]
    '''
    if atmosphere == 'air':
        return 3.454e12 * np.exp(-6133 / Tz)
    elif atmosphere == 'CO2':
        return 611.0 * np.exp((2.83e6 * M_WATER / R) * (1.0 / 273.16 - 1.0 / Tz))
    else:
        raise ValueError(f"Unknown atmosphere '{atmosphere}'; use 'air' or 'CO2'.")
    ### end saturation_vapor_pressure ###
    #####################################


def saturation_vapor_density(Tz, atmosphere='air'):
    '''
    Saturation vapor density over ice [kg m-3], from saturation_vapor_pressure (ideal gas).

    :param Tz: temperature [K]
    :param atmosphere: 'air' or 'CO2'
    '''
    return saturation_vapor_pressure(Tz, atmosphere) * M_WATER / (R * Tz)
    ### end saturation_vapor_density ###
    ####################################


def firn_tortuosity_inverse(rho, b=1.3):
    '''
    Inverse of the tortuosity of firn, 1/tau = 1 - b (rho/rho_i)^2 for densities below rho_i / sqrt(b), and
    0 above (Johnsen et al., 2000; b = 1.3 gives 804 kg m-3). The parameter was calibrated on observed
    isotope diffusion lengths. Used by the vapor transport and by the isotope diffusion, so that both have the
    same effective diffusivity of the vapor in the firn: D_air * porosity * (1/tau).

    :param rho: density [kg m-3]
    :param b: tortuosity parameter (configuration key vapor_tortuosity_b)
    '''
    rho = np.asarray(rho, dtype=float)
    return np.where(rho < RHO_I / np.sqrt(b), 1.0 - b * (rho / RHO_I) ** 2, 0.0)
    ### end firn_tortuosity_inverse ###
    ###################################


def vapor_diffusivity(Tz, P_atm, atmosphere='air', yosida_factor=None):
    '''
    Water vapor diffusivity within firn/snow pore space (before the porosity is applied).

    'air' (Earth): D_v = yosida_factor * 2.11e-5 * (T/273.15)**1.94 * (101325/P), the free-air value of
    isotopeDiffusion.py (Johnsen et al., 2000).
    'CO2' (Mars): D_v = yosida_factor * D_v_star * (T/273)**2 * (101300/P), derived below.

    Derivation of the CO2 form: Clow (1987) Eq. 10 gives the general gas-phase diffusivity
    form (Chapman and Cowling, 1952) for viscous flow; water vapor
    diffusivity is stated to follow "an expression analogous to Eq. (10)"
    with D_v_star, n, T_star, P_star as given below for a CO2 atmosphere.
    The yosida_factor (default 4.0) is Yosida's (1955) empirical finding
    that vapor diffuses ~4x faster within an actual snowpack than in the
    free gas, because ice grains act as vapor sources/sinks rather than
    passive obstacles -- this is a property of the snow/firn pore
    geometry, not of the background gas, so it is retained for both
    atmosphere options below.

    Verified: back-substituting this expression (with n=2) into Clow's
    Eqs. 12-13 (vapor heat flux, Clausius-Clapeyron) reproduces his Eq. 14
    for the vapor-diffusion thermal conductivity exactly, confirming the
    leading "4" there is this same Yosida factor and not an unrelated
    constant.

    AW00 (Arthern, Winebrenner, and Waddington, 2000) cite Clow (1987) as
    their source for D(T,P); this function is that source, made explicit.

    :param Tz: temperature [K]
    :param P_atm: total ambient gas pressure [Pa]
    :param atmosphere: 'air' (Earth, N2/O2; default) or 'CO2' (Mars)
    :param yosida_factor: snow-vs-free-gas enhancement (Yosida 1955). Default (None): 1 for 'air' (the free-air
        value, as isotopeDiffusion.py; no enhancement has been verified for Earth firn) and 4 for 'CO2'
        (Clow 1987 and AW00). Exposed for sensitivity tests (configuration key vapor_yosida_factor).

    :return: D_v [m2/s]
    '''
    if atmosphere == 'air':
        ### Earth: the free-air diffusivity of isotopeDiffusion.py (the same in both modules)
        yosida = 1.0 if yosida_factor is None else yosida_factor
        return yosida * D_V_STAR_AIR * (Tz / 273.15)**D_V_EXP_AIR * (101325.0 / P_atm)
    elif atmosphere == 'CO2':
        ### Mars: International Critical Tables (1929), via Clow (1987) Eq. 10
        yosida = 4.0 if yosida_factor is None else yosida_factor
        T_star = 273.0       # K
        P_star = 101300.0    # Pa (1013 mbar)
        n = 2.0              # Chapman and Cowling (1952)
        return yosida * D_V_STAR_CO2 * (Tz / T_star)**n * (P_star / P_atm)
    else:
        raise ValueError(f"Unknown atmosphere '{atmosphere}'; use 'air' or 'CO2'.")
    ### end vapor_diffusivity ###
    ##############################

def vapor_latent_conductivity(Tz, rho, config):
    '''
    Thermal conductivity that the vapor transport adds to the firn [W m-1 K-1]: the latent heat carried by the
    vapor, L * q, where the vapor flux q = -D * pore * d(rho_vs)/dz (vapor_mass_flux_divergence) follows the
    temperature gradient through rho_vs(T), so L * q = -K_v * dT/dz with

        K_v = L * D * pore * d(rho_vs)/dT,        d(rho_vs)/dT = rho_vs * (d ln p_sat/dT - 1/T)

    (D the diffusivity and pore the porosity times 1/tortuosity, as in the mass flux). The heat solvers add it to the
    conductivity (diffusion.py, when vapor_latent_heat and vaporTransport are on), which makes the latent heat
    implicit, and so stable for the thinnest layers; the same vapor that moves mass then moves its heat. It is
    evaluated per layer (the mass flux is evaluated on blocks, see vapor_mass_flux_divergence) and is not limited
    by the limiter of the mass flux. At 253 to 273 K it is 5 to 30% of the conductivity of firn at 150 kg m-3, 1 to
    6% at 300, 0.3 to 1.8% at 450, 0.1 to 0.4% at 600 and nothing from 804 kg m-3.

    The conductivity of firn that the model has by default (Calonne et al., 2019) is conduction through the ice and
    the air only (3-D image computations, no vapor), so adding K_v does not count the vapor twice there. Conductivities that
    are fits to measurements (Yen, Sturm, Van Dusen, ...) are apparent conductivities and may already include
    the vapor; the latent heat should not be switched on for them.

    :param Tz: temperature [K]
    :param rho: density [kg m-3]
    :param config: configuration (vapor_atmosphere, site_pressure, vapor_yosida_factor, vapor_tortuosity_b)
    :return: K_v [W m-1 K-1], same shape as Tz
    '''
    atmosphere = config.get('vapor_atmosphere', 'air')
    if atmosphere == 'air':
        L, dlnp_dT = LS_SUBLIMATION, 6133.0 / Tz**2
    else:
        L, dlnp_dT = 2.83e6, 2.83e6 * M_WATER / (R * Tz**2)   # CO2 (Mars): the constant latent heat of AW00, Eq. 6
    D = vapor_diffusivity(Tz, 100.0 * config.get('site_pressure', 1013.25), atmosphere=atmosphere,
                          yosida_factor=config.get('vapor_yosida_factor'))
    pore = np.maximum(1.0 - rho / RHO_I, 0.0) * firn_tortuosity_inverse(rho, config.get('vapor_tortuosity_b', 1.3))
    return L * D * pore * saturation_vapor_density(Tz, atmosphere) * (dlnp_dT - 1.0 / Tz)
    ### end vapor_latent_conductivity ###
    #####################################


def block_ids(dz, d_min):
    '''
    Group the layers into blocks of adjacent layers for the vapor flux: a block is at least d_min thick,
    except possibly the last one. Layers thicker than d_min are their own block; thin layers are grouped
    with their neighbours. Returns the block of each layer (0, 1, 2, ...; increasing downward).

    The layers are first put in bins of d_min by the depth of their top, and a block thinner than d_min/2 is
    added to the block below it (repeatedly), so that no thin block is left alone.

    :param dz: layer thickness [m]
    :param d_min: the thickness of a block [m]; 0 (or less) is one layer per block
    '''
    n = len(dz)
    if d_min <= 0 or n < 2:
        return np.arange(n)
    top = np.cumsum(dz) - dz
    bins = np.floor(top / d_min).astype(np.int64)
    ids = np.concatenate(([0], np.cumsum(np.diff(bins) > 0)))              # consecutive numbers: the layers of a bin are a block
    while True:
        thickness = np.bincount(ids, weights=dz)
        thin = np.flatnonzero(thickness[:-1] < 0.5 * d_min)                  # (the last block has no block below it)
        if len(thin) == 0:
            return ids
        below = np.arange(len(thickness))
        below[thin] = thin + 1                                               # a thin block goes with the one below it
        while True:                                                          # follow chains of thin blocks
            nxt = below[below]
            if np.array_equal(nxt, below):
                break
            below = nxt
        ids = np.concatenate(([0], np.cumsum(np.diff(below[ids]) > 0)))
    ### end block_ids ###
    #####################


def vapor_mass_flux_divergence(self, dt, max_loss=0.5, max_fill=0.5):
    '''
    Net mass flux into each layer from vapor diffusion, dmass_dt [kg m-2 s-1]
    (positive: the layer gains mass). The vapor flux is q = -D_eff * phi * (1/tau) * d(rho_v)/dz
    [kg m-2 s-1, positive downward] (Arthern et al. 2000, Eq. 5, with the effective pore space of the
    isotope diffusion: the porosity phi times the inverse tortuosity of Johnsen et al. 2000, which is
    zero from 804 kg m-3), and a layer gains the flux through its top face minus the flux through its
    bottom face.

    Reuses the same z_dummy/z_P finite-volume grid convention as heatDiff/refreezeDiff
    in diffusion.py, and the same face-conductance helper (_fv_geometry) used by every
    diffusion solver in solver.py -- but here Gamma_P = D_eff*phi/tau plays the role of
    a "conductivity" for a KNOWN field (rho_v(Tz)), not a state being solved for.

    The caller adds dt * dmass_dt to self.mass and recomputes the density at fixed layer
    thickness (see firn_density_nospin.py). The flux through the bottom is zero. At the top it is zero by default, which is
    consistent for Earth: the exchange of vapor with the atmosphere is the sublimation or deposition (the SUBLIM
    forcing, or the latent heat flux of the surface energy balance with a turbulent-flux method), which takes its mass
    from or gives it to the top layers (sublim.py), so a second flux through the surface would count it twice.

    THE SURFACE (configuration key vapor_surface_exchange). The pore space just beneath the surface is in contact with the
    atmosphere, and the vapor flux q(0,t) that enters the snowpack is part of the accumulation of AW00 (Arthern et al.,
    2000, Eq. 4: b_s = rho_s v(0,t) + q(0,t), the mass deposited on the surface plus the net vapor flux that enters; q(0)
    decides on Mars whether the pore space fills, their Figs. 9-10). AW00 take q(0) as the flux within the saturated snowpack
    immediately beneath the surface, from the temperature gradient there (Sec. 4.2). With the key true the top face has that flux:
    the flux of the top block with the gradient of the vapor density between the top two blocks,
    q(0) = -D pore_0 (rho_vs(T_1) - rho_vs(T_0)) / (z_P1 - z_P0), downward positive (the first interior face uses the same
    gradient with the mean pore space of the two blocks); it is limited like the other faces. (A flux between the top block and a
    surface saturated at the surface temperature would be zero: the top layer is at the surface temperature.) This mass comes
    from the atmosphere, so the mass of the column changes by dt q(0), and the accumulation of the model must be the mass deposited
    on the surface (b_s minus q(0)), not b_s. Do not combine it with sublimation or deposition of the same vapor (SUBLIM or the SEB
    turbulent fluxes).


    THE BLOCKS. This is a single explicit step, and the flux depends on the temperature gradient and not on
    the mass of a layer, so the relative loss of a layer grows as 1/thickness^2: with daily steps the
    near-surface grid has layers of 0.03 mm and 0.004 kg m-2 (days with almost no snowfall), which are
    smaller than a grain and which the flux empties in a step. A vapor molecule goes about sqrt(D dt) = 4 cm in a
    day, so the flux is calculated on blocks: adjacent layers are grouped (block_ids) in blocks of at
    least vapor_min_thickness (configuration, default 0.01 m), with the mean properties of their layers
    (temperature and density weighted by mass, the pore space by thickness). The flux goes between the blocks, and
    the mass change of a block is shared by its layers: in proportion to their mass for a loss, to their free
    pore space (rho_i * dz - mass) for a gain. Every layer of a block therefore loses the same fraction of
    its mass, and the model layers do not change (the merge function of merge.py changes them). Layers that are
    at least vapor_min_thickness thick are their own block, so this is the flux between layers for them.

    THE LIMITER. The flux between blocks is also limited, each face flux being scaled by its donor block
    (the block it leaves) so that no block loses more than max_loss of its mass in the step, and by its
    receiving block so that no block gains more than max_fill of its free pore space. Every face flux is
    scaled once, so the column sum stays zero (mass is conserved), the mass stays positive and the density
    stays below the ice density.

    :param self: FirnDensityNoSpin instance (needs self.z, self.dz, self.Tz, self.rho, self.mass and the configuration self.c,
        which may have vapor_atmosphere ('air', default, or 'CO2'), site_pressure [hPa, default 1013.25],
        vapor_yosida_factor (see vapor_diffusivity), vapor_tortuosity_b (see firn_tortuosity_inverse) and
        vapor_min_thickness [m, default 0.01; 0 is one layer per block])
    :param dt: time step [s]
    :param max_loss: the largest fraction of the mass that a block can lose in the step
    :param max_fill: the largest fraction of the free pore space that a block can fill in the step

    :return: dmass_dt [kg m-2 s-1], same shape as self.rho; sums to zero
    :return: number of blocks whose outflow or inflow was limited
    '''
    z, dz, T, rho, mass = self.z, self.dz, self.Tz, self.rho, self.mass
    config = getattr(self, 'c', {})
    atmosphere = config.get('vapor_atmosphere', 'air')

    ids = block_ids(dz, config.get('vapor_min_thickness', 0.01))
    nb = ids[-1] + 1
    if nb < 2:
        return np.zeros(len(dz)), 0
    total = lambda x: np.bincount(ids, weights=x, minlength=nb)   # the sum over the layers of each block
    safe = lambda x: np.where(x > 0, x, 1.0)

    room = np.maximum(RHO_I * dz - mass, 0.0)                     # free pore space of a layer, as mass
    H, M, P = total(dz), total(mass), total(room)
    T_b = np.where(M > 0, total(T * mass) / safe(M), total(T) / np.bincount(ids, minlength=nb))
    ### effective pore space for the vapor, as in the isotope diffusion: porosity times 1/tortuosity (Johnsen et al., 2000)
    pore = np.maximum(1.0 - rho / RHO_I, 0.0) * firn_tortuosity_inverse(rho, config.get('vapor_tortuosity_b', 1.3))
    pore_b = np.where(H > 0, total(pore * dz) / safe(H), 0.0)
    P_atm = 100.0 * config.get('site_pressure', 1013.25)           # hPa to Pa, as isotopeDiffusion.py
    Gamma_P = vapor_diffusivity(T_b, P_atm, atmosphere=atmosphere, yosida_factor=config.get('vapor_yosida_factor')) * pore_b
    rho_v = saturation_vapor_density(T_b, atmosphere)

    ### The grid of the blocks, as the convention of heatDiff/refreezeDiff in diffusion.py
    first = np.concatenate(([0], np.flatnonzero(np.diff(ids)) + 1))   # the first layer of each block
    z_edges = np.concatenate((z[first], [z[-1] + dz[-1]]))
    z_P = (z_edges[1:] + z_edges[:-1]) / 2
    dZ, deltaZ_u, deltaZ_d, Gamma_u, Gamma_d = _fv_geometry(z_edges, z_P, Gamma_P)

    ### Flux through the top face of each block [kg m-2 s-1, positive downward]: face i is
    ### the top of block i (face 0 is the surface, face nb the bottom of the column). Each
    ### interior face is computed once, from the block above it (its lower face).
    q_face = np.zeros(nb + 1)
    q_face[1:-1] = -(Gamma_d[:-1] / deltaZ_d[:-1]) * (rho_v[1:] - rho_v[:-1])
    ### q_face[-1] = 0: no flux through the bottom. The surface flux is zero too unless vapor_surface_exchange (see docstring)
    if config.get('vapor_surface_exchange', False):
        q_face[0] = -Gamma_P[0] * (rho_v[1] - rho_v[0]) / (z_P[1] - z_P[0])   # AW00: the flux just beneath the surface, from the gradient between the top two blocks

    ### Limit the outflow of every block (the flux through the bottom face leaves it if it is downward, and the flux
    ### through the top face leaves it if it is upward) and the inflow (the opposite)
    outflow = np.maximum(q_face[1:], 0.0) + np.maximum(-q_face[:-1], 0.0)
    inflow = np.maximum(q_face[:-1], 0.0) + np.maximum(-q_face[1:], 0.0)
    allowed_out = max_loss * M / dt
    allowed_in = max_fill * P / dt
    out_limited, in_limited = outflow > allowed_out, inflow > allowed_in
    scale_out, scale_in = np.ones_like(outflow), np.ones_like(inflow)
    scale_out[out_limited] = allowed_out[out_limited] / outflow[out_limited]
    scale_in[in_limited] = allowed_in[in_limited] / inflow[in_limited]
    face_scale = np.ones_like(q_face)   # a downward flux leaves the block above and enters the block below; an upward one the reverse
    face_scale[1:-1] = np.where(q_face[1:-1] > 0, np.minimum(scale_out[:-1], scale_in[1:]),
                                np.minimum(scale_out[1:], scale_in[:-1]))
    face_scale[0] = scale_in[0] if q_face[0] > 0 else scale_out[0]   # the surface flux enters the top block (downward) or leaves it
    q_face = q_face * face_scale
    change = q_face[:-1] - q_face[1:]   # net mass flux into each block [kg m-2 s-1]

    ### Share the change of each block among its layers: by mass for a loss, by free pore space for a gain
    loss = np.where((change < 0) & (M > 0), change / safe(M), 0.0)
    gain = np.where((change > 0) & (P > 0), change / safe(P), 0.0)
    return loss[ids] * mass + gain[ids] * room, int(np.sum(out_limited | in_limited))
    ### end vapor_mass_flux_divergence ###
    ######################################
