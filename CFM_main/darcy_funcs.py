# -*- coding: utf-8 -*-
"""
Functions required for the Darcy-type liquid water flow scheme

@author: Vincent
"""

import numpy as np
from constants import *

def hydrconducsat_Calonne(rad, rho):
    '''
    Saturated hydraulic conductivity, following Calonne et al., Eq. (6).

    Parameters
    ----------
    rad : ndarray
        Grain radius [m].
    rho : ndarray
        Node density [kg m-3].

    Returns
    -------
    bigksat : ndarray
        Saturated hydraulic conductivity [m s-1].
    '''
    mu = 0.001792   # dynamic viscosity of water at 273.15 K [kg m-1 s-1]
    bigksat = 3 * (rad)**2 * RHO_W_KGM * GRAVITY / mu * np.exp(-0.013 * rho)   # [m s-1]
    return bigksat

def vG_Yama(rad,rho,thetaeff):
    '''
    Pressure head and relative hydraulic conductivity computations from the
    van Genuchten (1980) model with the Yamaguchi et al. (2012)
    parameterisation.

    NOTE: This combines the same calculations performed separately by
    `vG_Yama_params` (van Genuchten shape parameters alpha, n, m),
    `phead_vG` (pressure head), and `krel_vG` (relative hydraulic
    conductivity) elsewhere in this file. `darcyscheme()` in melt.py calls
    those three functions individually rather than this combined version.
    This function is not currently called anywhere in melt.py, but is kept
    here in case it is used elsewhere in the broader CFM codebase -- verify
    before removing.

    Citation note: the pressure-head formula below is presented, in
    equivalent form, in both Hirashima et al. (2010), Eq. 9, and Hirashima
    et al. (2014), Eq. 3 (https://doi.org/10.1016/j.coldregions.2014.09.004).

    Parameters
    ----------
    rad : ndarray
        Grain radius [m].
    rho : ndarray
        Node density [kg m-3].
    thetaeff : ndarray
        Effective water saturation [-], 0-1.

    Returns
    -------
    head : ndarray
        Pressure head [m] (Hirashima et al. 2010, Eq. 9 / 2014, Eq. 3).
    bigkrel : ndarray
        Relative hydraulic conductivity [-] (Hirashima et al. 2010, Eq. 10).
    '''
    alpha = 4.4e6*(rho/(2*rad))**(-0.98) #.alpha parameter, Yamaguchi 2012 Eq.(6)
    n     = 1+2.7e-3*(rho/(2*rad))**0.61 #n parameter, Yamaguchi 2012 Eq. (7)
    m     = 1-(1/n) #m parameter, Yamaguchi 2012 p.7
    head  = 1/alpha * (thetaeff**(-1/m)-1)**(1/n) #head pressure, Hirashima 2010 (9) / 2014 (3)
    bigkrel = thetaeff**(1/2) * (1-(1-thetaeff**(1/m))**m)**2 # Hirashima 2010 (10)
    return(head,bigkrel)

def thetae_update(absfl, th_i, th_s, LWC, dz):
    '''
    Computes the effective water saturation that would result at two
    neighboring nodes if a given flux were transferred between them.

    Used by the iterative flux solvers (`flux_bisection`,
    `flux_newtonraphson`) to evaluate candidate flux guesses without
    committing them to the model state -- the upper node (index 0) is
    treated as losing `absfl` to outflow, and the lower node (index 1) as
    gaining `absfl` as inflow.

    Parameters
    ----------
    absfl : float
        Candidate water flux across the interface between the two nodes [m].
    th_i : ndarray, shape (2,)
        Irreducible volumetric water content at the two nodes [-].
    th_s : ndarray, shape (2,)
        Volumetric water content at saturation at the two nodes [-].
    LWC : ndarray, shape (2,)
        Current liquid water content at the two nodes [m].
    dz : ndarray, shape (2,)
        Thickness of the two nodes [m].

    Returns
    -------
    th_e : ndarray, shape (2,)
        Resulting effective water saturation at the two nodes [-], 0-1,
        clamped away from the exact bounds for numerical stability
        (Hirashima et al. 2010, Eq. 5).
    '''
    lw_in  = np.append(0, absfl)      # upper node: no inflow from this interface; lower node: gains absfl
    lw_out = np.append(absfl, 0)      # upper node: loses absfl; lower node: no outflow from this interface
    th_w   = (LWC + lw_in - lw_out) / dz   # resulting volumetric water content

    th_e = (th_w - th_i) / (th_s - th_i)   # effective water saturation, Hirashima et al. (2010), Eq. 5

    stab_e = 1e-9                           # stabilization bound
    th_e = np.maximum(stab_e, th_e)         # avoid non-positive effective saturation
    th_e = np.minimum(1 - stab_e, th_e)     # avoid effective saturation equal to 1

    return th_e

def thetaeff_equaliser(th_i2, th_s2, LWC2, dz2):
    '''
    Computes the flux, from node 0 to node 1, that would bring both nodes
    to equal effective saturation.

    This is a cheap closed-form estimate used as an initial guess for the
    iterative flux solvers (`flux_bisection`, `flux_newtonraphson`) when
    the flux has not been stable between sub-steps (see `darcyscheme()`).
    It is not the same equilibrium condition the solvers ultimately target
    (which additionally accounts for the gravitational head term, Hirashima
    et al. 2010 Eq. 20), but serves as a reasonable starting point for
    refinement.

    Derivation: setting th_e_after[0] == th_e_after[1], where
    th_e_after[0] = th_e0[0] - q / (dz[0]*(th_s[0]-th_i[0])) and
    th_e_after[1] = th_e0[1] + q / (dz[1]*(th_s[1]-th_i[1])), and solving
    for q, gives the expression implemented below.

    Parameters
    ----------
    th_i2 : ndarray, shape (2,)
        Irreducible volumetric water content at the two nodes [-].
    th_s2 : ndarray, shape (2,)
        Volumetric water content at saturation at the two nodes [-].
    LWC2 : ndarray, shape (2,)
        Current liquid water content at the two nodes [m].
    dz2 : ndarray, shape (2,)
        Thickness of the two nodes [m].

    Returns
    -------
    lwflux : float
        Flux from node 0 to node 1 that equalizes effective saturation
        between the two nodes [m].
    '''
    th_w  = LWC2 / dz2
    th_e0 = (th_w - th_i2) / (th_s2 - th_i2)   # effective water saturation, Hirashima et al. (2010), Eq. 5

    # Flux from index[0] to index[1] that ensures equal saturation between both nodes
    lwflux = ((dz2[0] * (th_s2[0] - th_i2[0]))**(-1)
              + (dz2[1] * (th_s2[1] - th_i2[1]))**(-1))**(-1) * (th_e0[0] - th_e0[1])
    return lwflux

def vG_Yama_params(rad, rho):
    '''
    Computes the van Genuchten (1980) shape parameters, following the
    grain-size- and density-dependent parameterization of Yamaguchi et al.
    (2012).

    Parameters
    ----------
    rad : ndarray
        Grain radius [m].
    rho : ndarray
        Node density [kg m-3].

    Returns
    -------
    alpha : ndarray
        Van Genuchten alpha parameter (Yamaguchi et al. 2012, Eq. 6).
    n : ndarray
        Van Genuchten n parameter (Yamaguchi et al. 2012, Eq. 7).
    m : ndarray
        Van Genuchten m parameter, m = 1 - 1/n (Yamaguchi et al. 2012, p.7).
    '''
    alpha = 4.4e6 * (rho / (2 * rad))**(-0.98)   # Yamaguchi et al. (2012), Eq. 6
    n     = 1 + 2.7e-3 * (rho / (2 * rad))**0.61  # Yamaguchi et al. (2012), Eq. 7
    m     = 1 - (1 / n)                            # Yamaguchi et al. (2012), p. 7
    return alpha, n, m

def phead_vG(alpha, n, m, thetaeff):
    '''
    Computes pressure (matric) head from van Genuchten shape parameters
    and effective water saturation.

    Citation note: this formula is presented, in equivalent form, in both
    Hirashima et al. (2010), Eq. 9, and Hirashima et al. (2014), Eq. 3
    (https://doi.org/10.1016/j.coldregions.2014.09.004).

    Parameters
    ----------
    alpha, n, m : ndarray
        Van Genuchten shape parameters (see `vG_Yama_params`).
    thetaeff : ndarray
        Effective water saturation [-], 0-1.

    Returns
    -------
    head : ndarray
        Pressure head [m].
    '''
    head = 1 / alpha * (thetaeff**(-1 / m) - 1)**(1 / n)   # Hirashima et al. (2014), Eq. 3
    return head

def krel_vG(m, thetaeff):
    '''
    Computes relative hydraulic conductivity from the van Genuchten m
    parameter and effective water saturation (van Genuchten-Mualem model).

    Parameters
    ----------
    m : ndarray
        Van Genuchten m parameter (see `vG_Yama_params`).
    thetaeff : ndarray
        Effective water saturation [-], 0-1.

    Returns
    -------
    bigkrel : ndarray
        Relative hydraulic conductivity [-], 0-1 (Hirashima et al. 2010,
        Eq. 10).
    '''
    bigkrel = thetaeff**(1 / 2) * (1 - (1 - thetaeff**(1 / m))**m)**2   # Hirashima et al. (2010), Eq. 10
    return bigkrel

def dfdg_derivative(th_sfull, th_ifull, th_efull, alphafull, nfull, mfull, dzfull):
    '''
    Computes the derivative of the equilibrium residual (f_eq) with respect
    to the water flux at an interface between two neighboring nodes.

    f_eq is defined by moving all terms of Eq. (20), Hirashima et al.
    (2010), to the right-hand side: f_eq = hd[upper] - hd[lower] - dltz,
    where hd is pressure head (`phead_vG`). Since increasing flux decreases
    effective saturation (and thus pressure head) at the upper node and
    increases it at the lower node, this derivative combines the head
    sensitivity to saturation at both nodes via the chain rule.

    Used by `flux_newtonraphson` to compute its Newton-Raphson update step.

    Implementation note: this function's parameters are suffixed "full" and
    are internally sliced into "upper" ([0:-1]) and "lower" ([1:]) subsets,
    suggesting it may have been designed to support computing this
    derivative across an entire column of interfaces at once. However, in
    the current codebase it is only ever called from `flux_newtonraphson`
    with inputs already reduced to a single interface's two nodes (i.e.
    2-element arrays), in which case the "upper"/"lower" slices each reduce
    to single-element arrays. The function is correct either way, but the
    apparent multi-interface vectorization capability is not currently
    exercised.

    Parameters
    ----------
    th_sfull, th_ifull, th_efull : ndarray
        Saturated, irreducible, and effective volumetric water content.
    alphafull, nfull, mfull : ndarray
        Van Genuchten shape parameters (see `vG_Yama_params`).
    dzfull : ndarray
        Node thickness [m].

    Returns
    -------
    dfdg : ndarray
        Derivative of f_eq with respect to flux at each interface.
    '''
    th_s, th_i, th_e = th_sfull[0:-1], th_ifull[0:-1], th_efull[0:-1]
    alpha, n, m, dz  = alphafull[0:-1], nfull[0:-1], mfull[0:-1], dzfull[0:-1]
    th_sd, th_id, th_ed = th_sfull[1:], th_ifull[1:], th_efull[1:]
    alphad, nd, md, dzd = alphafull[1:], nfull[1:], mfull[1:], dzfull[1:]

    dfdg = (
        1 / ((th_s - th_i) * alpha * n * m * dz) * th_e**(-1 * (1 + 1 / m)) * (th_e**(-1 / m) - 1)**((1 - n) / n)
        + 1 / ((th_sd - th_id) * alphad * nd * md * dzd) * th_ed**(-1 * (1 + 1 / md)) * (th_ed**(-1 / md) - 1)**((1 - nd) / nd)
    )
    return dfdg


def flux_bisection(gc, LWCav, glwcacm, th_i, th_s, lwc, dz, avG, nvG, mvG, eps_cvg):
    '''
    Bisection root-finder for the water flux (gc) between two neighboring
    nodes that brings their pressure heads into equilibrium.

    The equilibrium residual f_eq is defined by moving all terms of
    Eq. (20), Hirashima et al. (2010), to the right-hand side:
    f_eq = hd[0] - hd[1] - dltz. Bisection converges either on this
    numerical criterion, or on one of three physically-motivated early-exit
    conditions: the upper node has run dry, the lower node has reached
    saturation, or the guess has reached zero flux (the minimum physically
    allowed in this scheme).

    Called from `flux_newtonraphson` as a fallback when Newton-Raphson
    fails to converge or diverges.

    Parameters
    ----------
    gc : float
        Initial guess for the flux across the interface [m].
    LWCav : ndarray
        LWC available for flow, full column [m] (only LWCav[0], the upper
        node, is used to bound the guess).
    glwcacm : ndarray
        Remaining accommodation space, full column [m] (only glwcacm[1],
        the lower node, is used to bound the guess).
    th_i, th_s, lwc, dz : ndarray, shape (2,)
        Irreducible/saturated volumetric water content, current LWC, and
        thickness at the two nodes.
    avG, nvG, mvG : ndarray, shape (2,)
        Van Genuchten shape parameters at the two nodes.
    eps_cvg : float
        Convergence tolerance on the equilibrium residual [m].

    Returns
    -------
    gc : float
        Converged (or early-exited) flux guess [m].
    '''
    bisitmax = 100      # maximum bisection iterations
    bisit    = 0
    cvg_bis  = False
    dltz     = 1 / 2 * sum(dz)   # distance between the centres of the two nodes

    gth_e = thetae_update(gc, th_i, th_s, lwc, dz)
    ghd   = phead_vG(avG, nvG, mvG, gth_e)         # pressure head [m]
    f_eq  = ghd[0] - ghd[1] - dltz                  # Hirashima et al. (2010), Eq. 20, evaluated at this interface

    # Initialize bisection bounds
    g0 = 0.                                # lower bound
    g1 = min(LWCav[0], glwcacm[1])         # upper bound (limited by available LWC and accommodation space)

    while (not cvg_bis) and (bisit < bisitmax):
        gprev0 = np.copy(gc)   # flux guess at the previous iteration

        if f_eq < 0:
            # hd[0] too low relative to equilibrium -> increase outgoing flux
            g0 = np.copy(gc)
            gc = (g1 + g0) / 2
        elif f_eq > 0:
            # hd[0] too high relative to equilibrium -> decrease outgoing flux
            g1 = np.copy(gc)
            gc = (g1 + g0) / 2

        gth_e = thetae_update(gc, th_i, th_s, lwc, dz)
        ghd   = phead_vG(avG, nvG, mvG, gth_e)
        f_eq  = ghd[0] - ghd[1] - dltz

        if f_eq < 0 and gth_e[0] < 1e-8:
            # Equilibrium would require more outflow, but the upper node has already run dry
            cvg_bis = True
        elif f_eq < 0 and gth_e[1] > 0.95:
            # Equilibrium would require more outflow, but the lower node is already saturated
            cvg_bis = True
        elif f_eq > 0 and gc <= 1e-6:
            # Equilibrium would require less outflow, but the guess is already near the minimum (zero) flux
            gc = 0.
            cvg_bis = True

        if abs(f_eq) < eps_cvg or abs(gc - gprev0) < 1e-6:
            cvg_bis = True   # flux estimate has converged

        bisit += 1
        if bisit == bisitmax:
            print('Maximum iteration number reached in bisection algorithm')

    return gc

def flux_newtonraphson(gc, LWCav, glwcacm, th_i, th_s, lwc, dz, avG, nvG, mvG, eps_cvg):
    '''
    Newton-Raphson root-finder for the water flux (gc) between two
    neighboring nodes that brings their pressure heads into equilibrium.

    Uses the same equilibrium residual f_eq as `flux_bisection` (Eq. 20,
    Hirashima et al. 2010) and the same three physically-motivated
    early-exit conditions. If a Newton-Raphson step makes the residual
    worse, or the analytic derivative becomes very large (risking an
    unstable step size), the algorithm falls back to `flux_bisection`
    using the pre-step guess as its starting point.

    Parameters
    ----------
    gc : float
        Initial guess for the flux across the interface [m].
    LWCav, glwcacm : ndarray
        Full-column available LWC and remaining accommodation space [m]
        (passed through to `flux_bisection` if a fallback is triggered).
    th_i, th_s, lwc, dz : ndarray, shape (2,)
        Irreducible/saturated volumetric water content, current LWC, and
        thickness at the two nodes.
    avG, nvG, mvG : ndarray, shape (2,)
        Van Genuchten shape parameters at the two nodes.
    eps_cvg : float
        Convergence tolerance on the equilibrium residual [m].

    Returns
    -------
    gc : float
        Converged (or early-exited, or bisection-fallback) flux guess [m].
    '''
    nritmax = 20      # maximum Newton-Raphson iterations
    nrit    = 0
    cvg_nr  = False
    dltz    = 1 / 2 * sum(dz)   # distance between the centres of the two nodes

    gth_e = thetae_update(gc, th_i, th_s, lwc, dz)
    ghd   = phead_vG(avG, nvG, mvG, gth_e)
    f_eq  = ghd[0] - ghd[1] - dltz   # Hirashima et al. (2010), Eq. 20, evaluated at this interface

    while (not cvg_nr) and (nrit < nritmax):
        gprev0 = np.copy(gc)     # flux guess at the previous iteration
        fprev0 = np.copy(f_eq)   # residual at the previous iteration

        if f_eq < 0 and gth_e[0] < 1e-8:
            # Equilibrium would require more outflow, but the upper node has already run dry
            cvg_nr = True
        elif f_eq < 0 and gth_e[1] > 0.95:
            # Equilibrium would require more outflow, but the lower node is already saturated
            cvg_nr = True
        elif f_eq > 0 and gc <= 1e-6:
            # Equilibrium would require less outflow, but the guess is already near the minimum (zero) flux
            gc = 0.
            cvg_nr = True
        else:
            # Newton-Raphson step
            dfdg = dfdg_derivative(th_s, th_i, gth_e, avG, nvG, mvG, dz)
            deltaglw = -1 * fprev0 / dfdg
            gc = gprev0 + deltaglw

            gth_e = thetae_update(gc, th_i, th_s, lwc, dz)
            ghd   = phead_vG(avG, nvG, mvG, gth_e)
            f_eq  = ghd[0] - ghd[1] - dltz

            if abs(f_eq) > abs(fprev0) or abs(dfdg) > 1e6:
                # Newton-Raphson diverged or the step size risks instability: fall back to bisection
                gc = flux_bisection(gprev0, LWCav, glwcacm, th_i, th_s, lwc, dz, avG, nvG, mvG, eps_cvg)
                cvg_nr = True

        nrit += 1
        if abs(f_eq) < eps_cvg or abs(gc - gprev0) < 1e-6:
            cvg_nr = True   # flux estimate has converged

        if nrit == nritmax:
            print('Maximum iteration number reached in Newton-Raphson algorithm')

    return gc

def runoffZuoOerlemans(dt, slope, lwcexcess, inds):
    '''
    Computes lateral runoff following the Zuo & Oerlemans (1996)
    parameterization.

    NOTE: an identical formula is currently duplicated inline within
    bucket() in melt.py, rather than calling this shared function.
    Follow-up refactor opportunity: have bucket() call this function
    instead, to avoid maintaining the same formula in two places.

    See melt.py's bucket() documentation for a caveat regarding the
    dimensional interpretation of `slope` (best-supported inference:
    dimensionless rise/run, not explicitly defined in the source paper),
    and a note that this is a forward-Euler discretization of the loss
    term in Zuo & Oerlemans (1996) Eq. 21, not a direct implementation of
    the full equation (which also includes a production term not relevant
    here).

    Parameters
    ----------
    dt : float
        Time step duration [s].
    slope : float
        Surface slope [m/m], dimensionless (see note above).
    lwcexcess : ndarray
        Liquid water content in excess of irreducible water content [m].
    inds : ndarray
        Indices of the firn column where runoff should be computed.

    Returns
    -------
    rfout : ndarray
        Runoff at each node [m] (zero outside `inds`).
    '''
    c1zuo = 1.5 * 24 * 3600    # Zuo & Oerlemans (1996) constant [s]
    c2zuo = 25. * 24 * 3600    # Zuo & Oerlemans (1996) constant [s]
    c3zuo = 140.                 # Zuo & Oerlemans (1996) constant [-]
    tstar = c1zuo + c2zuo * np.exp(-1 * c3zuo * slope)   # Eq. 22, Zuo & Oerlemans (1996) [s]

    rfout = np.zeros(len(lwcexcess))
    rfout[inds] = dt * lwcexcess[inds] / tstar   # Euler discretization of the loss term, Eq. 21 [m]
    return rfout


def runoffDarcy(dt, slope, bigk, inds):
    '''
    Computes lateral runoff via a simplified Darcy's-law formulation,
    treating the local surface slope as a proxy driving gradient (since
    this is a single-column model with no neighboring column against which
    to compute an actual lateral head gradient), and assuming horizontally
    homogeneous hydraulic properties.

    This is the lateral runoff function actually used by darcyscheme() in
    melt.py (runoffZuoOerlemans is available as an alternative but is
    currently commented out there).

    Parameters
    ----------
    dt : float
        Time step duration [s].
    slope : float
        Surface slope [m/m], dimensionless, used here as a proxy driving
        gradient for lateral flow.
    bigk : ndarray
        Hydraulic conductivity [m s-1].
    inds : ndarray
        Indices of the firn column where runoff should be computed.

    Returns
    -------
    rfout : ndarray
        Runoff at each node [m] (zero outside `inds`).
    '''
    rfout = np.zeros(len(bigk))
    rfout[inds] = dt * bigk[inds] * slope   # total outgoing lateral Darcy flux
    return rfout





