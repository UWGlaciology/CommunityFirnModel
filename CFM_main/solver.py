#!/usr/bin/env python
'''
solver.py

Numerical solvers for the diffusion equation used by the CFM, covering
standard (no phase change) diffusion (heat, isotopes, gas) and
refreezing/melt scenarios where liquid water is present and latent
heat must be accounted for.

All solvers use a 1-D finite volume (FV) discretization following
Patankar (1980), Numerical Heat Transfer and Fluid Flow. Layers (control
volumes) are defined by z_edges (boundaries) and Z_P (centers).

Conventions:
    - Temperature: functions dealing with liquid water (transient_solve_enthalpy,
      transient_solve_ahc, transient_solve_decp) use degrees C, with the
      fusion/freezing point at T=0. transient_solve_TR uses Kelvin.
    - th_solid, th_liquid: MASS per total volume [kg/m3] (not volume fraction).
    - Gamma_P: thermal conductivity [W/m/K].
    - dt: time step [s]; dZ: layer thickness [m].

Liquid water content is INDEPENDENT STATE, not a function of temperature:
    This deserves emphasis because two solvers in this module have been
    broken by getting it wrong, in the same way, for the same reason.

    The published methods these solvers draw on -- apparent heat capacity
    (Dall'Amico 2010), and the NCZ enthalpy method (Tubini et al. 2021) --
    were developed for frozen soil and for pure water. In both of those
    settings, liquid content is a single-valued function of temperature: in
    soil via the soil freezing characteristic curve (SFCC), and in pure
    water via a narrow linearized ramp at the fusion point. The enthalpy
    function h(T) is then well defined, and the liquid/solid partition can
    be recovered from temperature alone.

    Firn does not work this way. Liquid water content is set by surface melt
    and percolation (see melt.py), not by the local temperature. A firn
    layer at T = 0 may hold anywhere between zero and its full mass in
    liquid, and two layers at identical temperature may have entirely
    different LWC. There is no single-valued h(T) that represents this.

    The practical consequence, when this is overlooked: any recovery of the
    form th_liquid = m * f(T), with m the TOTAL H2O mass and f -> 1 as
    T -> 0, assigns nearly all of a layer's mass to liquid whenever that
    layer sits microdegrees below fusion -- which is exactly where wet firn
    sits. th_solid then collapses toward zero, and since diffusion.py's
    refreezeDiff sets self.rho = th_solid and self.mass = rho * dz, the
    error propagates as division by ~zero in melt.py (rho = mass/dz,
    cold_content/(CP_I*mass)), then as NaN into SEB.py, which fails with an
    opaque "truth value of an empty array is ambiguous" from an empty
    root-finding result. The proximate crash is several modules away from
    the cause.

    The fix, applied in both transient_solve_ncz and transient_solve_ahc
    (2026-08): scale the latent-heat term and the liquid/solid partition by
    th_liquid_old -- the water actually present at the start of the step --
    rather than by total mass m. A dry layer then has no latent-heat
    contribution and cannot gain liquid through diffusion, and refreezing
    is bounded by the water that was there to begin with.

    A consequence of that fix, easy to miss and expensive to debug: scaling
    the latent heat by th_liquid_old means h(T) is NO LONGER A FIXED FUNCTION
    OF STATE. It is re-anchored at the start of every time step, because the
    scale factor w = th_liquid_old changes from step to step. The old and new
    enthalpies within a single step are therefore evaluated on the same curve,
    but consecutive steps use different curves. Any expression that compares
    an enthalpy across the step boundary must be checked for consistency.

    The instance of this that actually occurred (fixed 2026-08, in
    transient_solve_ncz): the old-state enthalpy was computed as
    h_of_T_only(T_old, m, w), which places the layer's liquid content at
    w * frac(T_old) with frac = (T_old + eps)/eps < 1. But w IS the liquid
    actually present at the start of the step -- not some fraction of it. The
    old state was thus credited with less latent heat than it held, so less
    energy had to be removed to reach any given new temperature, and freezing
    came out too cheap. On the Neumann benchmark the front advanced as t^0.60
    instead of t^0.50 and overshot the analytical solution by 2x, exceeding
    even the pure-Stefan bound obtained by ignoring the warm water entirely.
    The fix is to evaluate the old enthalpy with wet layers snapped to T = 0,
    i.e. at the top of the ramp where a layer holding liquid by definition
    sits. Sensible heat across the window    

    The other two solvers avoid the trap by different means, worth knowing
    when choosing among them: transient_solve_enthalpy tracks th_liquid as
    explicit state through its Picard iteration and clamps spurious growth
    (reporting the amount clawed back via claw_mushy/claw_dry), while
    transient_solve_decp rebalances latent heat explicitly after each
    diffusion sub-step, refreezing only water that exists.

    Early warning sign in any new solver here: refreezeDiff's "liquid
    gained in N layers" diagnostic firing, especially with "initially dry"
    counts above zero. Diffusion cannot create liquid water.

Live solver functions (called from diffusion.py's refreezeDiff / heatDiff):
    - transient_solve_TR: standard diffusion, no phase change (heat/isotope/gas).
    - transient_solve_ncz: refreezing via the nested Newton-Casulli-Zanolli
      enthalpy method (Tubini et al. 2021). RECOMMENDED -- most accurate of
      the four on the one case where an analytical solution exists, and the
      only one with a convergence proof at arbitrary time step.
    - transient_solve_enthalpy: refreezing via an enthalpy formulation with
      Picard iteration, a large effective heat capacity in mushy layers, and
      overshoot clamping. Reliable but ~5x less accurate than NCZ.
    - transient_solve_decp: refreezing by operator splitting (diffuse, then
      rebalance latent heat), sub-stepped internally. Converges to the same
      accuracy as transient_solve_enthalpy, but only at iters >> the default 10.
    - transient_solve_ahc: refreezing via an Apparent Heat Capacity
      formulation. RETAINED FOR COMPARISON ONLY -- fails the analytical
      benchmark at every smearing width tested, and frequently fails to
      converge. See its docstring before using it for anything.

Relative accuracy (2026-08), maximum freezing-front position error on the
    analytical Neumann problem of Tubini et al. (2021) Sect. 4.1 -- a
    semi-infinite pure-water column, initially +5 C, with a -5 C Dirichlet
    surface, run 100 days. All at dz = 0.01 m, dt = 3600 s; exact front
    position at 100 days is 0.6805 m::

        transient_solve_ncz        0.00273 m   (0.00018 m at dz=1mm, dt=60s)
        transient_solve_enthalpy   0.01376 m
        transient_solve_decp       0.01380 m   at iters=1000
                                   0.02777 m   at iters=10 (default)
        transient_solve_ahc        0.775 m     best case, W=0.1, and
                                               1071 of 2400 steps failed

    Caveat: this benchmark is pure water, where liquid content IS a
    single-valued function of temperature -- exactly the condition that firn
    violates (see below). It exercises the solvers' treatment of latent heat
    and their convergence behaviour, but it does not test the firn-specific
    modification described in the next section. It is evidence, not proof,
    that the ordering above carries over to firn. The benchmark driver lives
    in neumann_test.py.

    NOTE: these four refreezing solvers are intentionally kept side-by-side
    for comparison during development; see [add reference/notes doc here]
    for evaluation of their relative accuracy/conservation properties.

Deprecated:
    - transient_solve_EN: superseded by transient_solve_enthalpy. Retained
      for reference only; not called elsewhere in the CFM.

Helper functions:
    - solver(): generic tridiagonal matrix solve (LAPACK dgtsv).
    - enthalpy_of() / invert_enthalpy(): enthalpy <-> (T, th_solid, th_liquid)
      conversions used by transient_solve_enthalpy.
'''

import numpy as np
np.set_printoptions(precision=4)
# from scipy import interpolate
import scipy.integrate
from scipy.sparse import spdiags
import scipy.sparse.linalg as splin
from constants import *
import sys
from scipy.linalg import lapack

_DIAG = {'iii': [], 'solver': [], 'energy_resid': [], 'count': [],
         'n_warm': [], 'n_badiag': []}


def _diag_log(name, iii, energy_resid, count, n_warm=0, n_badiag=0):
    '''
    Append per-call diagnostics to a module-level list for post-run
    analysis. Deliberately unconditional and cheap (a few list appends);
    memory is ~100 bytes per call, so a 34,000-call DYE-2 run costs a few
    MB. Call _diag_write() at the end of a run to dump it. NOT thread safe
    and NOT reset between runs -- call _diag_reset() if driving multiple
    runs in one process.
    '''
    _DIAG['iii'].append(iii)
    _DIAG['solver'].append(name)
    _DIAG['energy_resid'].append(energy_resid)
    _DIAG['count'].append(count)
    _DIAG['n_warm'].append(n_warm)
    _DIAG['n_badiag'].append(n_badiag)


def _diag_reset():
    for k in _DIAG:
        _DIAG[k] = []


def _diag_write(path):
    import csv
    keys = list(_DIAG)
    with open(path, 'w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(keys)
        wr.writerows(zip(*[_DIAG[k] for k in keys]))
    print(f'wrote {len(_DIAG["iii"])} diagnostic rows to {path}')

def solver(a_U, a_D, a_P, b, use_dgtsv=True):
    '''
    Solve the tridiagonal linear system arising from the FV discretization.

    Given the FV coefficients (upper/lower neighbor coupling and diagonal)
    and RHS vector, solves for the updated field using either LAPACK's
    dedicated tridiagonal solver (dgtsv, faster; default) or a general
    sparse solver (scipy.sparse spsolve, slower but more general -- useful
    as an independent check if dgtsv results are ever in doubt).

    :param a_U: upper-neighbor coefficient at each layer [same units as a_P]
    :param a_D: lower-neighbor coefficient at each layer
    :param a_P: diagonal coefficient at each layer
    :param b: right-hand-side vector at each layer
    :param use_dgtsv: if True (default), use LAPACK dgtsv; if False, use
        scipy.sparse.linalg.spsolve

    :return: phi_t, solution vector (updated field at each layer)
    '''

    if use_dgtsv:
        dl = np.ascontiguousarray(a_U[1:])
        d  = np.ascontiguousarray(-a_P)
        du = np.ascontiguousarray(a_D[:-1])
        rhs = np.ascontiguousarray(-b)
        _, _, _, phi_t, _ = lapack.dgtsv(dl, d, du, rhs)
    else:
        nz = np.size(b)
        diags = np.append([a_U, -a_P], [a_D], axis=0)
        cols = np.array([1, 0, -1])
        big_A = spdiags(diags, cols, nz, nz, format='csc').T
        rhs = -b
        phi_t = splin.spsolve(big_A, rhs)

    return phi_t
    ###################
    #### end solver ###
    ###################

def _fv_geometry(z_edges, Z_P, Gamma_P):
    '''
    Compute finite volume geometry and face-conductance terms shared by
    all diffusion/refreezing solvers in this module.

    Extracted from what was previously duplicated inline in each solver
    function (transient_solve_TR, transient_solve_EN, transient_solve_enthalpy,
    transient_solve_ahc, transient_solve_decp). Computes layer thicknesses,
    distances to upper/lower neighbor centers, and harmonic-mean face
    conductances (Patankar 1980, eq. 4.9) between each layer and its
    neighbors.

    :param z_edges: layer (control volume) edges [m], length nz_P+1
    :param Z_P: layer (control volume) centers [m], length nz_P
    :param Gamma_P: conductivity/diffusivity at layer centers [W/m/K or
        equivalent], length nz_P

    :return: (dZ, deltaZ_u, deltaZ_d, Gamma_u, Gamma_d)
        dZ: layer thickness [m], length nz_P
        deltaZ_u: distance from each layer center to its upper neighbor's
            center [m] (first layer uses distance to its own edge/neighbor
            per boundary convention -- see source)
        deltaZ_d: distance from each layer center to its lower neighbor's
            center [m] (last layer analogous to deltaZ_u)
        Gamma_u: harmonic-mean face conductance at the upper face of each
            layer [W/m2/K when divided by deltaZ_u in the calling function]
        Gamma_d: harmonic-mean face conductance at the lower face of each
            layer

    NOTE: includes a np.maximum(..., 1e-12) floor on Gamma_P/Gamma_U/Gamma_D
        before dividing, to avoid division by zero. This guard was present
        in transient_solve_enthalpy/_ahc/_decp's inline versions but NOT in
        transient_solve_TR/transient_solve_EN's original inline versions --
        using this shared helper for those two functions will introduce a
        (harmless in normal cases) small numerical difference; see cleanup
        notes.

    Reference: Patankar (1980), Numerical Heat Transfer and Fluid Flow, eq. 4.9.
    '''

    dZ = np.diff(z_edges)
    Z_P_diff = np.diff(Z_P)

    deltaZ_u = np.zeros_like(Z_P)
    deltaZ_u[0] = Z_P_diff[0]
    deltaZ_u[1:] = Z_P_diff

    deltaZ_d = np.zeros_like(Z_P)
    deltaZ_d[:-1] = Z_P_diff
    deltaZ_d[-1] = Z_P_diff[-1]

    f_u = 1 - (Z_P - z_edges[:-1]) / deltaZ_u
    f_d = 1 - (z_edges[1:] - Z_P) / deltaZ_d

    Gamma_U = np.zeros_like(Gamma_P)
    Gamma_U[0] = Gamma_P[0]
    Gamma_U[1:] = Gamma_P[:-1]

    Gamma_D = np.zeros_like(Gamma_P)
    Gamma_D[:-1] = Gamma_P[1:]
    Gamma_D[-1] = Gamma_P[-1]

    Gamma_u = 1.0 / ((1 - f_u) / np.maximum(Gamma_P, 1e-12) + f_u / np.maximum(Gamma_U, 1e-12))
    Gamma_d = 1.0 / ((1 - f_d) / np.maximum(Gamma_P, 1e-12) + f_d / np.maximum(Gamma_D, 1e-12))

    return dZ, deltaZ_u, deltaZ_d, Gamma_u, Gamma_d

    ########################
    ### end _fv_geometry ###
    ########################

def _apply_bcs(a_U, a_D, a_P, b, deltaZ_u, deltaZ_d,
                bc_u_0, bc_type_u, bc_d_0, bc_type_d):
    '''
    Apply upper and lower boundary conditions to the FV coefficient arrays.

    Extracted from what was previously duplicated inline in each solver
    function. Supports Dirichlet (fixed value) or specified-gradient
    boundary conditions independently at each end of the domain, e.g. to
    allow a zero-flux (no-flux) surface boundary condition for testing,
    in addition to the previously-hardcoded fixed-temperature surface BC.

    Modifies a_U, a_D, a_P, b in place (and returns them for convenience).

    :param a_U: upper-neighbor coefficient at each layer (modified in place)
    :param a_D: lower-neighbor coefficient at each layer (modified in place)
    :param a_P: diagonal coefficient at each layer (modified in place)
    :param b: right-hand-side vector at each layer (modified in place)
    :param deltaZ_u: distance to upper neighbor [m] (see _fv_geometry);
        used to scale the lower boundary's gradient BC
    :param deltaZ_d: distance to lower neighbor [m] (see _fv_geometry);
        used to scale the upper boundary's gradient BC
    :param bc_u_0: upper boundary value (fixed value [K or deg C] if
        bc_type_u==1, or gradient [value/m] if bc_type_u==2)
    :param bc_type_u: upper boundary type; 1 = Dirichlet (fixed value),
        2 = specified gradient (e.g., 0.0 for a no-flux surface condition)
    :param bc_d_0: lower boundary value (same convention as bc_u_0)
    :param bc_type_d: lower boundary type (same convention as bc_type_u)

    :return: (a_U, a_D, a_P, b), same arrays as passed in (mutated)

    :raises ValueError: if bc_type_u or bc_type_d is not 1 or 2

    NOTE: bc_type_u/bc_type_d use "specified gradient" in units of
        value/meter (e.g., K/m), NOT a physical heat flux in W/m2. If you
        need a true physical flux BC (e.g., a specified W/m2 surface flux),
        you would need to convert flux -> gradient by dividing by the local
        conductivity (Gamma_P[0] or Gamma_P[-1]) before calling this function.
    '''

    # Upper boundary
    a_P[0] = 1
    if bc_type_u == 1:
        a_U[0] = 0
        a_D[0] = 0
        b[0]   = bc_u_0
    elif bc_type_u == 2:
        a_U[0] = 0
        a_D[0] = 1
        b[0]   = deltaZ_d[0] * bc_u_0
    else:
        raise ValueError(f"bc_type_u must be 1 or 2, got {bc_type_u}")

    # Lower boundary
    a_P[-1] = 1
    if bc_type_d == 1:
        a_U[-1] = 0
        a_D[-1] = 0
        b[-1]   = bc_d_0
    elif bc_type_d == 2:
        a_U[-1] = 1
        a_D[-1] = 0
        b[-1]   = deltaZ_u[-1] * bc_d_0
    else:
        raise ValueError(f"bc_type_d must be 1 or 2, got {bc_type_d}")

    return a_U, a_D, a_P, b
    ########################
    ### end _apply_bcs ###
    ########################


def transient_solve_TR(z_edges, Z_P, dt, Gamma_P, phi_0, nz_P, phi_s, c_vol, airdict=None):
    '''
    Standard transient 1-D diffusion solver (finite volume), no phase change.

    Used for heat, isotope, and firn-air diffusion where no liquid water is
    present (contrast with the refreezing solvers in this module, which
    handle latent heat from liquid water). If airdict is provided, includes
    additional physics for gas diffusion (gravitational fractionation,
    thermal fractionation, advection due to porosity change); otherwise
    solves plain diffusion with a volumetric heat-capacity source term
    (via c_vol).

    Grid / FV convention:
        z_edges, Z_P follow the same convention as _fv_geometry() in this
        module, though this function currently computes its own local
        version of that geometry inline rather than calling the helper.

    Temperature convention: Kelvin (contrast with the refreezing solvers,
        which use deg C with fusion at T=0).

    :param z_edges: layer (control volume) edges [m], length nz_P+1
    :param Z_P: layer (control volume) centers [m], length nz_P
    :param dt: time step [s]
    :param Gamma_P: diffusivity/conductivity at layer centers (thermal
        conductivity [W/m/K] for heat; diffusivity for isotope/gas)
    :param phi_0: field profile at start of step (temperature [K], isotope
        ratio, or gas concentration depending on use case)
    :param nz_P: number of layer centers (length of Z_P)
    :param phi_s: surface (upper boundary) value, used as Dirichlet BC
    :param c_vol: volumetric heat capacity [J/m3/K] (rho * cp), used to build
        the transient term a_P_0; only used in the non-gas (airdict is None) branch
    :param airdict: optional dict of gas-diffusion-specific parameters
        (por_op, d_eddy, gravity, thermal, deltaM, Tz, omega, dz, rho, z_co).
        If provided, adds gravitational/thermal fractionation and advection
        terms appropriate for firn-air diffusion. See inline comments for
        required keys.

    :return: phi_t (updated field profile), or (phi_t, w_p) if airdict is
        provided, where w_p is the advection velocity at layer centers [m/s]

    Reference: Patankar (1980), Numerical Heat Transfer and Fluid Flow.
    '''

    phi_t = phi_0
    phi_t_old = phi_t.copy()

    dZ, deltaZ_u, deltaZ_d, Gamma_u, Gamma_d = _fv_geometry(z_edges, Z_P, Gamma_P)
        
    #######################################
    # this part is for gas diffusion, which takes a bit more physics
    if airdict!=None:
        Gamma_Po    = Gamma_P * airdict['por_op'] #This is the diffusivity times the open porosity.

        Gamma_U     = np.append(Gamma_Po[0], Gamma_Po[0: -1] )
        Gamma_D     = np.append(Gamma_Po[1:], Gamma_Po[-1])
        Gamma_u     =  1 / ((1 - f_u) / Gamma_Po + f_u / Gamma_U) #Patankar Eq. 4.11
        Gamma_d     =  1 / ((1 - f_d) / Gamma_Po + f_d / Gamma_D)

        d_eddy_P    = airdict['d_eddy'] * airdict['por_op']
        d_eddy_U    = np.append(d_eddy_P[0], d_eddy_P[0:-1] )
        d_eddy_D    = np.append(d_eddy_P[1:], d_eddy_P[-1])
        d_eddy_u    =  1/ ( (1 - f_u)/d_eddy_P + f_u/d_eddy_U )
        d_eddy_d    =  1/ ( (1 - f_d)/d_eddy_P + f_d/d_eddy_D )

        if airdict['gravity']=="off" and airdict['thermal']=="off":
            S_C_0   = 0.0

        elif airdict['gravity']=='on' and airdict['thermal']=='off':
            S_C_0   = (-Gamma_d + Gamma_u) * (airdict['deltaM'] * GRAVITY / (R * airdict['Tz'])) / airdict['dz'] #S_C is independent source term in Patankar

        elif airdict['gravity']=='on' and airdict['thermal']=='on':
            # dTdz    = np.gradient(airdict['Tz'])/airdict['dz']
            dTdz    = np.gradient(airdict['Tz'], Z_P)
            Gamma_del = (Gamma_d-Gamma_u)
            # Gamma_del[Gamma_del<0]=1e-65
            S_C_0   = Gamma_del * (-(airdict['deltaM'] * GRAVITY / (R * airdict['Tz'])) + (airdict['omega'] * dTdz)) / airdict['dz'] # should thermal still work in LIZ? if so use d_eddy+diffu
            # S_C_0[S_C_0<0]=1.e-40
        else:
            print('Error at in solver.py at 119')
            sys.exit()

        S_C         = S_C_0 * phi_t
        b_0         = S_C * dZ

        rho_edges = np.interp(z_edges,Z_P,airdict['rho'])

        w_edges = w(airdict, z_edges, rho_edges, Z_P, dZ) # advection term (upward relative motion due to porosity changing)

        w_p = np.interp(Z_P,z_edges,w_edges) # Units m/s
        w_edges[z_edges>airdict['z_co']] = 0.0
        w_u = w_edges[0:-1]
        w_d = w_edges[1:]

        D_u = ((Gamma_u+d_eddy_u) / deltaZ_u) # Units m/s
        D_d = ((Gamma_d+d_eddy_d) / deltaZ_d)

        F_u =  w_u * airdict['por_op'] # Units m/s
        F_d =  w_d * airdict['por_op']

        P_u = F_u / D_u
        P_d = F_d / D_d

        op_ind              = np.where(z_edges<=airdict['z_co'])[0] #indices of all nodes wiht open porosity (shallower than CO)
        op_ind2             = np.where(z_edges<=airdict['z_co']+20)[0] # a bit deeper
        co_ind              = op_ind[-1]

        a_U = D_u * A( P_u ) + F_upwind(  F_u )
        a_D = D_d * A( P_d ) + F_upwind( -F_d )

        a_P_0 = airdict['por_op'] * dZ / dt
    #######################################
    ### end gas physics portion ###########
    #######################################

    #######################################
    else: # just for heat, enthalpy, isotope diffusion

        S_C = 0
        S_C = S_C * np.ones(nz_P)

        D_u = (Gamma_u / deltaZ_u)
        D_d = (Gamma_d / deltaZ_d)

        b_0 = S_C * dZ # first term of Patankar eq. 4.41d

        a_U = D_u # Patankar eq. 4.41a,b
        a_D = D_d # Patankar eq. 4.41a,b

        # a_P_0 = dZ / dt
        # a_P_0 = tot_rho * dZ / dt #  (old)
        a_P_0 = c_vol * dZ / dt # (new) Patankar eq. 4.41c
        # a_P_0 = RHO_I * c_firn * dZ / dt
    #######################################

    S_P     = 0.0
    a_P     = a_U + a_D + a_P_0 - S_P*dZ

    b       = b_0 + a_P_0 * phi_t #Patankar 4.41d

    #######################################
    ### Boundary conditions:
    ### type 1 is a specified value, type 2 is a specified gradient
    ### (units for gradient are degrees/meter)
    ### need to pay attention to surface boundary for gas
    
    bc_u_0    = phi_s      # or 0.0 with bc_type_u=2 for a no-flux surface test
    bc_type_u = 1
    
    bc_d_0    = 0
    bc_type_d = 2

    a_U, a_D, a_P, b = _apply_bcs(a_U, a_D, a_P, b, deltaZ_u, deltaZ_d,
                                    bc_u_0, bc_type_u, bc_d_0, bc_type_d)

    phi_t = solver(a_U, a_D, a_P, b)

    a_P = a_U + a_D + a_P_0

    if airdict!=None:
        return phi_t, w_p
    else:
        return phi_t

###################################
### end transient_solve_TR ########
###################################

############################################# 
### solvers for melt water freezing below ###
#############################################
def enthalpy_of(T_C, th_solid, th_liquid):
    '''
    Compute volumetric enthalpy from temperature and phase composition.

    Enthalpy is referenced to fusion at T=0 deg C: a fully-solid layer at
    T=0 has zero latent contribution, and enthalpy increases with both
    sensible heat (temperature away from 0) and latent heat (liquid
    fraction present). Air is treated as thermally inert (no contribution).

    :param T_C: temperature [deg C], fusion at T=0
    :param th_solid: solid (ice) mass per total volume [kg/m3]
    :param th_liquid: liquid water mass per total volume [kg/m3]

    :return: volumetric enthalpy [J/m3]
    '''
    return th_solid * CP_I * T_C + th_liquid * (CP_W * T_C + LF_I)
    #######################
    ### end enthalpy_of ###
    #######################

def _energy_residual(T_new, th_solid, th_liquid, T_old_C, th_solid_old,
                     th_liquid_old, dZ, Gamma_u, deltaZ_u, dt):
    '''
    Interior-cell energy conservation residual [J/m2], common to all four
    refreezing solvers so that their values are directly comparable.

    Uses enthalpy_of() for both states -- referenced to h = 0 for
    fully-solid firn at T = 0 -- rather than any solver's internal
    enthalpy function. transient_solve_ncz in particular works internally
    with h_of_T_only(), whose reference is h = 0 at T = -eps and whose
    latent term is scaled by th_liquid_old; residuals computed with that
    function are NOT comparable across solvers.

    Cell 0 is excluded: it carries a Dirichlet boundary condition, which is
    an unaccounted source/sink, so surface energy is not conserved by
    construction. Only the 0|1 face flux enters. The lower boundary is
    assumed zero-gradient (no flux).

    :return: dH_interior - flux_in [J/m2]; zero for a perfectly
        conservative scheme, up to the backward-Euler flux approximation.
    '''
    aU_top   = Gamma_u[1] / deltaZ_u[1]
    flux_top = aU_top * (T_new[0] - T_new[1]) * dt
    H_old = enthalpy_of(T_old_C, th_solid_old, th_liquid_old)
    H_new = enthalpy_of(T_new, th_solid, th_liquid)
    dH_interior = np.sum((H_new[1:] - H_old[1:]) * dZ[1:])
    return dH_interior - flux_top
    ############################
    ### end _energy_residual ###
    ############################    

def invert_enthalpy(Hhat, th_solid, th_liquid):
    '''
    Invert volumetric enthalpy to recover temperature and phase partition.

    Given total volumetric enthalpy and total H2O mass (th_solid + th_liquid,
    conserved), determines whether the layer is fully frozen (T<0), fully
    liquid (T>0), or mushy/isothermal at the fusion point (T=0, partial
    liquid fraction determined by available latent enthalpy).

    :param Hhat: volumetric enthalpy [J/m3] (see enthalpy_of)
    :param th_solid: solid (ice) mass per total volume [kg/m3] (used only to
        compute total mass m = th_solid + th_liquid; the solid/liquid split
        is recomputed from scratch based on Hhat, not adjusted incrementally)
    :param th_liquid: liquid water mass per total volume [kg/m3] (see above)

    :return: (T_C, th_solid_new, th_liquid_new)
        T_C: temperature [deg C], fusion at T=0
        th_solid_new: updated solid mass per total volume [kg/m3]
        th_liquid_new: updated liquid mass per total volume [kg/m3]
        (th_solid_new + th_liquid_new == th_solid + th_liquid, mass conserved)

    NOTE: if a layer has zero total mass (m=0) and Hhat==0 exactly, this
        function will divide by zero (0/0 -> NaN) in the "cold" branch
        rather than raising an error. Not expected in practice for real
        firn layers, but no explicit guard exists; see cleanup notes for
        a suggested fix (clamp m with np.maximum before dividing).
    '''

    m = th_solid + th_liquid                      # total H2O mass/vol [kg/m3]
    latent = m * LF_I                           # enthalpy to melt all of it
    T   = np.zeros_like(Hhat)
    solid = np.empty_like(Hhat)
    liquid = np.empty_like(Hhat)

    cold = Hhat <= 0.0
    warm = Hhat >= latent
    mush = ~cold & ~warm

    # cold: all ice, T<0
    solid[cold] = m[cold]; liquid[cold] = 0.0
    T[cold] = Hhat[cold] / (m[cold] * CP_I)

    # warm: all liquid, T>0
    liquid[warm] = m[warm]; solid[warm] = 0.0
    T[warm] = (Hhat[warm] - latent[warm]) / (m[warm] * CP_W)

    # mush: T=0, split by latent content
    liquid[mush] = Hhat[mush] / LF_I
    solid[mush] = m[mush] - liquid[mush]
    return T, solid, liquid
    ###########################
    ### end invert_enthalpy ###
    ###########################

def transient_solve_enthalpy(z_edges, Z_P, dt, Gamma_P, T_old_C, th_liquid_old, th_solid_old, iii,
                              big_slope=1e13, max_iter=200, tol=1e-8, T_eps=0.0, *,
                        bc_u=None, bc_type_u=1, bc_d=0.0, bc_type_d=2):
    '''
    Enthalpy-based refreezing solver.

    Solves heat diffusion with phase change using a fully-implicit enthalpy
    formulation with Picard iteration. Latent heat is handled via a locally
    large effective heat capacity ("big_slope") in mushy (partially liquid)
    layers, followed by an explicit enthalpy inversion each iteration to
    recover consistent temperature and solid/liquid mass. Includes an
    overshoot "clamp" correction that prevents diffusion from spuriously
    creating liquid water in layers that were dry, or increasing liquid
    beyond what energy balance allows, and tallies the amount clamped.

    Grid / FV convention:
        z_edges, Z_P follow the same convention as _fv_geometry() in this module.

    :param z_edges: layer (control volume) edges [m], length nz_P+1
    :param Z_P: layer (control volume) centers [m], length nz_P
    :param dt: time step [s]
    :param Gamma_P: thermal conductivity at layer centers [W/m/K]
    :param T_old_C: temperature profile at start of step [deg C], fusion at T=0
    :param th_liquid_old: liquid water mass per total volume at start of step [kg/m3]
    :param th_solid_old: solid (ice) mass per total volume at start of step [kg/m3]
    :param big_slope: effective dH/dT [J/m3/K] assigned to mushy layers to pin
        them near T=0 during iteration. Sensitivity-tested on a DYE-2, Greenland
        firn run (2026-07): 1e6 fails to converge; 1e9 (old default) converges
        but with occasional clamp corrections (claw_mushy/claw_dry); 1e12-1e13
        converge cleanly with no clamping; 1e15 converges but iterates slowly
        (possible mild ill-conditioning at very high values). Default updated
        to 1e13 based on this test; not exhaustively validated across other
        sites/conditions.
    :param max_iter: maximum Picard iterations (default 200).
    :param tol: convergence tolerance on max temperature change between
        Picard iterations [deg C]
    :param T_eps: legacy snap width [deg C]. If > 0, mushy layers whose
        post-solve temperature falls within +/- T_eps of fusion are set to
        exactly 0 before the enthalpy update. DEFAULT 0.0 (disabled).
        Enabling it injects up to big_slope * T_eps [J/m3] of spurious
        enthalpy per mushy layer per iteration, because the perturbed
        temperature is immediately multiplied by dHdT = big_slope. Retained
        only to reproduce pre-2026-08 behaviour; invert_enthalpy already
        returns T = 0 exactly for mushy layers, so the snap has no effect
        on the returned state and only corrupts Hhat_new.

    :return:
        dict with keys::

            TzC_return: updated temperature profile [deg C]
            th_solid: updated solid mass per total volume [kg/m3]
            th_liquid: updated liquid mass per total volume [kg/m3]
            count: number of Picard iterations used
            claw_mushy: liquid mass clamped back (spurious growth in already-mushy
                layers) this step [kg/m2]
            claw_dry: liquid mass clamped back (spurious creation in previously-dry
                layers) this step [kg/m2]
            energy_resid: interior-cell energy conservation residual [J/m2]

    Reference(s): Voller and Swaminathan (1991), eqs. 31-32;
        Voller, Swaminathan, and Thomas (1990), eq. 61.

    VALIDATION (2026-08): pure-water Neumann benchmark, dz=0.01 m, dt=3600 s:
        front-position error 0.01376 m against 0.00273 m for
        transient_solve_ncz on the same grid. Converges reliably (no failed
        steps), but is ~5x less accurate. Some of that gap is conductivity
        feedback rather than the scheme itself: big_slope pins mushy layers
        to T ~ 0 +/- rounding, which flips a temperature-keyed lambda(T)
        between water and ice values (a 3.5x jump) at the front. With
        conductivity keyed on liquid fraction instead, this solver and
        transient_solve_ncz agree to all printed digits (0.00882 m each).
    '''

    dZ, deltaZ_u, deltaZ_d, Gamma_u, Gamma_d = _fv_geometry(z_edges, Z_P, Gamma_P)

    a_U_base = Gamma_u / deltaZ_u
    a_D_base = Gamma_d / deltaZ_d

    a_U_work = np.empty_like(a_U_base)
    a_D_work = np.empty_like(a_D_base)

    beta = dZ / dt
    
    Hhat_old = enthalpy_of(T_old_C, th_solid_old, th_liquid_old)
    
    T = T_old_C.copy()
    th_solid = th_solid_old.copy()
    th_liquid = th_liquid_old.copy()
    
    dHdT = np.empty_like(T_old_C, dtype=float)

    count = 0
    # --- clamp backstops with tally ---
    claw_mushy = 0.0
    claw_dry   = 0.0
    
    change_history = []
    converged = True
    for i_time in range(max_iter): # Testing indicates that this should never need this many iterations

        Hhat_m = enthalpy_of(T, th_solid, th_liquid)

        m = th_solid + th_liquid
        wet = th_liquid > 0.0
        mushy = wet & (th_solid > 0.0)

        dHdT[:] = m * CP_I
        dHdT[wet] = (m * CP_W)[wet]
        dHdT[mushy] = big_slope
        np.maximum(dHdT, 1.0, out=dHdT)
        
        # _apply_bcs mutates in place, so the base arrays must be shielded
        a_U_work[:] = a_U_base
        a_D_work[:] = a_D_base
        a_U, a_D = a_U_work, a_D_work

        a_P = a_U + a_D + beta * dHdT
        b   = beta * (dHdT * T - Hhat_m + Hhat_old)

        ###############

        ### Boundary conditions:
        ### type 1 is a specified value, type 2 is a specified gradient
        ### (units for gradient are degrees/meter)
        bc_u_0    = T_old_C[0] if bc_u is None else bc_u
        bc_d_0    = bc_d

        a_U, a_D, a_P, b = _apply_bcs(a_U, a_D, a_P, b, deltaZ_u, deltaZ_d,
                                        bc_u_0, bc_type_u, bc_d_0, bc_type_d)

        #####
        T_new = solver(a_U, a_D, a_P, b) #sensible enthalpy. 0 for layers at freezing (have LWC), negative for dry layers
        #####
        
        ### The crux is to adjust liquid fraction and temperture field based on solution
        ### Note previous ways of solving in dev branch and previous releases.

        ### fix for spurious melting
        ### mushy layers must sit at fusion T; kill sub-tolerance residual
        # T_new[mushy & (np.abs(T_new) < T_eps)] = 0.0     # T_eps = 1e-8
        if T_eps > 0.0:
            T_new[mushy & (np.abs(T_new) < T_eps)] = 0.0

        Hhat_new = Hhat_m + dHdT * (T_new - T)           # dHdT = big_slope = 1e13

        T_cons, th_solid_new, th_liquid_new = invert_enthalpy(Hhat_new, th_solid, th_liquid)

        ###  "clamp" fix
        ### --- no spurious melting in layers that were already mushy ---
        ### diffusion cannot ADD liquid to a 0C layer; cap at the pre-iteration value
        grew = (th_liquid_new > th_liquid_old) & (th_liquid_old > 0.0)
        excess = np.where(grew, th_liquid_new - th_liquid_old, 0.0)
        claw_mushy += np.sum(excess * dZ)          # kg/m2 removed this step
        th_liquid_new -= excess
        th_solid_new  += excess
        T_cons[grew]   = 0.0

        dry_grew = (th_liquid_new > 0.0) & (th_liquid_old == 0.0)
        excess_d = np.where(dry_grew, th_liquid_new, 0.0)   # all of it is spurious
        claw_dry += np.sum(excess_d * dZ)          # kg/m2 removed this step
        th_liquid_new -= excess_d
        th_solid_new  += excess_d
        T_cons[dry_grew] = 0.0
        
        change = np.max(np.abs(T_cons - T))
        change_history.append(change)
        count += 1
        T, th_solid, th_liquid = T_cons, th_solid_new, th_liquid_new
        if change < tol:
            break

        ### END ITERATION LOOP
        ######################
    
    else:
        converged = False
        change_history_arr = np.array(change_history)

        # crude oscillation check: are we alternating up/down rather than
        # monotonically (or near-monotonically) decreasing?
        diffs = np.diff(change_history_arr)
        sign_flips = np.sum(np.diff(np.sign(diffs)) != 0)
        frac_flips = sign_flips / len(diffs) if len(diffs) > 0 else 0.0

        # crude plateau check: did the last N iterations stop improving much?
        N = 20
        tail = change_history_arr[-N:]
        tail_improvement = tail[0] - tail[-1]
        tail_ratio = tail_improvement / tail[0] if tail[0] != 0 else 0.0

        print(f"WARNING: enthalpy solver did not converge in {max_iter} iterations (iii={iii})")
        print(f"  final change: {change_history_arr[-1]:.3e} (tol={tol:.1e})")
        print(f"  first change: {change_history_arr[0]:.3e}")
        print(f"  min change over run: {change_history_arr.min():.3e}")
        print(f"  fraction of sign flips in diffs (oscillation indicator): {frac_flips:.2f}")
        print(f"  last {N} iters: improved by {tail_ratio*100:.1f}% (plateau if near 0)")

        # optional: dump the whole history for manual inspection
        # print(f"  change_history: {change_history_arr}")

    # --- interior-cell energy residual (same diagnostic as AHC) ---
    energy_resid = _energy_residual(T, th_solid, th_liquid, T_old_C,
                                    th_solid_old, th_liquid_old,
                                    dZ, Gamma_u, deltaZ_u, dt)

    _diag_log('enthalpy', iii, energy_resid, count,
              n_warm=int(np.sum(T > 0.0)))
    
    return dict(TzC_return=T, th_solid=th_solid, th_liquid=th_liquid,
                count=count, claw_mushy=claw_mushy, claw_dry=claw_dry,energy_resid=energy_resid)

####################################
### end transient_solve_enthalpy ###
####################################

def transient_solve_ahc(z_edges, Z_P, dt, Gamma_P, T_old_C,
                        th_liquid_old, th_solid_old, iii, W=1.0,
                        max_iter=50, tol=1e-10, *,
                        bc_u=None, bc_type_u=1, bc_d=0.0, bc_type_d=2):
    '''
    Apparent Heat Capacity (AHC) refreezing solver.

    Solves heat diffusion with phase change by folding latent heat into an
    "apparent" heat capacity, smeared over a fixed cold-side temperature
    window [-W, 0] rather than tracked via an explicit enthalpy inversion
    (contrast with transient_solve_enthalpy). Liquid/solid mass fractions
    are inferred after convergence from the liquid-fraction function f_l(T),
    not tracked as independent state during iteration. Uses the same FV
    stencil, face conductances, and boundary conditions as the other
    refreezing solvers in this module.

    Grid / FV convention:
        z_edges, Z_P follow the same convention as _fv_geometry() in this module.

    :param z_edges: layer (control volume) edges [m], length nz_P+1
    :param Z_P: layer (control volume) centers [m], length nz_P
    :param dt: time step [s]
    :param Gamma_P: thermal conductivity at layer centers [W/m/K]
    :param T_old_C: temperature profile at start of step [deg C], fusion at T=0
    :param th_liquid_old: liquid water mass per total volume at start of step [kg/m3]
    :param th_solid_old: solid (ice) mass per total volume at start of step [kg/m3]
    :param W: temperature window [K] over which the liquid water present in
        a layer refreezes (T in [-W, 0]). Smears only the latent heat of
        pre-existing liquid, not of total layer mass.

        DO NOT REDUCE W EXPECTING BETTER ACCURACY. Swept on the pure-water
        Neumann benchmark (2026-08, dz=0.01 m, dt=3600 s, 2400 steps):
        W=1.0 -> 1.094 m error, 373 non-converged steps; W=0.1 -> 0.775 m,
        1071 failed; W=0.01 -> 1.675 m, 978 failed; W=1e-3 -> 1.755 m;
        W=1e-4 -> 1.765 m. For reference the analytical front is 0.6805 m and
        transient_solve_ncz achieves 0.0027 m on the same grid, so AHC's best
        case is ~300x worse and is unreliable at every W tested.

        Mechanism of the small-W failure: as W shrinks, a cell can cool
        clean past the window [-W, 0] within a single time step, so
        dfl_dT is never sampled inside the ramp and the latent heat is
        simply never applied. The errors at W <= 0.01 saturate at the value
        obtained with NO latent heat at all. This is the discrete
        chain-rule failure described in Tubini et al. (2021) Sect. 2: AHC
        is analytically equivalent to the enthalpy formulation but not
        discretely equivalent, absent enforcement of Roe's condition, which
        this implementation does not do. Note that the front still grows as
        t^0.52 in the failed cases -- a latent-heat-free diffusion front is
        also self-similar, so sqrt(t) scaling does NOT detect this error.

        Recommendation: prefer transient_solve_ncz. This solver is retained
        for comparison, not for production use.
    :param max_iter: maximum Picard iterations (default 50)
    :param tol: convergence tolerance on max temperature change between
        Picard iterations [deg C]

    :return:
        dict with keys::

            TzC_return: updated temperature profile [deg C]
            th_solid: updated solid mass per total volume [kg/m3], computed as
                m - th_liquid where m = th_solid_old + th_liquid_old. Because
                th_liquid <= th_liquid_old always (see below), th_solid is
                guaranteed >= th_solid_old: this solver can only refreeze
                liquid, never melt ice.
            th_liquid: updated liquid mass per total volume [kg/m3], computed
                as th_liquid_old * f_l(T_final), i.e. the fraction of the
                layer's PRE-EXISTING liquid that remains unfrozen at the final
                temperature. Not tracked as independent state during iteration
                -- it is reconstructed from T_final and th_liquid_old after
                convergence. A layer with th_liquid_old == 0 therefore returns
                th_liquid == 0 exactly, regardless of T_final.
                Total H2O mass is conserved: th_solid + th_liquid ==
                th_solid_old + th_liquid_old, exactly, by construction.
            count: number of Picard iterations used
            energy_resid: interior-cell energy conservation residual [J/m2]

    NOTE: unlike transient_solve_enthalpy, this function does not return
        claw_mushy/claw_dry diagnostics (no overshoot clamping is applied).
        Callers (e.g., refreezeDiff) should not assume these keys are present.

    NOTE (fixed 2026-08): the liquid/solid partition previously scaled the
        liquid-fraction function f_l by TOTAL H2O mass m, i.e.
        th_liquid = m * f_l. Because f_l -> 1 as T -> 0 from below, a layer
        only microdegrees below fusion was assigned essentially all of its
        mass as liquid, driving th_solid (and hence rho and self.mass)
        toward zero and producing division-by-zero NaNs downstream in
        melt.py (rho = mass/dz, cold_content/(CP_I*mass)) and ultimately in
        SEB.py. The latent-heat term and the partition are now scaled by
        th_liquid_old -- the water actually present in the layer at the
        start of the step -- so a dry layer has no latent contribution and
        cannot gain liquid by diffusion. This mirrors the same correction
        applied in transient_solve_ncz; see that function's docstring for
        the fuller rationale.
    '''

    dZ, deltaZ_u, deltaZ_d, Gamma_u, Gamma_d = _fv_geometry(z_edges, Z_P, Gamma_P)    

    beta = dZ / dt

    a_U = Gamma_u / deltaZ_u
    a_D = Gamma_d / deltaZ_d

    m = th_solid_old + th_liquid_old        # total H2O mass/vol; conserved
    w = th_liquid_old.copy()        # add near "m = th_solid_old + th_liquid_old"
    T = T_old_C.copy()

    change_history = []
    converged = True
    
    for count in range(max_iter):                    
        ### apparent capacity evaluated at current iterate (Picard)
        f_l    = np.clip(1.0 + T / W, 0.0, 1.0)
        dfl_dT = np.where((T > -W) & (T < 0.0), 1.0 / W, 0.0)
        liq    = w * f_l                          # liquid mass/vol present
        C_app  = np.maximum((m - liq) * CP_I + liq * CP_W
                            + w * LF_I * dfl_dT, 1.0)

        a_P = a_U + a_D + beta * C_app
        b   = beta * C_app * T_old_C

        bc_u_0    = T_old_C[0] if bc_u is None else bc_u
        bc_d_0    = bc_d
        # bc_u_0    = T_old_C[0]
        # bc_type_u = 1
        # bc_d_0    = 0.0
        # bc_type_d = 2

        a_U, a_D, a_P, b = _apply_bcs(a_U.copy(), a_D.copy(), a_P, b, deltaZ_u, deltaZ_d,
                                        bc_u_0, bc_type_u, bc_d_0, bc_type_d)

        T_new = solver(a_U, a_D, a_P, b)

        change = np.max(np.abs(T_new - T))
        change_history.append(change)
        T = T_new
        if change < tol:                             # 8-space: loop body
            break                                    # 12-space: under the if

    else:
        converged = False
        # print(f"WARNING: AHC solver did not converge in {max_iter} iterations (iii={iii if 'iii' in locals() else '?'})")

    # infer phase from final temperature via f_l (AHC analog of th_solid/th_liquid)
    f_l    = np.clip(1.0 + T / W, 0.0, 1.0)
    th_liquid = w * f_l
    th_solid  = m - th_liquid

    # --- rigorous energy residual on interior cells (1..n-1) ---
    # top-interface conductance the solver actually used [W/m2/K]
    energy_resid = _energy_residual(T, th_solid, th_liquid, T_old_C,
                                    th_solid_old, th_liquid_old,
                                    dZ, Gamma_u, deltaZ_u, dt)

    _diag_log('ahc', iii, energy_resid, count,
              n_warm=int(np.sum(T > 0.0)))

    return dict(TzC_return=T, th_solid=th_solid, th_liquid=th_liquid,
            count=count, energy_resid=energy_resid, converged=converged, change=change)
###############################
### end transient_solve_ahc ###
###############################

def transient_solve_decp(z_edges, Z_P, dt, Gamma_P, T_old_C,
                         th_liquid_old, th_solid_old, iii, iters=10, *,
                        bc_u=None, bc_type_u=1, bc_d=0.0, bc_type_d=2, rebalance_top=True):
    '''
    Decoupled Energy-Conservation Parametrization (DECP) refreezing solver.

    Solves heat diffusion with phase change using operator splitting: each
    sub-step first solves sensible-heat-only FV diffusion (no latent heat
    term), then explicitly rebalances latent heat by refreezing liquid water
    in any layer left below fusion temperature after diffusion (partially or
    fully, depending on available "cold content" vs. latent heat required).
    Contrast with transient_solve_enthalpy (implicit, coupled via Picard
    iteration) and transient_solve_ahc (implicit, latent heat folded directly
    into an apparent heat capacity). Sub-stepped `iters` times per call at
    dt_sub = dt/iters to control operator-splitting error.

    Grid / FV convention:
        z_edges, Z_P follow the same convention as _fv_geometry() in this module.

    :param z_edges: layer (control volume) edges [m], length nz_P+1
    :param Z_P: layer (control volume) centers [m], length nz_P
    :param dt: total time step [s] (subdivided internally into `iters` sub-steps)
    :param Gamma_P: thermal conductivity at layer centers [W/m/K]
    :param T_old_C: temperature profile at start of step [deg C], fusion at T=0
    :param th_liquid_old: liquid water mass per total volume at start of step [kg/m3]
    :param th_solid_old: solid (ice) mass per total volume at start of step [kg/m3]
    :param iters: number of internal sub-steps (dt_sub = dt/iters).
        Sensitivity tested on the pure-water Neumann benchmark (2026-08):
        front-position error 0.0757 m at iters=1, 0.0278 at 10 (the default),
        0.0158 at 100, 0.0138 at 1000. Monotonic and roughly first order in
        dt_sub, converging to the same 0.0138 m floor reached by
        transient_solve_enthalpy -- the two schemes agree to 4 decimal places
        once splitting error is removed, as expected since both use the same
        stencil with explicit latent accounting. The default iters=10 is
        therefore ~2x its own converged error, at 10 tridiagonal solves per
        call. transient_solve_ncz reaches 0.00273 m on the same grid in ~2
        solves per call, i.e. more accurate and cheaper.
    :param rebalance_top: if True (default), layer 0 participates in the
        latent-heat rebalance. Correct for CFM, where layer 0 is a real firn
        layer with mass and LWC whose water must be allowed to refreeze even
        though its temperature is Dirichlet-pinned by SEB.py. Set False only
        for idealized tests in which node 0 represents a massless boundary.
        NOTE that when True, the rebalance can shift T[0] away from the
        imposed boundary value on the final sub-step (earlier sub-steps
        re-pin it), so the returned surface temperature may disagree with
        what SEB.py specified by up to the latent adjustment.

    NOTE (fixed 2026-08): the sensible heat capacity used in the diffusion
        half-step and the cold content used in the latent rebalance must be
        the SAME quantity. They were not: diffusion used a mass-weighted
        capacity while the rebalance used m * CP_I, so for liquid water the
        rebalance credited only CP_I/CP_W = 50% of the cold content that
        diffusion had actually created. Only half the expected mass refroze,
        the front lagged ~20%, and no amount of sub-stepping fixed it because
        the inconsistency is per-sub-step. Both now use c_vol, recomputed
        each sub-step as (m - th_liquid) * CP_I + th_liquid * CP_W. Any future
        change to one must be mirrored in the other.

    :return:
        dict with keys::

            TzC_return: updated temperature profile [deg C]
            th_solid: updated solid mass per total volume [kg/m3]
            th_liquid: updated liquid mass per total volume [kg/m3]
            count: number of sub-steps used (always equals `iters`; retained for
                return-signature consistency with the other solvers, which use
                `count` to report iterations actually needed for convergence)
            energy_resid: interior-cell energy conservation residual [J/m2],
                accumulated across all sub-steps

    NOTE: unlike transient_solve_enthalpy, this function does not return
        claw_mushy/claw_dry diagnostics (no overshoot clamping is applied
        in the same sense -- refreezing here is handled explicitly by the
        latent-heat rebalance step rather than by clamping a diffusion
        overshoot). Callers (e.g., refreezeDiff) should not assume these
        keys are present.
    '''

    dZ, deltaZ_u, deltaZ_d, Gamma_u, Gamma_d = _fv_geometry(z_edges, Z_P, Gamma_P)
    
    a_U = Gamma_u / deltaZ_u
    a_D = Gamma_d / deltaZ_d

    dt_sub = dt / iters
    beta   = dZ / dt_sub

    m = th_solid_old + th_liquid_old            # total H2O mass/vol; conserved
    Hhat_old = enthalpy_of(T_old_C, th_solid_old, th_liquid_old)

    T         = T_old_C.copy()
    th_solid  = th_solid_old.copy()
    th_liquid = th_liquid_old.copy()

    flux_top = 0.0                      # accumulate 0|1 face flux across sub-steps
    aU_top   = Gamma_u[1] / deltaZ_u[1] # 0|1 face conductance [W/m2/K]

    # c_vol = m * CP_I                    # sensible vol. heat cap [J/m3/K] (decoupled: no latent)
    
    # a_P0  = c_vol * beta                # a_P_0 per sub-step
    converged = True

    for count in range(iters):
        # ---- (1) diffusion sub-step (sensible only) ----
        # c_vol = np.where(th_liquid > 0.0, m * CP_W, m * CP_I)
        c_vol = (m - th_liquid) * CP_I + th_liquid * CP_W
        a_P0  = c_vol * beta
        aU = a_U.copy(); aD = a_D.copy()
        aP = aU + aD + a_P0
        b  = a_P0 * T                   # backward Euler, previous sub-step T

        bc_u_0    = T_old_C[0] if bc_u is None else bc_u
        bc_d_0    = bc_d
        # bc_u_0    = T_old_C[0]
        # bc_type_u = 1
        # bc_d_0    = 0.0
        # bc_type_d = 2

        aU, aD, aP, b = _apply_bcs(aU, aD, aP, b, deltaZ_u, deltaZ_d,
                                    bc_u_0, bc_type_u, bc_d_0, bc_type_d)

        T = solver(aU, aD, aP, b)
        flux_top += aU_top * (T[0] - T[1]) * dt_sub   # J/m2, this sub-step

        # ---- (2) latent-heat rebalance: refreeze where wet & cold ----
        wetcold      = (T < 0.0) & (th_liquid > 0.0)
        if not rebalance_top:
            wetcold[0] = False      # benchmark only: node 0 is a pure
                                    # Dirichlet boundary with no real mass.
                                    # In CFM, layer 0 IS a firn layer and
                                    # its liquid must be allowed to refreeze.
        cold_content = c_vol * (0.0 - T)   # same c_vol as diffusion half --
                                           # a mismatch here silently halves
                                           # the refreezing rate (fixed 2026-08)
        heattofreeze = th_liquid * LF_I

        partial = wetcold & (cold_content <  heattofreeze)
        full    = wetcold & (cold_content >= heattofreeze)

        # partial: not enough cold content -> T rises to 0, some liquid freezes
        refroze = cold_content / LF_I           # kg/m3
        th_liquid[partial] -= refroze[partial]
        th_solid[partial]  += refroze[partial]
        T[partial] = 0.0

        # full: enough cold content -> all liquid freezes, leftover cold warms ice
        th_solid[full]  += th_liquid[full]
        leftover = cold_content - heattofreeze  # J/m3, >=0 in full case
        T[full]  = -leftover[full] / (m[full] * CP_I)
        th_liquid[full] = 0.0

        ### clamp any layer above fusion back to 0 C
        # T[T > 0.0] = 0.0

    # --- interior-cell energy residual (same diagnostic as _enthalpy/_ahc) ---
    # aU_top   = Gamma_u[1] / deltaZ_u[1]             # 0|1 face conductance [W/m2/K]
    # flux_top = aU_top * (T[0] - T[1]) * dt          # J/m2 into interior (full dt)
    
    # --- interior-cell energy residual (flux accumulated over sub-steps) ---
    energy_resid_sub = (np.sum((enthalpy_of(T, th_solid, th_liquid)[1:]
                               - Hhat_old[1:]) * dZ[1:]) - flux_top)
    energy_resid = _energy_residual(T, th_solid, th_liquid, T_old_C,
                                    th_solid_old, th_liquid_old,
                                    dZ, Gamma_u, deltaZ_u, dt)

    _diag_log('decp', iii, energy_resid_sub, count,
              n_warm=int(np.sum(T > 0.0)))

    return dict(TzC_return=T, th_solid=th_solid, th_liquid=th_liquid, converged=converged,
        count=count, energy_resid=energy_resid, energy_resid_sub=energy_resid_sub)
################################
### end transient_solve_decp ###      
################################

################################
### Tubini method (ncz) ########
################################

def h_of_T_only(T_C, m, w, eps=1e-4, h_at_0=None):    
    '''
    Pure-temperature enthalpy h(T) [J/m3] for a firn layer.

    Adapted from Tubini et al. (2021) Eq. E1, with a key modification:
    the latent heat term is scaled by the liquid water actually present
    in the layer at the start of the step (w), NOT by total H2O mass (m).
    Tubini's formulation assumes liquid content is a single-valued
    function of T (valid for an SFCC or pure water); in firn, LWC is an
    independent state set by melt/percolation, so a dry layer must have
    no latent spike at all.

    Reference: h = 0 at T = -eps.

    :param T_C: temperature [deg C], fusion at T=0
    :param m: total H2O mass per volume [kg/m3] (sensible heat carrier)
    :param w: liquid water mass per volume at start of step [kg/m3]
        (latent heat carrier)
    :param eps: window [deg C] over which w refreezes
    :return: h(T) [J/m3]
    :param h_at_0: optional precomputed m*CP_I*eps + w*LF_I [J/m3]; pass
        from the caller to avoid recomputing it every Newton iteration.
    '''
    
    if h_at_0 is None:
        h_at_0 = m * CP_I * eps + w * LF_I

    cold  = T_C < -eps
    mushy = (T_C >= -eps) & (T_C < 0.0)
    warm  = T_C >= 0.0

    h = np.empty_like(T_C)
    h[cold] = m[cold] * CP_I * (T_C[cold] + eps)

    frac = (T_C[mushy] + eps) / eps
    h[mushy] = m[mushy] * CP_I * eps * frac + w[mushy] * LF_I * frac

    h[warm] = h_at_0[warm] + m[warm] * CP_W * T_C[warm]
    return h


def apparent_heat_capacity(T_C, m, w, eps=1e-4, Ca_mushy=None):
    '''Ca(T) = dh/dT, piecewise constant. Latent spike scales with w.'''
    if Ca_mushy is None:
        Ca_mushy = m * CP_I + w * LF_I / eps

    cold  = T_C < -eps
    mushy = (T_C >= -eps) & (T_C < 0.0)
    warm  = T_C >= 0.0

    Ca = np.empty_like(T_C)
    Ca[cold]  = m[cold] * CP_I
    Ca[mushy] = Ca_mushy[mushy]
    Ca[warm]  = m[warm] * CP_W
    return Ca


def _ca_max(m, w, eps):
    '''Peak of Ca(T); see earlier note on why the max is required.'''
    return np.maximum(m * CP_I + w * LF_I / eps, m * CP_W)


def p_of_T(T_C, m, w, eps=1e-4, Ca_max=None, Ca_mushy=None):
    '''p(T), Eq. 18.'''
    if Ca_max is None:
        Ca_max = _ca_max(m, w, eps)
    return np.where(T_C < 0.0,
                    apparent_heat_capacity(T_C, m, w, eps, Ca_mushy),
                    Ca_max)


def q_of_T(T_C, m, w, eps=1e-4, Ca_max=None):
    '''q(T), Eq. 18.'''
    if Ca_max is None:
        Ca_max = _ca_max(m, w, eps)
    return np.where(T_C < 0.0, 0.0, Ca_max - m * CP_W)


def h1_of_T(T_C, m, w, eps=1e-4, Ca_max=None, h_at_0=None):
    '''h1(T), Eq. 19.'''
    if Ca_max is None:
        Ca_max = _ca_max(m, w, eps)
    if h_at_0 is None:
        h_at_0 = m * CP_I * eps + w * LF_I
    return np.where(T_C < 0.0,
                    h_of_T_only(T_C, m, w, eps, h_at_0),
                    h_at_0 + Ca_max * T_C)


def h2_of_T(T_C, m, w, eps=1e-4, Ca_max=None, h_at_0=None):
    '''h2(T), Eq. 19.'''
    if Ca_max is None:
        Ca_max = _ca_max(m, w, eps)
    if h_at_0 is None:
        h_at_0 = m * CP_I * eps + w * LF_I
    return (h1_of_T(T_C, m, w, eps, Ca_max, h_at_0)
            - h_of_T_only(T_C, m, w, eps, h_at_0))

# def transient_solve_ncz(z_edges, Z_P, dt, Gamma_P, T_old_C, th_liquid_old,
#                         th_solid_old, iii, eps=1e-4, max_outer=20,
#                         max_inner=20, tol=1e-8):

def transient_solve_ncz(z_edges, Z_P, dt, Gamma_P, T_old_C, th_liquid_old,
                        th_solid_old, iii, eps=1e-4, max_outer=20,
                        max_inner=20, tol=1e-8, tol_E=1e-10, *,
                        bc_u=None, bc_type_u=1, bc_d=0.0, bc_type_d=2):

    '''
    Nested Newton-Casulli-Zanolli (NCZ) enthalpy solver for refreezing.

    Solves heat diffusion with phase change using the enthalpy formulation
    (Tubini et al. 2021, Eq. 1) discretized implicitly on the FV grid, and
    solves the resulting non-linear system with the nested Newton algorithm
    of Casulli and Zanolli (2010). The enthalpy function h(T) is split via
    a Jordan decomposition into two non-negative, non-decreasing pieces
    h1, h2 (Eq. 19), to which Newton's method is applied in a nested
    outer/inner iteration (Eqs. 21-26). This avoids the non-monotonic
    apparent heat capacity that causes Picard and Newton-Raphson to stall
    or cycle indefinitely (paper Fig. 5), and carries a proof of
    convergence for any time step size.

    Contrast with the other refreezing solvers in this module:
    transient_solve_enthalpy (Picard + big_slope + overshoot clamping),
    transient_solve_ahc (apparent heat capacity, Picard), and
    transient_solve_decp (operator-split). See the module header.

    Grid / FV convention:
        z_edges, Z_P follow the same convention as _fv_geometry().

    IMPORTANT -- deviation from the paper's enthalpy function:
        Tubini et al. define h(T) assuming liquid water content is a
        single-valued function of temperature (via an SFCC for soil, or
        via Eq. E1 for pure water). That assumption does not hold for
        firn, where LWC is an independent state variable set by surface
        melt and percolation: a firn layer at T=0 may hold anywhere from
        zero to its full mass in liquid. Accordingly, the latent heat term
        here is scaled by the liquid water actually present at the start
        of the step (th_liquid_old), NOT by total H2O mass. A dry layer
        therefore has no latent-heat spike and is solved as pure sensible
        heat. Without this modification the T -> liquid-fraction inversion
        assigns large liquid fractions to layers only microdegrees below
        freezing, driving th_solid (and hence rho) toward zero and
        producing NaNs downstream in melt.py. The modified Ca(T) still
        satisfies the paper's requirements C1 and C2, so the convergence
        guarantee is preserved.

    :param z_edges: layer (control volume) edges [m], length nz_P+1
    :param Z_P: layer (control volume) centers [m], length nz_P
    :param dt: time step [s]
    :param Gamma_P: thermal conductivity at layer centers [W/m/K]
    :param T_old_C: temperature at start of step [deg C], fusion at T=0
    :param th_liquid_old: liquid water mass per total volume at start of
        step [kg/m3]; also serves as the latent-heat carrier (see above)
    :param th_solid_old: solid (ice) mass per total volume at start of
        step [kg/m3]
    :param iii: current model time step index (used in warning messages)
    :param eps: temperature window [deg C] over which the liquid water
        present in a layer refreezes; analogous to the paper's epsilon
        (their Eq. E2, set to 1e-4 in their Sect. 4.1 tests). Unlike
        big_slope in transient_solve_enthalpy, this enters the exact
        Jordan decomposition rather than being an ad hoc stiffness.
        Not yet sensitivity-tested across sites/conditions.
    :param max_outer: cap on outer (h2-linearization) iterations
    :param max_inner: cap on inner (h1-linearization) iterations
    :param tol: convergence tolerance, applied to both loops. The inner loop
        uses it directly on max abs(dT) [deg C]. The outer loop accepts EITHER
        max abs(dT) < tol OR an energy-based test, resid_vol * dt < tol *
        h_scale, where resid_vol is the Eq. 23 residual divided by dZ
        [J/m3/s] and h_scale = max(h_at_0) [J/m3] is the layer enthalpy at
        fusion. Dividing by dZ makes the test grid-independent; the raw
        outer_resid scales with dZ and so would tighten or loosen with the
        mesh.

        The energy test is not a fallback -- it is the operative criterion.
        On the Neumann benchmark it fired on all 525600 outer-loop exits
        across nine (dz, dt) combinations (dz 0.01/0.005/0.001 m, dt
        3600/300/60 s) and the temperature test fired on none: when
        Ca ~ w * LF_I / eps (order 1e12 at eps=1e-4), a 1e-8 deg C
        tolerance is unreachable because the enthalpy balance is satisfied
        long before the temperature iterate stops moving at that level.
        Front-position errors are identical to the temperature-only version
        (0.00018-0.00373 m over the nine cases), but outer_count falls to
        2-3 (mean 2.00-2.83) from up to max_outer=20, and total inner
        iterations from ~41 to ~5. The extra passes were chasing a
        criterion carrying no information about the solution.

        outer_count was never 1 in any of those calls, so the outer
        (h2-linearization) loop is always exercised at least twice and the
        Jordan decomposition is not short-circuited by an early energy exit.
        Worth re-checking on any new problem: a first-pass exit would reduce
        this solver to a single linearization and forfeit the convergence
        guarantee.

        Tubini et al. likewise use an energy-based tolerance, rescaled by
        RHO_W_KGM * LF_I. Scaling by the layer's own h_at_0 instead
        generalizes to firn, where w varies by orders of magnitude between
        layers and a pure-water scale would be far too loose in dry firn.
        NOT yet tested on a firn run -- see outstanding DYE-2 re-run.

    :return:
        dict with keys::

            TzC_return: updated temperature [deg C]
            th_solid: updated solid mass per total volume [kg/m3]
            th_liquid: updated liquid mass per total volume [kg/m3]
                (th_solid + th_liquid == th_solid_old + th_liquid_old)
            count: total inner iterations summed over all outer iterations
            outer_count: number of outer iterations used
            energy_resid: interior-cell energy residual [J/m2]

    NOTE: does not return claw_mushy/claw_dry (no overshoot clamping is
        needed or applied). Callers must not assume those keys exist.

    NOTE: energy_resid is computed with this function's h_of_T_only, whose
        reference point is h=0 at T=-eps. The other solvers use
        enthalpy_of(), referenced to h=0 for fully-solid firn at T=0.
        The two are NOT directly comparable across solvers.

    VALIDATION (2026-08): verified against the analytical Neumann solution
            of the paper's Sect. 4.1 / Table E1 (semi-infinite pure-water column,
            T0 = +5 C, Dirichlet surface at -5 C, 100 days, front position by
            interpolation of the liquid-fraction profile). Maximum front-position
            error over nine (dz, dt) combinations -- dz in (0.01, 0.005, 0.001) m,
            dt in (3600, 300, 60) s -- was 0.00018-0.00373 m, against the paper's
            reported 0.00153-0.00905 m for the same method. Front growth exponent
            0.500-0.512 (exact: 0.5). No non-convergence in any of the nine runs,
            including dt = 3600 s, consistent with the paper's claim of
            convergence at arbitrary time step.

            Error is dominated by dt, not dz: at fixed dt = 3600 s, refining dz
            from 0.01 to 0.001 m made the error slightly WORSE (0.00273 ->
            0.00373 m), while refining dt at fixed dz improved it by an order of
            magnitude. Refining the grid without also refining the time step is
            not useful here. Note also that these are maxima over the full 100
            days and are dominated by early times, when the front spans only a
            few cells; final-time error at the finest grid was 1e-4 m.

            For comparison on the identical grid (dz=0.01, dt=3600):
            transient_solve_ncz 0.00273 m; transient_solve_enthalpy 0.01376 m;
            transient_solve_decp 0.01380 m (but only with iters >= 1000; 0.02777 m
            at its default iters=10); transient_solve_ahc 0.775 m at best (W=0.1),
            with 1071 of 2400 steps failing to converge.

    CAVEAT -- this solver can return T > 0 without melting ice:
        Earlier versions clamped the outer iterate to T <= 0 on every pass,
        which made T > 0 unreachable. That clamp was wrong: the paper's
        Algorithm 1 (line 11) constrains only the INITIAL guess, and clamping
        every accepted iterate destroys any genuinely warm layer -- it
        annihilated the +5 C column in the Neumann benchmark. The initial
        guess is now set to min(T_old, -eps), strictly below T* = 0 as the
        paper requires (at exactly T = 0 the q branch flips and the inner
        loop oscillates between branches), and accepted iterates are left
        alone.

        Consequence: a layer driven above fusion now reports T > 0, but the
        phase partition is refreeze-only (th_liquid = w * frac_unfrozen,
        bounded above by w), so its ice is NOT melted. Energy is conserved
        within this function's own enthalpy function, but the state is
        physically inconsistent: ice coexisting with T > 0. For CFM firn runs
        this should not arise -- the surface is Dirichlet-pinned by SEB.py at
        T <= 0 and interior layers have no warming source -- but it is not
        prevented, and melt.py is not known to handle it. Worth a diagnostic
        count of (T_final > 0) on any new site or forcing.

    CAVEAT -- non-positive diagonal at exactly T = 0:
        p - q is positive whenever both are evaluated at the same
        temperature, but p is evaluated at the inner iterate while q is
        frozen at the outer iterate. If the outer iterate sits at exactly
        T = 0 (where q_of_T takes its warm branch, q ~ w * LF_I / eps) and
        the inner iterate goes cold, p - q is strongly negative and the
        matrix loses positive definiteness. T = 0 exactly is COMMON in CFM,
        since melt.py sets wet layers there. Not observed to cause failure,
        but a np.any(a_P <= 0) check is left in the inner loop deliberately;
        do not remove it without testing on a wet firn run.

    NOTE on eps: the pure-water benchmark above used eps = 1e-4 C, the value
        used in the paper's Sect. 4.1. Sensitivity of firn results to eps has
        NOT been tested. Note that eps enters an exact Jordan decomposition
        here, unlike big_slope in transient_solve_enthalpy or W in
        transient_solve_ahc, both of which are ad hoc stiffness parameters --
        AHC in particular degrades catastrophically as W is reduced (see its
        docstring), whereas NCZ has no such failure mode.

    References:
        Tubini, N., Gruber, S., and Rigon, R. (2021). A method for solving
            heat transfer with phase change in ice or soil that allows for
            large time steps while guaranteeing energy conservation.
            The Cryosphere, 15, 2541-2568.
            https://doi.org/10.5194/tc-15-2541-2021
        Casulli, V. and Zanolli, P. (2010). A nested Newton-type algorithm
            for finite volume methods solving Richards' equation in mixed
            form. SIAM J. Sci. Comput., 32, 2255-2273.

    IMPLEMENTATION NOTE -- the inner diagonal is (A + P - Q), not (A + P + Q):
            The paper's Eqs. 22 and 25 print the correction operator as (A + Q),
            which reads as a plus sign on the frozen outer-iterate capacity. That
            reading is wrong for this discretization, and it is wrong in a way
            that converges quietly rather than failing loudly.

            With the plus form, collecting terms at the fixed point (where the
            inner and outer iterates coincide) leaves a spurious residual term
            2 * beta * q * T. The iteration still satisfies its own linearized
            system to machine precision, so dT -> 0 and both loops report
            success -- but the converged state solves the wrong equation. Because
            q is enormous inside the phase-change window (q ~ w * LF_I / eps,
            order 1e12 for pure water at eps=1e-4) while the conduction terms
            are order 1e2, the spurious term dominates by ~5 orders and pins
            every layer at or above fusion to T ~ 0. In firn this is invisible:
            wet layers sit at T = 0 anyway, and q = 0 wherever T < 0, so the
            error vanishes identically. It was found only by running the paper's
            own pure-water Neumann benchmark, where the +5 C initial column
            collapsed to 0 C on the first step.

            Diagnostic signature, if this is ever reintroduced: the outer
            residual (Eq. 23) settles at exactly 2 * beta * q * T while
            outer_change falls below tol. Verified numerically -- observed
            outer_resid 5.815e+01 against a predicted 2 * 2.78e-6 * 3.337e12 *
            3.1368e-6 = 5.82e+01.

            With the minus form the fixed point reduces to the intended
            conduction-plus-enthalpy balance, and outer_resid becomes a genuine
            convergence measure rather than a diagnostic-only quantity.
    '''

    dZ, deltaZ_u, deltaZ_d, Gamma_u, Gamma_d = _fv_geometry(z_edges, Z_P, Gamma_P)

    m = th_solid_old + th_liquid_old        # total H2O mass/vol, conserved
    w = th_liquid_old.copy()                # latent-heat carrier, fixed over step
    beta = dZ / dt

    # These depend only on m, w, eps -- all fixed for the whole time step.
    # Hoisted out of the Newton loops to avoid recomputing every iteration.
    Ca_mushy = m * CP_I + w * LF_I / eps
    Ca_max   = np.maximum(Ca_mushy, m * CP_W)
    h_at_0   = m * CP_I * eps + w * LF_I
    # Convergence scale for the energy-based exit test [J/m3]. Fixed for
    # the step, like h_at_0 itself.
    h_scale = max(np.max(h_at_0), 1.0)

    a_U = Gamma_u / deltaZ_u
    a_D = Gamma_d / deltaZ_d
    
    # a_P_diff = a_U + a_D

    T_ref = T_old_C.copy()
    wet = (w > 0.0) & (T_old_C > -eps) & (T_old_C < 0.0)
    T_ref[wet] = 0.0
    Hhat_old = h_of_T_only(T_ref, m, w, eps, h_at_0)

    bc_u_0    = T_old_C[0] if bc_u is None else bc_u
    bc_d_0    = bc_d
    # bc_u_0    = T_old_C[0]
    # bc_type_u = 1
    # bc_d_0    = 0.0
    # bc_type_d = 2

    T_outer = np.minimum(T_old_C.copy(), -eps)   # ensure T^0 <= T*=0, per paper's requirement
    outer_count = 0
    total_inner_count = 0
    inner_maxed = False
    outer_maxed = False
    n_badiag = 0

    for k in range(max_outer):

        Q_prev = q_of_T(T_outer, m, w, eps, Ca_max)
        h2_prev = h2_of_T(T_outer, m, w, eps, Ca_max, h_at_0)

        T_inner = T_outer.copy()
        inner_count = 0

        for l in range(max_inner):

            P_prev = p_of_T(T_inner, m, w, eps, Ca_max, Ca_mushy)
            h1_prev = h1_of_T(T_inner, m, w, eps, Ca_max, h_at_0)

            a_P = a_U + a_D + beta * (P_prev - Q_prev)
            if np.any(a_P <= 0.0):
                n_badiag += 1
                print(f"WARNING: non-positive diagonal, iii={iii} k={k} l={l}")

            b   = beta * (Hhat_old + h2_prev - Q_prev * T_outer
                        - h1_prev + P_prev * T_inner)

            a_U_bc, a_D_bc, a_P_bc, b_bc = _apply_bcs(
                a_U.copy(), a_D.copy(), a_P, b, deltaZ_u, deltaZ_d,
                bc_u_0, bc_type_u, bc_d_0, bc_type_d)

            T_new = solver(a_U_bc, a_D_bc, a_P_bc, b_bc)

            change = np.max(np.abs(T_new - T_inner))

            T_inner = T_new
            inner_count += 1
            if change < tol:
                break

        else:
            inner_maxed = True

        total_inner_count += inner_count
        T_outer_new = T_inner

        # Eq. 23 outer residual, interior rows only (BC rows are overwritten
        # by _apply_bcs and do not satisfy the interior energy equation).
        # Flux terms written as differences to avoid cancellation error.
        h_new = h_of_T_only(T_outer_new, m, w, eps, h_at_0)
        resid = (beta[1:-1] * (h_new[1:-1] - Hhat_old[1:-1])
                 + a_U[1:-1] * (T_outer_new[1:-1] - T_outer_new[:-2])
                 + a_D[1:-1] * (T_outer_new[1:-1] - T_outer_new[2:]))
        outer_resid = np.max(np.abs(resid))   # [J/m3/s], diagnostic only

        outer_change = np.max(np.abs(T_outer_new - T_outer))

        # T_outer = np.minimum(T_outer_new, 0.0)
        T_outer = T_outer_new

        outer_count += 1
        # Volumetric residual [J/m3/s]; dZ divided out so the criterion is
        # grid-independent, unlike outer_resid, which scales with dZ.
        resid_vol = np.max(np.abs(resid / dZ[1:-1]))
        hit_T = outer_change < tol
        hit_E = resid_vol * dt < tol_E * h_scale
        if hit_T or hit_E:
            if not hit_T:
                _diag_log('ncz_energy_exit', iii, resid_vol * dt,
                          outer_count)
            break

    else:
        outer_maxed = True
        print(f"WARNING: NCZ outer loop hit max_outer={max_outer} "
              f"(iii={iii}); final outer_change={outer_change:.3e} "
              f"(tol={tol:.1e}), resid_vol*dt={resid_vol*dt:.3e} "
              f"(tol_E={tol_E*h_scale:.3e}), "
              f"total inner iters={total_inner_count}, "
              f"inner loop also maxed: {inner_maxed}")

    # recover final phase partition from T via h_of_T_only's branches
    T_final = T_outer

    # Liquid remaining is a fraction of the water that was present at the
    # start of the step (w), never of total mass m -- a dry layer stays dry.
    frac_unfrozen = np.clip((T_final + eps) / eps, 0.0, 1.0)
    th_liquid = w * frac_unfrozen
    th_solid  = m - th_liquid

    energy_resid = _energy_residual(T_final, th_solid, th_liquid, T_old_C,
                                    th_solid_old, th_liquid_old,
                                    dZ, Gamma_u, deltaZ_u, dt)

    n_warm = int(np.sum((T_final > 1e-12) & (th_solid > 0.0)))
    if n_warm > 0:
        print(f"NCZ: {n_warm} layers above fusion (iii={iii}), "
              f"max={T_final.max():.3e} C")

    _diag_log('ncz', iii, energy_resid, total_inner_count,
              n_warm=int(np.sum(T_final > 0.0)), n_badiag=n_badiag)

    return dict(TzC_return=T_final, th_solid=th_solid, th_liquid=th_liquid,
                count=total_inner_count, outer_count=outer_count, converged=not (outer_maxed or inner_maxed),
                outer_change=outer_change, outer_resid=outer_resid,
                energy_resid=energy_resid)

################################
### end Tubini method ##########
################################

def transient_solve_EN_old(z_edges, Z_P, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s,
                            mix_rho, c_vol, LWC, mass_sol, dz, ICT, rho_firn,
                            iii=0, max_iter=200):
    '''
    Legacy enthalpy-based refreezing solver (renamed from transient_solve_EN).

    Retained specifically for testing/comparison against transient_solve_enthalpy
    and the other current refreezing solvers, and to exactly reproduce
    MC_2609's enthalpyDiff/transient_solve_EN numerics; not used by default,
    but reachable via config key "meltwater_solver": "legacy" (see
    enthalpyDiff_old in diffusion.py for the corresponding legacy call site).
    # claude, 26/09/10: previously unreachable from any config option and
    # would have crashed on call (enthalpyDiff_old had a NameError on tot_rho,
    # and wasn't importing this function) -- both fixed in diffusion.py, and
    # this function is now wired up via firn_density_nospin.py's
    # meltwater_solver dispatch (26/09/11: was the now-removed "LWC_heat" key).

    Uses Picard iteration with an explicit
    liquid-fraction overshoot correction (0.6 relaxation factor), similar in
    spirit to transient_solve_enthalpy but with a different internal
    bookkeeping approach (tracks g_liq/g_sol volume fractions and H_tot
    directly, rather than using enthalpy_of()/invert_enthalpy()).

    Grid / FV convention:
        z_edges, Z_P follow the same convention as _fv_geometry() in this module.

    :param z_edges: layer (control volume) edges [m], length nz_P+1
    :param Z_P: layer (control volume) centers [m], length nz_P
    :param nt: NOTE currently unused inside this function (iteration count is
        controlled solely by max_iter); retained for call-site compatibility
        with enthalpyDiff_old, which still computes nt based on LWC presence
    :param dt: time step [s]
    :param Gamma_P: thermal conductivity at layer centers [W/m/K]
    :param phi_0: temperature profile at start of step [deg C], fusion at T=0
    :param nz_P: NOTE currently unused inside this function
    :param nz_fv: NOTE currently unused inside this function
    :param phi_s: NOTE currently unused inside this function (upper BC value
        is instead taken from phi_t_old[0] internally each iteration)
    :param mix_rho: NOTE currently unused inside this function
    :param c_vol: NOTE currently unused inside this function; a local value
        (c_vol1 = RHO_I * CP_I) is computed and used instead, ignoring this
        argument entirely
    :param LWC: liquid water volume [m3] per layer at start of step
    :param mass_sol: solid (ice) mass [kg] per layer at start of step
    :param dz: layer thickness [m] (note: distinct from dZ computed
        internally via _fv_geometry from z_edges; should be equivalent in
        practice but is a separate input here)
    :param ICT: NOTE currently unused inside this function (deprecated
        "Iteration Count Threshold", per caller's inline comment)
    :param rho_firn: NOTE currently unused inside this function
    :param iii: current model time step index; used only in the
        non-convergence warning message
    :param max_iter: maximum Picard iterations (default 200)

    :return:
        tuple (phi_t_out, g_liq, count, iterdiff, g_sol)::

            phi_t_out: updated temperature [deg C], clamped to <=0 (fusion) in
                layers with liquid present
            g_liq: updated liquid volume fraction (of ice+liquid volume, porosity
                ignored), unitless
            count: number of Picard iterations used
            iterdiff: difference in total g_liq between the final two iterations
                prior to break (diagnostic only; not otherwise used downstream)
            g_sol: updated solid/ice volume fraction, unitless

    Reference: Voller and Swaminathan (1991), eqs. 31-32;
        Voller, Swaminathan, and Thomas (1990), eq. 61.

    NOTE: H_lat (used in cond2's threshold check) is computed once before
        the iteration loop and never updated per-iteration, while a
        per-iteration H_lat_iter is computed but not used in any condition.
        This may be intentional or may be a latent inconsistency inherited
        from earlier development -- flagged here for awareness, not
        corrected, since the intended physics wasn't confirmed.
    '''

    # phi_t = phi_0.copy()
    LWC_old = LWC.copy()
    phi_in = phi_0.copy()

    vol_S       = mass_sol / RHO_I     # volume_Solid, i.e. volume of the ice (solid) portion of each control volume
    vol_SL      = vol_S + LWC    # volume of solid and liquid in each control volume
    mass_liq    = LWC * RHO_W_KGM  # mass of liquid water in each control
    mass_tot    = mass_liq + mass_sol # total mass (solid +liquid) of each control
    rho_liq_eff = RHO_W_KGM*dz #new 10/31
    g_liq       = LWC / dz    #  use liquid volume fraction of total volume of the control, which will net us the enthalpy/volume
    g_sol       = vol_S / dz #unitless

    H_L_liq = RHO_W_KGM*LF_I #volumetric latent enthalpy [J/m3]

    phi_t = phi_0 # phi_t is just the temperature

    phi_t_old = phi_t.copy() # initial temperature
    g_liq_old = g_liq.copy() # initial liquid fraction
    g_sol_old = g_sol.copy() # initial solid fraction

    itercheck = 0.9
    count = 0

    ### Big H stands for latent enthalpy, little h is sensible enthalpy

    ### H_tot is the sum of latent and sensible enthalpy
    H_lat      = H_L_liq*g_liq # Latent enthalpy for each layer
    H_lat_old  = H_lat.copy()
    H_tot         = phi_t * g_sol * RHO_I * CP_I + H_lat  #total enthalpy, voller 1990b eq 4a
    H_tot_old     = H_tot.copy()
    h_old         = phi_t * g_sol * RHO_I * CP_I
    h_updated     = h_old.copy()

    update_gsol = True

    dZ, deltaZ_u, deltaZ_d, Gamma_u, Gamma_d = _fv_geometry(z_edges, Z_P, Gamma_P)
    ### Gamma has units J/s/m/K (W/m/K)

    for i_time in range(max_iter): # Testing indicates that this should never need this many iterations

        H_tot_iter  = H_tot.copy()
        phi_iter    = phi_t.copy()
        g_liq_iter  = g_liq.copy() #unitless
        g_sol_iter  = g_sol.copy()
        H_lat_iter  = H_lat.copy()
        h_iter      = h_updated.copy()

        #### version with working dt ####
        ### S_C is independent
        S_P = np.zeros_like(Gamma_P)
        S_C = H_L_liq  * (g_liq_old - g_liq_iter) # J/m3/s = W/m3 Latent heat as source term.
        # S_C = H_L_liq  * (g_liq_iter) # J/m3/s = W/m3 Latent heat as source term.
        # S_C = RHO_I * CP_I * phi_iter * (g_sol_old - g_sol_iter) + (RHO_W_KGM * CP_W * phi_iter + H_L_liq)  * (g_liq_old - g_liq_iter)

        D_u = (Gamma_u / deltaZ_u) # [W/m2/K]
        D_d = (Gamma_d / deltaZ_d)

        a_U   = D_u        #* dt # [W/m2/K]
        a_D   = D_d        #* dt # [W/m2/K]

        c_vol1 = RHO_I * CP_I # gets multiplied by g_sol below

        a_P_0 = c_vol1 * dZ / dt # [W/m2/K] (new) Patankar eq. 4.41c, this is b_p in Voller (1990; Eq. 30)

        if update_gsol:
            a_P   = a_U + a_D + a_P_0 * g_sol_iter - S_P * dZ #* dt # check the multiply on the S_P
        else:
            a_P   = a_U + a_D + a_P_0 * g_sol_old - S_P * dZ #* dt # check the multiply on the S_P

        b_0   = S_C * dZ/dt #* dt # [W/m2]
        b     = b_0 + a_P_0 * g_sol_old * phi_t_old # By this phi_t_old has to be in K

        # print(f'b: {b}')

        ###############

        ### Boundary conditions:
        ### type 1 is a specified value, type 2 is a specified gradient
        ### (units for gradient are degrees/meter)
        bc_u_0    = phi_t_old[0]
        bc_type_u = 1
        bc_d_0    = 0
        bc_type_d = 2

        a_U, a_D, a_P, b = _apply_bcs(a_U, a_D, a_P, b, deltaZ_u, deltaZ_d,
                                        bc_u_0, bc_type_u, bc_d_0, bc_type_d)
        #####
        
        phi_t = solver(a_U, a_D, a_P, b) #sensible enthalpy. 0 for layers at freezing (have LWC), negative for dry layers
        # print(f'phi_t: {phi_t}')
        #####
        
        '''
        ####
        The crux is to adjust liquid fraction and temperture field based on solution
        Calculations are (partially) based on the fact that the freezing temp is 0,
        which means that if H_tot<0 there is no liquid and you can calculate temperature, and if H_tot>0 there is liquid and T is 0.

        Note previous ways of solving in dev branch and previous releases.
        ###
        '''

        ### The best way to solve: calculate new g_liq, and apply overshoot correction on g_liq
        ### Break loop if iteration (prior to overshoot) is the same as previous solution

        h_updated = phi_t * CP_I * RHO_I * g_sol_old # updated sensible enthalpy after solver. g_sol is volume_solid/dz
        delta_h = h_updated - h_old # change in sensible enthalpy, relative to initial (not iteration)

        ### Figure out what delta_h and g_liq should be based on different conditions
        cond0 = ((delta_h>0) & (LWC_old>0)) # Layers where there was sensible enthalpy increased and there is liquid water
        delta_h[cond0] = 0 # sensible enthalpy should not increase if LWC present (should either stay at 0C or cool down)

        ndh = -1*delta_h #negative delta_h (which makes it positive), makes corrections below easy

        ### everything refreezes if the calculated change in enthalpy is greater than the latent enthalpy
        cond1 = ((ndh>=H_lat_old) & (g_liq_old>0)) #layers where energy change is larger than the needed to refreeze, and where there is water
        H_tot[cond1] = (delta_h[cond1] + H_tot_old[cond1]) # total enthalpy in those layers is change+total, should be net negative
        g_liq[cond1] = 0 #no liquid left

        ### partial refreezing if the delta_h is less than the latent enthalpy
        cond2 = ((ndh<H_lat) & (g_liq_old>0))
        H_tot[cond2] = 0
        g_liq[cond2] = (H_lat_old[cond2] + delta_h[cond2])/H_L_liq #remaining liquid

        ### Make sure that there is no liquid in layers that did not have liquid at start
        cond3 = (g_liq_old<=0)
        g_liq[cond3] = 0
        H_tot[cond3] = h_updated[cond3]

        ### if this iteration gave the same solution as the last iteration, break
        iterdiff = (np.sum(g_liq_iter) - np.sum(g_liq)) # Deprecated? Used to use to calulate time to break loop
        if ((np.allclose(phi_iter,phi_t,rtol=1e-4,atol=1e-3))):
            break
        elif ((np.allclose(g_liq_iter,g_liq,rtol=1e-4,atol=1e-4))):
            break

        ### otherwise apply overshoot correction to g_liq and iterate again.
        delta_g_liq = g_liq - g_liq_iter # change in liquid fraction. Should always be negative.
        g_liq = g_liq_iter + 0.6*delta_g_liq
        g_liq[g_liq<0]=0
        g_liq[g_liq>1]=1
        ################

        g_sol = g_sol - delta_g_liq * (RHO_W_KGM / RHO_I)  # convert liquid-fraction change to solid-fraction change

        ## Now update temperatures after liquid corrections
        phi_t[H_tot>=0] = 0 # H_tot>0 means liquid present, T=0
        phi_t[H_tot<0] = H_tot[H_tot<0] / (CP_I * RHO_I * g_sol_old[H_tot<0]) # H_tot<0 means no liquid; all enthalpy is sensible

        phi_t[g_liq>0] = 0
        #############

        count += 1

        ### END ITERATION LOOP
        ######################
    else:
        print(f"WARNING: enthalpy solver did not converge in {max_iter} iterations (iii={iii if 'iii' in locals() else '?'})")

    phi_t_out = phi_t
    phi_t_out[g_liq>0] = 0
    phi_t_out[(phi_t_out>0)] = 0

    return phi_t_out, g_liq, count, iterdiff,g_sol

###################################
### end transient_solve_EN ########
###################################

#################################
### Functions below are for firn air
### Works, but consider to be in beta

def w(airdict, z_edges, rho_edges, Z_P, dZ):
    '''
    Function for downward advection of air and also calculates total air content.
    '''
    if airdict['advection_type']=='Darcy':
        por_op_edges=np.interp(z_edges,airdict['z'],airdict['por_op'])
        T_edges = np.interp(z_edges,airdict['z'],airdict['Tz'])
        p_star = por_op_edges * np.exp(M_AIR *GRAVITY*z_edges/(R*T_edges))
        dPdz = np.gradient(airdict['air_pressure'],airdict['z'])
        dPdz_edges=np.interp(z_edges,airdict['z'],dPdz)

        # perm = 10.0**(-7.29) * por_op_edges**3.71 # Adolph and Albert, 2014, eq. 5, units m^2
        # perm = 10.0**(-7.7) * por_op_edges**3.4 #Freitag, 2002
        perm = 10.0**(-7.7) * p_star**3.4 #Freitag, 2002
        visc = 1.5e-5 #kg m^-1 s^-1, dynamic viscosity, source?
        flux = -1.0 * perm / visc * dPdz_edges # units m/s
        # w_ad = flux / airdict['dt']  / por_op_edges # where did I get this?
        w_ad = flux / p_star / airdict['dt']
        # w_ad = flux / por_op_edges / airdict['dt']

    elif airdict['advection_type']=='Christo':
        por_tot_edges       = np.interp(z_edges,Z_P,airdict['por_tot'])
        por_cl_edges        = np.interp(z_edges,Z_P,airdict['por_cl'])
        por_op_edges        = np.interp(z_edges,Z_P,airdict['por_op'])
        w_firn_edges        = np.interp(z_edges,Z_P,airdict['w_firn']) # units m/s
        T_edges             = np.interp(z_edges,Z_P,airdict['Tz'])
        p_star              = por_op_edges * np.exp(M_AIR *GRAVITY*z_edges/(R*T_edges))
        dscl                = np.gradient(por_cl_edges,z_edges)
        C                   = np.exp(M_AIR*GRAVITY*z_edges/(R*T_edges))

        op_ind              = np.where(z_edges<=airdict['z_co'])[0] #indices of all nodes wiht open porosity (shallower than CO)
        op_ind2             = np.where(z_edges<=airdict['z_co']+20)[0] # a bit deeper
        co_ind              = op_ind[-1]
        cl_ind1             = np.where(z_edges>airdict['z_co'])[0] #closed indices
        cl_ind              = np.intersect1d(cl_ind1,op_ind2)

        # print('depth co_ind',z_edges[co_ind])

        Xi                  = np.zeros((len(op_ind2),len(op_ind2)))
        Xi_up               = por_op_edges[op_ind2]/np.reshape(por_op_edges[op_ind2], (-1,1))
        Xi_down             = (1 + np.log( np.reshape(w_firn_edges[op_ind2], (-1,1))/ w_firn_edges[op_ind2] ))
        Xi                  = Xi_up / Xi_down # Equation 5.10 in Christo's thesis; Xi[i,j] is the pressure increase (ratio) for bubbles at depth[i] that were trapped at depth[j]

        integral_matrix     = (Xi.T*dscl[op_ind2]*C[op_ind2]).T
        integral_matrix_sum = integral_matrix.sum(axis=1)

        p_ratio_t           = np.zeros_like(op_ind2)
        p_ratio             = np.zeros_like(z_edges)
        p_ratio[op_ind]         = integral_matrix_sum[op_ind]   #5.11
        p_ratio[cl_ind]         = p_ratio[co_ind]*Xi[cl_ind, co_ind] # 5.12
        p_ratio[cl_ind[-1]+1:]  = p_ratio[cl_ind[-1]]

        flux                = w_firn_edges[co_ind-1] * p_ratio[co_ind-1] * por_cl_edges[co_ind-1]

        velocity            = np.minimum(w_firn_edges ,((flux + 1e-10 - w_firn_edges * p_ratio * por_cl_edges) / ((por_op_edges + 1e-10 * C))))

        # velocity            = (flux + 1e-10 - w_firn_edges * p_ratio * por_cl_edges) / ((por_op_edges + 1e-10 * C))
        # velocity = flux / p_star# / airdict['dt']

        # w_ad = velocity
        w_ad              = (velocity - w_firn_edges)

        # w_ad[w_ad>0] = 0

        # w_ad[co_ind:+1] = 0

        # veldiff = velocity

    elif airdict['advection_type']=='zero':
        w_ad = np.zeros_like(rho_edges)

    return w_ad


def A(P):
    '''Power-law scheme, Patankar eq. 5.34'''
    A = np.maximum( (1 - 0.1 * np.abs( P ) )**5, np.zeros(np.size(P) ) )
    return A

def F_upwind(F):
    ''' Upwinding scheme '''
    F_upwind = np.maximum( F, 0 )
    return F_upwind


# def w(z_edges,rho_edges,por_op,T,p_a,por_tot,por_cl,Z_P,dz, w_firn): # Function for downward advection of air and also calculates total air content.

#     por_tot_edges=np.interp(z_edges,Z_P,por_tot)
#     por_cl_edges=np.interp(z_edges,Z_P,por_cl)
#     por_op_edges=np.interp(z_edges,Z_P,por_op)
#     teller_co=np.argmax(por_cl_edges)
#     # w_firn_edges=Accu*rho_i/rho_edges #Check this - is there a better way?
#     w_firn_edges=np.interp(z_edges,Z_P,w_firn)

#     # if ad_method=='ice_vel':
#     #     w_ad=w_firn_edges
#     #     trapped = 0.0
#     #     bubble_pres = np.zeros_like(z_edges)


#     ### Christo's Method from his thesis (chapter 5). This (maybe) could be vectorized to speed it up.

#     bubble_pres = np.zeros_like(z_edges)
#     # print(len(np.diff(por_cl)))
#     # print(len(dz))
#     # dscl = np.append(0, np.diff(por_cl)/dz)
#     dscl = np.append(0, np.gradient(por_cl,dz))
#     T_edges = np.interp(z_edges,Z_P,T)
#     C=np.exp(M_AIR*GRAVITY*z_edges/(R*T_edges))
#     strain = np.gradient(np.log(w_firn),dz)
#     s=por_op_edges+por_cl_edges

#     for teller1 in range (0,teller_co+1):
#         integral = np.zeros(teller1+1)
#         integral2 = np.zeros(teller1+1)

#         for teller2 in range(0,teller1+1):
#             # integral[teller2] = dscl[teller2]*C[teller2]*(s[teller2]/s[teller1])/(1+scipy.integrate.trapz(strain[teller2:teller1+1],dz)) #need to get this indexing correct 6/19/14: I think it is fine.
#             integral[teller2] = dscl[teller2]*C[teller2]*(s[teller2]/s[teller1])/(1+scipy.integrate.trapz(strain[teller2:teller1+1],z_edges[teller2:teller1+1])) #need to get this indexing correct 6/19/14: I think it is fine.
#             if dscl[teller2]==0:
#                 dscl[teller2]=1e-14
#             integral2[teller2] = dscl[teller2]

#         bubble_pres[teller1] = (np.mean(dz)*np.sum(integral))/(np.mean(dz)*np.sum(integral2))

#     bubble_pres[teller_co+1:] = bubble_pres[teller_co]*(s[teller_co]/s[teller_co+1:])/(w_firn_edges[teller_co+1:]/w_firn_edges[teller_co])

#     bubble_pres[0] = 1
#     #print 'bubble pressure = %s' % bubble_pres

#     flux= w_firn_edges[teller_co]*bubble_pres[teller_co]*por_cl[teller_co]

#     velocity = np.minimum(w_firn_edges ,((flux+(1e-10)-w_firn_edges*bubble_pres*por_cl_edges)/((por_op_edges+1e-10)*C)))
#     #velocity = velocity * 2
#     w_ad=velocity


#     return w_ad #, bubble_pres
