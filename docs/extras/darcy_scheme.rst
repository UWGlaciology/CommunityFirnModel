The Darcy percolation scheme
------------------------------

Overview
========

The Darcy scheme (``darcyscheme()`` in ``melt.py``) is a physically detailed
alternative to the bucket scheme (:doc:`bucket_scheme`) for routing surface
meltwater and rainfall through a firn column, selected with
**liquid** = ``darcy`` (see :doc:`../running/json`). Rather than treating
each layer as a discrete reservoir with a storage capacity, it solves for
actual water flux between adjacent layers using an unsaturated Darcy-flow
formulation, following Hirashima et al. (2010), with van Genuchten
unsaturated-flow parameters and an upwind/downwind conductivity scheme
(Szymkiewicz, 2009) at layer interfaces.

Because Darcy flow can require much finer time resolution than a single
outer model time step to remain numerically stable, this scheme runs its
own internal adaptive sub-stepping loop within each outer model time step --
expect it to be noticeably slower than the bucket scheme.

Originally coded by Vincent Verjans. Supporting functions (hydraulic
conductivity, van Genuchten parameterizations, flux root-finders, lateral
runoff) live in ``darcy_funcs.py``.

Algorithm outline
==================

For each outer model time step, ``darcyscheme()`` advances an internal timer
from 0 to ``self.dt[iii]`` in adaptive sub-steps (``dtsub``, bounded by
``dtmin``/``dtmax``). At each sub-step, in order:

1. **Surface melting and regridding** -- mass is supplied at a constant rate
   (``meltflux_mass``) over the sub-step; nodes are consumed from the
   surface downward and the grid is regridded, analogous to
   ``bucket()``'s surface-melting logic.

2. **Ice lens detection** -- impermeable barriers are (re-)identified only
   when the density profile has changed since the last sub-step (a caching
   optimization), using the same ``_find_impermeable_layers`` helper as
   ``bucket()`` (see :doc:`bucket_scheme`), with
   ``force_bottom_impermeable=False``.

3. **Density capping for flow calculations** -- node density is capped at
   **rholim_dcy** for the purposes of computing flow properties only,
   giving ice lenses a small residual permeability rather than an abrupt
   wall of zero porosity. Any leftover LWC in these density-modified nodes
   is forced to run off at the end of the routine.

4. **Surface liquid input distribution** -- melt and rain input for this
   sub-step is distributed downward into available accommodation space,
   based on physical pore space (not combined refreeze+retention capacity
   as in ``bucket()``).

5. **Effective saturation and hydraulic properties** -- van Genuchten shape
   parameters, saturated and relative hydraulic conductivity, and pressure
   head are computed for each node.

6. **Flow-eligibility masking** -- nodes that cannot have downward outflow
   are identified: dry nodes, impermeable nodes, nodes directly above a
   barrier, and nodes in a fully-saturated column above a barrier (a
   perched water table).

7. **Iterative flux solving** -- for each remaining candidate interface,
   the limiting flux ("qlim", Hirashima et al. 2010 Eq. 20) is solved for
   via bisection or Newton-Raphson, using an analytic or warm-started
   initial guess. See "The equilibrium-flux-solving system" below.

8. **Flux application and lateral runoff** -- actual water transfer this
   sub-step is computed via exponential relaxation toward the limiting
   flux (Eq. 23) and applied to update LWC; lateral runoff is optionally
   computed for water with no downward outflow path.

9. **Refreezing** -- applied based on each node's cold content and
   available pore volume.

10. **Adaptive time-stepping** -- the sub-step size grows or shrinks based
    on how close the resolved flux came to the limiting flux, then is
    clamped to stay within ``[dtmin, dtmax]`` and to not overshoot the
    outer time step.

After the sub-stepping loop completes, water-balance and dry-cold-firn
sanity checks are performed once, followed by a runtime-logging check.

Physical & mathematical formulation
=====================================

Darcy flux
~~~~~~~~~~

Water flux between adjacent nodes follows Darcy's law: flux is proportional
to hydraulic conductivity times the gradient of total head (pressure head
plus gravitational head). In code::

    q0 = bigkedg * (dhdz + 1)

where ``dhdz`` is the pressure-head gradient between node centers, and the
``+ 1`` term represents the gravitational contribution (water flows
downward under gravity even absent a pressure-head gradient).

Hydraulic conductivity
~~~~~~~~~~~~~~~~~~~~~~~

Conductivity at each node is the product of:

- **Saturated conductivity** (``bigk_s``, via ``hydrconducsat_Calonne``) --
  a function of grain size and density.
- **Relative conductivity** (``bigk_r``, via ``krel_vG``) -- a van
  Genuchten-Mualem relative permeability curve, a function of current
  effective saturation.

Upwind/downwind interface conductivity (Szymkiewicz, 2009)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At each interface between two nodes, the conductivity used to compute flux
(``bigkedg``) is chosen based on the *sign* of the total head gradient: if
flow is driven downward, the upper node's conductivity is used; if the
gradient would drive flow upward against gravity, the lower node's
conductivity is used instead. This avoids numerical instability that a
naive central-conductivity averaging scheme would introduce.

Effective saturation and van Genuchten parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Volumetric water content is rescaled onto a 0-1 effective-saturation scale
(``theta_e``) between the irreducible content (``theta_i``) and saturated
content (``theta_s``), following Hirashima et al. (2010) Eq. 5. Van
Genuchten shape parameters (``avG``, ``nvG``, ``mvG``, via
``vG_Yama_params``) govern the pressure head (``phead_vG``) and relative
conductivity (``krel_vG``) as functions of this effective saturation.

The limiting flux ("qlim") and iterative solving
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than applying the instantaneous Darcy flux (``q0``) directly over
each sub-step (which could overshoot equilibrium and cause unphysical
oscillation), the scheme first solves for a *limiting* flux (``qlim``,
i.e. ``glw`` in code) at each interface -- the flux that would bring the
pressure heads on both sides of the interface into equilibrium (Hirashima
et al. 2010, Eq. 20). This is solved iteratively via bisection or
Newton-Raphson, using either an analytic first guess
(``thetaeff_equaliser``) or a warm-started guess from the previous
sub-step's flux if flux has been stable.

The *actual* flux applied over the sub-step (``qstep``) is then computed
via an exponential relaxation toward this limiting flux (Hirashima et al.
2010, Eq. 23)::

    qstep = qlim * (1 - exp(-dtsub * q0 / qlim))

This asymptotically approaches ``qlim`` as ``q0 * dtsub`` grows, without
ever overshooting it -- a numerically stable alternative to either a naive
``q0 * dtsub`` estimate or jumping directly to ``qlim``.

The equilibrium-flux-solving system
=====================================

Several functions in ``darcy_funcs.py`` work together to solve, at each
interface between two nodes, for the water flux that brings both nodes'
pressure heads into equilibrium (accounting for gravity) -- this is the
"qlim" concept from Hirashima et al. (2010), Eq. 20.

**The equilibrium residual**, referred to throughout as ``f_eq``, is
defined by moving all terms of Eq. 20 to one side::

    f_eq = hd[upper] - hd[lower] - dltz

where ``hd`` is pressure head (``phead_vG``) and ``dltz`` is the distance
between the two nodes' centers. At equilibrium, ``f_eq == 0``.

Solving process, as orchestrated by ``darcyscheme()``:

1. An initial guess for the flux is made -- either analytically
   (``thetaeff_equaliser``, which finds the flux that equalizes
   *saturation* rather than *head*, a cheap approximation) or via a
   warm-start from the previous sub-step's flux, if flux has been stable.
2. That guess is evaluated: ``thetae_update`` computes the resulting
   saturation at both nodes if the guess were applied, then ``phead_vG``
   converts that to pressure head, and ``f_eq`` is computed.
3. If the initial guess's residual is large, ``flux_bisection`` is used to
   refine it (robust, guaranteed to converge, but slower). If the residual
   is already small, ``flux_newtonraphson`` is used instead (faster
   convergence near the root).
4. ``flux_newtonraphson`` internally uses ``dfdg_derivative`` to compute
   the analytic derivative of ``f_eq`` with respect to flux, taking Newton
   steps until convergence. If a step makes the residual worse, or the
   derivative becomes very large (risking numerical instability), it
   abandons Newton-Raphson and falls back to ``flux_bisection``.
5. Both root-finders share three physically-motivated early-exit
   conditions, checked before attempting further refinement:

   - The upper node has already run dry (``theta_e < 1e-8``) -- no more
     water can be extracted from it.
   - The lower node is already saturated (``theta_e > 0.95``) -- no more
     water can be pushed into it.
   - The guess has reached the minimum physically allowed flux (~0) --
     this scheme does not represent upward/capillary-rise flow.

This system is deliberately conservative about not overshooting either
node's physical bounds (dryness, saturation) or the zero-flux floor, at the
cost of sometimes converging to an early-exit "good enough" answer rather
than a numerically exact root.

The density-capping trick (``rholim_dcy``)
=============================================

A key design choice in this scheme: rather than treating an ice lens (a
node at or above **RhoImp_darcy**, e.g. 873 kg/m\ :sup:`3`) as having
exactly zero porosity -- which would make Darcy flow calculations
numerically unstable (division by zero, infinite pressure head) -- node
density is *artificially capped* at **rholim_dcy** (e.g. 910 kg/m\
:sup:`3`) for the purposes of computing flow properties only. This gives
ice lenses a small but nonzero residual porosity and conductivity, avoiding
a hard numerical wall while still representing them as much less permeable
than typical firn.

The nodes whose density was modified this way are tracked (``rhomod``), and
at the very end of the sub-stepping loop, any LWC still present in those
nodes is forced to run off -- since their *true* density does not actually
support long-term water storage; the density relaxation was a numerical
convenience for flow calculations, not a physical storage allowance.

**RhoImp_darcy** must always be set below **rholim_dcy** -- the scheme
checks this at startup and exits with an error if violated.

Configuration options
======================

These are the Darcy-scheme-relevant json keys; see :doc:`../running/json`
for types and defaults. ``ColeouLesaffre`` and ``IrrVal`` (shared with the
bucket scheme, see :doc:`bucket_scheme`) are also read by this scheme.

.. list-table::
   :header-rows: 1

   * - Key
     - Description
   * - **dtmin_darcy**
     - Minimum internal sub-step size.
   * - **dtmax_darcy**
     - Maximum internal sub-step size.
   * - **dtsub_darcy**
     - Initial internal sub-step size (adapted iteratively thereafter).
   * - **RhoImp_darcy**
     - Density threshold above which a node is part of an impermeable ice
       lens. Must be less than **rholim_dcy**.
   * - **ThickImp_darcy**
     - Minimum ice-lens thickness required for impermeability.
   * - **lat_runoff_darcy**
     - If true, compute lateral runoff for nodes with available water but
       no downward outflow path.
   * - **slope_proxy_darcy**
     - Surface slope used in the lateral-runoff parameterization. As with
       ``bucket()``'s **Slope**, best understood as rise/run, inferred
       from context.
   * - **rholim_dcy**
     - Density cap used only for Darcy flow calculations, giving ice
       lenses residual permeability. Must exceed **RhoImp_darcy**.
   * - **eps_cvg_darcy**
     - Convergence tolerance for the equilibrium-head solve (default
       ``0.1e-3``). It was read but overwritten by the hardcoded value in
       earlier versions; it now takes effect.
   * - **runoff_method_darcy**
     - Selects which lateral-runoff parameterization to use (only
       relevant if **lat_runoff_darcy** is true): the Zuo & Oerlemans
       (1996) time-scale formulation (``'ZuoOerlemans'``, via
       ``runoffZuoOerlemans`` in ``darcy_funcs.py``) or a simplified
       Darcy's-law formulation using slope as a proxy gradient
       (``'Darcy'``, via ``runoffDarcy``). An unrecognized value causes
       the model to exit.

Key variables reference
=========================

.. list-table::
   :header-rows: 1

   * - Variable
     - Meaning
   * - ``timer``, ``timetot``, ``dtsub``
     - Internal sub-stepping clock [s]: current elapsed time, total
       duration to cover (``self.dt[iii]``), and current adaptive
       sub-step size.
   * - ``meltflux_mass``, ``rainflux_vol``
     - Melt and rain input, expressed as constant rates over the outer
       time step [kg/s], [m w.e./s] -- unlike ``bucket()``, which treats
       melt as a single instantaneous pulse.
   * - ``ncv``
     - Number of nodes (control volumes) -- analogous to ``bucket()``'s
       ``nnd``.
   * - ``rhodcy``
     - Node density capped at **rholim_dcy**, used only for
       flow-property calculations.
   * - ``rhomod``
     - Indices of nodes whose density was modified (``rhodcy != self.rho``)
       -- any LWC remaining in these nodes is forced to run off at the end
       of the routine.
   * - ``imp``, ``imptop``
     - Indices of impermeable nodes, and the top node of each impermeable
       lens (see ``_find_impermeable_layers``, :doc:`bucket_scheme`).
   * - ``theta_w``, ``theta_s``, ``theta_i``, ``theta_e``
     - Volumetric water content: current, at saturation, irreducible, and
       effective (rescaled 0-1) saturation.
   * - ``LWCav``
     - LWC available for Darcy flow (beyond irreducible content).
   * - ``LWCacm``
     - Remaining accommodation space per node (beyond current LWC, up to
       saturation).
   * - ``i_input_bottom``
     - Deepest node receiving distributed surface liquid input this
       sub-step (analogous to ``bucket()``'s variable of the same name).
   * - ``runoff0``
     - Runoff from surface liquid input that could not be accommodated
       (e.g. an impermeable lens at or near the surface).
   * - ``avG``, ``nvG``, ``mvG``
     - Van Genuchten shape parameters (alpha, n, m), grain-size- and
       density-dependent.
   * - ``bigk_s``, ``bigk_r``, ``bigk``
     - Saturated, relative, and actual hydraulic conductivity [m/s].
   * - ``bigk_d``
     - Hydraulic conductivity staggered down by one node (the "downwind"
       value), used in the upwind/downwind interface selection.
   * - ``hd``, ``dhdz``
     - Pressure head [m] and its gradient between node centers.
   * - ``bigkedg``
     - Conductivity actually used at each interface, selected via the
       Szymkiewicz (2009) upwind/downwind rule.
   * - ``indsdry``, ``i0dry``
     - Dry nodes, and the depth below which all nodes are dry.
   * - ``imp_d``, ``satnofl``, ``satnofl_d``
     - Nodes blocked from outflow: directly above an impermeable barrier,
       in a saturated no-outflow stack above a barrier, and directly
       above such a stack.
   * - ``i00``, ``i11``
     - Nodes that must have zero outflow, and candidate nodes that may
       have outflow > 0 (passed to the flux solver).
   * - ``glw``
     - The limiting flux ("qlim") at each interface, solved for
       iteratively.
   * - ``glwc``, ``glwcacm``, ``gtheta_e``, ``ghd``
     - Working (hypothetical) LWC, accommodation space, effective
       saturation, and pressure head, updated as ``glw`` is resolved
       interface-by-interface, bottom-up.
   * - ``glwflux_d1``, ``glwflux_d2``
     - Resolved flux from the previous sub-step and two sub-steps ago --
       used to warm-start the next guess if flux has been stable.
   * - ``qlim``, ``q0``, ``qstep``
     - Limiting flux, instantaneous Darcy flux, and actual flux applied
       this sub-step.
   * - ``lwcin``, ``lwcout``
     - Per-node inflow and outflow this sub-step.
   * - ``rfout``, ``runoff1``
     - Lateral runoff per node, and its sub-step total.
   * - ``refr_tot``, ``runofftot``
     - Running totals of refrozen water and runoff across the entire
       outer time step (all sub-steps).

Function reference (``darcy_funcs.py``)
=========================================

.. list-table::
   :header-rows: 1

   * - Function
     - Purpose
   * - ``hydrconducsat_Calonne(rad, rho)``
     - Saturated hydraulic conductivity (Calonne et al., Eq. 6).
   * - ``vG_Yama_params(rad, rho)``
     - Van Genuchten shape parameters (alpha, n, m), grain-size- and
       density-dependent (Yamaguchi et al. 2012, Eqs. 6-7).
   * - ``phead_vG(alpha, n, m, thetaeff)``
     - Pressure (matric) head from van Genuchten parameters and
       effective saturation (Hirashima et al. 2010 Eq. 9 / 2014 Eq. 3).
   * - ``krel_vG(m, thetaeff)``
     - Relative hydraulic conductivity (van Genuchten-Mualem; Hirashima
       et al. 2010, Eq. 10).
   * - ``thetae_update(absfl, th_i, th_s, LWC, dz)``
     - Effective saturation resulting from a candidate flux between two
       nodes.
   * - ``thetaeff_equaliser(th_i2, th_s2, LWC2, dz2)``
     - Closed-form flux estimate that equalizes saturation between two
       nodes; used as an initial guess.
   * - ``dfdg_derivative(...)``
     - Analytic derivative of the equilibrium residual with respect to
       flux, for Newton-Raphson.
   * - ``flux_bisection(...)``
     - Bisection root-finder for the equilibrium-limiting flux ("qlim")
       at an interface.
   * - ``flux_newtonraphson(...)``
     - Newton-Raphson root-finder for the equilibrium-limiting flux;
       falls back to bisection if divergent.
   * - ``runoffZuoOerlemans(dt, slope, lwcexcess, inds)``
     - Lateral runoff via Zuo & Oerlemans (1996). Used when
       **runoff_method_darcy** is ``'ZuoOerlemans'``. An equivalent
       formula is duplicated inline in ``bucket()`` rather than calling
       this shared implementation -- a follow-up refactor opportunity,
       left as-is for now.
   * - ``runoffDarcy(dt, slope, bigk, inds)``
     - Lateral runoff via a simplified Darcy's-law formulation using
       slope as a proxy gradient. Used when **runoff_method_darcy** is
       ``'Darcy'`` (the default).
   * - ``vG_Yama(rad, rho, thetaeff)``
     - Combined pressure head + relative conductivity (van Genuchten /
       Yamaguchi et al. 2012). **Not currently called** anywhere in the
       CFM -- superseded by the three functions above, used
       individually. Kept in case it is used elsewhere; verify before
       removing.

Comparison to the bucket scheme
=================================

.. list-table::
   :header-rows: 1

   * - Aspect
     - ``bucket()``
     - ``darcyscheme()``
   * - Time resolution
     - Single instantaneous update per outer time step
     - Adaptive internal sub-stepping within each outer time step
   * - Water movement
     - Discrete "storage capacity" per node; excess passed downward
     - Continuous flux driven by pressure-head gradients and hydraulic
       conductivity
   * - Melt input
     - Single pulse at the start of the time step
     - Constant rate spread across the time step
   * - Ice lens permeability
     - Fully impermeable (binary)
     - Small residual permeability via density capping (**rholim_dcy**)
   * - Irreducible water
     - Coléou & Lesaffre (1998) or fixed fraction
     - Same options, shared config keys
   * - Computational cost
     - Lower (single pass)
     - Higher (iterative flux solving at every sub-step)

Known limitations & caveats
=============================

- **Shared irreducible-water configuration**: **ColeouLesaffre** and
  **IrrVal** use the *same* config keys as ``bucket()``, on the reasoning
  that a user switching between schemes would likely want a consistent
  physical choice for irreducible water content. If this scheme is used
  side-by-side with ``bucket()`` in a way that requires independent
  settings, this sharing would need to be revisited. One consequence:
  enabling **ColeouLesaffre** (the default) changes the Darcy scheme's
  irreducible water content from the old hardcoded 0.02 to a
  density-dependent value roughly 80% higher at typical firn densities --
  the single largest numerical change when adopting this release's
  defaults, since the Darcy scheme previously ignored this key entirely.

- **Tight water-balance tolerance**: the final sanity check uses a
  ``1e-12`` tolerance, tighter than ``bucket()``'s ``1e-3``/``1e-5``.
  Given the many sub-steps typically run per outer time step, this may be
  prone to frequent false-positive warnings from accumulated
  floating-point error; this has not yet been revisited in practice.

- **Gravitational constant**: ``hydrconducsat_Calonne`` sources gravity
  from ``constants.py`` (``GRAVITY = 9.8``) rather than the literal
  ``9.81`` used elsewhere historically. This introduces a ~0.1% change in
  computed saturated hydraulic conductivity, for consistency with the
  constant used elsewhere in the CFM.

- **``vG_Yama`` appears unused** -- see the Function Reference table
  above.

- **``dfdg_derivative``'s "full" parameter naming**: this function's
  parameters are suffixed "full" and are internally sliced into
  upper/lower subsets, suggesting a possible original design intent to
  support computing the derivative across an entire column of interfaces
  at once. In the current codebase, it is only ever called with inputs
  already reduced to a single interface's two nodes. The function is
  correct either way; this is noted for future maintainers, not flagged
  as a bug.

- **Citation discrepancy (resolved, not a bug)**: ``phead_vG``'s
  pressure-head formula cites "Hirashima 2014, Eq. 3," while the identical
  formula in ``vG_Yama`` cites "Hirashima 2010, Eq. 9." Both papers present
  the same/equivalent van Genuchten pressure-head relation at different
  equation numbers -- confirmed as two distinct, correct citations, not a
  copy-paste typo.

References
==========

- Calonne, N. et al. -- saturated hydraulic conductivity parameterization.
- Hirashima, H. et al. (2010). Numerical modeling of liquid water movement
  through layered snow based on new measurements of the water retention
  curve. *Cold Regions Science and Technology*, 64(2), 94-103.
- Hirashima, H. et al. (2014). https://doi.org/10.1016/j.coldregions.2014.09.004
- Szymkiewicz, R. (2009) -- upwind/downwind conductivity selection.
- Yamaguchi, S. et al. (2010, 2012) -- van Genuchten parameter
  parameterization.
- Zuo, Z. and Oerlemans, J. (1996). Modelling albedo and specific balance
  of the Greenland ice sheet: calculations for the Søndre Strømfjord
  transect. *Journal of Glaciology*, 42(141), 305-317.
