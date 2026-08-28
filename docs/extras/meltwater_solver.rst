Choosing a meltwater solver
----------------------------

When **MELT** is on, the CFM has to solve heat diffusion with phase change: liquid water refreezing releases latent heat, which couples the temperature and liquid-water-content fields. The CFM offers four interchangeable numerical schemes for this problem, selected with the **meltwater_solver** json key (see :doc:`../running/json`). All four solve the same physics; they differ in numerical formulation, accuracy, and robustness.

Solver choice involves a genuine trade-off, and ranking the four on a single
number is misleading. Two properties matter, and they are independent:

**Freezing-front accuracy** -- how well the scheme places a refreezing front in
space and time. Relevant if you care about where a refrozen lens forms within a
melt season.

**Energy conservation** -- whether the scheme creates or destroys energy over
many steps. Relevant if you care about cumulative firn temperature, firn air
content, or multi-decadal SMB.

The orderings differ. ``decp`` is the least accurate of the three working
solvers on the front, and among the best at conserving energy. ``ncz`` is the
reverse. Choose on the axis that matters for your application.

Quick recommendation
=====================

For a new run, use ``ncz`` if front accuracy is your priority, or leave the
default ``enthalpy`` if cumulative energy conservation is. Do not use ``ahc``.

The default is ``enthalpy``. After a 2026-08 correction it conserves energy to
round-off, and most of its apparent accuracy disadvantage relative to ``ncz``
is a conductivity-coupling artifact rather than a property of the scheme.

.. warning::

   The ``enthalpy`` solver changed in 2026-08: a redundant temperature-snapping
   step that injected spurious energy was removed. Runs made before that change
   are **not** bit-reproducible with current CFM, and refreezing totals differ.
   Front-position accuracy is unaffected. If you need to reproduce the old
   behavior, ``diffusion.py`` retains it as ``enthalpyDiff_old`` (calling
   ``transient_solve_EN_old``); it is not exposed via ``meltwater_solver`` and
   must be called directly in place of ``refreezeDiff`` in
   ``firn_density_nospin.py``.

The solvers
============

``ncz``
  Nested Newton-Casulli-Zanolli enthalpy method (Tubini et al., 2021). Solves
  the enthalpy formulation directly using a nested Newton algorithm, which
  avoids the non-monotonic apparent heat capacity that causes simpler iteration
  schemes to stall or cycle. Most accurate of the four on the analytical
  benchmark, and the only one with a convergence proof independent of time step
  size -- verified here at time steps from 60 s to one day with no failures.

  Its energy conservation is less good than ``enthalpy`` or ``decp`` by several
  orders of magnitude, because the liquid/solid split is recovered from
  temperature across a finite window ``eps``: a layer just below fusion retains
  a small spurious amount of liquid. The resulting error is systematically
  signed -- ``ncz`` consistently under-freezes -- so it accumulates over long
  runs rather than cancelling.

  This is controllable. Narrowing ``eps`` reduces the conservation error in
  direct proportion, with no measurable change in front accuracy and no
  stability or convergence penalty. Consider ``eps`` smaller than its 1e-4
  default for multi-decadal runs.

``enthalpy``
  Enthalpy formulation with Picard iteration, a large effective heat capacity in
  mushy layers, and explicit overshoot clamping. This is the scheme the CFM used
  before the other three were added, and it remains the default.

  Conserves energy to round-off. Its front-position error is roughly 5x larger
  than ``ncz`` under temperature-keyed conductivity, but most of that gap is
  conductivity feedback, not the phase-change treatment: with conductivity keyed
  on liquid fraction the two schemes agree to all printed digits.

``decp``
  Decoupled/operator-split scheme: each sub-step diffuses heat with no latent
  term, then explicitly refreezes liquid water in any layer left below fusion
  temperature.

  Conserves energy to round-off, at any ``iters``, because its latent-heat step
  is explicit bookkeeping that moves energy between reservoirs and cannot leak.
  Its error is entirely in the *timing* of that transfer, which is what
  sub-stepping controls. Front-position error falls monotonically with ``iters``
  toward a floor equal to ``enthalpy``'s -- the two are the same discrete
  solution once splitting error is removed. The default ``iters=10`` sits at
  roughly twice its own converged error while performing ten tridiagonal solves
  per call; either ``iters=1`` (cheap, larger splitting error) or ``iters>=100``
  (converged) is a more coherent choice.

``ahc``
  Apparent heat capacity method: folds latent heat into an effective heat
  capacity smeared over a fixed temperature window ``W``.

  **Not recommended.** It fails the analytical benchmark at every window width
  tested, frequently fails to converge, and destroys energy at a rate comparable
  to the refreezing signal itself. Retained for comparison only.

  Narrowing ``W`` makes it worse, not better: a layer can cool clean through the
  window within one time step, so the latent-heat term is never sampled and
  never applied, and errors saturate at the value obtained by ignoring latent
  heat entirely. Note that the front still grows as t^0.52 in these failed
  cases, so confirming square-root growth does **not** verify that a scheme is
  handling latent heat at all.

Why keep four solvers around
==============================

These schemes were originally developed and validated for frozen soil (apparent heat capacity) or pure water (the NCZ enthalpy method), where liquid water content is a single-valued function of temperature. That assumption does not hold in firn: liquid water content is set independently by surface melt and percolation, so a firn layer sitting at 0 C may hold anywhere from none to all of its mass as liquid. Adapting each method to firn required a firn-specific correction -- scaling the latent-heat term by the liquid water actually present, rather than by total mass -- described in detail in the ``solver.py`` module docstring, along with the failure mode that occurs when it's missed.

Keeping all four side by side let us check that correction, and the solvers' relative accuracy and conservation properties, against each other rather than trusting any one implementation blind.

Benchmark and caveat
======================

The accuracy comparison above comes from the analytical Neumann problem in Tubini et al. (2021), Sect. 4.1: a semi-infinite pure-water column, initially +5 C, with a -5 C Dirichlet surface, run 100 days, where the freezing-front position has a closed-form solution. At dz = 0.01 m, dt = 3600 s (exact front position 0.6805 m after 100 days):

- ``ncz``: 0.00273 m (0.00018 m at dz=1mm, dt=60s)
- ``enthalpy``: 0.01376 m
- ``decp``: 0.01380 m at iters=1000; 0.02777 m at the default iters=10
- ``ahc``: 0.775 m, best case (W=0.1); 1071 of 2400 steps failed to converge

This benchmark is pure water, where liquid content *is* a single-valued function of temperature -- exactly the condition firn violates (see above). It is a fair test of each scheme's treatment of latent heat and its convergence behavior, but it does not test the firn-specific correction. Treat the ordering above as evidence that it carries over to firn, not proof.

Full implementation details, caveats, and references for each solver are documented in the ``transient_solve_*`` function docstrings in ``solver.py``.

Reference:
Tubini, N., Gruber, S., and Rigon, R. (2021). A method for solving heat transfer with phase change in ice or soil that allows for large time steps while guaranteeing energy conservation. *The Cryosphere*, 15, 2541-2568. https://doi.org/10.5194/tc-15-2541-2021
