Choosing a meltwater solver
----------------------------

When **MELT** is on, the CFM has to solve heat diffusion with phase change: liquid water refreezing releases latent heat, which couples the temperature and liquid-water-content fields. The CFM offers four interchangeable numerical schemes for this problem, selected with the **meltwater_solver** json key (see :doc:`../running/json`). All four solve the same physics; they differ in numerical formulation, accuracy, and robustness.

Quick recommendation
=====================

If you are starting a new run and are not trying to reproduce older CFM output, use ``ncz``. It is the most accurate of the four on the one case where an analytical solution exists, and it is the only one of the four with a convergence proof at arbitrary time step.

The default remains ``enthalpy`` for backward compatibility -- it reproduces the CFM's prior refreezing behavior, so existing configs that don't set **meltwater_solver** are unaffected by the addition of the other three schemes.

The solvers
============

``ncz``
  Nested Newton-Casulli-Zanolli enthalpy method (Tubini et al., 2021). Solves the enthalpy formulation directly using a nested Newton algorithm, which avoids the non-monotonic apparent heat capacity that causes simpler iteration schemes to stall or cycle. Most accurate of the four on the analytical benchmark below, and the only one with a convergence proof independent of time step size.

``enthalpy``
  Enthalpy formulation with Picard iteration, a large effective heat capacity in mushy layers, and explicit overshoot clamping. Reliable and well-tested (this is the scheme the CFM used before the other three were added), but roughly 5x less accurate than ``ncz`` on the benchmark below.

``decp``
  Decoupled/operator-split scheme: each sub-step diffuses heat with no latent term, then explicitly refreezes liquid water in any layer left below fusion temperature. Converges to the same accuracy as ``enthalpy``, but only once the internal sub-stepping (``iters``) is pushed well above its default of 10 -- meaning it needs to be noticeably more expensive than ``enthalpy`` to match it.

``ahc``
  Apparent heat capacity method: folds latent heat into an effective heat capacity smeared over a fixed temperature window. Retained for comparison only. It fails the analytical benchmark at every window width tested, and frequently fails to converge outright. Prefer one of the other three unless you specifically need this method for comparison.

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
