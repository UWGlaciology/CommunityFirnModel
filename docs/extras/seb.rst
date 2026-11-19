.. _seb-page:

The surface energy balance (SEB)
--------------------------------

.. DRAFT NOTES FOR MAX (RST comments are not rendered; search this file for "VERIFY"):
   Written from SEB.py, SEB_turbulent.py, RCMpkl_to_spin.makeSpinFiles (SEB branch), firn_density_nospin.py (the SEB call
   and what is done with its result), diffusion.py (the surface boundary of the heat diffusion), the changelog entries for
   2.0.0 / 3.0.0 / 4.0.0, and the json.rst entries of the SEB keys. Everything below describes what the code does;
   things I could not check against a source are marked VERIFY.

Overview
========

With **SEB** set to true, the CFM calculates the surface temperature and the surface melt from the energy balance
of the snow surface, in the time-stepping loop. Without it, the surface temperature is the input temperature
(``TSKIN``) and the melt is an input too (``SMELT``, for example from a degree-day model, see :doc:`../running/inputs`).
The module is ``SEB.py`` (class ``SurfaceEnergyBudget``), with the turbulent-flux calculations in ``SEB_turbulent.py``.

Each model time step it does the following:

1. It adds up the energy fluxes at the surface: net shortwave radiation, downwelling longwave radiation, the sensible
   and latent turbulent heat fluxes, the heat carried by rain, and the heat conducted from the firn below.
2. It finds the surface temperature at which the energy that goes into the top layer of the firn is balanced by the
   outgoing longwave radiation and by the warming of that layer.
3. If that temperature would be above the melting point, the surface is at 273.15 K instead, and the surplus energy,
   after the cold content of the layer has been used up, melts snow.
4. The temperature of the top layer is set to the surface temperature, the melt is passed to the meltwater scheme
   (see :doc:`bucket_scheme`) in place of a melt forcing, and the surface temperature is the boundary condition for
   the heat diffusion.

The SEB module was introduced in version 2.0.0 and has been developed since (a sub-time-step in 3.0.0, turbulent
fluxes in 4.0.0; see the changelog). **It is still in beta** (see "Status and limitations" below).

Setting it up
=============

* ``"SEB": true`` in the configuration file, with **MELT** true. (If **MELT** is false, the model prints a warning
  and sets it to true, because the SEB computes melt.) The SEB is only available with the forcing as a
  dataframe or as the dictionary ``climateTS`` (see :doc:`../running/inputs`), not with csv files, and only with
  ``"timesetup": "exact"``; the model stops with a message otherwise. **InputFileNamemelt** is not needed (the melt is
  computed), and a ``SMELT`` column in the dataframe is dropped.
* The forcing needs these columns in the dataframe: ``T2m`` (2-m air temperature), ``TSKIN``, ``ALBEDO`` (only with ``albedo_method`` ``prescribed``), ``SW_d``,
  ``LW_d``, ``QH``, ``QL`` (with ``SEB_turb_method`` ``prescribed``; see below), ``RAIN`` and ``BDOT``, and ``SUBLIM`` if
  **SUBLIM** is true. For the other turbulent-flux methods, add ``WS``, ``QV2m`` and ``PS``. Units and signs are in the table of
  :doc:`../running/inputs`; **every column of the dataframe must be one of the recognized names**.
* ``makeSpinFiles`` (see :doc:`../running/running`) resamples the mass fluxes and temperature to the model time step
  and returns, in addition to ``climateTS``, the dictionary ``SEBfluxes``, which holds the energy fluxes at the
  time resolution of the original dataframe. Give both to ``FirnDensityNoSpin``. The ``run_CFM_example.py`` script
  uses the SEB with the example DYE-2 forcing (MERRA-2 energy fluxes, daily).
* The albedo comes from the method that **albedo_method** selects (see :ref:`seb-albedo`). With ``prescribed`` (the default), missing (NaN) values of the
  albedo are replaced by the mean albedo (MERRA-2 gives no albedo when there is no shortwave radiation), and
  **albedo_factor** multiplies the albedo; it is a simple way to tune the melt.
* In SEB mode, the temperature forcing of the model before the energy balance is solved, from which the mean annual
  temperature :math:`T_m` used by the densification models comes, is ``T2m`` and not ``TSKIN``.
* The SEB is run through the whole series, including the repeated reference climate interval of the spin-up (see
  :doc:`../running/spinup`). The spin-up therefore includes the melt and refreezing of the reference climate.


The energy balance
==================

The energy available to the top layer is the sum of the fluxes that are positive towards the surface:

.. math::

    E_{net} = SW_d\,(1-\alpha) + LW_d - \varepsilon_s \sigma T_s^4 + Q_H + Q_L + Q_{rain} + G ,

with :math:`SW_d` and :math:`LW_d` the downwelling shortwave and longwave radiation, :math:`\alpha` the albedo,
:math:`\varepsilon_s = 0.98` the emissivity of snow and :math:`\sigma` the Stefan-Boltzmann constant (so
:math:`\varepsilon_s \sigma T_s^4` is the upwelling longwave radiation; the downwelling radiation is absorbed
completely), :math:`Q_H` and :math:`Q_L` the sensible and latent heat fluxes, :math:`Q_{rain}` the heat of rain,
and :math:`G` the heat conducted from the firn below. All are in W m\ :sup:`-2`. The radiation is absorbed at the
surface; the energy balance has no penetration of shortwave radiation into the snow.

**Turbulent fluxes** (:math:`Q_H`, :math:`Q_L`): positive when they heat the surface, so a sensible heat flux is
positive when the air is warmer than the surface and a latent heat flux is positive when the air is moister than
the surface (condensation or deposition). How they are calculated is described below.

**Heat of rain.** Rain falling on a colder snow surface cools to the melting point and releases heat (Hock and
Holmgren, 2005, their Eq. 19):

.. math::

    Q_{rain} = c_w\,R\,\left(T_{air} - T_{melt}\right) ,

with :math:`R` the rain rate in kg m\ :sup:`-2` s\ :sup:`-1` and :math:`c_w = 4180` J kg\ :sup:`-1` K\ :sup:`-1`. The rain is assumed to be at the air
temperature and the surface at the melting point, and the term is zero if the air temperature (``T2m``) is at or below
the melting point.

**Ground heat flux.** :math:`G` is conduction between the top of the column and the first node at or below 1 m
depth,

.. math::

    G = K\,\frac{T_{1\,m} - T_s}{z_{1\,m}} ,

with :math:`K = K_{ice}\,(\rho/\rho_i)^{2 - 0.5\,\rho/\rho_i}` and :math:`K_{ice} = 9.828\exp(-0.0057\,T)` W m\ :sup:`-1` K\ :sup:`-1`
(the conductivity of ice; Cuffey and Paterson, 2010, Eq. 9.2, after Yen, 1981, Eq. 33), where
:math:`T` and :math:`\rho` are the mean temperature and density of the column above that node. This expression is the
SEB's own; it is not the conductivity of the heat-diffusion module (**conductivity**). (When the turbulent fluxes
are prescribed, the analytic solution described below needs :math:`G` to be independent of :math:`T_s`; there the temperature of the
uppermost node at the beginning of the step is used for :math:`T_s` in :math:`G`.)

.. VERIFY: K_ice is confirmed (Yen 1981, Eq. 33 is in ~/papers/Yen_ThermalSnow_1981.pdf; Cuffey and Paterson Eq. 9.2 per Max).
   The density factor (rho/rho_i)^(2 - 0.5 rho/rho_i) is NOT the Yen (1981) rho^1.885 form that the Yen options of
   **conductivity** use (diffusion.py), and SEB.py gives no source for it. Where does it come from?

The top layer and the surface temperature
=========================================

The energy balance is applied to a single "top layer": the nodes from the surface down to, and including, the first
node whose depth is at least **SEB_TL_thick** (default 0.05 m). The layer has mass :math:`m` and a mean temperature
:math:`T_{TL}` (weighted by mass), and it takes the energy of the surface in an implicit (backward Euler)
step of length :math:`\Delta t`, in which its heat storage is :math:`c_i\,m\,(T_s - T_{TL})/\Delta t`:

.. math::

    c_i\,m\,\frac{T_s - T_{TL}}{\Delta t} = E_{net}(T_s) .

* When :math:`Q_H` and :math:`Q_L` do not depend on :math:`T_s` (``prescribed``), this is a quartic equation in :math:`T_s`
  (from the :math:`T_s^4` term) that is solved analytically (the "fast quartic solver", ``FQS``).
* When they do (``bulk_richardson`` and ``monin_obukhov``), it is solved with Newton-Raphson iterations, with the fluxes
  recalculated at each iterate.

If the solution is at or above 273.15 K, the surface temperature is set to 273.15 K and the energy that is left, after
the layer has been warmed to the melting point (its **cold content**), is used for melting:

.. math::

    M = \max\left(0,\ \frac{\left[E_{net}(273.15\,\mathrm{K})\right]\Delta t - c_i\,m\,(273.15 - T_{TL})}{L_f}\right) ,

the melt :math:`M` in kg m\ :sup:`-2` per step (with :math:`c_i = 2097` J kg\ :sup:`-1` K\ :sup:`-1`
and :math:`L_f = 333\,500` J kg\ :sup:`-1`).

**Thickening the top layer.** A layer cannot melt more than its own mass. If :math:`M` is more than :math:`m`, the layer
is made thicker: the result is the thinnest layer (in whole nodes) that is at least as heavy as the melt it produces. A thicker
layer has more cold content and heat capacity, so its melt can only be smaller and its mass larger, which is why the
solution can be found by bisection. If the melt is more than the mass of the whole column, the whole column is used.

**What is done with the result.** The temperature of every node of the top layer is set to :math:`T_s`. The surface node
temperature is the fixed temperature at the top of the heat diffusion, so below the top layer the energy goes into the
firn by the usual heat diffusion (see the **heatDiff** and **conductivity** keys). New snow added at the surface in
the step has the temperature :math:`\min(T_s, 273.15\,\mathrm{K})`. The melt, converted to m ice equivalent per
year, replaces the melt forcing: it is what the meltwater scheme uses, and it is written to the ``climate`` output and
to the forcing record ``CFMforcing.hdf5``. The surface temperature is the ``Ts`` of the ``climate`` output.

The turbulent fluxes
====================

**SEB_turb_method** selects how :math:`Q_H` and :math:`Q_L` are obtained:

``prescribed`` (default)
    The fluxes are taken from the forcing (``QH`` and ``QL``, for example from MERRA-2). They do not respond to the
    surface temperature that the CFM calculates. This was the only choice before version 4.0.0.

    ``QH`` and ``QL`` must be positive towards the surface (into the snow), like all the other fluxes, **in every
    way of giving the forcing to the model** (before version 4.1.0 the module multiplied them by -1 if ``climateTS``
    was given without ``SEBfluxes``; it no longer does). Forcing products differ in their sign convention: MERRA-2
    fluxes are positive upward, so they must be multiplied by -1 (the example scripts do this when they make the dataframe),
    while ERA5 fluxes have the opposite sign to MERRA-2 and are positive towards the surface. The table in
    :ref:`sign-conventions` lists the conventions of the variables of both products.

    The module checks the sign of ``QH`` (only with this method), and prints a warning if it looks positive upward, for
    example when the multiplication by -1 was forgotten. The check is described in :ref:`sign-conventions`; it is a safety net
    and does not replace knowing the convention of your forcing.

``bulk_richardson``
    A bulk-aerodynamic calculation from the wind speed :math:`U`, the air temperature, the humidity and the
    pressure, with the stability correction as a function of the bulk Richardson number (after Louis, 1979, as used for
    glaciers by Braithwaite, 1995). With the transfer coefficients :math:`C_{h}` and :math:`C_{q}` for neutral conditions
    (from the roughness lengths and the measurement heights) multiplied by a stability factor :math:`F(Ri_b)`,

    .. math::

        Q_H = \rho_a c_{p}\,C_h F\,U\,(T_a - T_s) , \qquad Q_L = \rho_a L\,C_q F\,U\,(q_a - q_s(T_s)) ,

    where :math:`q_s` is the saturation specific humidity over ice at the surface temperature, and
    :math:`L` is the latent heat of sublimation (2.834 :math:`\times` 10\ :sup:`6` J kg\ :sup:`-1`) or, if the surface is at
    the melting point, of vaporization (2.501 :math:`\times` 10\ :sup:`6`). The wind speed has a lower limit of 0.1 m s\ :sup:`-1`.

``monin_obukhov``
    Monin-Obukhov similarity theory, with an iteration on the Obukhov length (up to 20 iterations). The stability
    functions are those of Paulson (1970) for unstable conditions and Beljaars and Holtslag (1991) for stable
    conditions.

Both calculated methods need the forcing variables ``WS`` (wind speed), ``QV2m`` (specific humidity) and ``PS``
(surface pressure), with ``T2m`` as the air temperature. The model stops at the start of the run, with a message that
names them, if they are missing. They use the roughness lengths and the measurement heights set by:

* **SEB_z_ref_wind**, **SEB_z_ref_T**, **SEB_z_ref_q**: the heights of the wind, air temperature and humidity forcing
  (default 2 m, which is right for the 2-m MERRA-2 fields).
* **SEB_z0m**, **SEB_z0h**, **SEB_z0q**: the roughness lengths for momentum, heat and humidity (default
  :math:`10^{-3}` m), or the forcing columns ``z0m``, ``z0h``, ``z0q`` (which take priority).
* **SEB_roughness_mode**: ``fixed`` uses the values above; ``parameterized`` derives the roughness lengths for heat and
  humidity from :math:`z_{0m}` and the friction velocity, after Andreas (1987).

The turbulent fluxes are sensitive to the roughness lengths. The entry of **SEB_roughness_mode** in :doc:`../running/json`
gives the results of tests at DYE-2 against the fluxes of MERRA-2, with advice on the choice. Computing the fluxes
from the surface temperature of the CFM, and not taking them from MERRA-2, gave about 16% less melt at DYE-2 for 1980-1999
(see the changelog for version 4.0.0).

.. _seb-albedo:

The albedo
----------

The shortwave radiation that the surface absorbs is :math:`SW_d (1 - \alpha)`, with :math:`\alpha` the albedo. The module ``albedo.py`` gives
the albedo to the SEB, and the key **albedo_method** selects how it is found. The methods are:

* ``prescribed``: the albedo of the forcing (the column ``ALBEDO``), with the treatment of missing values and **albedo_factor**
  described above.
* ``gardner_sharp``: the physically based parameterization of Gardner and Sharp (2010) (their equations 7 to 12), see below.
* ``greuell_konzelmann``: the albedo from the density of the surface snow of Greuell and Konzelmann (1994), see below.

The model asks for the albedo once for each model time step (an albedo for each SEB sub-step, if the SEB forcing is finer than the model step;
see :ref:`seb-substeps`), from the state of the column at the start of the step. An albedo that depends on the melt or on the surface temperature
therefore uses the values of the last step. The albedo that was used (the mean over the sub-steps) is the output ``albedo`` (see :doc:`../running/outputs`),
which is written whenever **SEB** is true: with a parameterization it is the way to see what the model did.

**Gardner and Sharp (2010).** The broadband albedo of snow is a function of the specific surface area :math:`\hat{S}` (cm\ :sup:`2` g\ :sup:`-1`) of the snow, the
concentration :math:`c` (ppmw) of light-absorbing carbon, the cosine :math:`u` of the solar zenith angle and the cloud optical thickness :math:`\tau`:

.. math::

    \alpha_{\hat{S}} = 1.48 - \hat{S}^{-0.07}, \qquad
    d\alpha_c = \max\left(0.04 - \alpha_{\hat{S}},\; \frac{-c^{0.55}}{0.16 + 0.6\hat{S}^{0.5} + 1.8\,c^{0.6}\hat{S}^{-0.25}}\right), \qquad
    \alpha_c = \alpha_{\hat{S}} + d\alpha_c ,

.. math::

    u' = 0.64\,x + (1-x)\,u, \quad x = \min\left[\left(\frac{\tau}{3u}\right)^{0.5}, 1\right], \qquad
    d\alpha_{u'} = 0.53\,\alpha_{\hat{S}}\,(1-\alpha_c)\,(1-u')^{1.2}, \qquad
    d\alpha_\tau = \frac{0.1\,\tau\,\alpha_c^{1.3}}{(1 + 1.5\tau)^{\alpha_{\hat{S}}}} ,

.. math::

    \alpha = \alpha_{\hat{S}} + d\alpha_c + d\alpha_{u'} + d\alpha_\tau .

It was fitted to a radiative transfer model of snow (Se = 0.01, maximum error 0.06) for :math:`\hat{S}` of 0.07 to 1300 cm\ :sup:`2` g\ :sup:`-1`,
:math:`c` up to 2 ppmw, solar zenith angles up to 85 degrees and :math:`\tau` up to 30. The finite-snow-depth correction of the paper (its equation 13, for a layer of snow over ice) is not included.
How the CFM finds the inputs:

* :math:`\hat{S}` is :math:`3/(\rho_i r_e)` with :math:`\rho_i = 910` kg m\ :sup:`-3` and :math:`r_e` the grain radius, from :math:`\sqrt{r^2}`: **r2s0** (a constant) without the grain growth, and the
  grain radius of each node with it (**physGrain**). It is the mean of :math:`\hat{S}` of the nodes down to **albedo_grain_depth** (default 0.05 m, the top layer of the SEB), weighted by mass.
  Be aware of two things. The grain growth of the CFM is slow: at DYE-2 the radius of the nodes of the top 5 cm grows by about 0.5% (from 0.1000 to 0.1004 mm), so **r2s0** is the value that sets the albedo,
  and the albedo of the method is that of pure, new snow of that radius (0.81 for 0.1 mm), which is higher than the albedo of MERRA-2 (0.79 on average at DYE-2, with the 2-year test run giving about a quarter of the melt).
  Use a larger **r2s0**, or **albedo_carbon**, to tune it. It is also not checked that the grain radius of the CFM is the optical effective radius of the paper.
* The solar zenith angle is calculated from the latitude (**albedo_latitude**, or **lat_val**), the longitude (**albedo_longitude**, or **lon_val**; it only matters for steps shorter than a day) and the time of the forcing.
  The albedo of a step is the mean over the step, weighted by :math:`u` (the shortwave radiation of a clear sky). A step is the interval that starts at its time, as for a daily mean of a dataframe.
* :math:`c` is **albedo_carbon** (default 0) and :math:`\tau` is the column ``TAU`` of the forcing if there is one, else **albedo_cloud_tau** (default 0, clear sky). The cloud optical thickness is not estimated from the radiation.

**Greuell and Konzelmann (1994).** The albedo changes linearly with the density :math:`\rho_{ss}` of the surface snow, from the albedo :math:`\alpha_s` of fresh snow at the density :math:`\rho_s` to the albedo :math:`\alpha_i` of ice at :math:`\rho_i`:

.. math::

    \alpha = \alpha_i + (\rho_{ss} - \rho_i)\,\frac{\alpha_s - \alpha_i}{\rho_s - \rho_i} .

:math:`\alpha_s` and :math:`\alpha_i` are **albedo_snow** and **albedo_ice**. They have no default, because they are tuning constants that depend on the application. :math:`\rho_s` and :math:`\rho_i` are **albedo_rho_snow** (default 300 kg m\ :sup:`-3`)
and **albedo_rho_ice** (default 910 kg m\ :sup:`-3`). The density of the surface is that of the top node, and it is limited to the range from :math:`\rho_s` to :math:`\rho_i`, so that the albedo is between :math:`\alpha_i` and :math:`\alpha_s`.
The cloud dependence of Greuell and Konzelmann (1994) is not included. With a constant surface density (**rhos0**) the albedo is nearly constant, since the top node is new snow at each step with accumulation.

To add a parameterization, write a class in ``albedo.py`` that inherits from ``AlbedoModel`` and calculates the albedo of the steps in ``step()``, and add it to
``ALBEDO_METHODS``. The docstring of the module has the details: what ``step()`` gets (the model state, the forcing) and the other files that list the
options (the validation of the configuration, the example configuration files, and the documentation).

.. _seb-substeps:

Sub-daily forcing (sub-steps)
=============================

The energy fluxes can have a finer time resolution than the model time step; for example, 4-hourly fluxes with a
daily model step capture the diurnal cycle, which matters because melt is a threshold process. ``makeSpinFiles``
keeps the energy fluxes at the resolution of the dataframe (in ``SEBfluxes``), and the ratio of the model step to that
resolution is ``SEBfluxes['dtRATIO']`` (it should be a whole number; it is 1 when they are the same, as with the daily example
forcing). The solver then divides each model step into ``dtRATIO`` sub-steps:

* in each sub-step, the surface temperature and melt are found as above with that sub-step's fluxes. The temperature of
  the top layer is carried from one sub-step to the next, starting from :math:`T_{TL}`;
* the melt of the sub-steps is summed over the model step;
* the surface temperature returned for the model step is the mean over the sub-steps, or 273.15 K if there was melt
  in any of them, and the top layer is set to that temperature.

The conductivity :math:`K` and the temperature at 1 m depth are taken once per model step. The heat of the rain is
included in every sub-step.

Status and limitations
======================

* **The SEB is in beta**, and the turbulent-flux methods are new in version 4.0.0. The stability functions, the
  expression for the saturation humidity over ice (``q_sat_ice``) and the roughness-length regression are marked ``VERIFY``
  in ``SEB_turbulent.py``: the coefficients have not all been checked against the literature. Check them before
  relying on the calculated fluxes.
* The surface is a single lumped layer of at least **SEB_TL_thick**; the results depend on this thickness (a thinner layer
  warms and melts more readily). The temperature of the whole layer is set to the surface temperature at the end of the step.
* The sub-stepped solver for ``prescribed`` fluxes does not thicken the top layer when the melt is more than its mass; the
  others do. The melt it finds is passed to the meltwater scheme, which removes that mass from the column, so on a
  heavy-melt day it can remove more than a thin top layer holds. For the solvers that do thicken, the thickening changes the melt
  by about 1% on such days (27 of 600 random heavy-melt test days, and no change in the 20-year DYE-2 runs; see the changelog
  for version 4.0.0).
* With ``albedo_method`` ``prescribed`` the albedo is an input (``ALBEDO``); the model does not calculate it, so the effect of melt or of the accumulation of
  new snow on the albedo is only in the forcing. For this reason **IceMeltFactor** is not used with the SEB.
* **Sublimation and the latent heat flux.** With prescribed fluxes, sublimation and deposition are taken from the forcing
  (``SUBLIM``), and the forcing ``QL`` is used in the energy balance; with one product (for example ``QL`` and the evaporation
  of MERRA-2) the two are consistent, which the model checks (see :ref:`sign-conventions`). With ``bulk_richardson`` or
  ``monin_obukhov``, :math:`Q_L` is calculated from the surface temperature of the CFM, and so is the mass that goes with
  it: at every time step the SEB returns :math:`Q_L/L\,\Delta t` (kg m\ :sup:`-2`), with :math:`L` the same latent heat as in
  the flux (sublimation below the melting point, vaporization at it; sub-daily forcing is summed over the sub-steps, each with
  its own surface temperature). A negative value is sublimation and is removed from the surface like any sublimation; a
  positive value is deposition and is added to the accumulation. The forcing ``SUBLIM`` is then only the first estimate used
  for ``iceout`` and the grid set-up (use zero or leave it out if you have no estimate, but note that ``iceout`` is then too
  high by the mean sublimation, which the end-of-run mass-consistent ``iceout`` diagnostic reports). The ``climate`` output
  shows the sublimation that was used. **SUBLIM** has to be true for this; the configuration checker warns otherwise.
* Rain heat assumes that the surface is at the melting point (see above), and the rain has no effect on the albedo.
* ``SEB_fqs`` and ``SEB_fqs_subdt`` were renamed ``solve_Ts`` and ``solve_Ts_subdt`` in version 4.0.0.

References
==========

* Andreas, E. L. (1987), A theory for the scalar roughness and the scalar transfer coefficients over snow and sea ice, *Boundary-Layer Meteorol.*, 38(1-2), 159-184.
* Beljaars, A. C. M., and Holtslag, A. A. M. (1991), Flux parameterization over land surfaces for atmospheric models, *J. Appl. Meteorol.*, 30(3), 327-341.
* Braithwaite, R. J. (1995), Aerodynamic stability and turbulent sensible-heat flux over a melting ice surface, the Greenland ice sheet, *J. Glaciol.*, 41(139), 562-571.
* Cuffey, K. M., and Paterson, W. S. B. (2010), *The Physics of Glaciers*, 4th ed., Butterworth-Heinemann, Oxford.
* Hock, R., and Holmgren, B. (2005), A distributed surface energy-balance model for complex topography and its application to Storglaciären, Sweden, *J. Glaciol.*, 51(172), 25-36.
* Louis, J.-F. (1979), A parametric model of vertical eddy fluxes in the atmosphere, *Boundary-Layer Meteorol.*, 17, 187-202.
* Paulson, C. A. (1970), The mathematical representation of wind speed and temperature profiles in the unstable atmospheric surface layer, *J. Appl. Meteorol.*, 9, 857-861.
* Yen, Y.-C. (1981), Review of thermal properties of snow, ice and sea ice, *CRREL Report 81-10*, U.S. Army Cold Regions Research and Engineering Laboratory, Hanover, NH.

.. VERIFY (references): Andreas (1987), Beljaars and Holtslag (1991), Braithwaite (1995) and Hock and Holmgren (2005) were
   checked against the reference list / text of Hock_JGlac_2005.pdf (the Hock and Holmgren paper, Eq. 19 is the rain heat).
   Louis (1979) and Paulson (1970) are NOT in ~/papers and are from memory (the code names them but gives no citation):
   Max will check them. Cuffey and Paterson (2010) is cited from Max's pointer (Eq. 9.2) and the diffusion.py docstring; I did not
   open the book. Yen (1981): title/report number from the PDF header. The CFM paper (Stevens et al., 2020) does not describe the SEB module.
