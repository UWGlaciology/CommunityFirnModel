.. _vapor-page:

Vapor transport in the pore space (experimental)
------------------------------------------------

.. note::
   **Experimental (beta).** The vapor transport is off by default (``vaporTransport`` false) and, when it is off, the
   model results are identical to those of a model without the vapor code. It conserves mass and its energy
   consistency has been checked, but it has **not been validated against measured firn density, temperature or
   isotope profiles**. Do not use its results for science without checking them. The status is listed at the end of
   this page.

Overview
========

Water vapor in the pore space of firn moves from the warmer to the colder layers: the saturation vapor density rises
steeply with temperature, so a temperature gradient is a vapor density gradient, and the vapor diffuses along it. In the
firn that it enters it condenses (the pore space is saturated), and where it leaves, ice sublimates. The result is a
redistribution of mass between layers, without compaction. It is the process that moves water isotopes in firn
(see :doc:`../files/isotopeDiffusion`), and on Mars, where the temperature gradients are large and the pressure is low, it
dominates the densification of the near-surface firn (Arthern et al., 2000).

The module is ``vapor.py``. It is called from the time loop of ``firn_density_nospin.py`` after the densification step.

What it does
============

**A mass change at fixed volume.** The layers of the model have a mass and a thickness; the density is their ratio.
Densification (the ``drho_dt`` of the :doc:`densification models <densification>`) changes the volume at fixed mass:
the layer gets thinner. Vapor transport does the opposite: it changes the mass at fixed volume, because the vapor fills
or empties the existing pores. The thickness of a layer does not change, so the vapor is not a densification term and
is never added to ``drho_dt``.

**The flux.** Between two layers the downward mass flux of vapor is

.. math::

   q = -D\,\phi\,\frac{1}{\tau}\,\frac{\partial \rho_{vs}}{\partial z}

with :math:`\rho_{vs}(T)` the saturation vapor density over ice, :math:`D` the diffusivity of water vapor in the
gas of the pores, :math:`\phi = 1 - \rho/\rho_i` the porosity, and :math:`1/\tau` the inverse tortuosity of the firn. The
form is that of Whillans and Grootes (1985) and Arthern et al. (2000, Eq. 5). The pieces:

* The saturation vapor pressure over ice is :math:`p_{sat} = 3.454\times10^{12}\exp(-6133/T)` Pa in air (the one of the isotope
  diffusion, Johnsen et al., 2000), and Eq. 6 of Arthern et al. (2000) for the CO\ :sub:`2` atmosphere of Mars (611 Pa at the
  triple point, latent heat of sublimation 2.83 :math:`\times` 10\ :sup:`6` J kg\ :sup:`-1`); :math:`\rho_{vs} = p_{sat} M / (R T)`.
* The diffusivity in air is :math:`D = 2.11\times10^{-5}\,(T/273.15)^{1.94}\,(101325/P)` m\ :sup:`2` s\ :sup:`-1` times the
  factor **vapor_yosida_factor** (default 1 for air), with :math:`P` the pressure **site_pressure**. For CO\ :sub:`2`
  it is the expression of Clow (1987) with a temperature exponent of 2 and the Yosida factor of 4 (the vapor diffuses
  faster in snow than in free gas, Yosida, 1955).
* The inverse tortuosity is :math:`1/\tau = 1 - b(\rho/\rho_i)^2` (Johnsen et al., 2000), which is zero from
  :math:`\rho_i/\sqrt{b}` = 804 kg m\ :sup:`-3` (**vapor_tortuosity_b**, :math:`b` = 1.3): no vapor moves through firn that is
  denser than that.

The ingredients are shared with the isotope diffusion: ``isotopeDiffusion.py`` takes the saturation pressure, the
diffusivity and the tortuosity from ``vapor.py`` (its results are unchanged, bit for bit).

**Conservation.** The flux through each face between two layers is calculated once and used by both neighbours, and the
flux through the top and the bottom of the column is zero (unless **vapor_surface_exchange**, below), so the mass of the
column is conserved to rounding. A layer cannot gain more than its free pore space or lose more than it has.

**Time stepping.** The vapor step is explicit and uses the temperature of the last step. A vapor molecule moves
about :math:`\sqrt{D\,\Delta t}` = 4 cm in a day, which is much more than the thickness of the thin near-surface layers
that daily forcing produces (a fraction of a millimeter on days without snowfall), and the flux does not depend on the mass of
a layer. Two measures keep this stable:

* *Blocks.* Adjacent layers are grouped into blocks that are at least **vapor_min_thickness** thick (default 1 cm; 0 is one
  layer per block). The flux is calculated between the blocks (temperature weighted by mass, pore space by thickness), and the mass
  change of a block is shared by its layers: in proportion to their mass for a loss, to their free pore space for a gain, so
  all layers of a block lose the same fraction. The layers of the model are not changed (unlike with **merging**).
* *A limiter.* The flux of each face is scaled so that no block loses more than half of its mass or fills more than half of
  its free pore space in a time step. A note is printed once if the limiter acts. At DYE-2 and Summit it did not act in the
  daily runs of 1980 to 2025 with the default block thickness.

Latent heat
===========

Sublimation absorbs heat and condensation releases it, so the vapor flux also carries heat, :math:`L q`. Because
:math:`q` follows the temperature gradient, this is a thermal conductivity,

.. math::

   K_v = L\,D\,\phi\,\frac{1}{\tau}\,\frac{d\rho_{vs}}{dT},

which is added to the conductivity of the firn in the heat diffusion (``heatDiff`` and ``refreezeDiff`` in ``diffusion.py``)
when **vapor_latent_heat** is true. It is solved implicitly with the rest of the heat equation, so it is stable even in thin
layers, and it is the same heat that goes with the mass that the vapor module moves. Arthern et al. (2000, Eq. B3) do the same
(with the vapor term of Clow, 1987).

:math:`K_v` is 5 to 30% of the conductivity of firn at 150 kg m\ :sup:`-3` and 253 to 273 K, 1 to 6% at 300 kg m\ :sup:`-3`,
0.3 to 1.8% at 450, 0.1 to 0.4% at 600, and zero from 804 kg m\ :sup:`-3`: it matters in the surface snow near the melting point.

**Do not count the vapor twice.** The default conductivity (**conductivity** ``Calonne2019``, Calonne et al., 2019) comes from
3-D image computations with air and ice only (conduction, no vapor), so adding :math:`K_v` is correct. Conductivities
that are fits to measurements (``Yen_*``, ``Sturm``, ``VanDusen``, ``Anderson``, ``Jiawen``, ...) are apparent conductivities
that may already include the vapor; do not use **vapor_latent_heat** with them.

Checked: summing the heat that :math:`K_v` carries and :math:`L` times the vapor mass flux over every face and time step of a
daily run, the two agree to 0.5 to 2.8% (net) at DYE-2 and 0 to 1.2% at Summit, with the largest difference in the top meter
where the layers are thin (the block average).

How to use it
=============

.. code-block:: json

   "vaporTransport": true,
   "vapor_latent_heat": true

It is used with any ``physRho`` model; with **no_densification** the mass is frozen too. The keys are described in :doc:`../running/json` (**vaporTransport**, **vapor_atmosphere**,
**vapor_yosida_factor**, **vapor_tortuosity_b**, **vapor_min_thickness**, **vapor_latent_heat**, **vapor_surface_exchange**).
Use **vapor_latent_heat** only with the ``Calonne2019`` conductivity (see above). If **isoDiff** is also true the model
prints a warning: the isotope diffusion does not include the vapor transport (the mass that the vapor moves carries no
isotopic composition).

How big is it?
==============

Daily steps, MERRA-2 forcing 1980 to 2025, example setups (reference climate interval spin-up), vapor on minus off, after 1981:

.. list-table::
   :header-rows: 1

   * - Site
     - Vapor mass only
     - Mass and latent heat (compared with the mass only)
   * - DYE-2 (Greenland, melt)
     - firn air content +0.14 m (+1.8%); density -1 to -2 kg m\ :sup:`-3` at 1 to 100 m; runoff -0.7%
     - firn air content +0.39 m more (+3.6%); density -3 to -5.5 kg m\ :sup:`-3` more at 5 to 100 m; temperature +0.6 K at 1 m, -0.3 K at 50 m
   * - Summit (Greenland, no melt)
     - firn air content -0.003 m; density changes under 0.1 kg m\ :sup:`-3` below 1 m
     - firn air content -0.50 m (-1.35%); density +4 to +5 kg m\ :sup:`-3` at 10 to 100 m; temperature +0.15 K

At both sites the latent heat matters more than the mass transport itself, through the temperature and the temperature
dependence of the densification. Monthly steps give about half the DYE-2 mass effect, because the monthly temperature
gradients are smoother. These are model sensitivities, not validated results.

Mars
====

With ``"vapor_atmosphere": "CO2"`` and **site_pressure** set to the surface pressure of Mars (hPa) the module uses the saturation
vapor pressure of Arthern et al. (2000) and the diffusivity of Clow (1987). The check ``mars_aw00_check.py`` (in the root
of the repository; ``tests/test_vapor.py`` runs it) evaluates their critical accumulation rate,
:math:`b_s^{crit} = \bar q(0)/(1-\rho_s/\rho_i)` (Eq. 22; the accumulation that is needed to keep a porous snowpack), with
the functions of ``vapor.py``, for a mean temperature of 175 K, a seasonal amplitude of 30 K and a surface density of
360 kg m\ :sup:`-3`:

.. list-table::
   :header-rows: 1

   * - Diurnal amplitude
     - Arthern et al. (2000), read from Figs. 9 and 10
     - ``vapor.py`` (their porosity / with the Johnsen pore space), 600 Pa
   * - 0 K
     - about 0.5 kg m\ :sup:`-2` Mars year\ :sup:`-1`
     - 0.75 / 0.60
   * - 10 K
     - about 5
     - 7.6 / 6.1

The pressure at the pole is not given in that paper (600 Pa is assumed; the flux is proportional to 1/P, and about 900 Pa
would match). The diurnal cycle increases the critical rate by a factor of 10 in both. A column with the CFM heat solver and
**vapor_surface_exchange** gives the flux of their Eq. 21 to 1% in a snow half space.

**vapor_surface_exchange** (default false) gives the top of the column the flux just beneath the surface, :math:`q(0,t)`, which in
Arthern et al. (2000) is part of the accumulation, :math:`b_s = \rho_s v(0,t) + q(0,t)`. The mass comes from the atmosphere,
so the mass of the column changes and the accumulation of the model must be the mass deposited on the surface,
:math:`b_s - q(0)`. Do not use it together with sublimation or deposition of the same vapor (**SUBLIM**, or the turbulent fluxes
of the **SEB**): on Earth that exchange is already in the model.

What is still missing for Mars (gravity, the densification of Arthern et al. (2000), forcing without an SEB, a steady-state
solver for the density profile) is listed in ``mars_todo.md`` in the root of the repository.

Status and limitations
======================

* **Not validated against data.** The checks are internal: mass conservation, the energy consistency of the latent heat, the
  Mars analytic result. A comparison with measured density, temperature or isotope profiles, or with another model, has not
  been done.
* **Earth constants.** The Yosida factor of 1 for air (no enhancement of the diffusivity in snow) is the free-air value of the
  isotope diffusion and not verified for bulk vapor transport; a factor of 4, as used for Mars, would multiply the flux by 4.
  The pore space uses the tortuosity of Johnsen et al. (2000), which was calibrated on isotope diffusion lengths; Arthern et al.
  (2000) use the porosity alone (20% more flux at 360 kg m\ :sup:`-3`, and the vapor goes on above 804 kg m\ :sup:`-3`).
* **Isotopes** are not carried by the vapor (a warning is printed).
* **The surface is closed on Earth.** The exchange with the atmosphere is the sublimation and deposition of the model (``SUBLIM``
  or the latent heat flux of the SEB), so a second flux through the surface would count it twice. The vapor that the pore
  space exchanges with the air beyond that (a rough bound is 1 kg m\ :sup:`-2` a\ :sup:`-1` at 250 K, a few percent of the sublimation at Summit
  or DYE-2) is not included.
* **Resolution.** The flux is calculated on blocks of at least **vapor_min_thickness**; the sensitivity of the results to
  that value was not mapped. A thin layer of snow over ice gives a top flux (**vapor_surface_exchange**) that is too small by
  about a factor of two, because the gradient between layer centers includes half of the ice cell.
* **Two sites.** It has only been run at DYE-2 and Summit with MERRA-2 forcing, and with the bucket scheme.
* **Knudsen diffusion** (pores smaller than about 10 micrometer) is not included (Arthern et al., 2000, Sec. 5).

References
==========

.. VERIFY: the page ranges of Johnsen et al. (2000), Clow (1987) and Yosida (1955) are from memory; check them against the papers.

* Arthern, R. J., Winebrenner, D. P., and Waddington, E. D. (2000). Densification of water vapor-saturated firn under Martian conditions. *Icarus*, 144, 367-386.
* Calonne, N., Milliancourt, L., Burr, A., Philip, A., Martin, C. L., Flin, F., and Geindreau, C. (2019). Thermal conductivity of snow, firn, and porous ice from 3-D image-based computations. *Geophys. Res. Lett.*, 46, 13079-13089.
* Clow, G. D. (1987). Generation of liquid water on Mars through the melting of a dusty snowpack. *Icarus*, 72, 95-127.
* Johnsen, S. J., Clausen, H. B., Cuffey, K. M., Hoffmann, G., Schwander, J., and Creyts, T. (2000). Diffusion of stable isotopes in polar firn and ice: the isotope effect in firn diffusion. In *Physics of Ice Core Records*, 121-140.
* Whillans, I. M., and Grootes, P. M. (1985). Isotopic diffusion in cold snow and firn. *J. Geophys. Res.*, 90, 3910-3918.
* Yosida, Z. (1955). Physical studies on deposited snow. I. Thermal properties. *Contributions from the Institute of Low Temperature Science*, 7, 19-74.
