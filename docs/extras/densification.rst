.. _densification-page:

Firn densification models (physRho)
-----------------------------------

.. DRAFT NOTES FOR MAX (RST comments are not rendered; search this file for "VERIFY"):
   Written from physics.py, the CFM v1.0 paper (Stevens et al., 2020, Sect. 2.2) and the original papers
   in ~/papers. Every coefficient listed below was compared with the code; where a paper could be
   searched, the coefficient was compared with the paper too (HL, Barnola, LZ11/15, Helsen, LIG, KM,
   Brils, GSFC2020/Medley 2022, Breant, Crocus/van Kampenhout). Items I could not confirm are marked VERIFY.

Overview
========

The key **physRho** selects the firn densification model, i.e. the equation that gives the rate of density
change, :math:`d\rho/dt`, of each model layer (node). The CFM contains the densification models listed
below, coded as described in their original publications, with the corrections and changes noted in
each entry (see Stevens et al., 2020, Section 2.2). Which one to use depends on the site, the climate
and the question; the models can differ substantially from each other, in particular outside the climate
conditions they were calibrated for (see Lundin et al., 2017, for an intercomparison of several of them).

All models are used in the same way. At every time step, the CFM calls the method of the class
``FirnPhysics`` (``physics.py``) that has the name given by **physRho**. The method receives the state
of the firn column (density, temperature, overburden stress, mean accumulation rate, and so on) and returns
the densification rate, in kg m\ :sup:`-3` s\ :sup:`-1`, and a viscosity. The viscosity (the output ``viscosity``) is that of a
Newtonian viscous fluid with the same densification rate, :math:`(1/\rho)\,d\rho/dt = \sigma / (2\eta)`, i.e.
:math:`\eta = \rho\,\sigma / (2\,d\rho/dt)`, with :math:`\sigma` the overburden stress. It is a diagnostic: it is not used
to calculate the densification. (Before version 4.1.0 most models returned this value divided by :math:`\sigma^2`.) The density is then advanced with an explicit scheme,

.. math::

    \rho_{new} = \rho_{old} + \frac{d\rho}{dt}\,dt .

The effects of horizontal strain (**strain_softening**, **horizontal_divergence**) are applied after this
(see :doc:`../running/json`), and **no_densification** switches the densification off.

The initial density profile of a run is generated with the analytic Herron and Langway (1980) solution
whatever **physRho** is (see :doc:`../running/spinup`).

Common features of the models
=============================

**Accumulation rate or stress.** Most models use the accumulation rate as a proxy for the overburden
stress: for constant accumulation, the stress at depth :math:`z` is :math:`\sigma = \dot b\,g\,\tau(z)`, with
:math:`\tau` the age of the firn. A few models are formulated with the stress itself (marked "stress" in the
table). For the models that use the accumulation rate, **bdot_type** selects which one:

* ``mean`` (the default): the mean accumulation rate over the lifetime of each parcel of firn,

  .. math::

      \bar{\dot b}(z,t) = \frac{1}{\tau_{z,t}} \int_{t-\tau_{z,t}}^{t} \dot b(t')\,dt' ,

  after Li and Zwally (2011, 2015). It is used because a model that depends on the instantaneous
  accumulation rate predicts no densification when the accumulation rate is zero. In steady state the mean
  and the instantaneous rate are the same.
* ``instant``: the accumulation rate of the current time step. Several models print a warning that this is
  not recommended, and two do not work with it (see "Limitations" below).
* ``stress``: listed as an option of **bdot_type**, but still in progress (see :doc:`../running/json`).

The accumulation rates in the input and in the configuration are always in m ice equivalent per year; each
model converts them internally to the units of its original publication (m w.e. per year for the Herron and
Langway and Li and Zwally families, kg m\ :sup:`-2` per year for the Arthern family).

**Temperature.** Models differ in whether they use the temperature :math:`T` of each parcel of firn, the mean
annual temperature :math:`T_m`, or both. In the CFM, :math:`T_m` is a 10-year running mean of the surface temperature, updated
every time step (it is the mean of the whole record if the running mean cannot be calculated). The densification models of the Arthern family use :math:`T_m` in the grain-growth term of the
Arrhenius factor, as in the latest form of Kuipers Munneke et al. (2015) (Steger et al., 2017).

**Stages.** Densification is usually divided into stages by density. Unless noted, stage 1 is
:math:`\rho \le 550` kg m\ :sup:`-3` (``RHO_1``) and stage 2 is above it. The models differ in how they treat
the transition and in what they do above about 800 kg m\ :sup:`-3`.

**Wet firn.** When **MELT** is true, the densification models still describe the dry compaction. The
density change by refreezing is handled by the meltwater scheme (see :doc:`bucket_scheme`). Of the models
here only Crocus has a dependence on the liquid water content.

Summary
=======

.. list-table::
   :header-rows: 1
   :widths: 17 24 9 50

   * - **physRho**
     - Reference
     - Driven by
     - Developed for, and remarks
   * - ``HLdynamic``
     - Herron and Langway (1980)
     - :math:`\bar{\dot b}`
     - Ice-core sites (17 cores). Dynamic form, in which the temperature is that of each parcel. Original model of the family.
   * - ``HLSigfus``
     - Herron and Langway (1980), Johnsen's stress form
     - :math:`\bar{\dot b}`, stress
     - Stage 1 as ``HLdynamic``; stage 2 written in terms of the overburden stress.
   * - ``Barnola1991``
     - Barnola et al. (1991)
     - :math:`\bar{\dot b}`, stress
     - Ice-core gas-age calculations. Stage 1 as Herron and Langway; creep law above.
   * - ``Goujon2003``
     - Arnaud et al. (2000); Goujon et al. (2003)
     - stress
     - Ice-core gas-age reconstructions (Antarctica, Greenland). Grain-scale physics; the stage 1 to 2 transition density depends on :math:`T_m`.
   * - ``Breant2017``
     - Bréant et al. (2017)
     - stress
     - Goujon model with three activation energies and a fixed transition density; deglaciation. Results made with versions before 4.1.0 are not valid.
   * - ``Li2004``
     - Li and Zwally (2004)
     - :math:`\bar{\dot b}`
     - Shallow firn. The vapor-transfer term is not coded.
   * - ``Li2011``
     - Li and Zwally (2011)
     - :math:`\bar{\dot b}`
     - Tuned for Greenland. Not suitable for wet firn.
   * - ``Li2015``
     - Li and Zwally (2015)
     - :math:`\bar{\dot b}`
     - Tuned for Antarctica. Not suitable for wet firn. Fails with ``instant``.
   * - ``Helsen2008``
     - Helsen et al. (2008)
     - :math:`\bar{\dot b}`
     - Antarctica. Same form as Li and Zwally, with one coefficient for both stages.
   * - ``Arthern2010S``
     - Arthern et al. (2010), steady-state form
     - :math:`\bar{\dot b}`
     - Based on compaction measurements in Antarctica. Basis of the models below it.
   * - ``Arthern2010T``
     - Arthern et al. (2010), transient form
     - stress
     - Needs **physGrain** true. Sensitive to the surface grain size; not in common use.
   * - ``Ligtenberg2011``
     - Ligtenberg et al. (2011)
     - :math:`\bar{\dot b}`
     - Antarctica (IMAU-FDM, RACMO2 subsurface scheme).
   * - ``Simonsen2013``
     - Simonsen et al. (2013)
     - :math:`\bar{\dot b}`
     - Greenland, tuned with airborne radar. The factors :math:`f_0`, :math:`f_1` are site-dependent.
   * - ``KuipersMunneke2015``
     - Kuipers Munneke et al. (2015)
     - :math:`\bar{\dot b}`
     - Greenland (IMAU-FDM).
   * - ``Brils2022``
     - Brils et al. (2022)
     - :math:`\bar{\dot b}`
     - Greenland (IMAU-FDM v1.2G).
   * - ``Veldhuijsen2023``
     - Veldhuijsen et al. (2023)
     - :math:`\bar{\dot b}`
     - Antarctica (IMAU-FDM v1.2A). In development; see "Limitations".
   * - ``GSFC2020``
     - Medley et al. (2022)
     - :math:`\bar{\dot b}`
     - Greenland and Antarctica, one calibration against 256 density profiles with MERRA-2 forcing. Only ``bdot_type`` ``mean``.
   * - ``Crocus``
     - Vionnet et al. (2012); van Kampenhout et al. (2017)
     - stress
     - Viscous snow model for seasonal and polar snow, with a dependence on liquid water. Used here as a hybrid above 450 kg m\ :sup:`-3`.
   * - ``Morris2014``
     - Morris and Wingham (2014)
     - stress (stage 1)
     - Greenland, stage 1 only; Herron and Langway for stage 2. In development.

The models
==========

Herron and Langway (1980) and relatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many later models share two assumptions of Herron and Langway (1980): the densification rate is linear in
the porosity (:math:`d\rho/dt = c\,(\rho_i - \rho)`), and it has an Arrhenius dependence on temperature,
with :math:`c = k \exp(-Q/RT)\,\dot b^{a}`. For HL,

.. math::

    c = \begin{cases}
        11\,\exp(-10160/RT)\,\dot b & \rho \le 550 \\
        575\,\exp(-21400/RT)\,\dot b^{0.5} & \rho > 550
    \end{cases}

with :math:`\dot b` in m w.e. per year and :math:`Q` in J mol\ :sup:`-1`. The units of :math:`k` differ between the two
stages because the exponent of :math:`\dot b` differs. The activation energies are lower than those of most
later models, so HL is less sensitive to temperature variations within a year.

``HLdynamic``
    The "dynamic" form: the equations above, in which :math:`T` is the temperature of each parcel (the original
    model used the mean annual temperature) and :math:`\dot b` is the mean or the instantaneous accumulation rate.
    With a steady firn-air run (**FirnAir**, ``runtype`` ``steady``, ``instant`` accumulation) the steady
    temperature is used. The CFM contains three formulations of the Herron and Langway model (Lundin et al.,
    2017), which agree in steady state but differ in transient runs: this dynamic form, the stress-based form
    (``HLSigfus``), and the analytic solution, which is not a choice of **physRho** but generates the initial
    profile of every run.

``HLSigfus``
    The Herron and Langway model with a modification suggested by reviewer Sigfús J. Johnsen, in which the
    overburden stress replaces the accumulation rate in stage 2 (Herron and Langway, 1980; Lundin et al., 2017).
    A thin transition zone between 550 and 551 kg m\ :sup:`-3` uses the dynamic stage-2 form, because the
    stress-based rate would otherwise become unbounded at the start of stage 2.

``Barnola1991``
    Barnola et al. (1991), developed for ice-core gas-age calculations. Stage 1 is Herron and Langway
    (with :math:`a=1`). Above 550 kg m\ :sup:`-3` the densification rate is a power-law creep,

    .. math::

        \frac{d\rho}{dt} = \rho\,A_0 \exp\!\left(-\frac{Q}{RT}\right) f\,\sigma^{n} ,

    with :math:`A_0 = 2.54\times10^{4}` MPa\ :sup:`-3` s\ :sup:`-1`, :math:`Q = 60` kJ mol\ :sup:`-1`, :math:`n = 3`,
    and the effective stress :math:`\sigma` taken to be the overburden stress. For
    :math:`550 \le \rho \le 800` kg m\ :sup:`-3`, :math:`f` is an empirical polynomial of the density that
    matches the rate and its derivative at the stage 1 to 2 transition and at close-off
    (:math:`\alpha=-37.455`, :math:`\beta=99.743`, :math:`\delta=-95.027`, :math:`\gamma=30.673`, exponent of 10).
    Above 800 kg m\ :sup:`-3`, :math:`f` is that of Pimienta and Duval (1987), as used by Barnola et al. (1991).

``Goujon2003``
    The Goujon et al. (2003) model, which adds heat diffusion to the grain-scale model of Arnaud et al. (2000),
    developed for ice-core gas-age reconstructions in Antarctica and Greenland. It is stress-based. Stage 1
    follows grain-boundary sliding (Alley, 1987) and does not depend on temperature; stage 2 follows sintering
    theory (Arzt, 1982) with activation energy 60 kJ mol\ :sup:`-1`; above a relative density of 0.9 the effective stress
    includes the pressure of the gas in closed bubbles, and above 0.98 a further term applies. The transition
    relative density is :math:`D_0 = 0.00226\,T_m + 0.03` (:math:`T_m` in K), limited to 0.59, and the
    ice density depends on temperature. A scaling factor :math:`\gamma` in stage 1 is iterated at every time step so
    that the densification rate is continuous at the transition (the implementation follows the
    suggestion of Buizert et al., 2015, to place the transition at :math:`D_0 + \epsilon`). Stage 2 uses a
    grain contact area and coordination number from the original paper. The CFM keeps the state of this
    iteration between time steps.

``Breant2017``
    The Goujon model with a different creep parameter, introduced to improve the temperature sensitivity at low
    temperatures (Bréant et al., 2017). The creep parameter of Goujon is replaced by a sum of three
    thermally activated mechanisms,

    .. math::

        A = A_0 \left( a_1 e^{-Q_1/RT} + a_2 e^{-Q_2/RT} + a_3 e^{-Q_3/RT} \right) ,

    with :math:`Q_1 = 110`, :math:`Q_2 = 75`, :math:`Q_3 = 1.5` kJ mol\ :sup:`-1` and
    :math:`a_1 = 1.05\times10^{9}`, :math:`a_2 = 1400`, :math:`a_3 = 6\times10^{-15}` (their Table 1), and a fixed relative transition density
    :math:`D_0 = 0.56`. According to Bréant et al. (2017) the equivalent activation energy stays between 54 and
    61 kJ mol\ :sup:`-1`, close to the 60 of Goujon et al. (2003), so only moderate changes in the behavior are expected. The
    **Before version 4.1.0 the third term of Eq. 7 was coded with** :math:`a_1`
    **instead of** :math:`a_3`, **which made the densification rate about 10**\ :sup:`10` **times too large; results made with this option
    in earlier versions are not valid.**

Li and Zwally family
~~~~~~~~~~~~~~~~~~~~

These models (Li and Zwally, 2002, 2004, 2011, 2015; Helsen et al., 2008) were developed to predict surface
elevation changes caused by accumulation and temperature variability and by firn compaction. They share
the form

.. math::

    \frac{d\rho}{dt} = \beta\; 8.36\,(273.15 - T)^{-2.061}\;\bar{\dot b}\;(\rho_i - \rho) ,

with :math:`T` the temperature of the parcel in K and :math:`\bar{\dot b}` in m w.e. per year. They differ in
the coefficient :math:`\beta`. The rate becomes very large near the melting point (it is infinite at
273.15 K), so these models are not suited to simulations of wet firn (Stevens et al., 2020).

``Li2004``
    :math:`\beta = 139.21 - 0.542\,T_m` (:math:`T_m` in K), one value for both stages. The vapor-transfer term
    of Li and Zwally (2004) is not coded.

.. VERIFY (Li2004): the coefficients 139.21 and 0.542 are not in the text I could extract from Li_Zwally_2004.pdf;
   a comment in physics.py says they come from Arthern et al. (2010), Eq. 2, "not sure where Rob got that".
   Please confirm the source of beta for Li2004 (or say how it should be described).

``Li2011``
    Tuned for Greenland. For stage 1, :math:`\beta_1 = -9.788 + 8.996\,\dot b_m - 0.6165\,T_{m,C}`, and for stage 2,
    :math:`\beta_2 = \beta_1 / (-2.0178 + 8.4043\,\dot b_m - 0.0932\,T_{m,C})`, with :math:`T_{m,C}` the mean annual
    temperature in :math:`^\circ`\ C and :math:`\dot b_m` the long-term mean accumulation rate of the site (in the CFM,
    the mean over the column). With ``instant``, the coefficient is calculated from the accumulation rate of the
    time step.

``Li2015``
    Tuned for Antarctica: :math:`\beta_1 = -1.218 - 0.403\,T_{m,C}` and
    :math:`\beta_2 = \beta_1\,(0.792 - 1.080\,\dot b_m + 0.00465\,T_{m,C})`. Does not run with ``instant``
    (see "Limitations").

``Helsen2008``
    Antarctica. The form of Li and Zwally with a single coefficient for both stages,
    :math:`\beta = 76.138 - 0.28965\,T_m` (:math:`T_m` in K).

Arthern family
~~~~~~~~~~~~~~

Arthern et al. (2010) derived their model from compaction-rate measurements at sites in Antarctica; it was the
first model based on measured compaction rates and not on rates derived with Sorge's law. They identified two
processes, diffusion of water molecules through the lattice (activation energy :math:`E_c = 60` kJ mol\ :sup:`-1`) and
grain growth (:math:`E_g = 42.4` kJ mol\ :sup:`-1`), whose combination gives a lower effective activation energy.

``Arthern2010T``
    The transient model of the appendix of the paper,
    :math:`d\rho/dt = k_c\,(\rho_i - \rho)\exp(-E_c/RT)\,\sigma / r^2`, with :math:`r^2` the grain radius squared
    (:math:`k_c = 9.2\times10^{-9}` for stage 1 and :math:`3.7\times10^{-9}` for stage 2), from the evolving stress and grain
    size. It needs **physGrain** to be true (it prints a message and does not run otherwise), is sensitive to
    the surface grain size (**r2s0**), which is poorly constrained over the range of climates of the ice
    sheets, and is not in common use (Stevens et al., 2020).

``Arthern2010S``
    The steady-state model of the main text of the paper, in the form of Herron and Langway with

    .. math::

        c_0 = 0.07\,\dot b\,g\,\exp\!\left(-\frac{E_c}{RT} + \frac{E_g}{R T_m}\right) \;(\rho \le 550), \qquad
        c_1 = 0.03\,\dot b\,g\,\exp\!\left(-\frac{E_c}{RT} + \frac{E_g}{R T_m}\right) \;(\rho > 550) ,

    with :math:`\dot b` in kg m\ :sup:`-2` per year. This is the basis of the next models. They differ from it in how
    :math:`c_0` and :math:`c_1` are scaled or calibrated:

``Ligtenberg2011``
    The model of IMAU-FDM for Antarctica (the subsurface scheme of RACMO2). :math:`c_0` and :math:`c_1` are multiplied by
    :math:`M_0 = 1.435 - 0.151\ln\dot b` and :math:`M_1 = 2.366 - 0.293\ln\dot b` (:math:`\dot b` in kg m\ :sup:`-2` per year),
    with a lower limit of 0.25.

``KuipersMunneke2015``
    As Ligtenberg et al. (2011), refitted for Greenland: :math:`M_0 = 1.042 - 0.0916\ln\dot b`,
    :math:`M_1 = 1.734 - 0.2039\ln\dot b` (lower limit 0.25). The grain-growth term uses :math:`T_m`, as
    modified by Steger et al. (2017).

``Brils2022``
    The model of IMAU-FDM v1.2G for Greenland, with refitted coefficients:
    :math:`M_0 = 0.6688 + 0.0048\ln\dot b` and :math:`M_1 = 1.7465 - 0.2045\ln\dot b` (lower limit 0.25).
    Brils et al. (2022) also changed the fresh-snow density and the thermal conductivity; in the CFM these are
    the options ``Brils22`` of **srho_type** and ``Calonne2019`` of **conductivity**.

``Veldhuijsen2023``
    The model of IMAU-FDM v1.2A for Antarctica (Veldhuijsen et al., 2023). It is in development, and the code
    prints a message to that effect. **The coefficients in the code are those of Kuipers Munneke et al. (2015); the
    coefficients published by Veldhuijsen et al. (2023) are not implemented yet.** See "Limitations".

``Simonsen2013``
    Greenland. :math:`c_0` is multiplied by :math:`f_0` and :math:`c_1` by :math:`f_1\,\gamma`, with
    :math:`\gamma = 61.7\,\dot b^{-0.5} \exp(-3800/(R T_m))`. Simonsen et al. (2013) tuned the factors to each site
    with a Monte Carlo inversion of airborne radar data; the CFM uses the defaults :math:`f_0 = 0.8` and
    :math:`f_1 = 1.25` (S. Simonsen, personal communication, 2015), which should be tuned for your site.

``GSFC2020``
    The dry-snow densification scheme calibrated by Medley et al. (2022) for the GSFC firn densification model
    (GSFC-FDMv1.2.1), with the Arthern form and calibrated exponents of the accumulation rate and activation
    energies,

    .. math::

        c_0 = 0.07\,\dot b^{\alpha_0}\,g\,\exp\!\left(-\frac{E_{c0}}{RT} + \frac{E_g}{R T_m}\right), \qquad
        c_1 = 0.03\,\dot b^{\alpha_1}\,g\,\exp\!\left(-\frac{E_{c1}}{RT} + \frac{E_g}{R T_m}\right) ,

    with :math:`\alpha_0 = 0.91`, :math:`\alpha_1 = 0.644`, :math:`E_{c0} = 59\,500` J mol\ :sup:`-1`,
    :math:`E_{c1} = 56\,870` J mol\ :sup:`-1`. The parameters were calibrated against 256 firn density profiles from
    Greenland and Antarctica with a single calibration and with MERRA-2 forcing and a reference climate interval
    (see :doc:`../running/spinup`); Medley et al. (2022) stress that the forcing of a run should be consistent with
    the forcing used in the calibration, because biases in the forcing bias the coefficients. It works only with
    ``bdot_type`` ``mean`` (the model stops with ``instant``). The name refers to the year of the first version
    (Medley et al., 2020, the discussion paper); the coefficients in the code are those of the final paper.

.. VERIFY (GSFC2020): the docstring in physics.py says "Taken from Medley et al. (2020)"; the coefficients
   (alpha 0.91/0.644, Ec 59500/56870 J/mol) match Eq. 18 of Medley et al. (2022) exactly.

Other models
~~~~~~~~~~~~

``Crocus``
    The viscous constitutive law of the Crocus snow model (Vionnet et al., 2012; Brun et al., 1992), which was
    developed for seasonal snow:

    .. math::

        \frac{d\rho}{dt} = \frac{\rho\,\sigma}{\eta}, \qquad
        \eta = f_1 f_2\,\eta_0\,\frac{\rho}{c_\eta}\,\exp\!\left(a_\eta (T_{melt} - T) + b_\eta \rho\right) ,

    with :math:`\eta_0 = 7.62237\times10^{6}` kg s\ :sup:`-1` m\ :sup:`-1`, :math:`a_\eta = 0.1` K\ :sup:`-1`,
    :math:`b_\eta = 0.023` m\ :sup:`3` kg\ :sup:`-1`, :math:`f_2 = 4` (angular grains), and
    :math:`c_\eta = 358` kg m\ :sup:`-3` (van Kampenhout et al., 2017; Vionnet et al. give 250). The factor
    :math:`f_1 = (1 + 60\,\mathrm{LWC}/\Delta z)^{-1}` depends on the liquid water content of the layer
    (:math:`f_1 = 1` for dry firn). **The CFM uses the model as a hybrid: above 450 kg m**\ :sup:`-3` **the Crocus rate is replaced by the
    stage-2 form of Kuipers Munneke et al. (2015)**, because the Crocus law, which was developed for seasonal snow, does not densify
    firn fast enough at higher densities. The hybrid is switched on by a setting in the code (``CHybrid`` in ``physics.py``),
    which is true; it is not a configuration key.

.. VERIFY (Crocus): Max confirmed that the hybrid is intended. Only the origin of f1 = 1/(1 + 60 LWC/dz) remains
   undocumented ("From Vincent" in the code; the CFM v1.0 paper used f1 = 1).

``Morris2014``
    Morris and Wingham (2014) derived a compaction equation for dry firn below 550 kg m\ :sup:`-3` from
    neutron-probe density logs in boreholes at sites on the Greenland ice sheet:

    .. math::

        \dot\epsilon = \frac{k_0^*}{\rho_w g}\left(\frac{\rho_i - \rho}{\rho}\right)\left(1 - M_0 m\right)
        \frac{1}{H(\tau)}\exp\!\left(-\frac{E_\alpha}{RT}\right)\sigma ,

    where :math:`H(\tau)` is a temperature history function (the CFM keeps it as the variable ``Hx``, and switches
    the temperature-history output on for this model), :math:`m` is the normalized deviation of the stage-1
    density profile from a quadratic fit, and :math:`M_0 = 3.3`. The activation energy :math:`E_\alpha` is the key
    **QMorris** (preferred value of the paper, 110 kJ mol\ :sup:`-1`; 60 and 200 were also tested). The constant
    :math:`k_0^*` is not the value of the paper, which was found to be in error (E. Morris, personal
    communication, 2019; see :doc:`MWfix`). The CFM calculates the corrected value from the revised Table 2
    of Stevens et al. (2020), as a regression on :math:`Q` and :math:`T_m`. Morris and Wingham give no stage 2; the CFM uses
    the dynamic Herron and Langway equation above 550 kg m\ :sup:`-3`. The code prints a message that the model is still
    under development.

Options in the code that cannot be selected
============================================

``physics.py`` contains more methods than the validated options of **physRho**. The configuration check
(see :doc:`../running/json`) rejects them:

* ``ArthernWingham1998``: the mechanistic pressure-sintering model of Arthern and Wingham (1998), with
  grain-boundary sliding, dislocation creep, and boundary and lattice diffusion. It is experimental and has not
  been validated; the dislocation-creep term in particular should be checked against the paper.
* ``Max2018``, ``Max2018b``, ``MaxSP``: experimental, based on firn compaction data of the FirnCover campaign.
* ``Yamazaki1993``: the viscous model for wet snow of Yamazaki et al. (1993). It is not a choice of **physRho**: it
  is used for the top of the column (density below 300 kg m\ :sup:`-3`) when **stage_zero** is true. This is in
  development.
* ``Marshall1999``: a placeholder for the model of Marshall et al. (1999), for the densification of snow
  during rain; it does nothing yet.

Limitations
===========

Things to know about the current implementation:

* ``Veldhuijsen2023`` uses the coefficients of ``KuipersMunneke2015`` (see above).
* ``Li2015`` with ``bdot_type`` ``instant`` stops with an error (a variable is used before it is defined).
  ``GSFC2020`` stops with ``instant``, and ``Arthern2010T`` needs **physGrain**.

Adding a model
==============

A densification model is a method of the class ``FirnPhysics`` in ``physics.py``. The model is called with
the dictionary ``PhysParams``, whose entries are available as attributes of ``self``: among others ``iii`` (time
step), ``steps``, ``gridLen``, ``bdotSec`` and ``bdot_mean`` (accumulation rates), ``bdot_type``, ``Tz`` (temperature
of each node), ``T_mean``, ``rho``, ``mass``, ``sigma`` (overburden stress), ``age``, ``z``, ``dz``, ``dt``,
``r2`` (grain radius squared), ``LWC``, and the strain rates ``eps_xx`` and ``eps_yy``. If your model needs
other variables, add them to the dictionary in ``firn_density_nospin.py`` (and in ``firn_density_spin.py``). The method must return
the dictionary ``self.RD`` with the entry ``drho_dt`` (kg m\ :sup:`-3` s\ :sup:`-1`, positive when the density increases),
or, for a model whose constitutive law already depends on the horizontal strain rates, the entry ``eps_zz``
(see ``get_drho_dt`` in ``physics.py``), and the entry ``viscosity``. To make the name a valid choice of
**physRho**, add it to the allowed values in ``config_validation.py`` and to the key's description in
:doc:`../running/json`.

References
==========

.. VERIFY: all entries below were taken from the reference list of Stevens et al. (2020), from the papers in
   ~/papers, or from the physics.py docstrings, except where marked in the comments above.

* Alley, R. B. (1987), Firn densification by grain-boundary sliding: a first model, *Journal de Physique Colloques*, 48, C1-249.
* Arnaud, L., Barnola, J.-M., and Duval, P. (2000), Physical modeling of the densification of snow/firn and ice in the upper part of polar ice sheets, in *Physics of Ice Core Records*, Hokkaido University Press, 285-305.
* Arthern, R. J., and Wingham, D. J. (1998), The natural fluctuations of firn densification and their effect on the geodetic determination of ice sheet mass balance, *Climatic Change*, 40, 605-624.
* Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R. (2010), In situ measurements of Antarctic snow compaction compared with predictions of models, *J. Geophys. Res.*, 115, F03011, doi:10.1029/2009JF001306.
* Arzt, E. (1982), The influence of an increasing particle coordination on the densification of spherical powders, *Acta Metall. Mater.*, 30, 1883-1890.
* Barnola, J.-M., Pimienta, P., Raynaud, D., and Korotkevich, Y. S. (1991), CO2-climate relationship as deduced from the Vostok ice core: a re-examination based on new measurements and on a re-evaluation of the air dating, *Tellus B*, 43, 83-90.
* Bréant, C., Martinerie, P., Orsi, A., Arnaud, L., and Landais, A. (2017), Modelling firn thickness evolution during the last deglaciation: constraints on sensitivity to temperature and impurities, *Clim. Past*, 13, 833-853, doi:10.5194/cp-13-833-2017.
* Brils, M., Kuipers Munneke, P., van de Berg, W. J., and van den Broeke, M. (2022), Improved representation of the contemporary Greenland ice sheet firn layer by IMAU-FDM v1.2G, *Geosci. Model Dev.*, 15, 7121-7138, doi:10.5194/gmd-15-7121-2022.
* Brun, E., David, P., Sudul, M., and Brunot, G. (1992), A numerical model to simulate snow-cover stratigraphy for operational avalanche forecasting, *J. Glaciol.*, 38, 13-22.
* Buizert, C., et al. (2015), The WAIS Divide deep ice core WD2014 chronology - Part 1: Methane synchronization (68-31 ka BP) and the gas age-ice age difference, *Clim. Past*, 11, 153-173, doi:10.5194/cp-11-153-2015.
* Goujon, C., Barnola, J.-M., and Ritz, C. (2003), Modeling the densification of polar firn including heat diffusion: Application to close-off characteristics and gas isotopic fractionation for Antarctica and Greenland sites, *J. Geophys. Res.*, 108, 4792, doi:10.1029/2002JD003319.
* Helsen, M. M., van den Broeke, M. R., van de Wal, R. S. W., van de Berg, W. J., van Meijgaard, E., Davis, C. H., Li, Y., and Goodwin, I. (2008), Elevation changes in Antarctica mainly determined by accumulation variability, *Science*, 320, 1626-1629, doi:10.1126/science.1153894.
* Herron, M. M., and Langway, C. C. (1980), Firn densification: an empirical model, *J. Glaciol.*, 25(93), 373-385, doi:10.3189/S0022143000015239.
* Kuipers Munneke, P., Ligtenberg, S. R. M., Noël, B. P. Y., Howat, I. M., Box, J. E., Mosley-Thompson, E., McConnell, J. R., Steffen, K., Harper, J. T., Das, S. B., and van den Broeke, M. R. (2015), Elevation change of the Greenland Ice Sheet due to surface mass balance and firn processes, 1960-2014, *The Cryosphere*, 9, 2009-2025, doi:10.5194/tc-9-2009-2015.
* Li, J., and Zwally, H. J. (2002), Modeled seasonal variations of firn density induced by steady-state surface air-temperature cycle, *Ann. Glaciol.*, 34, 299-302.
* Li, J., and Zwally, H. J. (2004), Modeling the density variation in the shallow firn layer, *Ann. Glaciol.*, 38, 309-313.
* Li, J., and Zwally, H. J. (2011), Modeling of firn compaction for estimating ice-sheet mass change from observed ice-sheet elevation change, *Ann. Glaciol.*, 52, 1-7, doi:10.3189/172756411799096321.
* Li, J., and Zwally, H. J. (2015), Response times of ice-sheet surface heights to changes in the rate of Antarctic firn compaction caused by accumulation and temperature variations, *J. Glaciol.*, 61, 1037-1047, doi:10.3189/2015JoG14J182.
* Ligtenberg, S. R. M., Helsen, M. M., and van den Broeke, M. R. (2011), An improved semi-empirical model for the densification of Antarctic firn, *The Cryosphere*, 5, 809-819, doi:10.5194/tc-5-809-2011.
* Lundin, J. M. D., et al. (2017), Firn Model Intercomparison Experiment (FirnMICE), *J. Glaciol.*, 63, 401-422, doi:10.1017/jog.2016.114.
* Marshall, H. P., Conway, H., and Rasmussen, L. A. (1999), Snow densification during rain, *Cold Reg. Sci. Technol.*, 30, 35-41.
* Medley, B., Neumann, T. A., Zwally, H. J., Smith, B. E., and Stevens, C. M. (2022), Simulations of firn processes over the Greenland and Antarctic ice sheets: 1980-2021, *The Cryosphere*, 16, 3971-4011, doi:10.5194/tc-16-3971-2022.
* Morris, E. M., and Wingham, D. J. (2014), Densification of polar snow: Measurements, modeling, and implications for altimetry, *J. Geophys. Res. Earth Surf.*, 119, 349-365, doi:10.1002/2013JF002898.
* Pimienta, P., and Duval, P. (1987), Rate controlling processes in the creep of polar glacier ice, *J. Phys. Colloques*, 48, C1-243-C1-248, doi:10.1051/jphyscol:1987134.
* Simonsen, S. B., Stenseng, L., Adalgeirsdóttir, G., Fausto, R. S., Hvidberg, C. S., and Lucas-Picher, P. (2013), Assessing a multilayered dynamic firn-compaction model for Greenland with ASIRAS radar measurements, *J. Glaciol.*, 59, 545-558, doi:10.3189/2013JoG12J158.
* Steger, C. R., et al. (2017), Firn meltwater retention on the Greenland ice sheet: a model comparison, *Front. Earth Sci.*, 5, 3, doi:10.3389/feart.2017.00003.
* Stevens, C. M., Verjans, V., Lundin, J. M. D., Kahle, E. C., Horlings, A. N., Horlings, B. I., and Waddington, E. D. (2020), The Community Firn Model (CFM) v1.0, *Geosci. Model Dev.*, 13, 4355-4377, doi:10.5194/gmd-13-4355-2020.
* van Kampenhout, L., Lenaerts, J. T. M., Lipscomb, W. H., Sacks, W. J., Lawrence, D. M., Slater, A. G., and van den Broeke, M. R. (2017), Improving the representation of polar snow and firn in the Community Earth System Model, *J. Adv. Model. Earth Syst.*, 9, 2583-2600, doi:10.1002/2017MS000988.
* Veldhuijsen, S. B. M., van de Berg, W. J., Brils, M., Kuipers Munneke, P., and van den Broeke, M. R. (2023), Characteristics of the 1979-2020 Antarctic firn layer simulated with IMAU-FDM v1.2A, *The Cryosphere*, 17, 1675-1696, doi:10.5194/tc-17-1675-2023.
* Vionnet, V., Brun, E., Morin, S., Boone, A., Faroux, S., Le Moigne, P., Martin, E., and Willemet, J.-M. (2012), The detailed snowpack scheme Crocus and its implementation in SURFEX v7.2, *Geosci. Model Dev.*, 5, 773-791, doi:10.5194/gmd-5-773-2012.
* Yamazaki, T., Kondo, J., Sakuraoka, T., and Nakamura, T. (1993), A one-dimensional model of the evolution of snow-cover characteristics, *Ann. Glaciol.*, 18, 22-26, doi:10.3189/S0260305500011204.

.. VERIFY (references): checked against the PDFs in ~/papers, the reference list of Stevens et al. (2020), or (Pimienta and Duval 1987,
   Yamazaki et al. 1993) the HAL / Cambridge pages that Max gave. Not opened: Alley (1987) (the PDF has no text layer; details are from the
   Stevens et al. list). Note that the CFM v1.0 paper cites "Pimienta (1987)"; I used the HAL record Max gave (Pimienta and Duval, 1987).
