.. _json-page:

**************************************
The .json-formatted configuration file
**************************************

The CFM uses a .json-formatted file to configure individual model runs. JSON (JavaScript Object Notation) is a data-interchange file format. It consists of a number of names, each associated with a value. Values can be strings, Booleans, integers, floats, or arrays. Comments are not allowed, but can be added by considering the comment as a name/value pair. For the CFM, it provides a file format that is both easy to read and easy to alter in order to specify parameters for a particular model run. The configuration file is passed to the CFM, and the name/value pairs are read by the model and incorporated into the model run. The file format is editable in any text editor, and the name/value pairs are given by name: value, and different name/value pairs are separated by commas.

The specific names that are in the configuration .json file for the CFM are as follows. If any of the name/value pairs are missing, the model will generally return a message that that that name/value pair is missing and will use a default instead. For some name/value pairs the model run will fail. Note that in the .json file true/false are lowercase, but in the .py files they are True/False (first letter capitalized). The model automatically converts this. :math:`\rho_{s}` is the surface.

.json keys
~~~~~~~~~~

InputFileFolder
---------------
  Directory where the input csv files are located (usually a subdirectory of the directory that contains main.py, but user can specify an absolute paths as well.) Use '' if the input files are in the same directory as main.py.
      
  :type: ``string``
  :example: ``inputdata``


InputFileNameXXXX
-----------------
  The names of the input files for temperature, accumulation/smb, water isotopes, surface density, and melt. See 'Inputs for the CFM' section for more details.

  :type: ``string``
  :Example: ``example_XXXX.csv``

resultsFolder
-------------
  Folder in which results are stored.

  :type: ``string``
  :Example: ``example_results``

initfirnFile
------------
  File containing initial conditions if you are using firn measurements/data (e.g. temperature, density) to begin the model run. See 'Inputs for the CFM' section for more details.

  :type: ``string``
  :Example: ``example_firndata.csv``

initprofile
-----------
  Whether or not the CFM should use the initfirnFile to generate an initial condition.

  :type: ``boolean``
  :default: ``False``

input_type
----------
  (New in version 1.1.0)
  Specify what type of inputs you want to use - .csv (historic behavior) or pandas dataframe that is stored in a pickle.
  
  :type: ``string``
  :default: ``csv``
  :options: ``csv``,``dataframe``

DFresample
----------
  (New in version 1.1.0)
  Specify the resolution you want for your model run, which will be the resample interval for the dataframe (this only has functionality when input_type is ``dataframe``)
  See https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Timedelta.html

  :type: ``pandas Timedelta (string)``
  :Example: ``1D``

DFfile
------
  The filename of the pickle containing the climate dataframe.

  :type: ``string``
  :example: ``example.pkl``

physRho
-------
  The firn-densification physics to use for the model run.

  :type: ``string``
  :Options: ``HLdynamic``, ``HLSigfus``, ``Li2004``, ``Li2011``, ``Helsen2008``, ``Arthern2010S``, ``Arthern2010T``, ``Li2015``, ``Simonsen2013``, ``Goujon2003``, ``Barnola1991``, ``Morris2014``, ``KuipersMunneke2015``, ``Brils2022``, ``Veldhuijsen2023``, ``GSFC2020``, ``Crocus``, ``Ligtenberg2011``

MELT
----
  Whether or not to include meltwater percolation physics in the model run.

  :type: ``boolean``
  :default: ``False``

ReehCorrectedT
--------------
  Only used when **MELT** is True. If True, applies the Reeh (1991, 2008) latent-heat temperature correction during spin up: the mean firn temperature is warmed to account for the latent heat released by refreezing meltwater, using the superimposed-ice rate (capped at 0.6 of the annual accumulation, following Reeh's PMAX). The correction raised temperature is :math:`T + 26.6\,\textrm{SIR}` (bounded at 273.15 K). If **MELT** is False, the model exits with a warning.

  :type: ``boolean``
  :default: ``False``

FirnAir
-------
  Whether or not to run the firn air module with the model run.

  :type: ``boolean``
  :default: ``false``

AirConfigName
-------------
  Name of the .json configuration files that contains the parameters for the firn air module.

  :type: ``string``
  :default: ``AirConfig.json``

TWriteInt
---------
  How often to write the results to file, relative to the time-step size, i.e. 1 will write at every time step, 10 will write at every 10th time step, etc.

  :type: ``int``
  :default: ``1``

TWriteStart
-----------
  The time at which to start saving model results. The time is model time, so must correspond to the time in the input forcing files.

  :type: ``float``
  :Example: If your model run is from 1958 to 2018, but you only want outputs from 2000 onwards, 'TWriteStart' should be 2000.

int_type
--------
  How to interpolate from the input file times to the model time. Use linear e.g. if you have sparse ice core climate data. Use nearest e.g. if you have monthly climate data and want to take monthly time steps (the model time may not be exactly the same as the input time).

  :type: ``string``
  :options: ``nearest``, ``linear``

SeasonalTcycle
--------------
  Whether or not to add a seasonal temperature cycle (on top of the forcing data). Use this only if you are using sub-annual time steps and your forcing data does not have a seasonal cycle already. Usually this would be if your forcing data is annual (or coarser resolution).

  :type: ``boolean``
  :default: ``false``

SeasonalThemi
-------------
  If 'SeasonalTCycle' is True, specify which hemisphere you are modeling to get the summer/winter timing correct.

  :type: ``string``
  :options: ``north``, ``south``

coreless
--------
  Only used when **SeasonalTcycle** is True and **SeasonalThemi** is ``south``. If True, adds a "coreless winter" to the seasonal temperature cycle: a second harmonic (following Orsi) is superimposed on the annual cosine, giving the flat mid-winter temperature plateau characteristic of the Antarctic interior. Has no effect in the northern hemisphere.

  :type: ``boolean``
  :default: ``false``

TAmp
----
  If 'SeasonalTCycle' is True, specify the amplitude of the cycle.

  :type: ``float``
  :default: ``10``
  :units: :math:`K`

physGrain
---------
  Whether or not to track grain size evolution. Must be True for Arthern2010S physics.

  :type: ``boolean``
  :default: ``false``

calcGrainSize
-------------
  True uses a parameterization to get a surface grain-size at each time step, and False uses a set grain size at the surface. 

  :type: ``boolean``
  :default: ``false``

GrGrowPhysics
-------------
  Which equation to use to calculate grain size evolution.

  :type: ``string``
  :options: ``Arthern``,``Katsushima``

heatDiff
--------
  Whether or not to include heat diffusion.

  :type: ``boolean``
  :default: ``True``

conductivity
------------
  Which parameterization for heat conductivity to use.

  :type: ``string``
  :options: ``Schwander``,``Yen_fixed``,``Yen_var``,``Anderson``,``Yen_b``,``Sturm``,``VanDusen``,``Schwerdtfeger``,``Riche``,``Jiawen``,``mix``,``Calonne2011``,``Calonne2019``

variable_srho
-------------
  Whether to vary the surface density through time. False uses a constant density.

  :type: ``boolean``
  :default: ``False``

srho_type
---------
  If variable_srho is true, how to vary the surface density through time. 'userinput' uses a csv file with surface density though time (must be specified with InputFileNamerho); 'param' uses a parametization; 'noise' adds noise at each time step to the value specified by **rhos0**.

  :type: ``string``
  :options: ``userinput``,``param``,``noise``

rhos0
-----
  Surface density at each time step if using a constant surface density, or the mean value if **variable_srho** is true and **srho_type** is 'noise'.

  :type: ``float``
  :default: ``350.0``
  :units: :math:`\textrm{kg m}^{-3}`

r2s0
----
  Surface grain size at each time step if **calcGrainSize** is false.

  :type: ``float``
  :default: ``1e-8``
  :units: :math:`\textrm{mm}^{2}`

AutoSpinUpTime
--------------
  Calculate the spin up time automatically based on the input accumulation rate and specified model domain depth; should be long enough to refresh the entire firn column during spin up.

  :type: ``boolean``
  :default: ``false``

yearSpin
--------
  How many years to spin up for. Only used when **AutoSpinUpTime** is False; if **AutoSpinUpTime** is True, the spin-up length is calculated automatically and this value is ignored. The spin up should be long enough to refresh the entire firn column at least once (i.e. long enough for a parcel deposited at the surface to advect to the bottom of the domain).

  :type: ``float``
  :units: years

stpsPerYearSpin
---------------
  **DEPRECATED** How many time steps per year to take during spin up. Previously the CFM gave the option to have different values for spin up and main run; now spin up uses **stpsPerYear**.

stpsPerYear
-----------
  How many time steps per year to take. E.g. 12 will make the model take monthly time steps, 1 will give annual time stepping. Take care to coordinate this value with your input files and the 'int_type'.

  :type: ``float``

H
---
  Thickness of the ice sheet in meters. This is a bit confusing. Probably keep it at 3000 or so. That would mean the surface of the firn is 3000 m above the bed.

  :type: ``float``
  :default: ``3000``
  :units: :math:`\textrm{m}`

HbaseSpin
---------
  The elevation of the bottom of the model domain above the bed. So, if you want to model to 250 m depth, and H is 3000, HbaseSpin will be 2750. Likewise, if you wanted to model just the top 50 m of firn, HbaseSpin will be 2950 (assuming H is 3000). This is an initial value at the start of the spin up. The depth of the model domain will change due to the fact the model is Lagrangian with a fixed number of nodes; e.g. if the accumulation rate increases, each node will be thicker, and the base of the domain will be deeper.

  :type: ``float``
  :units: :math:`\textrm{m}`

D_surf
------
  The CFM features a generic layer tracker called *D_con*; it can be used for a number of things. This is the value to assign a new layer at the surface at each time step.

  :type: ``float``
  :default: ``1``

bdot_type
---------
  The type of accumulation rate to use for the densification physics. ‘Instant’ is the instantaneous value (i.e. at that time step) of accumulation, ‘mean’ is the mean accumulation over the lifetime of a parcel of firn. (‘Stress’ is in progress and will use the stress directly).

  :type: ``string``
  :default: ``mean``
  :options: ``mean``,``instant``,``stress``

grid_outputs
------------
  Whether or not to put the outputs on a regular grid (i.e. evenly spaced vs. the internal variable grid)

  :type: ``boolean``
  :default: ``True``

grid_output_res
---------------
  If grid_output is ``True``, this is the spacing of the grid nodes in meters.

  :type: ``float``
  :default: ``0.1``

isoDiff
-------
  Whether or not to include water isotope diffusion in the model run.

  :type: ``boolean``
  :default: ``False``

iso
---
  If isoDiff is true, which isotopes to model. ``NoDiffusion`` will include the isotopes but does not diffuse them at each time step, to allow analysis of the effects of advection and compaction alone (it uses the d18O forcing).

  :type: ``list of strings``
  :default: ``["d18O", "dD"]``
  :options: ``d18O``, ``dD``, ``NoDiffusion``

spacewriteint
-------------
  **DEPRECATED / NON-FUNCTIONAL.** Was intended to set the spatial resolution interval saved to results (1 = every node, 2 = every other, etc.), but is not read by the current code. To reduce output size, use **grid_outputs** / **grid_output_res** or **truncate_outputs** instead.

  :type: ``int``
  :default: ``1``

horizontal_divergence
---------------------
  Whether to include the effect of horizontal divergence (from dynamic ice-sheet/glacier flow) on the firn. When True, the mass in each layer is rescaled at each time step according to the horizontal strain rate, which thins (divergence) or thickens (convergence) the column. The strain-rate forcing is supplied via **InputFileNameStrain**.

  (Replaces the older ``strain`` key, which is still accepted and automatically converted.)

  :type: ``boolean``
  :default: ``False``

strain_softening
----------------
  Whether to include strain softening in the stage-2 (power-law creep) densification regime. When True, densification rates are scaled by the effective horizontal strain rate. Requires a strain-rate forcing via **InputFileNameStrain**.

  :type: ``boolean``
  :default: ``False``

residual_strain
---------------
  Regularization threshold for strain-softening calculations. Vertical strain rates with magnitude below this value are set to zero to avoid singularities.

  :type: ``float``
  :default: ``2e-4``
  :units: :math:`\textrm{a}^{-1}`

tuning_bias_correction
----------------------
  If True, applies a bias correction to the strain-softening scheme to account for the strain-softening signal already implicitly captured by the tuned Herron-Langway densification model. Only relevant when **strain_softening** is True.

  :type: ``boolean``
  :default: ``False``

outputs
-------
  Which outputs to save.

  :type: ``list of strings``
  :example: ``["density", "depth"]``
  :options: ``density``, ``depth``, ``temperature``, ``age``, ``Dcon``, ``bdot_mean``, ``climate``, ``compaction``, ``grainsize``, ``temp_Hx``, ``isotopes``, ``BCO``, ``DIP``, ``DIPc``, ``LWC``, ``PLWC_mem``, ``viscosity``, ``runoff``, ``refrozen``, ``meltoutputs``, ``gasses``

resultsFileName
---------------
  Name of the .hdf5 file that results are saved in.

  :type: ``string``
  :default: ``CFMresults.hdf5``

spinFileName
------------
  Name of the .hdf5 file that the spin up results are saved in.

  :type: ``string``
  :default: ``CFMspin.hdf5``

doublegrid
----------
  Whether or not to use the feature that keeps a high-resolution grid near the surface and a lower-resolution grid at greater depth.

  :type: ``boolean``
  :default: ``false``

nodestocombine
--------------
  If **doublegrid** is True, this is how many nodes are combined into a single node at the high/low resolution boundary. So, if it is 50, at every 50th time steps 50 nodes will be combined into a single node.

  :type: ``int``
  :default: 50

multnodestocombine
------------------
  If **doublegrid** is True, this is how many nodes are combined into a single node at the boundary between the low and very low resolution grid. For example, if nodestocombine is 50, multnodes will combine 'multnodestocombine' of those 50-node thick layers into a single node.

  :type: ``int``
  :default: 6

grid1bottom
-----------
  If **doublegrid** is True, the depth (m) at which the high-resolution grid nodes are combined.

  :type: ``float``
  :default: 10

grid2bottom
-----------
  If **doublegrid** is True, the depth (m) at which the low-resolution grid nodes are combined to make the very-low resolution grid.

  :type: ``float``
  :default: 20

spinup_climate_type
-------------------
  What climate to use for the spin up. 'initial' uses the very first value in the input .csv files and 'mean' uses the mean of the values in those files.

  :type: ``string``
  :options: ``initial``, ``mean``

manual_climate
--------------
  Manually specify the background climate (long-term means). This is useful if you are doing a very short model run, in which the input csv files may not be representative of the long-term climate.

  :type: ``boolean``
  :default: ``false``

deepT
-----
  If manual_climate is true, this is the long term site temperature (the temperature that would be measured at the bottom of a borehole).

  :type: ``float``
  :units: :math:`\textrm{K}`

bdot_long
---------
  If manual_climate is true, this is the long-term mean accumulation rate.

  :type: ``float``
  :units: :math:`\textrm{m ice eq. a}^{-1}`

manual_iceout
-------------
  Allows the user to specify the ice that is effectively removed from the bottom of the firn due to ice sheet thinning from ice flow. In steady state, iceout is the same as the long-term ice equivalent accumulation rate (and that is what is used if manual_iceout is false).

  :type: ``boolean``
  :default: ``false``

iceout
------
  If manual_iceout is True, this is the value.

  :type: ``float``
  :units: :math:`\textrm{m ice eq. a}^{-1}`

QMorris
-------
  The Morris and Wingham (2014) model allows for different activation energies; specify it here.

  :type: ``float``
  :default: ``110.0e3``
  :units: :math:`\textrm{kJ mol}^{-1}`

timesetup
---------
  How to set up the time step size. 'Exact' uses the input files to find the times at which a time step occurs and the corresponding time-step size *dt*; 'interp' uses a uniform *dt* and interpolates the input data onto the timeline that the model generates with uniform time steps. 'retmip' is a specialty for the RETMIP experiment and may not be fully functional.

  :type: ``string``
  :options: ``exact``, ``interp``, ``retmip``

liquid
------
  If **MELT** is true, which percolation scheme to use. ``bucket`` is the standard single-bucket scheme (``bucketVV`` / ``percolation_bucket`` are accepted aliases for the Verjans variant); ``darcy`` solves flow with Darcy's law; ``prefsnowpack`` and ``resingledomain`` are the preferential-flow and single-domain Richards-equation snowpack schemes (in development).

  :type: ``string``
  :options: ``bucket``, ``darcy``, ``resingledomain``, ``prefsnowpack``

meltwater_solver
----------------
  If **MELT** is true, which numerical scheme is used to solve heat diffusion with meltwater refreezing (i.e., the coupled temperature / latent-heat problem for wet firn). All of the schemes solve the same physics but differ in their numerical formulation and conservation properties:

  - ``enthalpy`` -- solves for enthalpy directly (default).
  - ``ahc`` -- apparent-heat-capacity method.
  - ``decp`` -- decoupled (operator-split) scheme.
  - ``ncz`` -- Newton-based scheme (Jordan-style).
  - ``legacy`` -- the enthalpy solver as it stood before mid-July 2026, retained to reproduce older results. Provided for backward comparison only; use one of the schemes above for new work.

  If the key is omitted, the model defaults to ``enthalpy``. For the four current schemes, refreezing is dispatched through ``refreezeDiff`` in ``diffusion.py``, which calls the corresponding ``transient_solve_*`` function in ``solver.py``; ``legacy`` instead calls ``enthalpyDiff_old`` (``diffusion.py``) and ``transient_solve_EN_old`` (``solver.py``), bypassing ``refreezeDiff``. Per-time-step solver diagnostics are written to ``diagnostics_<solver>.csv`` at the end of the run. See :doc:`../extras/meltwater_solver` for guidance on which scheme to choose.

  This key replaces the former **LWC_heat**, which selected the same thing and has been removed.

  :type: ``string``
  :default: ``enthalpy``
  :options: ``enthalpy``, ``ahc``, ``decp``, ``ncz``, ``legacy``

merging
-------
  If a model volume gets too thin, merge it with another. Needed for numerical stability with melt schemes.

  :type: ``boolean``
  :default: ``false``

merge_min
---------
  If merging is true, the thickness threshold at which merging should occur.

  :type: ``float``
  :default: ``1e-4``


manualT
-------
  Option to use manual temperature measurements, e.g. from a thermistor string.

  :type: ``boolean``
  :default: ``false``

no_densification
----------------
  Option to set densification to false (perhaps you are simulating temperature diffusion in a core in a lab)

  :type: ``boolean``
  :default: ``false``

rad_pen
-------
  Option to turn on radiation penetration module.

  :type: ``boolean``
  :default: ``false``

site_pressure
-------------
  Set the pressure at the site, which can affect isotope diffusion.

  :type: ``float``
  :default: ``1013.25``

output_bits
-----------
  Set the bits for the outputs.

  :type: ``string``
  :default:  ``float32``

spinUpdate
----------
  Specify if you want to update the spin file at some date.

  :type: ``boolean``
  :default: ``false``

spinUpdateDate
--------------
  Specify the date at which to update the spin file. should correspond to the start of your reference climate interval.

  :type: ``float``
  :default: ``1980.0``

DIPhorizon
----------
  Depth horizon at which to calculate DIP/FAC (because the bottom of the domain varies a bit).

  :type: ``float``
  :default: ``100.0``

NewSpin
-------
  Whether or not to perform a new spin up (if the spin file exists already.)

  :type: ``boolean``
  :default: ``false``

Input file names
~~~~~~~~~~~~~~~~~

The CFM reads each climate/boundary forcing from its own file (when **input_type** is ``csv``) or from a column of the climate dataframe (when **input_type** is ``dataframe``). The ``InputFileName*`` keys give the file name (relative to **InputFileFolder**) for each forcing. Only the files needed for the enabled physics are required.

InputFileNameTemp
-----------------
  Surface (skin) temperature forcing time series. Required.

  :type: ``string``
  :Example: ``example_TSKIN.csv``

InputFileNamebdot
-----------------
  Accumulation-rate (surface mass balance) forcing time series. Required.

  :type: ``string``
  :Example: ``example_BDOT.csv``

InputFileNamemelt
-----------------
  Surface-melt forcing time series. Required when **MELT** is True and **SEB** is False (when **SEB** is True, melt is computed internally).

  :type: ``string``
  :Example: ``example_SMELT.csv``

InputFileNameRain
-----------------
  Rain forcing time series. Required when **RAIN** is True.

  :type: ``string``
  :Example: ``example_RAIN.csv``

InputFileNameSublim
-------------------
  Sublimation/deposition forcing time series. Used when **SUBLIM** is True; if omitted, sublimation is inferred from negative values of the accumulation forcing.

  :type: ``string``
  :Example: ``example_SUBLIM.csv``

InputFileNameIso
----------------
  Surface water-isotope forcing time series. Required when **isoDiff** is True.

  :type: ``string``
  :Example: ``example_ISOTOPE.csv``

InputFileNamerho
----------------
  Surface-density forcing time series. Required when **variable_srho** is True and **srho_type** is ``userinput``.

  :type: ``string``
  :Example: ``example_RHOS.csv``

InputFileNameStrain
-------------------
  Horizontal strain-rate forcing, used when **horizontal_divergence** or **strain_softening** is True. The file may contain one column (divergence), two columns (the two principal strain rates), or three columns (:math:`\dot\epsilon_{xx}`, :math:`\dot\epsilon_{yy}`, :math:`\dot\epsilon_{xy}`). Replaces the deprecated ``InputFileNamedudx`` (still accepted and auto-converted).

  :type: ``string``
  :Example: ``example_STRAIN.csv``

ManualTFilename
---------------
  File containing a 2-D (depth × time) temperature field, e.g. from a thermistor string. Used only when **manualT** is True. The first row is the time vector, the first column is the depth vector, and the remaining entries are the temperature matrix.

  :type: ``string``

forcingFileName
---------------
  Name of the .hdf5 file into which the climate forcing is written alongside the model output.

  :type: ``string``
  :default: ``CFMforcing.hdf5``

Surface energy balance (SEB)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SEB
---
  Whether to run the surface-energy-balance module, which computes surface temperature and melt from energy fluxes rather than using a prescribed melt forcing. When True, melt is calculated internally within the time-stepping loop and **InputFileNamemelt** is not needed.

  :type: ``boolean``
  :default: ``false``

SEB_TL_thick
------------
  Thickness of the surface "top layer" over which the surface energy fluxes are integrated to compute melt. Smaller values make the surface temperature/melt response more sensitive. Only used when **SEB** is True.

  :type: ``float``
  :default: ``0.05``
  :units: :math:`\textrm{m}`

albedo_factor
-------------
  Scaling factor applied to the input albedo (albedo is multiplied by this value). Use for tuning surface reflectivity. Only used when **SEB** is True.

  :type: ``float``
  :default: ``1``

Melt, liquid water, and runoff
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These keys configure the bucket percolation scheme (**liquid** = ``bucket``).

RAIN
----
  Whether to include rain as a liquid-water input at the surface. Requires **InputFileNameRain**.

  :type: ``boolean``
  :default: ``false``

ColeouLesaffre
--------------
  How to set the irreducible water content (the fraction of pore space that retains water against gravity). If True, uses the Coléou and Lesaffre (1998) parameterization; if False, uses the constant value **IrrVal**.

  :type: ``boolean``
  :default: ``true``

IrrVal
------
  Irreducible water content as a fraction of available pore space. Used only when **ColeouLesaffre** is False.

  :type: ``float``
  :default: ``0.02``

RhoImp
------
  Density at or above which a layer is treated as an impermeable ice lens that blocks percolation.

  :type: ``float``
  :default: ``830``
  :units: :math:`\textrm{kg m}^{-3}`

ThickImp
--------
  Minimum thickness for an ice lens (density :math:`\geq` **RhoImp**) to be treated as impermeable. Set to 0 to make every ice layer impermeable. Used only when **DownToIce** is False.

  :type: ``float``
  :default: ``0.1``
  :units: :math:`\textrm{m}`

DownToIce
---------
  If True, meltwater bypasses all ice lenses and percolates down to the depth where density reaches **RhoImp**; if False, water stops at the first impermeable ice lens (see **ThickImp**).

  :type: ``boolean``
  :default: ``false``

Ponding
-------
  If True, meltwater that is blocked by an impermeable barrier ponds in the layers above it rather than running off.

  :type: ``boolean``
  :default: ``false``

DirectRunoff
------------
  Fraction of blocked (ponded) water that runs off immediately rather than ponding. Only used when **Ponding** is True.

  :type: ``float``
  :default: ``0.0``

RunoffZuoOerlemans
------------------
  If True, lateral runoff of ponded water is computed with the Zuo and Oerlemans (1996) parameterization. Only used when **Ponding** is True.

  :type: ``boolean``
  :default: ``false``

Slope
-----
  Surface slope (rise/run) used in the Zuo and Oerlemans (1996) lateral-runoff calculation. Only used when **RunoffZuoOerlemans** is True.

  :type: ``float``
  :default: ``0.1``

keep_firnthickness
------------------
  Controls grid behavior when melt removes surface layers. If True, the domain thickness is maintained by adding nodes at the base; if False, the original layer thicknesses are kept.

  :type: ``boolean``
  :default: ``false``

LWC_heat
--------
  **Deprecated -- replaced by meltwater_solver.** This key formerly selected the method for handling the latent heat of liquid water during heat diffusion, duplicating the job of **meltwater_solver** (see above). Use **meltwater_solver** instead; the ``highC`` / ``Teff`` / ``LWCcorr`` branches have been retired.

  For backward compatibility, a config that still contains ``LWC_heat`` is handled as follows, with a warning printed either way:

  - ``LWC_heat`` present and **meltwater_solver** absent: the run uses ``meltwater_solver = legacy``. Such a config predates the rename, and therefore also predates the 2026-08 change to the ``enthalpy`` solver, so falling back to ``legacy`` reproduces the numerics the config was written against rather than silently substituting a different scheme. If the ``LWC_heat`` value was one of the retired options, an additional warning notes that ``legacy`` is the old ``enthalpy`` scheme and not a like-for-like replacement.
  - Both keys present: ``LWC_heat`` is ignored and **meltwater_solver** governs.

  Remove ``LWC_heat`` from your .json to silence the warnings. The related internal flags ``LWCheat``, ``LWCcorr_subdt``, and ``correct_therm_prop`` are likewise inactive and kept only for backward compatibility.

Sublimation
~~~~~~~~~~~

SUBLIM
------
  Whether to include sublimation/deposition. When False, negative values in the accumulation forcing (interpreted as sublimation) are set to zero.

  :type: ``boolean``
  :default: ``true``

bdm_sublim
----------
  Whether to include sublimation when computing the mean accumulation rate used for spin up. If True, sublimation is subtracted from accumulation; if False, accumulation alone is used.

  :type: ``boolean``
  :default: ``true``

Densification and thermal options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MQ
--
  Activation energy for the Morris and Wingham (2014) densification model. Only used when **physRho** is ``Morris2014``.

  :type: ``float``
  :default: ``60``
  :units: :math:`\textrm{kJ mol}^{-1}`

THist
-----
  Whether to track the temperature history of each layer (required by physics such as Morris2014). Set automatically to True when **physRho** is ``Morris2014``.

  :type: ``boolean``
  :default: ``false``

stage_zero
----------
  **In development.** Enables a "stage-zero" (fresh-snow) densification stage below the transition density **s_zero_rho**, using the parameterization selected by **snow_model**, before the usual firn densification physics take over.

  :type: ``boolean``
  :default: ``false``

snow_model
----------
  Which stage-zero (snow) densification parameterization to use. Only used when **stage_zero** is True.

  :type: ``string``
  :default: ``Yamazaki1993``

s_zero_rho
----------
  Transition density separating the stage-zero snow regime from firn densification. Only used when **stage_zero** is True.

  :type: ``float``
  :default: ``200``
  :units: :math:`\textrm{kg m}^{-3}`

iceblock
--------
  If True, the firn column is initialized at a uniform density (**iceblock_rho**) instead of using the Herron-Langway analytic spin-up profile. Useful e.g. for simulating a solid ice block.

  :type: ``boolean``
  :default: ``false``

iceblock_rho
------------
  Uniform initial density used when **iceblock** is True.

  :type: ``float``
  :default: ``917``
  :units: :math:`\textrm{kg m}^{-3}`

Outputs and run metadata
~~~~~~~~~~~~~~~~~~~~~~~~~

truncate_outputs
----------------
  If True, reduces output file size by writing only every 5th depth node for the large fields (density, temperature, LWC, age). If False, every node is written.

  :type: ``boolean``
  :default: ``false``

runID
-----
  Optional numeric identifier used to label a run, convenient for batch runs and parameter studies.

  :type: ``int`` or ``float``
  :default: ``-9999``

lat_val
-------
  Latitude of the site being modeled. Used for bookkeeping/metadata and for site lookups in the example scripts.

  :type: ``float``
  :units: degrees

lon_val
-------
  Longitude of the site being modeled. Used for bookkeeping/metadata and for site lookups in the example scripts.

  :type: ``float``
  :units: degrees

lat_int / lon_int
-----------------
  Requested latitude/longitude used to select the nearest climate grid point when building forcing from a regional climate model (see the example run scripts). **lat_val** / **lon_val** hold the actual coordinates of the selected grid point.

  :type: ``float``
  :units: degrees


















