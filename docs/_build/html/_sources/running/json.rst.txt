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
  :Options: ``HLdynamic``, ``HLSigfus``, ``Li2011``, ``Helsen2008``, ``Arthern2010S``, ``Arthern2010T``, ``Li2015``, ``Goujon2003``, ``Barnola1991``, ``Morris2014``, ``KuipersMunneke2015``, ``Crocus``, ``Ligtenberg2011``

MELT
----
  Whether or not to include meltwater percolation physics in the model run.

  :type: ``boolean``
  :default: ``False``

ReehCorrectedT
--------------
  If melt is enabled, (NEED TO FILL IN WHAT THIS DOES)

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
  If 'SeasonalTCycle' is True, add the coreless winter. ADD MORE INFO HERE

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
  How many years to spin up for. COULD EXPAND ON THIS

  :type: ``float``

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
  If isoDiff is true, which isotopes to model. 'NoDiffusion' will include the isotopes but does not diffuse them at each time step to allow analysis of the effects of advection and compaction alone (it uses the d18O forcing).

  :type: ``list of strings``
  :default: ``["18", "D", "NoDiffusion"]``
  :options: ``18``, ``D``, ``NoDiffusion``

spacewriteint
-------------
  NOT WORKING CURRENTLY. Spatial resolution interval to save to results. 1 is every node; 2 is every other, etc.

  :type: ``int``
  :default: ``1``

strain
------
  Whether or not to include layer thinning due to horizontal strain from dynamic ice sheet/glacier flow.

  :type: ``boolean``
  :default: ``False``

du_dx
-----
  If strain is true, this is the horizontal strain rate. Future work will allow this to vary in time. NEED TO CHECK UNITS ARE CORRECT.

  :type: ``float``
  :default: ``1e-5``
  :units: :math:`\textrm{m a}^{-1}`

outputs
-------
  Which outputs to save.

  :type: ``list of strings``
  :example: ``["density", "depth"]``
  :options: ``density``, ``depth``, ``temperature``, ``age``, ``dcon``, ``bdot_mean``, ``climate``, ``compaction``, ``grainsize``, ``temp_Hx``, ``isotopes``, ``BCO``, ``LIZ``, ``DIP``, ``LWC``, ``gasses``

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
  If **MELT** is true, which percolation scheme to use.

  :type: ``string``
  :options: ``percolation_bucket``, ``bucketVV``, ``resingledomain``, ``prefsnowpack``

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


















