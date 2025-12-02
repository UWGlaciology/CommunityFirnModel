# CFM Change Log
All notable changes to the Community Firn Model should be documented in this file. Contributors to the CFM who are unfamiliar with changelogs should review the notes at the end of this document.

TL;DR: Write down the changes that you made to the the model in this document and update the version number here, in CITATION.cff, and in main.py, then update main branch on github.

General update protocol:
- First, get the staging branch to match main branch.
- Then, while on staging branch checkout files from whatever branch you want (e.g., git checkout dev melt.py)
- While on staging, test functionality. At minimum, python main.py example.json should work.
- If there are new keys for the json, update: *example.json*, *example_df.json*, *example_csv.json*, and *run_CFM_example_notebook.ipynb*

### Tests to do prior to release:
1) Example run with csv input:
    a) set "input_type": "csv" and "resultsFolder": "CFMoutput_example/csv" in example.json
    b) run python main.py example.json -n (make sure to do -n!)
2) Example run with df input:
    a) set "input_type": "dataframe" and "resultsFolder": "CFMoutput_example/df" in example.json
    b) run python main.py example.json -n (make sure to do -n!)

Then, switch to main branch and merge from staging: 
git checkout main
git merge staging

Then, double check that changelog is updated on main.

Then:
git commit -a -m "updating to vX.Y.Z. Details in changelog."
git push
git tag -a vX.Y.Z -m "CFM version vX.Y.Z"
git push origin vX.Y.Z

Then, on github do a release, which will trigger an updated DOI. 

## Current Version
3.1.0

## Full Documentation

https://communityfirnmodel.readthedocs.io/en/latest/

## Work in progress and known issues

- *Issues* 
	- If data is not written at each time step, dH that is saved/written to file does not represent the change since the last write. 
	- The dH output does not sum to zero when the model reaches or approaches steady state.
	- The siteClimate_from_RCM.py file is a mess; if you want to use it the best thing to do is to get in touch with me (maxstev@umd.edu) and I can help you get it running for your needs.


- *Work in progress*
	- If data is not written at each time step, dH that is saved/written to file does not represent the change since the last write. 
	- Testing percolation modules from Vincent Verjans that solves Richard's equation and includes a dual-domain approach to handle preferential flow (these modules are included but may not work properly yet)
	- Melt will likely be changed to have its own class
	- Documentation for the CFM
	- Goujon physics work, but could possibly be implemented more elegantly (it would be nice to avoid globals)
	- Not exactly in progress, but at some point adding a log file that gets saved in the results folder would be a good idea.
	- I am working on adding additional physics to simulate near-surface snow compaction ("stage zero compaction")
	- I am working on adding turbulent flux calculations to the SEB module


## [3.1.0] 2025-09-23
### Notes
- This release has a number of small feature updates and bug fixes.
- Much of the CFM-related work in the last several months has focused on developing scripts to do gridded runs of the CFM over the ice sheets on HPC clusters. This is a work in progress, but those scripts can be found in this repository: https://github.com/maximusjstevens/ATL_masschange
- Note that those scripts include python scripts and slurm batch scripts (.j files) to leverage HPC resources. If you have use for these scripts and have questions, please email me at maxstev@umd.edu

### New
- *RCMpkl_to_spin.py, firn_density_spin.py* These scripts both use mean annual accumulation rate (bdot_mean) to calculate the initial density profile and the grid structure. In areas of the ice sheet where sublimation is greater than the accumulation, the bdot_mean would be calculated as negative, and thusly the CFM would fail. There is a new key in the .json config called "bdm_sublim" (bdot mean sublimation). When bdm_sublim is True, sublimation is included in the calculation of bdot_mean. When False, sublimation is excluded. In this case (False), the grid is initialized using the mean ansnowfall snowfall, which will create a grid with thicker layers. During the model run, the total domain thickness will decrease substantially. The net effect of this is still not fully tested, but potentially could lead to a scenario in which the spin up (if automatically calculated in RCMpkl_to_spin.py) will not be long enough. Likely the best option in this case is to use the new "iceblock" feature.
- *firn_density_spin.py* There is a new option to just initialize the firn column with a constant density. This is controlled with the "iceblock" key in the configuration. The default value is 917 kg/m3, but this can also be customized using the key "iceblock_rho". This is recommended in high-melt areas or where sublimation is greater than snowfall. This sometimes fails when there are very thin layers - it seems that the numerics of a large density contrast in combination with thin layers causes a numerical in the SEB calculation.
- *firn_density_nospin.py* There is a new feature to add a densification scheme from snow model to low-density snow on top of the firn. I am calling this "stage zero densification", as the classic (e.g., Herron and Langway) literature defines stage one densification as being from the surface to the 550 kg/m3 density horizon. There are three new keys in the .json config file to enable this: *"stage_zero"* (True/False) is whether or not to use this feature, *"snow_model"* is the name of the densification model to use (default Yamazaki1993; same formatting and functionality as the present *physRho* key), and *"s_zero_rho"* is the stage zero/stage one transition density (default 200, but this will likely change as this model feature is further developed).
- *writer.py* New feature ("truncate_outputs") that will write only a subset of the matrix outputs ('rho','Tz','LWC','age') to save memory. The subsetting is along the time dimension, so if truncate_outputs is true, it only writes those outputs at every 5th time step. Vector outputs (e.g., DIP) are still saved at every time step.

### Changed
- *RCMpkl_to_spin.py* now foregoes the creation of a pandas dataframe (in "option 3") and puts the climate forcing data straight into a dictionary, which saves a bit of memory.
- *SEB.py* small changes to implementation of FQS solver to improve cross-platform compatibility
- *diffusion.py* small edit to the way 'z_dummy' is calculated, which is used to calculate finite volume centers.
- *diffusion.py* changed variable names from z_P_vec and z_edges_vec to z_P and z_edges, which makes it easier to track them into solver.py.
- *firn_density_spin.py* when using climate forcing from a dictionary (i.e., climateTS), we now assign the time written to the spinup file (CFMspin.hdf5) to match the start of climate forcing.
- *siteClimate_from_RCM.py* Added code to load MAR data, including dealing with varible "MSK" in MAR data and unit calculation for mass fluxes in MAR (bottom of MAR section in code - double check the MAR data being used!)
- *siteClimate_from_RCM.py* changed defaults for SEB and melt in function "getClimate" to True
- *solver.py* moved code that sets up grid in enthalpy solver outside of the iteration loop
- *writer.py* changed code to no longer write model forcings in main results. Now, forcing data is written to a separate forcing file (CFMforcing.hdf5, or similar).
- *firn_density_nospin.py, firn_density_spin.py* There is a small change in how these scripts call the densification physics class. Previously, the syntax was to call a particular densification scheme directly with its name, e.g.: *FirnPhysics(PhysParams).HL_dynamic*. To implement this, all possible options had to be coded to be stored in a dictionary, which was verbose in the code and made it difficult to add new physics options. Now, the code uses *getattr* as a generic way to instantiate the class. As a result, new densification physics schemes can now be easily added by just adding code to *physics.py*. The name of the physics scheme in *physics.py* must now match the name specified in the config file (previously, these were different, e.g. HLdynamic vs. HL_dynamic, for some reason).

### Fixed
- *firn_density_nospin.py, firn_density_spin.py* fixed bug in code the converts temperature to Kelvin to avoid (rare) situation where first temperature was above zero but units provided were celsius. (This happed in a few cases using mountain glacier air temperature data that started in melt season.) Now, code checks if the mean of the temperature time series is negative, rather than the first value of the time series.
- *physics.py* small bug fix to viscosity calculation in Li physics
- *regrid.py* fixed an issue in which under high melt or sublimation scenarios there would not be enough layers in "grid1" to split, and the length of the grid would reduce and cause a model failure.


## [3.0.0] 2024-10-15
### Notes
- This is a major release and will not necessarily be backwards compatible with version 2 depending on how the CFM scripts are being called. 
- The biggest changes are: (1) an improved SEB module, which uses a 'sub' time step to calculate melt at a higher temporal resolution; and (2) the spin routine now saves the model state periodically so that CFM can be restarted if a run fails without having to start again from the beginning.
- There are now two example json files - one for running with dataframe inputs, and one for running with csv inputs. 
- There is a new .py script and a jupyter notebook that are examples/templates for running the CFM without first editing a .json and then calling from the command line - they are a bit more "all in one" methods to configure and then run.

### New
- *SEB.py* includes a "sub" time step for the SEB calculation. So, melt and temperature are calculated at a higher temporal resolution (e.g., 4 hourly) to capture diurnal signal, but main timestepping loop in CFM stays the same (e.g., 1 day). The melt from the sub steps is summed, and the temperature is averaged. *RCMpkl_to_spin.py* now includes code to set up the fluxes for those sub time steps.
- *writer.py, firn_density_nospin.py* model forcing now gets written to its own hdf5 file at end of run, sublim added to forcing that gets written.  
- *firn_density_nospin.py* added feature that spinupdate will periodically write model state into spin file (e.g., every 100 years). If the model shuts off it can be restarted from its last write during spin up. All of the periodic writes are saved, so it is recommended that the perodicity is not too frequent.
- *solver.py* new tridiagonal solver from lapack. Speeds solver significantly. Thanks to Ben Smith (UW APL).
- Added *run_CFM_example.py* and *run_CFM_example_notebook.ipynb*, which can be used to do some example runs, or they can be edited to help you set up your own run.
- *json config* added 'TL_thick': thickness of the top layer to consider for SEB calculations.

### Changed
- *diffusion.py/solver.py* Previously, the finite volume solver centered the volumes on the z coordinates, (e.g., z = 0, 0.1, 0.2 would have dz=0.1, and the edges of the volumes for the solver were 0.05, 0.15, 0.25, etc.). It makes more intuitive sense for the layers and finite volumes to be the same - so a finite volume should extend from z=0 to z=0.1 (i.e., those are the edges), and the volume centers are at 0.05, etc. Testing indicates that this makes negligible difference in the modeled temperature. (have not thoroughly tested with LWC). This means that in the model outputs, if e.g. depth = 0, 0.1, 0.2, ..., and density = 350, 375, 400, ..., the firn layer from 0 to 0.1m has density of 350, and the layer from 0.1 to 0.2 has density 375, etc.
- *firn_density_nospin.py* In concert with the SEB updates, __init__ now takes SEBfluxes as in input, which is a dictionary similar to climateTS. It contains the energy fluxes at the time resolution for the "sub" time step. If SEBfluxes are not provided but SEB is on, the energy fluxes need to be in climateTS and will be at the same resolution as the other climateTS variables.
- *firn_density_nospin.py* added code to track mass conservation of liquid input/output
- *ModelOutputs.py, firn_density_nospin.py* "climate" key now has ,bdot,Ts,snowmelt,rain,sublim
- *siteClimate_from_RCM.py* getClimate() now returns df_CLIM (as before), and also the spin date start and end dates.
- *reader.py* read initial condition from most recent restart (model state) in spin file.
- *siteClimate_from_RCM.py* added functionality to get SEB fluxes from RCM/GCM data
- *siteClimate_from_RCM.py* now returns a dictionary rather than a dataframe; the dictionary includes spin date start and end
- *RCMpkl_to_spin.py* now includes calcSEB, a general solver to calculate skin temperature and melt flux based on energy inputs, and includes FQS solver to get the melt from those fluxes

### Fixed
- *firn_density_nospin.py* fixed an issue with gridtrack. self.gridtrack is now set to None if doublegrid is not turned on. merge module now uses gridtrack, which fixes an issue in which the tracking scheme would not update when layers were merged, thereby mis-assigning grid track numbers and potentially leading to incorrect regridding with the doublegrid scheme.
- *firn_density_nospin.py* fix to get Morris physics more user friendly (set THist = True automatically)
- *firn_density_nospin.py* fix issue to turn rain off if not in climateTS
- *firn_density_nospin.py* fixed issue where spinupdate would continually update the spinfile if using it over and over again.
- *firn_density_nospin.py* fixed issue that was incorrectly adjusting grid at timesteps that had high melt and 0 snow accumulation.
- *melt.py* fixed issue that would cause melt to fail if grainsize was turned off
- *melt.py* fixed issue with impermeable layers; automatically set the bottom layer of the grid to be impermeable regardless of density
- *merge.py* fixed issue that would cause melt to fail if grainsize was turned off
- *merge.py* added gridtrack to merge routine


## [2.3.2] 2024-03-04
### Notes
- This minor release fixes an issue in 2.3.0 - I did not manage to push the correct time step changes as described in the 2.3.0 release notes.

### Fixed
- Issue that dt and modeltime were not the same length vectors. 

## [2.3.1] 2023-11-02
### Notes
- Changes for 2.3.0 did not all merge correctly.

## [2.3.0] 2023-11-02
### Notes
- This version fixes two issues: when using the 'spinupdate' feature, the code would take the first time step off of each subsequent run when spinupdate was enabled. Now it does not. There was also an issue that merge.py did not correctly change the gridtrack variable when merging layers.
- I also made a change to how CFM handles time steps when timesetup is set to 'exact'. The forcing data has a decimal date, and differencing the input times gives a vector (dt) of time-step size. Because length of dt is one less than length of input time, CFM previously used input_time[1:]. Now, CFM appends a value (the mean of dt) to the start of dt. The upshot is where previously the values in the results represented the values in the firn at the start of the time interval, now they represent the values at the end. E.g., if you are using daily time steps, previously the results were the state of the firn at the start of that day, and now they are the state of firn at the end of the day. The reason I made this change is that now the forcing vectors line up exactly with dt, whereas before they were offset by one timestep.

### Fixed
- *merge.py* Merge.py now includes code to corretly reassign gridtrack values to each layer when it runs. The error only occured in isolated scenarios when melt or sublimation would leave the top layer very thin.
- *firn_density_nospin.py* The code has been changed to only trigger a spinup file update if it is not the first timestep of the model run. This leaves a typical workflow of: run full model with spinup; at the end of the spinup spinup file (i.e. the restart) gets updated; any subsequent runs will use that restart but not update it. Previously, if 'spinupdate' was true, the restart in the subsequent runs would get updated during the first time step, thereby cutting off the first time step from the previous run. 

## [2.2.0] 2023-06-26
### Notes
- This version fixes an issue in the SEB module where new snow that was added was too warm. It also updates the enthalpy solving routine. There are numerous small changes and fixes, detailed below.

### Fixed
- *SEB.py, firn_density_nospin.py* There was an bug when running using the SEB module that new snowfall was always at the melting temperature (which led to very warm firn). It now correctly sets the new snow temperature to be the minium of the freezing temperature and the air temperature

- *SEB.py* Fixed an issue to set the albedo to some value if it is NaN in the input (MERRA-2 albedo is assinged as NaN when there is no SW flux). Even with no SW flux the NaNs caused an issue in the SEB calculation.

- *RCMpkl_to_spin.py* Fixed an issue that the domain depth was always being set to be the depth of the 916 density horizon, rather than being configurable.

- *melt.py* Fixed an issue for timesteps where there was rain input but no melting, in which the grid was still being altered (removing top layers, adding new layers to the bottom to keep number of layers constant). Now rain-only timesteps do not alter the grid structure (just adds mass).

- *sublim.py* Fixed an issue where the sublimated mass was an numpy array with a single value rather than a float, which caused issues elsewhere when using that value.

- *firn_density_nospin.py* Fixed an issue where surface temperature was concatenated (effectively advecting the temperature profile) rather than reset during timesteps without accumulation. 

### Changed
- *SEB.py* There is a new solver built into the code (FQS, or fast quartic solver) to find the surface temperature based on the previous surface temperature. The code (presently just in SEB.py, not as an option in the config file) uses a defined thickness (presently 1 cm) as the top layer from an SEB perspective; i.e., the energy is assumed to warm/cool that layer, and that layer's temperature is set to the new temperature at the end of the SEB routine (regardless of the resolution of the CFM grid). (Note that there are several other SEB solving schemes, but FQS is the most tested code. In theory they should all give the same value.)

- *firn_density_spin.py* The firn column is initialized as a solid-ice column is ReehCorrectedT is true in the .json. The effect of this should be removed during spin up, but in areas around the ELA the column was sometimes disappearing or becoming too thin.

- *melt.py, sublim.py* Previously, the CFM kept the number of grid layers constant when there was melt or sublimation by dividing the bottom layer into thinner layers. In places with a lot of melt this led to the column becoming very thin. Now, the model instead adds layers to the bottom that have the same properties (temperature, density) as the current bottom layer and the same thickness as was melted/sublimated away. Note that if the bottom of the domain is not close to the ice density this will affect the DIP calculation. This is toggled in the .json as 'keep_firnthickness'. 'keep_firnthickness' = False is the old behavior. 

- *example.json* Addition of 'keep_firnthickness' (see above).

- *solver.py* Another update to the enthalpy solver. Previously, the overshoot factor applied to each iteration was applied to the new enthalpy solution, and now it just applied to the liquid portion. Additionally, the iteration's break statement is now prior to the overshoot, which (seems to) solves an issue where the solution was bouncing between values on each iteration instead of converging.

### Removed
- *plotter.py* Plotter.py was written long ago and not updated. There are better ways of plotting the model outputs so I removed this file from the repository. Please let me know if you need a script or jupyter notebook to help with output processing/plotting and I can send you something.

## [2.1.0] 2023-05-31
### Notes
- This version fixes an issue with how DIP was saved and updates the isotope diffusion module to allow inputs from a dictionary (parallel structure to loading temperature, accumulation rate, etc.)

### Changed
- *firn_density_nospin.py, isotopeDiffusion.py, firn_density_spin.py* These files have been updated to import the isotope forcing data in the climate import dictionary (climateTS).

## [2.0.0] 2023-02-28
### Notes
- This is the first major version release (i.e., to 2.x.x) due to several large changes. The first is the addition of a surface energy balance module, SEB.py. The second is a major overhaul in the enthalpy solver.

- The master branch is being renamed to the main branch.

### New
- *SEB.py* The new surface energy balance module caculates surface temperature and melt volume based on albedo, shortwave and longwave fluxes, and turbulent fluxes (those fields come from an RCM or AWS.) Presently does not calculate the turbulent fluxes from humidity and wind but I will add that capability in the future. THIS MODULE SHOULD BE CONSIDERED TO BE IN BETA PRESENTLY. IF YOU ARE USING IT, I SUGGEST GETTING IN TOUCH WITH ME AND WE CAN DISCUSS DETAILS.

- *firn_density_nospin.py, solver.py, diffusion.py* There is now an option to solve heat/refreezing in different ways. The enthalpy method is still present, but there are three other options now to compare soling methods. These features are in beta and will be better described in a planned paper. Feel free to email me to ask for details.

### Changed
- *firn_density_nospin.py, solver.py, diffusion.py* The enthalpy solver has been updated, again. The crux of this problem is that solver has to iterate to converge on a solution. At the end of the iteration, the liquid fraction must be updated, which is non trivial because a layer can have some mass refreeze but still be at the freezing temperature. The upshot is that in testin the new solver produces warmer firn than the older solver. I'd be happy to chat details about this if you have questions. 

- *firn_density_nospin.py* Temperature history, THist, was previously only calculated when using Morris 2014 physics; now it can be calculated with any densification scheme. This is for e.g. allowing experiments looking at how the temperature history may have affected other properties like grain size.

- *physics.py* The new densications from the Utrecht group (Brils et al., 2022, Veldhuijsen et al., 2023) have been added.

## [1.1.11] 2022-12-13
### Notes
- This release fixes an issue in the estimated surface elevation change *dh*, which was not accounting properly for elevation change due to sublimation and melt processes (it was just considerine dh due to firn compaction and new snow acccumulation). The new code explicitly includes *dh_melt* and *dh_acc*, which are the elevation change (for that time step) due to melt and accumulation + sublimation. **The elevation change calculation should be considered to be in beta.** 

### Changed
- *firn_density_nospin.py, melt.py, sublim.py* firn_density_nospin's *update_dh* function is now:
*dH = (sdz_new - sdz_old) + dh_acc + dh_melt - (iceout \* t[iii])*
The sdz terms are the sum of the layer thicknesses before and after compaction (different is thus dh from firn compaction); dh_acc is the elevation change due to new snow accumulation minus the sublimated volume; dh_melt is the elevation decrease due to surface melt; iceout (rate of m ice e.q./year) times the time step size is the elevation change due to ice flow. Note that this assumes steady state, and that the ice flow is calculated using the spin up climate (iceout is set to be the mean ice-equivalent accumulation rate during the spin up), unless set explicitly in the config file. 

>> *dh_melt* is calculated in *melt.py*, and *dh_acc* is calculated using dh_sub, which is now returned by *sublim.py*.

## [1.1.10]
### Notes
- Version 1.1.10 is a minor update to fix an issue with the strain softening routine introduced in v1.1.9. (The code would throw an error if InputFileNamedudx was not in the .json).
- There is a small update to the documentation regarding BCO outputs.

## [1.1.9]
### Notes
- Version 1.1.9 is an update to add strain softening as described by Oraschewski and Grinsted (2021; https://doi.org/10.5194/tc-2021-240) 

## [1.1.8]
### Notes
- This is a minor update. The main point is to add example .csv files and update the input climate .pkl file. The example .csv files were made using the .pkl files, so the outputs should match.
- There is a readme.txt file in the CFMinput_example directory with some details about those forcing files.

### Fixed
- *isotopeDiffusion.py* This file had not been updated to deal with the 'updatedStartDate' feature, which it now does. Current fuctionality limits its use to using .csv files as inputs.

### Updated
- Documentation: I updated the model output documentation.

## [1.1.7] 2022-04-19
### Notes
- This update fixes an issue where sublimation was automatically turned off.
- Hopefully this is a short-lived release; the next release should include an improved enthaply solver to avoid the ICT threshold error and still run quickly. A surface energy balance module is also in development.

### Fixed
- *firn_density_nospin.py, sublim.py* firn_density_nospin had a line that set any accumulation below a threshold to 0. This prevented sublimation from occuring because sublimation flux was inferred from the assumption that any negative values in the accumulation input was sublimation. The code now explicitly takes sublimation inputs, either through a field in the input dictionary (climateTS in firn_density_nospin.py) or a .csv file ('InputFileNameSublim' in .json). The .json should now include a boolean key/value 'SUBLIM' (default True). If that key is not in .json the code will automatically set to be true. (The design is that sublimation can be turned off, but only deliberately.) The code retains the ability to infer sublimation from negative values of accumulation (it does this if (1) inputs come from climateTS but SUBLIM is not a field in climateTS; or (2) if 'InputFileNameSublim' is NOT in the .json.)
- Note that sublim.py still does not alter the temperature. The reason for this presently is that we assume that the surface (skin) temperature calcuated by the forcing RCM (e.g. TS in MERRA-2) already accounts for that energy balance, so we just set the surface temperature to equal the input skin temperature at that time step.
- *physics.py* The viscosity calculations would throw a divide by zero error for layers that are at the ice density. This is now fixed so that those layers will have viscosity equal to zero (which is of course not true, but can be dealt with in post processing.)


## [1.1.6] 2021-11-22
### Notes
- The enthalpy solver scheme's routine to correct the amount of LWC after each iteration was wrong, which caused the firn to densify too quickly when there was melt present.

### Fixed
- *solver.py* See above note. The g_liq correction at the end of the iteration loop is reverted to what was used in v1.1.2.  

## [1.1.5] 2021-11-17
### Notes 
- There was an issue with the previous release in the order of operations within the model when melt was enabled. Prior to CFMv1.1.3, within a time step the firn first densified, then the melt routine occured, then heat diffusion, and then addition of a new snow layer with some density and temperature. Within the melt scheme, layers would melt, and the uppermost layer after melt (ie. the new surface) would take on the temperature of the old surface (which would often be colder than the melting temperature), which was not realistic. I switched the code so that that layer (called the 'partial melt' layer in the CFM) would have temperature of 273.15. However, this caused the surface temperature to always become the melting temperature, and the cold from a new snow layer would not ever diffuse into the firn if the next time step included melt. The solution was to move the diffusion to the end of the time step, so now the process goes: densification, meltwater percolation and refreezing (due to cold content), addition of new layer, heat diffusion. (Also: admitedly in reality the skin temperature should be zero during a time step when there is melt, but in model world with larger time steps, e.g. daily, there can be melt but the mean temperature for the day is still below freezing.)

### Fixed
- *ModelOutputs.py* There was an issue with the gridding feature for liquid water content, which is fixed. (The interpolation was just a linear interpolation of LWC, but it needs to be calculated by linearly interpolating the cumulative sum and then differencing to get the mass correct.)

### Changed
- *firn_density_nospin.py* See above under notes. Heat diffusion now comes at the end of the time step loop.
- *melt.py* When ponding is turned on, layers that reached impermeable density during the refreeze process are now set to have zero available pore space for accomodating excess water. These layers could form above a volume of excess liquid. This change caused some nodes to have very small negative LWC, which was a rounding issue; now the code sets negative LWC to be zero.

## [1.1.4] 2021-11-12
### Notes 
- Addition of melt scheme parameters to .json and several small bug fixes.

### New
- *melt.py, example.json* I added the 'user choices' from melt.py to the .json config file. They are: 
	- `ColeouLesaffre` (True/False; whether or not to use the ColeouLesaffre parameterization for irreducible water content); 
	- `IrrVal` (float [%], default 0.02; Irreducible water content if ColeouLesaffre is false); 
	- `RhoImp` (float, [kg m^-3], default 830; impermeable density threshold); 
	- `DownToIce` (True/False; if True all water will permeate through any ice lenses above the firn/ice transition); 
	- `Ponding` (True/False; True allows water to saturate (pond) atop impermeable layers, e.g. to form an aquifer. Otherwise water above irreducible content becomes runoff); 
	- `DirectRunoff` (float; applicable if Ponding==True; fraction of excess LWC not considered for ponding but running off directly [between 0 and 1]); 
	- `RunoffZuoOerlemans` (True/False; applicable if Ponding==True; computes lateral runoff following Zuo and Oerlemans (1996) Eqs.(21,22)); 
	- `Slope` (float, used only if RunoffZuoOerlemans==True: slope value used in Zuo and Oerlemans (1996) Eq.(22)).

### Fixed
- *melt.py* Fixed issue in the ColeouLesaffre parameterization for the case when rho = rho_ice, which put a zero in a denominator
- *firn_density_nospin.py, firn_density_spin.py* Fixed issue with reader.py, which now returns 4 values (associated with new feature that saves model forcing to the output file)

## [1.1.3] 2021-11-10
### Notes
- There are a fair number of changes and fixes in this release, and admittedly I did a poor job of documenting them over the last few months as I worked on things. Most of the work deals with the meltwater bucket scheme and the associated enthalpy scheme.

### In Progress
- I have moved away from using main.py to run the CFM. Instead, I am using either (1) a script similar to main.py but that includes forcing file generation and model configuration right in that script (i.e. you edit json parameters in that script rather than in a .json file directly); or (2) a jupyter notebook to configure the run and then start the run. Please let me (Max) know if you want those scripts or notebooks prior to them making it to Github.

### Fixed
- *diffusion.py, solver.py* The enthalpy scheme was continuing to give me issues. The solver needs to iterate to converge on a solution. In order to save computing time, the iteration loop would terminate once one iteration was within some percentage threshold of the previous iteration - but that is not full convergence. It turns out that that method could lead to a mass conservation issue. There is a new variable called ICT (itercheck threshold) that determines how close one iteration needs to be to the previous one in order to terminate the loop. By default I set this to zero (full convergence), but you can change it to something like 1e-8 if you want to speed things up a little bit. (The metric I am using for convergence is the total LWC in the firn). If there are more than 100 iterations, ICT is set to 1e-8 (if initially zero), or multiplied by 10 (if ICT was initially >0). After 200 iterations, ICT is again multiplied by 10. My testing has indicated that this works, but please let me know if you are having issues with this. 
- *melt.py* There were a few issues with the bucket scheme. In particular, there is a loop to deal with distributing water based on available pore space and cold content, which would throw an error when there was available pore space at a shallower depth then where the water existed. Also, the code that allows ponding atop impermeable layers did not remove 'excess' water (more than could be accomodated with full saturation), which led to mass conservation issues. Finally, the temperature of the uppermost node where melt occurs (the PM, or partial melt, layer) did not previously have its temperature set to the melting temperature; it now does. 

### New
- *ModelOutputs.py, writer.py* the output file now by default includes the original forcing data (including spin up period) that was fed to the CFM. This facilitates reproducibility. The columns in the 'forcing' output variable are: [decimal time, skin temperature, smb, melt, rain].

## [1.1.2] 2021-06-16
### Notes
- The main work in this release is improving the enthalpy scheme for resolving heat diffusion when there is liquid water present. I cannot say for certain how 'wrong' the previous scheme was (I don't think it was too wrong!), but Vincent Verjans identified that with the new meltwater schemes water would refreeze more quickly than expected. The newest code is actually slower because the solver needs to iterate to find a solution; previously I used a set value of iterations, but now a while loop ensures that the iterations continue until convergence.

### Fixed
- *diffusion.py, solver.py* See note above regarding the enthalpy routine.

### Changed
- *ModelOutputs.py* GridOutputs is now compatible with melt functions.
- The outputs from running the melt module have been reduced to three: **LWC** (liquid water content, in each model node, at each timestep [m w.e.]); **refreeze** (total refreezing within the firn at each time step [m w.e.]); and **runoff** (total runoff from the firn at each timestep [m w.e.]). Note that this assumes a 1m x 1m firn column, so if you want total runoff for e.g. a MERRA-2 grid cell, you will need to multiply by the area of the grid cell. 

## [1.1.1] 2021-05-19
### Notes
- This release could be buggy; it is the first release of two new liquid water schemes developed by Vincent Verjans. My limited testing indicates that they are working, but I am still working on testing. Please let me know if you encounter errors.
- I am still working on making the pre-CFM workflow smooth; i.e. taking data directly from an RCM, creating a timeseries of the climate varibles, and passing those to the CFM (scripts: siteClimate_from_RCM.py, RCMpkl_to_spin.py). If you are using these (or interested), please let me know if you have suggestions on how to make this workflow easier.

### New
- *melt.py* Previously there were two bucket schemes (one written by Vincent Verjans, and one written by Max Stevens). Vincent wrote a new bucket scheme that took the best of each of those, so now there is one bucket scheme that we think works as accurately as a bucket scheme could. It has more options than the previous schemes, e.g. set the density and/or thickness of impermeable layers; set the irreducible saturation. At present, you need to dig into the code to change them, but these will be integrated into a .json in a future release.
- *melt.py* Vincent also wrote a meltwater scheme that uses Darcy flow. The details are on the CFM's readthedocs page. FYI: It is slow compared to the bucket scheme because it uses a sub-time step. 
- *darcy_funcs.py* This is a new file that contains functions needed to make the darcy scheme work.

### Fixed
- *firn_density_nospin.py* There are more changes than this (I have been very lazy with documentation recently! Sorry!), but the main fix is adding a bit of code to deal with the instance when the forcing data are passed in a dictionary (i.e. climateTS) but spinUpdate is false. This is implemented around line 185.

## [1.1.0] 2021-03-10
### Notes
- This is the first major change to the code structure that warrants going from 1.0 to 1.1. Previously, the spin up was done by calling firn_density_spin, getting the output, and then calling firn_density_nospin. This was a bit clunky, because (1) if you were writing your own API (i.e. not using main.py), you would have to call both of those; and (2) if you were using a variable climate spinup, firn_density_spin was just a formality to create an initial condition. 
- I have now wrapped firn_density_spin into firn_density_nospin. If the model needs to spin up or needs an initial condition, it is called; otherwise it is bypassed. THE UPSHOT IS THAT YOU ONLY NEED TO CREATE AN INSTANCE OF firn_density_nospin IN YOUR API. The behavior of the CFM should be the same as before.
- You can now input climate data from a dictionary into firn_density_nospin.

### New
- *RCMpkl_to_spin.py* This new script takes climate data that is packaged in a pandas dataframe and creates a time series (including spin up) that forces the CFM. Returns a dictionary.
- *siteClimate_from_RCM.py* This script takes climate data from an RCM or analysis product for a specified lat/lon and puts it in a dataframe, which in turn will feed into RCMpkl_to_spin.py.
- *firnbatch_generic.py* This script is an alternative API to using main.py - it is useful to batch a bunch of runs when you have a common (baseline) .json file and need to change something slightly for each run.

### Changed
- See Notes above for an overview.
- *main.py* Now it only calls firn_density_nospin. 
- *main.py* Changed so that it can now take a pickled pandas dataframe (specified in .json; key: 'input_type') as opposed to .csv files; calls new script 'RCMpkl_to_spin' to generate a spin up time series. That spin up series is put in a dictionary, which is passed to firn_density_nospin
- *firn_density_nospin.py* Now calls firn_density_spin if no spin file is found, of if NewSpin=True in the .json, or if -n is included as an argument when calling main.py. Otherwise spin up is skipped and uses the spin file for the initial condition.
- *firn_density_nospin.py* The 'stpsPerYear' parameter in the .json is no longer needed if 'timesetup' is set to exact. In this case, CFM looks at the input climate files and figures out the time step size and updates 'stpsPerYear'; for 'exact' runs, it is only used to calculate the number of nodes in the grid in firn_density_spin. If 'timestep' is 'interp', the model behavior is the same as before (you need to set stpsPerYear correctly.)

### New keys in the .json file
- *"NewSpin"* true/false - this is whether or not you want to redo the spin if a spinup file exists already. Default false.
- *"input_type"* 'dataframe' or 'csv' - if the inputs to CFM come from a csv (as has been the case until now) or from a dataframe. Default csv. This only matters at the main.py level. If you write your own API to firn_density_nospin you can pass a dictionary as the climateTS argument passed to firn_density_nospin (as is done with firnbatch_generic.py).
- *"DFresample"* pandas Timedelta (string) - the resampling frequency (e.g. '5D') for the climate input data. 
- *"DFfile"* filename as string - the filename, located in "InputFileFolder", of the dataframe to load if "input_type" is 'dataframe'.

## [1.0.9]
### Added
- *firn_density_nospin.py* (and outputs): The DIP output now includes an additional column (the last), which is the DIP to a specific depth horizon (DIPhz). This is specified in the .json with key "DIPhorizon". This feature is helpful because the bottom of the model domain at each time step can change through time, which results in inconsistencies in how the (total) DIP changes through time. The default value for DIPhorizon is 80% of the initial bottom of the domain (e.g. if the model domain was initially 200m, the DIPhorizon is 160m). If the firn thickness changes and the bottom of the domain becomes less than DIPhorizon, the DIPhz value output will be NaN. In this case, the density and depth outputs can be used to calculate DIP to any depth horizon. (This is the previous model behavior.) The very first value in the DIPhz column is the horizon depth. 

### Fixed
- *ModelOutputs.py* When using the grid outputs option, previous behavior was that values would be extrapolated using the same value for all points outside the interpolation. For example, if the output grid was 0 to 150 m, but the actual model domain only went to 140 m, the same values for density, temperature, etc. would be assigned to the outputs for all depths beyond 140m. Now they are filled in using NaN.

### Changed
- *physics.py* B. Medley at NASA GSFC re-ran her calibration, and the model parameters have changed slightly.
- *firn_density_nospin.py* The CFM would include 'dr2_dt' in the outputs associated with grain size. That variable was initialized, but not updated in the main code (time-stepping loop), which led to the saved file including a bix matrix of zeros. For now this has been disabled (i.e. there is no 'dr2_dt' at all.)

## [1.0.8]
### Notes
- there is only one quick fix in this release. 
- Related to this fix, I should note that I use 365.25 to calculate seconds in a year. So, when creating forcing files from climate data, which requires calculating mass fluxes in units of m/year, you should use 365.25 days per year. (And, this then gets divided out when you use 'exact' time stepping.) If you have any suggestions on how to improve this, please let me know! The reality is that it is a bit challenging because each year does not have the same number of seconds due to leap years.

### Fixed
- *firn_density_nospin.py* Using 'exact' time setup still used the 'stpsPerYear' in the .json file to figure out the specific amount of accumulation (and melt, and rain) at each time step. This led to small errors in the calculated accumulation rate (especially if the time step size was slightly different for each time step.) It now uses the delta time (dt) value to calculate those mass fluxes.


## [1.0.7]
### Notes
- This release is a bit preemptive, but pushing it now to fix issue with timing of temperature diffusion and new surface temperature in time loop.
- A few major changes in this release: 
- The way that the CFM saves its outputs and writes to file has been updated (including a new module, ModelOuputs.py)
- The regridding scheme has been updated by Vincent Verjans to add a third subgrid.

### Added
- *ModelOutputs.py, firn_density_nospin.py* ModelOutputs.py is a new module that takes the place of a bunch of code in firn_density_nospin.py. Previously, firn_density_nospin initialized arrays for each variable so be saved/written to file (e.g self.rho_out). ModelOutputs is a class, and it sets up a dictionary to contain all of the output arrays. These arrays are updated at each time step by calling ModelOutputs. The upshot is that firn_density_nospin is cleaned up a fair bit and it is easier to add new features that you might want to write to file.
- *diffusion.py* Added thermal conductivity parameterization from Calonne et al. (2019) (which is the new default).

### Changed
- *writer.py* With the new ModelOutputs module, writer.py now loops through the variables to write with a generic line of code, rather than unique bits of code for each variable.
- *regrid.py* regrid now splits the grid into 3 subgrids of different resolutions. This new scheme is very similar to the original regrid scheme. The difference is that it uses a coarser resolution below grid2, in a grid called grid22. There is still grid1 (high resolution), then grid2 at a coarse resolution and then grid22 at a very coarse resolution. Below grid22, there is the grid23 which is again the same resolution as grid2. As before, the grid3 provides the stock of layers to be removed at each accumulation event. To use this, you need to add 3 new keys to the .json: "multnodestocombine" specifies how many nodes in grid2 to combine to make grid22; "grid2bottom" specifies the depth of the grid2 to grid22 transition. More information about the doublegrid feature is in the CFM documentation.
- *writer.py, ModelOutputs.py* There is now an option to grid your outputs onto a common grid, which can save some space on your hard drive. Add 'grid_outputs' (True/False) to the .json to enable it, and add 'grid_output_res' to the resolution to specify the resolution (in meters) that you want the grid to be at (e.g. 0.1 will make a 10-cm grid.)
- *diffusion.py* Calonne 2019 is now the default thermal conductivity.

### Fixed
- *firn_density_nospin.py* Order of operations in the time-stepping loop has been changed (reverted) - the surface temperature was updated before diffusion was called, which resulted in getting slightly different answers when using doublegrid.
- *regrid.py* Previously, the CFM did not account for changes in self.gridtrack due to removal of nodes by melt/sublimation. This is fixed for the sublime and bucketVV schemes.
- *melt.py, bucketVV* When the melt amount was such that the pm_dz became extremely small, the model crashed. pm_dz is now forced to be at least 1e-6 m thick.
- melt.py, bucketVV* A previous update got rid of a for loop in bucketVV when calculating runoff. There was a small error that computed the runoff term as a vector instead of a scalar. 


## [1.0.6]  2020-11-04
### Notes
- This is a long overdue release, and there will likely be a number of changes that are not documented here. I am pretty confident that nothing should break from 1.0.5 to 1.0.6. 

- The CFM paper to cite has finally been published; it is at https://doi.org/10.5194/gmd-13-4355-2020


### Added
- *physics.py* The NASA Goddard densification model (Medley et al. 2020) has been added.
- *diffusion.py* There is now a function that allows the user to choose which parameterization for conductivity to use. There a quite a few! The .json configuration file needs a new key, 'conductivity', to utilize this functionality. Default is 'Anderson'.
- *firn_density_nospin.py* New functionality to use firn temperature observations (e.g. from a thermistor string)
- *firn_density_nospin.py* New functionality to turn off densification (in case e.g. you only want to model temperature)
- *radiation penetration* Still in progress!
- *firn_density_nospin.py* New functionality to set output bit size. Previous and current default is float32. Add key 'output_bits' to .json file to specify. Recommend keeping at float32 except in the case that you might be feeding the results into another model run.
- *firn_density_nospin.py, writer.py* New functionality to write to the spin file during the 'main' (non-spin) run, which allows you to run a long model run (e.g. effective spin up with a reference, variable climate). Then, future runs can reuse that effective spin up. See documentation for more information. Two new keys in the .json file for this: 'spinUpdate' (boolean) and 'spinUpdateDate' (float).

### Fixed
- *physics.py* Morris and Wingham (2014) physics has a bug (the units were wrong on the deltaE parameter)

### Changed
- *firn_density_nospin.py* the final step of the model run now gets written - prior behavior was that the time steps to write were by interval, e.g. every 3rd time step. New behavior will now write the final time step no matter what.



## [1.0.5] - 2019-11-05
### Notes
- This release includes quite a number of changes, which I admittedly did a poor job documenting. The biggest change is the directory structure (sorry if this causes things to break for you, but it simplifies things in the long haul.) *CFM_main* was moved out of directory *firnmodel* and is now just a subdirectory just below CommunityFirnModel. Directory *gasmodel* is gone (all that code was in the firn air module, so it was legacy code). Directory *firnmodel* is also gone.
 
- The other notable changes are: fixing Morris (2014) physics; changing the way the timestepping is set up and handled; and isotope diffusion now is its own class.

- I have begun to start documenting the CFM using readthedocs.io, so there is now a directory under *CommunityFirnModel* called docs which holds that documentation information. Please let me know if you have suggestions.

### Fixed
- *physics.py* Barnola physics had an error where the exponent *n* for zone 2 densification was always 3; it should be 1 (and now is) when the stress is less than 0.1 MPa.
- *physics.py* Fixed Morris (2014) physics based on correspondence with L. Morris (spring 2019)

### Added
- *isotopeDiffusion.py* Isotope diffusion is now in a class. Each isotope (dD and d18O) gets its own instance. This module took the code previously in diffusion.py (written by Emma Kahle). There are several bug fixes (notably b in the tortuosity factor, Johnsen 2000 eqn. 18), and diffusion length is now included as an output.
- *merge.py* Code to merge very thin layers into a thicker layer.
- *fcts_snowpackflow.py, prefflow_snowpack.py, re_snowpack.py* Scripts for Vincent Verjans' implemention of the reynolds equation and preferential flow (not yet fully tested for compatibility in master branch.)
- *sublim.py* Script to include sublimation as a process. 

### Changed
- *firn_density_nospin.py, firn_density_spin.py* Previously the CFM took time steps of size dt that were the same size. Now dt is a vector and can vary. There is new field 'timesetup' in the .json configuration file, which can be 'exact', 'interp', or 'retmip'. 'interp' is the old method; the model takes the start and end dates from the input files and the number of timesteps per year from the .json to put together (uniform dt) vector of modeltimes. 'exact' actually takes the series of decimal years (e.g. 2015.12, 2015.24, 2015.37) from the input files and uses the spacing between those dates to get dt. 'retmip' is specific to the retmip experiment - its functionality is not tested within the main framework and it may be removed in the future.


## [1.0.4] - 2019-07-17
### Notes
- This release includes several components: I fixed implementation of Morris and Wingham (2014) physics, changed how mean temperature is calculated, added a corrected/modified DH calculation, and changed how T10m is calculated. There may be a few things (e.g. random print statements) that are hanging around. I am still working on implementing V. Verjans melt schemes.

### Changed
- *physics.py, firn_density_spin.py, firn_density_nospin.py* A number of the models use the mean annual surface temperature in their calculations (e.g. the Arthern family, Li and Zwally). On short time scales, this is generally equal to the mean of the entire temperature input. In longer simulations it changes through time. self.T_mean is now calculated at the beginning of firn_density_nospin as the mean temperature of the last 10 years. This should be quite close to the value of T10m if there is no melt. I am still working on how to implement this; it currently uses a Hamming filter with pandas 'rolling' feature. 
- *diffusion.py* T10m previously used an interpolation to find the temperature at 10m exactly; it now is just the temperature at the first node at 10m depth or greater. It should be slightly faster.
- *physics.py* I fixed how the Morris and Wingham (2014) model was coded. There was an error in the original text. I fixed it based on correspondence with L. Morris in May 2019. Email me (maxstev@uw.edu) if you want more detail; it will be published in a forthcoming paper.
- *firn_density_nospin.py* DH is calculated based on compaction at each time step, the thickness of the new snow layer, and the thickness of ice (called *iceout* in the model) that is removed at each time step. If the bottom of the model domain is less than ice density, the CFM did not account for the additional compaction that occured between the bottom of the firn and the ice (917) density horizon. I added a 2 new fields, DHcorr and DHcorrC, that is added to the 'DIP' model output, which mirror the previous fields DH and DHc (elevation change since last time step and cumulative since begining of when output begins writing).

## [1.0.3] - 2019-05-13
### Notes
- This release fixes several bugs. (I am working on implementing V. Verjan's preferential flow melt routine; there are several things that may you may notice related to that in this release, e.g. a note about 'merging' and 'Reeh01' in the .json file. This will be fixed in version 1.0.4.)

### Changed
- *physics.py* HLSigfus was changed slightly so that boolean masks are created for the different zones, which makes the code easier to read. HLdynamic (instant version) was changed so that in the case of negative accumulation, A_instant is set to zero.
- *firn_density_nospin.py* Changed dH calculation to set to be zero at the first write time step.
- *firn_density_nospin.py* Added **comp_firn** as an additional variable written with DIP. It is the total firn compaction during the previous time step, not including the ice dynamics or new accumulation (i.e. the summed thickness change of all of the layers).
- *firn_density_nospin.py* Variable **crate_out** has been changed to **comp_out** to reflect that it is just compaction, not the rate.
- *writer.py* Output name has changed from 'compaction_rate' to 'compaction'
- *firn_density_nospin.py and firn_density_spin.py* Changed the 'SeasonalTCycle' feature to exit model run if turned on and hemisphere is not selected.

### Fixed
- *constants.py* Latent heat and specific heat were incorrectly converted to kJ units - fixed, and values were changed to be consistent with Cuffey and Paterson. Note that specific heat changes with temperature; it can be set to reflect that in diffusion.py (Thanks to Jiajin Chen and Gang Hai for pointing this out.)
- *firn_density_spin.py* the seasonalTcycle feature was hard coded to be only for Greenland (it had been updated to be modular in firn_density_nospin) - now 'spin' matches 'nospin'.
- *melt.py*, in def percolation_bucket: The compaction output was incorrect because the partial melt (PM) was given a new thickness (in dzn) of zero. Fixed so that the boxes the melt entirely get a value of zero in dzn but the PM box gets its new thickness after melt
- *general* Running the example, using example.json, did not work because bdot_type was set in it to instant, and the example smb file contains a negative value. firn_density_nospin.py now prints an error and exits if this happens and example.json is set to bdot_type = 'mean'.

### Added
- *physics.py, firn_density_nospin.py, writer.py* Added code to calculate the firn viscosity predicted by each model and optionally write that to file (Brita Horlings contribution).
- *firn_density_nospin.py* Added Li2015, which uses a different value for beta than Li and Zwally 2011.
- *firn_density_spin.py* Added 'manualclimate', which must also be added to .json file to use. It allows the user to set the long-term temperature and accumulation rate manually; helpful e.g. in the case of short model runs.
- *firn_density_nospin.py* Added 'manual_iceout', which must also be added to .json file to use. It allows the user to set the long-term ice-dynamic thinning manually, e.g. the vertical ice velocity at the bottom of the firn column; helpful e.g. if it is different than the mean accumulation rate for the modeled period.
- *example.json* several new fields were added.

### Removed
- *melt.py* Max's original bucket scheme, which did not work well, was removed.

## [1.0.2] - 2019-02-08
### Fixed
- *diffusion.py* The enthalpy method has been redone using a source-based method. The code is based on work from Voller and Swaminathan, 1990. The function in diffusion.py is called enthalpyDiff(). The prior code is still there (for now) and is called enthalpyDiff_old. 

### Added
- *melt.py*: An additional bucket scheme, coded by Vincent Verjans, has been added. It is meant to be used with standard heat diffusion. Testing is underway to use it with the enthalpy diffusion method. To use this method (for now) you need to comment/uncomment in firn_density_nospin.py.

### Changed
- *diffusion.py, solver.py*: Previously, Gamma_P in diffusion.py was defined as `Gamma_P = K_firn / (c_firn)`. This, however, assumes that c_firn is constant, but it is not (varies with temperature). Now, `Gamma_P = K_firn`. A new variable, c_vol, is added - it is the heat capacity: `c_vol = rho * c_firn`. This comes into solver.py where a_P_0 is defined: previously it was `a_P_0 = tot_rho * dZ / dt`, but now it is `a_P_0 = c_vol * dZ / dt`.

- *plotter.py*: Changed to be simpler and more generic.

## [1.0.1] - 2018-11-05
### Added

- User can now initialize the model using density and temperature measurements, rather than using the firn profile predicted by the spin-up climate. Two fields need to be added to the .json file. The .csv file with the input data is column-wise. The first column is depth, the second is density, and the third is temperature.
	- *example.json*: field **"initfirnFile"**, (filename.csv)
	- *example.json*: field **"initprofile"**, (true/false)

- There is now an option to write refrozen and runoff (*writer.py*)

- *example.json*: field **"SeasonalThemi"**, (north/south), which specifies the hemisphere that the model run is simulating. This only works in concert with with SeasonalTCycle: true
- *example.json*: field **"coreless"**, (true/false): if SeasonalTCycle: true and hemisphere is south, the user can choose to implement a 'coreless winter' seasonal temperature cycle.

- *All files:* Added #/usr/bin/env python at start
- *All files:* Converted all tabs to spaces (python standard)

### Fixed

- There was an issue that during time steps with no accumulation self.dzNew was not set to zero, leading to the model over-predicting surface-elevation changes when there are small time steps (i.e. many with zero accumulation). Thanks to Brooke Medley for pointing this out.

## Contributor notes
A changelog is meant to keep a record of changes that are made to software. Any changes made should be recorded here. When you are working on the code, I recommend keeping this file open and logging changes as you make them. You should (1) keep a version of changelog.md that is exclusive to your branch and (2) keep a running list at the top of what features you are working on. When you push your code to master, update this file (i.e. the version on your branch) to reflect what it should look like on the master branch.

If the changes are pushed to the master branch, the version should be updated. 

Changes should be categorized into *Added*, *Changed*, *Deprecated*, *Removed*, or *Fixed*. Because there are numerous files in the CFM that can be changed, the affected file should be listed in italics prior to the change when appropriate. E.g.

- *constants.py* Added specific heat of water

Changes to the .json files (most frequently new fields that are added to the configuration .json files) should be updated in the file example.json, as then recorded in this document specifying the addition/change, e.g.

- *example.json* Added field **initprofile**, (true/false) which specifies whether to use an initial condition from depth/density/temperature measurements, rather than a spin up climate

### Versioning
Any push to the master branch should get a new version, which should be changed in main.py as well as in this changelog and readme.md. Small changes warrant changing the version in the 3rd digit (i.e. 1.0.X), and larger changes warrant changing the 2nd digit (i.e. 1.X.0). I am not sure what a version 2 of the CFM would look like.

I used information from [KeepAChangeLog](https://keepachangelog.com/en/1.0.0/) for this document.


