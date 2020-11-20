# CFM Change Log
All notable changes to the Community Firn Model should be documented in this file. Contributors to the CFM who are unfamiliar with changelogs should review the notes at the end of this document.

TL;DR: Write down the changes that you made to the the model in this document and update the version number here and in main.py, then update master on github.

## Current Version
1.0.5

## Full Documentation

https://communityfirnmodel.readthedocs.io/en/latest/

## Work in progress and known issues

- *Issues* 
	- If data is not written at each time step, dH that is saved/written to file does not represent the change since the last write. 

- *Work in progress*
	- If data is not written at each time step, dH that is saved/written to file does not represent the change since the last write. 
	- Testing percolation modules from Vincent Verjans that solves Richard's equation and includes a dual-domain approach to handle preferential flow (these modules are included but may not work properly yet)
	- Melt will likely be changed to have its own class
	- Documentation for the CFM
	- Goujon physics work, but could possibly be implemented more elegantly (it would be nice to avoid globals)
	- Not exactly in progress, but at some point adding a log file that gets saved in the results folder would be a good idea.

## [1.0.6]  2020-11-04
### Notes
- This is a long overdue release, and there will likely be a number of changes that are not documented here. I am pretty confident that nothing should break from 1.0.5 to 1.0.6. 

- The CFM paper to cite has finally be published; it is at https://doi.org/10.5194/gmd-13-4355-2020


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

Changes to the .json files (most frequently new fields that are added to the configuration .json files) should be updated in the file generic.json, as then recorded in this document specifying the addition/change, e.g.

- *example.json* Added field **initprofile**, (true/false) which specifies whether to use an initial condition from depth/density/temperature measurements, rather than a spin up climate

### Versioning
Any push to the master branch should get a new version, which should be changed in main.py as well as in this changelog and readme.md. Small changes warrant changing the version in the 3rd digit (i.e. 1.0.X), and larger changes warrant changing the 2nd digit (i.e. 1.X.0). I am not sure what a version 2 of the CFM would look like.

I used information from [KeepAChangeLog](https://keepachangelog.com/en/1.0.0/) for this document.


