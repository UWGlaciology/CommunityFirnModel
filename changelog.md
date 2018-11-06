# CFM Change Log
All notable changes to the Community Firn Model should be documented in this file. Contributors to the CFM who are unfamiliar with changelogs should review the notes at the bottom.

## Current Version
1.0.1

## Work in progress and known issues

- *Issues* 
	- Enthalpy diffusion for meltwater percolation does not work as it should. Right now, the enthalpy diffusion works by setting the diffusivity of the nodes with liquid water very high. The correct way of doing this is to include a heat sink term. This should be resolved by the end of 2018.

- *Work in progress*
	- V. Verjans is working on a percolation module that solves Richard's equation and includes a dual-domain approach to handle preferential flow
	- Documentation for the CFM

## [1.0.1] - 2018-11-06
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

### Changed

- *firn_density_nospin.py* The order of grain growth is changed. Now, grain growth occurs independently of adding a new layer of surface grain size, which allows for grain growth when there is not accumulation. 

## Contributor notes
A changelog is meant to keep a record of changes that are made to software. Any changes made should be recorded here. When you are working on the code, I recommend keeping this file open and logging changes as you make them. You should (1) keep a version of changelog.md that is exclusive to your branch and (2) keep a running list at the top of what features you are working on. When you push your code to master, update this file (i.e. the version on your branch) to reflect what it should look like on the master branch.

If the changes are pushed to the master branch, the version should be updated. 

Changes should be categorized into *Added*, *Changed*, *Deprecated*, *Removed*, or *Fixed*. Because there are numerous files in the CFM that can be changed, the affected file should be listed in italics prior to the change when appropriate. E.g.

- *constants.py* Added specific heat of water

Changes to the .json files (most frequently new fields that are added to the configuration .json files) should be updated in the file generic.json, as then recorded in this document specifying the addition/change, e.g.

- *example.json* Added field **initprofile**, (true/false) which specifies whether to use an initial condition from depth/density/temperature measurements, rather than a spin up climate

### Versioning
Any push to the master branch should get a new version, which should be changed in main.py as well as in this changelog and readme.md. Small changes warrant changing the version in the 3rd digit, and larger changes warrant changing the 2nd digit. I am not sure what a version 2 of the CFM would look like.

I used information from [KeepAChangeLog](https://keepachangelog.com/en/1.0.0/) for this document.


