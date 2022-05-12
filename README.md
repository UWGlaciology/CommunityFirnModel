[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3585885.svg)](https://doi.org/10.5281/zenodo.3585885)

# The Community Firn Model

Welcome to the repository for the **Community Firn Model (CFM)**. 

The CFM is a comprehensive firn-model framework. It is designed to be modular, which allows the user to choose which physical processes in firn she or he would like to model. Its base function is to model firn densification and heat transfer in firn.

The CFM is meant to be useful for a variety of applications, including ice core science and estimations of ice-sheet surface-elevation change.

A few novel aspects of the CFM include:

- The user can choose which firn densification physics she/he would like to use for a model run.
- A firn-air model is optionally coupled to the firn-densification model, which allows studies of firn-air transport under non-steady-state firn conditions.
- Simulations of meltwater percolation, refreezing, and runoff. This feature currently only uses a basic bucket scheme. Work is underway to model meltwater percolation using Richard's equation and/or a dual domain scheme.
- Ability to simulate diffusion of water isotopes in firn.

Funding for this project has come from the National Science Foundation (NSF) and NASA.

## Documentation

Full documentation available here: https://communityfirnmodel.readthedocs.io/en/latest/ 

Documentation is still (always) a work in progress. Please let me know of glaring omissions!

## Installation

The CFM is coded in python. To use the software, you need to clone the repository to your computer.

## Running the CFM

The CFM is meant to be run from the command line or through an ipython session. All of the details for a model run are specifed in a .json file; the repository includes example.json. The values for each key in the .json file can be altered for a particular model run. The model is run from the command line using:

>>> python main.py example.json

More details can be found in the full documentation. 

## Dependencies

The CFM is coded in python 3. It is not backwards compatible with python 2.

We recommend the Anaconda python distribution, which includes all of the python packages needed to run the CFM. The specific packages you will need are: numpy, scipy, h5py, and pandas. If you want to plot your results, you will need matplotlib.

## Technical support and bug reports

The software does not include techinal support, but we are happy to help provide assistance if we can. Please submit questions and/or bug reports to Max Stevens: maxstev@uw.edu or maxstev@umd.edu

## Feature Requests

If there is a feature you would like to see (e.g. simulating different processes; writing different outputs), please let us us know and we can try to integrate that for you in the next release.

## Citing

If you use the CFM for your research, please cite your use. The correct citation is: https://doi.org/10.5194/gmd-13-4355-2020

If you are using the CFM, please consider sending a note (maxstev@umd.edu) to let us know that you are using it. We try to keep track of user numbers, which helps us keep the project going.

## License

MIT © C. Max Stevens

The CFM is open-source software, but we ask that you please cite your use:

Stevens, C., Waddington E.D., Conway, H., & Koutnik, M. (2018). Investigations of physical processes in polar firn through modeling and field measurements, ProQuest Dissertations and Theses.

## Contributors
C. Max Stevens
Vincent Verjans
Emma Kahle
Vasileios Gkinis
Brita Horlings
Annika Horlings
Brooke Medley
Jessica Lundin
Huong Vo
Ed Waddington
