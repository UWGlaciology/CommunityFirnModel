[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3585884.svg)](https://doi.org/10.5281/zenodo.3585884)

# The Community Firn Model

Welcome to the repository for the **Community Firn Model (CFM)**. 

The CFM is a comprehensive firn-model framework. It is designed to be modular, which allows the user to choose which physical processes in firn she or he would like to model. Its base function is to model firn densification and heat transfer in firn.

The CFM is meant to be useful for a variety of applications, including ice core science and estimations of ice-sheet surface-elevation change.

A few novel aspects of the CFM include:

- The user can choose which firn densification physics she/he would like to use for a model run.
- A firn-air model is optionally coupled to the firn-densification model, which allows studies of firn-air transport under non-steady-state firn conditions.
- Simulations of meltwater percolation, refreezing, and runoff. This feature currently only uses a basic bucket scheme. Work is underway to model meltwater percolation using Richard's equation and/or a dual domain scheme.
- Ability to simulate diffusion of water isotopes in firn.

Presently, this work is funded by NASA grant 80NSSC25K7216. Previous funding for this project has come from the National Science Foundation (NSF) and NASA.

## Documentation

Full documentation available here: https://communityfirnmodel.readthedocs.io/en/latest/ 

Documentation is still (always) a work in progress, and in the past few years I have been especially negligent of it. Please let me know of glaring omissions, and send me an email (maxstev at umd.edu) if you can't figure something out.

## Installation

The CFM is coded in python. To use the software, you need to clone the repository to your computer.

## Running the CFM

The CFM can be run from the command line using the main.py script. It is also relatively easy to create a separate script or jupyter notebook to configure and run the CFM; this can make it easier to do a large number of runs with similar parameters. 

All of the details for a model run are specifed in a .json file; the repository includes example.json. The values for each key in the .json file can be altered for a particular model run. For basic use, the model can be run from the command line using:

>>> python main.py example_df.json -n

or 

>>> python main.py example_csv.json -n

The df version forces the model using climate data stored in a pandas dataframe, while the csv version uses climate data from csv files (old/original behavior, which will continue to be supported, but is not recommended).

Starting with version 3.0.0, the repository also includes a script (run_CFM_example.py) and a jupyter notebook (run_CFM_example_notebook.ipynb) that allow the user to configure and run the CFM from the same script (or notebook). They work for example runs or can be edited for a specific use case. Instructions on running those are found in the files themselves.

More details can be found in the full documentation.

## Example Runs

Starting with version 3.0.0, there are two example .json files as well as a new .py script and notebook that can be used to run CFM examples. (See previous section, 'Running the CFM'.) Further instructions for using run_CFM_example.py are in the script itself. 

The CFMinput_example directory contains example forcing files to run the CFM. It includes .csv files and .pkl files, each of which contain a pandas dataframe with climate data. The source for these forcings is MERRA2 data for Summit, Greenland (72.5 N, -38.75 W) or DYE-2 Greenland (66.5 N, -46.25 W). 

The example_*.csv files were created using the data in the .pkl file; effectively I am running RCMpkl_to_spin.py, which returns a dictionary full of arrays of climate variables and saving those arrays as csv files. I do this in a jupyter notebook called CFM_create_examples.ipynb. I am not including that on the repository but would be happy to share.

I also create an artificial surface isotope record using the skin temperature. I use equations from Jouzel and Merlivat, 1984 to do so. These forcings only go into the .csv files.

The CFMoutput_example directory contains results generated when running the CFM using these forcings with the settings prescribed in example.json. The only thing I change in example.json between the two example runs is "input_type": "csv" or "dataframe"; and "resultsFolder": "CFMoutput_example/csv" or "CFMoutput_example/df". 

The CFM automatically moves the .json file for a particular run to the directory containing the results, so the exact .json files for the example runs are included in the CFMoutput_example/csv and CFMoutput_example/df subdirectories. You could move those .json files to CFM_main and then e.g. run CFM using >>>python main.py example_df.json. (Note that this will overwrite the results in those example results directories.)

## Dependencies

The CFM is coded in python 3. It is not backwards compatible with python 2.

We recommend the miniconda python distribution. Specific packages you will need are: numpy, scipy, h5py, and pandas, and xarray. If you want to plot your results, you will need matplotlib.

## Technical support and bug reports

I am happy to help provide technical assistance if I can. I appreciate feedback regarding what makes the CFM difficult to use and/or suggestions on how to improve ease of use. Additionally, I very much appreciate bug reports so I can fix them.

Please submit questions and/or bug reports to Max Stevens: maxstev@uw.edu or maxstev@umd.edu

## Feature Requests

If there is a feature you would like to see (e.g. simulating different processes; writing different outputs), please let us us know and we can try to integrate that for you in the next release.

## Contributing

We welcome contributions to the CFM. If you are interested in adding code, it is best to get in touch (maxstev@umd.edu) to discuss the contribution and then to submit a pull request.

## Citing

If you use the CFM for your research, please cite your use. The correct citation is: 

Stevens, C. M., Verjans, V., Lundin, J. M. D., Kahle, E. C., Horlings, A. N., Horlings, B. I., and Waddington, E. D.: The Community Firn Model (CFM) v1.0, Geosci. Model Dev., 13, 4355–4377, https://doi.org/10.5194/gmd-13-4355-2020, 2020.

If I provide substantial assistance (e.g., building custom scripts, advising on interpreting model results, etc.), I ask to be considered for co-authorship on any manuscripts coming from the work. Note this is not a requirement, and I leave it to the CFM user to decide whether any contribution I make is worthy of co-authorship.

If you are using the CFM, please consider sending a note (maxstev@umd.edu) to let us know that you are using it. We try to keep track of user numbers, which helps us keep the project going.

## License

MIT © C. Max Stevens

The CFM is open-source software, but we ask that you please cite your use.

## Contributors
C. Max Stevens
Vincent Verjans
Emma Kahle
Vasileios Gkinis
Ilyse (Brita) Horlings
Annika Horlings
Brooke Medley
Jessica Lundin
Huong Vo
Ed Waddington
Falk Orachewski
Ben Smith
