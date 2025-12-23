Running the CFM
===============

NOTE: For a more up-to-date list of changes to the CFM, please check the changelog which is in the CFM's root folder.

The CFM is mostly easily run from the command line, though it can alternatively be run from an integrated development environment (IDE) such as Spyder (included with the Anaconda Python distribution). To run the model, the user must set up a .json-formatted configuration file that specifies the settings for a particular model run. Here, we generically use example.json as the name of that configuration file.

The model is run by calling main.py and specifying the name of the configuration file:

.. code-block:: bash

	>>> python main.py config.json -n

If the results folder (specified in config.json) already exists and contains the results of a spin-up run (the spin-up file name is typically CFMspin.hdf5), the model will not run the spin-up routine again. If the user wishes to include the spin-up run, he/she should add the –n at the end of the above command to force the spin-up to run; if he/she wished to omit the spin up run, the –n should be omitted.

Running the example
===================

The CFM includes an example .json file (cleverly called example.json). To test that the CFM is running on your system properly, you need to make sure that you have a directory called CFMinput with the appropriate file inside of it (``MERRA2_CLIM_df_72.5_-38.75.pkl``).

Then run:

.. code-block:: bash

	>>> python main.py example.json

It takes ~60 seconds to run on my laptop. You can compare the outputs to the resuls provided in the 'example_results_compare' directory.
