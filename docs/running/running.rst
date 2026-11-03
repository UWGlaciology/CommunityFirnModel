.. _running-page:

Ways to run the CFM
===================

NOTE: For a more up-to-date list of changes to the CFM, please check the changelog which is in the CFM's root folder.

The model is a Python class, ``FirnDensityNoSpin`` (in ``firn_density_nospin.py``). To run it you create an instance with a configuration file (see :doc:`json`) and the climate forcing (see :doc:`inputs`), and call its method ``time_evolve()``. The spin-up (see :doc:`spinup`) is done by the same class. The scripts below differ in how much they do for you. Run them from the folder ``CFM_main`` of the repository, because the example input files and the configuration files are found relative to it. The CFM can be run from the command line or from an integrated development environment (IDE) such as Spyder.

Recommended: ``run_CFM_example.py`` and the notebook
----------------------------------------------------

``run_CFM_example.py`` and ``run_CFM_example_notebook.ipynb`` each prepare the forcing, configure the model and run it, in one place. The notebook is a guided, step-by-step introduction for a single site; the script is the way the CFM is run in batch mode (for example by a job on a cluster), and some of its code is there to illustrate a workflow, not to be used as it is. Both work with the example forcing in ``CFMinput_example``, which are daily MERRA-2 data for two sites in Greenland, DYE-2 (66.5 N, 46.25 W) and Summit (72.5 N, 38.75 W).

To run the script for DYE-2:

.. code-block:: bash

	>>> python run_CFM_example.py "66.5,-46.25"

(or ``"72.5,-38.75"`` for Summit; the quotes are needed, and the coordinates are latitude, longitude in degrees east). On a laptop the DYE-2 run takes about a minute and a half. The script does the following:

1. **Gets the forcing.** It reads the climate dataframe for the site from ``CFMinput_example`` (the default, ``climate_source = 'dataframe'``). The option ``'zarr'`` makes the dataframe from a MERRA-2 zarr store, which is not included. The format of the dataframe is described in :doc:`inputs`.
2. **Reads a configuration.** It loads ``example_df.json`` into a dictionary (see :doc:`json`; a ``.toml`` file works the same way), and changes the keys that are specific to the run, such as the densification physics (``GSFC2020``), the model time step (``DFresample``, here ``1d``), the surface energy balance (``SEB``, switched on), melt, and the surface density.
3. **Prepares the forcing and the spin-up.** It calls ``RCMpkl_to_spin.makeSpinFiles``, with the reference climate interval 1980-1985 (``sds`` and ``sde`` in the script), to resample the dataframe, repeat the reference climate interval ahead of the record, and return the forcing dictionary ``climateTS`` (see :doc:`spinup`). It also uses what the function returns to set the number of time steps per year, the depths at which the grid is coarsened, and the depth of the model domain in the configuration.
4. **Saves the configuration, creates the model and runs it.** The configuration is written as a ``.json`` file (``CFMconfig_<lat>_<lon>_GSFC2020_MERRA2_example.json``), and ``FirnDensityNoSpin(configuration, climateTS=climateTS, NewSpin=True, SEBfluxes=SEBfluxes)`` is run with ``time_evolve()``.
5. **Moves the configuration into the results folder**, so that it is stored with the results of the run.

The results go to ``CFMoutput_example/CFMresults_<lat>_<lon>_GSFC2020/``. The folder holds ``CFMresults.hdf5`` (the output; see :doc:`outputs`), ``CFMforcing.hdf5`` (the forcing that the model used), ``CFMspin.hdf5`` (the initial profile from ``FirnDensitySpin``), and the configuration. The part of the model time before 1980 is the spin-up.

Writing your own script
-----------------------

The scripts above are a template for your own. The steps are always the same: get the forcing into a dataframe, load and edit a configuration, call ``makeSpinFiles``, and run the model. A minimal script, with the temperature, accumulation, sublimation, rain and melt from the DYE-2 dataframe for 1980-1999, a five-year reference climate interval, no surface energy balance, and the melt taken from the ``SMELT`` column, is:

.. code-block:: python

    import json
    import pandas as pd
    from config_loader import load_config
    from firn_density_nospin import FirnDensityNoSpin
    import RCMpkl_to_spin as RCM

    # 1. Forcing: a dataframe with a DatetimeIndex. Mass fluxes are in kg m-2 per time step of the dataframe.
    df = pd.read_pickle('CFMinput_example/CFM_example_66.5_-46.25.pkl')
    df = df.loc['1980':'1999']

    # 2. Configuration: start from an example file and change what you need.
    c = load_config('example_df.json')
    c['resultsFolder'] = 'my_results'
    c['DFresample'] = '1d'           # model time step
    c['MELT'] = True
    c['SEB'] = False
    c['physRho'] = 'GSFC2020'

    # 3. Forcing for the model, with the spin-up: the reference climate interval (here 1980-01-01
    #    to 1984-12-31) is repeated ahead of the record. See the spin-up page.
    climateTS, stpsPerYear, depth_S1, depth_S2, grid_bottom, SEBfluxes = RCM.makeSpinFiles(
        df, timeres=c['DFresample'], spin_date_st=1980.0, spin_date_end=1984.0,
        melt=c['MELT'], desired_depth=None, SEB=c['SEB'], rho_bottom=900)
    climateTS['forcing_data_start'] = 1980.0         # first date of the record

    # 4. Settings that depend on the time step and on the depth of the domain.
    c['stpsPerYear'] = float('%.2f' % stpsPerYear)
    c['grid1bottom'] = float('%.1f' % depth_S1)
    c['grid2bottom'] = float('%.1f' % depth_S2)
    c['HbaseSpin'] = float('%.1f' % (c['H'] - grid_bottom))
    c['nodestocombine'] = 30
    c['multnodestocombine'] = 12
    c['keep_firnthickness'] = True

    with open('my_config.json', 'w') as f:
        json.dump(c, f, indent=4)

    # 5. Run.
    firn = FirnDensityNoSpin('my_config.json', climateTS=climateTS, NewSpin=True, SEBfluxes=SEBfluxes)
    firn.time_evolve()

This takes a bit more than a minute on a laptop. Notes on the choices:

* ``nodestocombine`` and ``multnodestocombine`` are the numbers of nodes that are merged at the two depths ``grid1bottom`` and ``grid2bottom`` (see :doc:`../extras/doublegrid`). They should be adjusted to the time step: with daily steps, 30 nodes make a layer of about one month of accumulation, and 12 of those make about a year.
* ``SEBfluxes`` holds the forcing of the surface energy balance at the original time resolution of the dataframe, because the energy balance can take shorter steps than the model; it is ``None`` when ``SEB`` is false.
* The dictionary keeps the forcing in the units that the model uses (m ice equivalent per year). If the forcing is not in a dataframe, build ``climateTS`` yourself (see :doc:`inputs`), or leave it out and give the model csv files.
* If your climate record is in a different time step or covers different years, change the dates of the reference climate interval, ``forcing_data_start``, and ``TWriteStart`` in the configuration (the time at which the output starts to be written; set it to the first date of the record to leave the spin-up out of the output).

Legacy: ``main.py``
-------------------

``main.py`` is the original way to run the CFM from the command line, and is still available. It reads a configuration file and runs the model, and with csv input it needs nothing else:

.. code-block:: bash

	>>> python main.py example_df.json -n

The configuration file can be in the ``.json`` or ``.toml`` format. ``main.py`` takes the forcing as set by **input_type**: ``dataframe`` reads the pickle given by **InputFileFolder** and **DFfile** and prepares the forcing with ``makeSpinFiles``, and ``csv`` makes the model read csv files (``example_csv.json`` is a configuration for the csv example). When it is finished, ``main.py`` copies the configuration into the results folder.

Two things to know about ``main.py``:

* **Spin-up.** With dataframe input it uses the default reference climate interval of ``makeSpinFiles``, **1980-1995**, and the depth of the domain (``H`` - ``HbaseSpin``) to set the number of repeats; they are not in the configuration. To use another interval, or to set the other options of ``makeSpinFiles``, use your own script (above).
* **The option -n.** If the results folder (**resultsFolder**) already contains a spin file (**spinFileName**, usually ``CFMspin.hdf5``), the model does not make a new initial profile and uses the existing one. Add ``-n`` to the command to make a new one (the key **NewSpin** in the configuration file does the same). With the reference-climate-interval spin-up, the repeated climate is run in every case, as it is part of the forcing, and the initial profile is cheap, so ``-n`` is the better choice.

The example runs
----------------

The folder ``CFMoutput_example`` of the repository holds the results of two runs made with ``main.py``: ``python main.py example_csv.json -n`` (csv forcing, results in the subfolder ``csv``) and ``python main.py example_df.json -n`` (dataframe forcing, results in ``df``). They can be used to check that the CFM gives the same results on your computer. Both run in a few seconds, since they use monthly time steps. Note that running them overwrites those results, unless you change **resultsFolder**.

The test suite (see :doc:`testing`) checks the model in a more systematic way.
