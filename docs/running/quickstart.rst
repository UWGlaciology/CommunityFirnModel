.. _quickstart-page:

Installation and quick start
============================

The CFM is a set of Python scripts. You get the code from GitHub, make a Python environment with the packages that it needs, and run the scripts from the folder ``CFM_main``. This page gets you from nothing to a first model run and a plot of the result in about ten minutes.

Install
-------

1. **Get the code.**

   .. code-block:: bash

       git clone https://github.com/UWGlaciology/CommunityFirnModel.git
       cd CommunityFirnModel

   (Or download the repository as a zip file from GitHub. Releases are tagged; ``git checkout <tag>`` gives you a released version, for example ``4.0.0``. The ``main`` branch is the latest release, and ``dev`` is the work in progress.)

2. **Make a Python environment.** The CFM needs Python 3.11 or newer (3.11 is the first version with the ``tomllib`` module that reads ``.toml`` configuration files) and the packages in ``requirements.txt``: numpy, scipy, pandas, h5py, xarray, netCDF4, matplotlib, psutil and pytest. With conda:

   .. code-block:: bash

       conda create -n cfm python=3.12
       conda activate cfm
       pip install -r requirements.txt

   Any other way of making an environment (``venv``, ``uv``, a system Python) works as well. The supported way to use the CFM is to run the scripts from the cloned folder ``CFM_main``. (``setup.py`` exists, but ``pip install .`` would only install the Python modules, without the example data and configuration files, and it does not list ``psutil`` as a requirement.)

3. **Check the installation (optional).** From the repository root:

   .. code-block:: bash

       python -m pytest -m "unit or scheme"

   runs the fast tests in a minute or so. See :doc:`testing` for the other groups.

The first run
-------------

Go to the folder with the scripts and run the example for DYE-2 in the Greenland percolation zone:

.. code-block:: bash

    cd CFM_main
    python run_CFM_example.py "66.5,-46.25"

(Use ``"72.5,-38.75"`` for the Summit example. Keep the quote marks.) It takes about two minutes on a laptop. The script reads the MERRA-2 forcing for the site from ``CFMinput_example``, sets up the spin-up, runs the model with the forcing from 1980 to the middle of 2025 with the surface energy balance and writes its results to a new folder, ``CFMoutput_example/CFMresults_66.5_-46.25_GSFC2020``, which contains:

* ``CFMresults.hdf5``: the model output (see :doc:`outputs`),
* ``CFMforcing.hdf5``: the forcing that the model used,
* ``CFMspin.hdf5``: the spin-up state,
* a copy of the configuration file, ``CFMconfig_66.5_-46.25_GSFC2020_MERRA2_example.json``.

If you prefer a step-by-step introduction, open the notebook ``run_CFM_example_notebook.ipynb`` in Jupyter instead. It does the same thing, with explanations. How the scripts work is in :doc:`running`.

A quicker run, for checking that everything works, is the older command-line interface (it takes a few seconds because it uses a monthly time step and no surface energy balance, and writes into ``CFMoutput_example/df``, where it replaces the stored example results):

.. code-block:: bash

    python main.py example_df.json -n

The option ``-n`` makes a new spin-up (see :doc:`running`).

Look at the result
------------------

The output is in HDF5 files that ``h5py`` reads. For example, a profile of density at the end of the run and the firn air content through time:

.. code-block:: python

    import h5py
    import numpy as np
    import matplotlib.pyplot as plt

    folder = 'CFMoutput_example/CFMresults_66.5_-46.25_GSFC2020'
    with h5py.File(f'{folder}/CFMresults.hdf5', 'r') as f:
        t   = f['density'][:, 0]       # time of each row; column 0 holds the date
        rho = f['density'][-1, 1:]     # density profile at the last time step
        z   = f['depth'][-1, 1:]       # depth of the nodes (m)
        fac = f['DIP'][:, 1]           # firn air content (m)

    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    ax[0].plot(rho, z); ax[0].invert_yaxis()
    ax[0].set_xlabel('density (kg m$^{-3}$)'); ax[0].set_ylabel('depth (m)')
    ax[1].plot(t[1:], fac[1:])             # skip row 0, the initial state
    ax[1].set_xlabel('year'); ax[1].set_ylabel('firn air content (m)')
    plt.tight_layout(); plt.show()

The layout of the file, with all of the outputs and the units, is in :doc:`outputs`.

Where to go next
----------------

* **Your own site or forcing.** Put the temperature, accumulation (and, if you have it, melt, rain, sublimation and the surface energy balance forcing) into a pandas dataframe in the format of :doc:`inputs`, which also gives the sign conventions, and copy ``run_CFM_example.py`` or the minimal script in :doc:`running`.
* **How the model starts.** :doc:`spinup` describes the spin-up with a reference climate interval and the restart files.
* **The settings.** :doc:`json` is the reference for every configuration key, in ``.json`` or ``.toml`` format. The densification equation is chosen with **physRho** (:doc:`../extras/densification`).
* **Problems.** The model prints messages about configuration problems before it starts (``config_validation.py``). If something does not work, check the changelog (``changelog.md``) for the known issues, and email Max Stevens (maxstev@umd.edu).
