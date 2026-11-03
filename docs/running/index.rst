Running the CFM
===============

New to the CFM? Start with :doc:`quickstart`. A CFM run needs two things: a configuration file, which sets the physics and the model settings (see :doc:`json`), and the climate forcing, which is usually a dataframe (see :doc:`inputs`). Before the period of interest the model is spun up (see :doc:`spinup`). The model is a Python class that you can start in three ways, described on the next page:

1. ``run_CFM_example.py`` or the notebook ``run_CFM_example_notebook.ipynb`` (recommended). They read the forcing, set up the spin-up, configure the model and run it, all in one place.
2. Your own script, which does the same with your forcing and your settings.
3. ``main.py`` with a configuration file. This is the original way to run the CFM and is being deprecated.

.. toctree::
    :maxdepth: 2

    quickstart.rst
    running.rst
    spinup.rst
    inputs.rst
    json.rst
    outputs.rst
    testing.rst
