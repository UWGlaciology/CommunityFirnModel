.. _testing-page:

Testing the CFM
===============

The CFM has a test suite in the folder ``tests/`` of the repository. It uses `pytest <https://pytest.org>`_.

Run the tests
-------------

Run the tests from the repository root, not from ``CFM_main/``:

.. code-block:: bash

	pip install -r requirements.txt
	python -m pytest

The file ``pytest.ini`` in the repository root sets the test folder and the markers.

The CFM modules in ``CFM_main/`` are flat scripts, not an installable package. The tests put ``CFM_main/`` on ``sys.path`` and run each test in that folder. You do not have to do this yourself.

Test groups
-----------

Each test has one marker. To select a group, use ``-m``:

.. code-block:: bash

	python -m pytest -m "unit or scheme"    # the fast tests
	python -m pytest -m conservation        # the short model runs

.. list-table::
   :header-rows: 1
   :widths: 20 60 20

   * - Marker
     - What it tests
     - Run time
   * - ``unit``
     - Single functions and equations
     - Seconds
   * - ``scheme``
     - One meltwater scheme on a synthetic column
     - Seconds
   * - ``conservation``
     - Mass and energy budgets over a short model run
     - Minutes

What the tests check
--------------------

A test file can contain tests of more than one group. The sections below list the tests by group.

Unit tests (``unit``)
~~~~~~~~~~~~~~~~~~~~~

A unit test calls one function with a small, specified input. Where a published equation or an analytic solution exists, the test compares the result with it.

In ``test_hl_analytic.py``, the test file calculates the reference values from the published equations, with literal constants. This test file does not call the CFM function under test to get the reference, and it does not import ``constants.py``. Thus an error in the code cannot also be in the reference. A change to the code that changes a result makes the test fail, so a change to the physics must also change the test. This makes an accidental change less likely.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Contents
   * - ``test_hl_analytic.py``
     - ``hl_analytic`` (the spin-up profile) against the Herron and Langway (1980) density and age equations, in zone 1 and zone 2. A node on the boundary between the two zones stays in the profile.
   * - ``test_enthalpy_solver.py``
     - The ``enthalpy`` solver with liquid water in a layer below 0 °C. The solver refreezes this water with its latent heat.
   * - ``test_melt_diagnostics.py``
     - The record of a failed meltwater check, with and without **melt_strict** (see :doc:`../extras/melt_checks`).
   * - ``test_fixtures.py``
     - The ``synthetic_column`` fixture gives the state arrays that the meltwater schemes read.

Scheme tests (``scheme``)
~~~~~~~~~~~~~~~~~~~~~~~~~

A scheme test calls one meltwater scheme once, on a synthetic column (see `Fixtures for new tests`_). It does not run the model.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Contents
   * - ``test_melt_schemes.py``
     - ``bucket`` and ``darcyscheme``: the melted thickness ``dh_melt`` agrees with the melt mass divided by the density.
   * - ``test_fixtures.py``
     - ``bucket`` closes its water budget on a synthetic column and passes all its checks in strict mode. ``darcyscheme`` runs on a synthetic column and passes its checks.

Conservation tests (``conservation``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A conservation test runs the model (spin-up and ``time_evolve()``) for two years of the example forcing (see `Fixtures for new tests`_). It checks that the run completes, and some tests also check the water budget.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Contents
   * - ``test_short_runs.py``
     - The model completes a run with each of these configurations: **MELT** true with **SUBLIM** true and with **SUBLIM** false; **darcy**, **prefsnowpack** and **resingledomain**, each with output from 1980 and from the first step.
   * - ``test_water_budget.py``
     - The water budget for each of the four **liquid** schemes, with the DYE-2 forcing. Melt plus rain over the run agrees with refreezing plus runoff plus the liquid water in the column at the end.
   * - ``test_fixtures.py``
     - A short run with **melt_strict** true completes.

Fixtures for new tests
----------------------

The file ``tests/conftest.py`` gives two fixtures. The helpers are in ``tests/cfm_testing.py``.

``synthetic_column``
  Makes a small firn column with the state arrays that the meltwater schemes read. The density, temperature and liquid water content can be a scalar, an array, or a function of depth. The configuration starts from ``example_df.json``. Call a scheme on the column directly:

  .. code-block:: python

	import melt

	def test_bucket_refreezes(synthetic_column):
	    column = synthetic_column(n=60, rho_profile=lambda z: 350.0 + 40.0 * z,
	                              T_profile=265.0, melt=0.02,
	                              config={'melt_strict': True})
	    melt.bucket(column, 0)
	    assert column.melt_diagnostics['refrozen'] > 0.0

``short_config``
  Makes a short model run from ``example_df.json``. It cuts the forcing to a few years and writes the output to a temporary folder. Keyword arguments replace .json keys. Call ``.run()`` to run the spin-up and ``time_evolve()``:

  .. code-block:: python

	def test_short_run_completes_in_strict_mode(short_config):
	    run = short_config(years=2, melt_strict=True)
	    run.run()
	    assert run.results_path.exists()

Continuous integration
----------------------

The GitHub Actions workflow ``.github/workflows/tests.yml`` runs the test suite on each push to ``main`` and on each pull request. It has two jobs, which run in parallel:

- **Fast tests**: ``python -m pytest -m "unit or scheme"``.
- **Full tests**: ``python -m pytest -m conservation``.

To run the workflow manually, use **Run workflow** in the **Actions** tab of the repository on GitHub.

**References**

*Herron, M. M., and Langway, C. C. (1980), Firn densification: an empirical model, J. Glaciol., 25(93), 373–385, doi:10.3189/S0022143000015239.*
