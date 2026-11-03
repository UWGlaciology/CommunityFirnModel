.. _inputs-page:

Inputs for the CFM
==================

The CFM is forced by surface-temperature and accumulation-rate boundary conditions. Additionally, the user can specify the surface melt, rain, sublimation, surface density, water isotopes, horizontal strain rate and (for the surface energy balance) radiation and turbulent fluxes. The forcing can be given to the model in three ways:

1. **A pandas dataframe** (recommended). All of the forcing is in the columns of one dataframe with a time index, usually stored in a pickle. This is how the example script and the example notebook work, and it is the format that the reference-climate-interval spin-up needs (see :doc:`spinup`).
2. **A dictionary of arrays** (``climateTS``) that you pass directly to the model. It is what ``RCMpkl_to_spin.makeSpinFiles`` makes from a dataframe, so you only make it yourself if you want to bypass the dataframe. The model is then not reading any forcing from file.
3. **csv files**, one for each variable. This is the original input format of the CFM. It is still supported and is the usual choice for long, coarse climate histories such as those from ice cores.

The configuration key **input_type** (``csv`` or ``dataframe``) selects between 1 and 3 when the model is started with ``main.py``. When you write your own script (see :doc:`running`), the forcing is passed as the ``climateTS`` argument of ``FirnDensityNoSpin``: if it is given, the model uses it, and if it is ``None`` the model reads csv files.

Dataframe input
---------------

The dataframe has a ``DatetimeIndex``. The time step of the dataframe can be anything (hourly and 3-hourly data are common), and it does not have to match the time step of the model. ``RCMpkl_to_spin.makeSpinFiles`` resamples it to the model time step that you give as **DFresample** (a pandas offset string, for example ``1d`` for daily or ``1ME`` for monthly). Mass fluxes are summed, temperatures are averaged, and (for the surface energy balance) energy fluxes and the other variables are averaged. The model time step is therefore a choice that you make when you prepare the forcing.

The columns, using the names the model expects, are:

.. list-table::
   :header-rows: 1
   :widths: 14 50 36

   * - Column
     - Meaning
     - Units in the dataframe
   * - ``TSKIN``
     - Surface (skin) temperature. Required.
     - K (or C; the model assumes C if the mean is below 0)
   * - ``BDOT``
     - Accumulation: the snowfall. Required.
     - kg m\ :sup:`-2` per dataframe time step
   * - ``SUBLIM``
     - Sublimation (negative) and deposition (positive). Used if **SUBLIM** is true.
     - kg m\ :sup:`-2` per dataframe time step
   * - ``RAIN``
     - Rain. Used if **MELT** and **RAIN** are true.
     - kg m\ :sup:`-2` per dataframe time step
   * - ``SMELT``
     - Surface melt, from a degree-day or other model. Used if **MELT** is true and **SEB** is false. (With **SEB** true, melt is calculated by the model, and this column is dropped.)
     - kg m\ :sup:`-2` per dataframe time step
   * - ``SRHO``
     - Surface density. Used if **variable_srho** is true and **srho_type** is ``userinput``.
     - kg m\ :sup:`-3`
   * - ``eps_divergence``, ``eps_1``, ``eps_2``, ``eps_xx``, ``eps_yy``, ``eps_xy``
     - Horizontal strain rates (optional; see **horizontal_divergence**, **strain_softening**).
     - a\ :sup:`-1`
   * - ``T2m``, ``ALBEDO``, ``SW_d``, ``LW_d``, ``QH``, ``QL``
     - Surface-energy-balance forcing (needed if **SEB** is true): 2 m air temperature, albedo, downwelling short- and long-wave radiation, sensible and latent heat flux.
     - K; unitless; W m\ :sup:`-2`
   * - ``WS``, ``QV2m``, ``PS``, ``LW_u``, ``z0m``, ``z0h``, ``z0q``
     - More SEB forcing: wind speed, specific humidity and surface pressure (needed for the ``bulk_richardson`` and ``monin_obukhov`` methods of **SEB_turb_method**; see that key for the details), and upwelling long-wave radiation and roughness lengths (both optional).
     - m s\ :sup:`-1`; kg kg\ :sup:`-1`; Pa; W m\ :sup:`-2`; m

A few details that are easy to get wrong:

* The **mass fluxes are in mass per time step of the dataframe**: if the dataframe is hourly, the accumulation is in kg m\ :sup:`-2` per hour (the mass of snow that fell in that hour), not a rate. ``makeSpinFiles`` converts them to the unit that the model uses, **m ice equivalent per year**, after it has resampled them.
* **Sign of sublimation:** the model expects sublimation to be negative (a mass loss from the surface) and deposition to be positive. Some reanalyses use the opposite sign: MERRA-2 evaporation (``EVAP``) is positive for sublimation, so it must be multiplied by -1 when the dataframe is made (``run_CFM_example.py`` does this, and the example dataframes already have the CFM sign). ``makeSpinFiles`` does not change the sign of a column named ``SUBLIM``. If the sign is wrong the model treats sublimation as deposition and the reverse, which changes the accumulation by the full sublimation flux (see the note on CFM 4.1.0 below).
* **Sign of the turbulent fluxes:** ``QH`` and ``QL`` are positive when they heat the surface (**into the snow**). MERRA-2 fluxes are positive upward, so they must be multiplied by -1 when the dataframe is made, and the example scripts do this. See "Sign conventions of the forcing" below; the model also warns if ``QH`` looks like it has the wrong sign.
* If **SEB** is true, **every column of the dataframe must be one of the names above**. A column with another name stops ``makeSpinFiles`` with a ``KeyError``. (Without the surface energy balance, columns that are not used are dropped.)
* Some older MERRA-2 column names are converted by ``makeSpinFiles`` when the surface energy balance is not used: ``TS`` is renamed ``TSKIN``, ``SUBLIM`` is made from ``EVAP`` **with the sign changed** (``SUBLIM = -EVAP``, because ``EVAP`` is the raw MERRA-2 variable, positive upward), ``BDOT`` is taken from ``PRECSNO`` and ``RAIN`` is ``PRECTOT`` - ``PRECSNO``. The pickle ``CFM_example_input_df.pkl``, which ``example_df.json`` and ``main.py`` use, has these names.
* The dataframe must cover the period you want to model, and the spin-up reference climate interval must be inside it. Its first date should be the start of a year (see :doc:`spinup`).

The example dataframes in the folder ``CFMinput_example`` (``CFM_example_66.5_-46.25.pkl`` for DYE-2 and ``CFM_example_72.5_-38.75.pkl`` for Summit, both in Greenland) are daily MERRA-2 data from 1980 (with ``QH``, ``QL`` and ``SUBLIM`` in the CFM sign convention) and have the columns ``T2m``, ``TSKIN``, ``SUBLIM``, ``QH``, ``QL``, ``SW_d``, ``LW_d``, ``RAIN``, ``BDOT``, ``ALBEDO`` and ``SMELT``. They are a template for your own dataframe. The comments at the top of ``run_CFM_example.py`` describe how they were made. The files ``CFM_example_<lat>_<lon>.csv`` in the same folder hold the same dataframes as text.

The configuration keys that go with it are **input_type**, **InputFileFolder**, **DFfile** (the name of the pickle, when using ``main.py``) and **DFresample**.

.. _sign-conventions:

Sign conventions of the forcing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Forcing products do not agree on the sign of the surface fluxes, and the CFM does not know which product your
dataframe comes from, so **the signs have to be right when the dataframe is made**. The CFM convention is that the
energy fluxes are positive **into the snow** (towards the surface) and that mass fluxes out of the snow are negative:

.. list-table::
   :header-rows: 1
   :widths: 14 22 32 32

   * - CFM column
     - CFM convention
     - MERRA-2
     - ERA5
   * - ``QH``, ``QL`` (sensible and latent heat flux)
     - positive **into the snow**; negative when the surface loses heat to the air
     - ``HFLUX``, ``EFLUX``: positive **upward**. **Multiply by -1.**
     - ``sshf``, ``slhf``: positive **downward** (towards the surface). Use as they are.
   * - ``SUBLIM`` (sublimation, deposition)
     - **negative = sublimation** (mass loss), positive = deposition
     - ``EVAP``: positive = evaporation or sublimation (upward). **Multiply by -1.**
     - ``e``: negative = evaporation (ECMWF defines processes from the surface to the atmosphere as negative). Use as it is.
   * - ``SW_d``, ``LW_d``
     - downwelling radiation, positive
     - ``SWGDN``, ``LWGAB``: positive downward. Use as they are.
     - ``ssrd``, ``strd``: positive downward. Use as they are.

Other products (regional climate models, station data) have their own conventions; read the documentation of
yours. In every case the question is the same: is the flux positive when it carries energy (or mass) towards the
surface?

.. VERIFY (ERA5): ECMWF's documentation page that I could fetch states the sign of evaporation (negative; Guideline 15 of the
   ERA5 data documentation) but not the heat fluxes. The statement that sshf and slhf are positive downward is the ECMWF
   convention for fluxes and agrees with Max's statement that ERA5 and MERRA-2 have opposite signs. Please check the exact
   wording in the ECMWF parameter database before publishing. I did not check the radiation parameters of ERA5 in detail
   (ssrd and strd are "downwards" by name).

.. note::

   **Change in CFM 4.1.0.** Before 4.1.0 the example forcing (the pickles and csv files in ``CFMinput_example``, including the
   legacy ``CFM_example_input_df.pkl`` and ``example_SUBLIM.csv``) held the raw MERRA-2 sign of ``SUBLIM`` (positive =
   sublimation), the opposite of what the model expects, and the sign was only corrected inside ``run_CFM_example.py`` and the notebook, after
   ``makeSpinFiles`` had set up the spin-up. The consequences were that (1) ``main.py`` runs with the example inputs treated sublimation as
   deposition and the reverse, and (2) the spin-up of the example script and notebook used the wrong accumulation to choose the number
   of reference-climate-interval repetitions and the depth of the grid (accumulation +8% at DYE-2 and -6% at Summit). For the shipped
   ``main.py`` examples the firn air content changes by 1.5 to 3%, the depth of the 830 kg m\ :sup:`-3` level by 2% and the density at
   10 m depth by 1.6%. The inputs, the scripts and the example results in ``CFMoutput_example`` have been corrected. If you made your own
   forcing by copying the old example, or ran ``main.py`` or ``firnbatch_generic.py`` with MERRA-2 data and did not multiply ``EVAP`` by -1,
   check the sign of your ``SUBLIM`` and rerun.

**A check for the sublimation sign.** The latent heat flux and the sublimation are the same process: sublimation (``SUBLIM`` < 0) removes
energy from the surface (``QL`` < 0), so in the CFM convention the two are positively correlated. If a dataframe has both columns,
``makeSpinFiles`` calculates their correlation and prints a warning if it is clearly negative (r < -0.5), which means that one of the
two columns (usually ``QL``, ``QH`` or ``SUBLIM``, whichever was multiplied by -1 without the others) has the wrong sign. The check
cannot be made if there is no ``QL`` column (for example without the surface energy balance), and it does not catch the case where both
have the wrong sign.

**A check for the sensible heat flux.** When the surface energy balance uses prescribed fluxes
(**SEB_turb_method** ``prescribed``), the model checks the sign of ``QH`` and prints a warning if it looks positive
upward. The check uses the physics: with the CFM convention, ``QH`` is positively correlated with the difference between the
air and the skin temperature (``T2m - TSKIN``), and over an ice sheet it is positive on average. If the correlation and the
mean are both clearly negative, the fluxes were probably not multiplied by -1. The check needs at least a year of
data, because in winter the correlation can have the opposite sign even for correct data; it was tested with the MERRA-2
example forcing for DYE-2 and Summit (it is not known how it behaves for other products or for Antarctica), it checks ``QH`` only, and it only
warns: a flux with the wrong sign that passes the check, or a ``QL`` with a different sign than ``QH``, is not caught.

The ``climateTS`` dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``makeSpinFiles`` returns the forcing as a dictionary, ``climateTS``, with the keys ``time`` (decimal date), ``TSKIN``, ``BDOT``, and, when they are in the dataframe, ``SUBLIM``, ``RAIN``, ``SMELT``, ``SRHO`` and the strain rates. The arrays are on the model time step. The mass fluxes are **rates in m ice equivalent per year**. The series already includes the repeated reference climate interval ahead of the record, and the key ``forcing_data_start`` is the first date of the record. The model reads this dictionary from the ``climateTS`` argument; this is how the example script and notebook run it. If you make the dictionary yourself, for example from a source other than a dataframe, use the same keys and units; ``time`` has to be in increasing order, and the spin-up has to be part of it if you want one.

csv input
---------

For csv input the forcing of each variable is in its own file, whose name is set in the configuration file (see the ``InputFileName*`` keys in :ref:`the configuration file <json-page>`), in the folder **InputFileFolder**. The files are csv-formatted. The first row is time (decimal date, i.e. 2015.3487) and the second row is the corresponding temperature, accumulation rate or other boundary-condition value at that time. Time must be going forward, i.e. the first column is a date some time ago and the last column is the most recent. (If the model is being forced with ice-core data, the user must be careful to ensure this is the case as ice-core data are often presented as years before present.) The times in the various input files do not need to be the same; they are interpolated onto a common axis (see **int_type**).

The units for temperature can be K or C. The CFM uses K, but it will change the temperature to K if you use C. The units for accumulation rate/surface mass balance (and for melt, rain and sublimation) are **m ice equivalent per year**, which is different from the dataframe input, where the fluxes are masses per time step. Sublimation is negative (a loss from the surface).

With csv input there is no reference climate interval: the model has no information on the climate before the start of the files, and the spin-up is the steady-state spin-up (see :doc:`spinup`). The surface energy balance needs the dataframe or ``climateTS`` input.
