.. _spinup-page:

Model spin-up
=============

A firn model needs an initial condition: the density, temperature and age of the firn column at the start of the run. If the column starts in a state that is not consistent with the climate that forces the model, the run begins with a transient that is an artifact of the initial condition and not a result. A spin-up is the part of the run that removes it: the model is run, usually for much longer than the time it takes to refresh the firn column (the time for a parcel of snow to reach the bottom of the model domain), so that the column is in equilibrium with the climate before the period of interest begins.

The CFM has two ways of doing this. They were developed for different applications, and which one is used depends on how the forcing is prepared, not on a single setting.

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * -
     - **Steady-state spin-up** (original)
     - **Reference-climate-interval (RCI) spin-up**
   * - Typical use
     - Ice-core studies; runs of 10s of thousands of years with a long, coarse climate history.
     - Recent, high-resolution runs (e.g. daily forcing from a regional climate model or reanalysis, for altimetry or mass-balance work), especially with melt.
   * - Forcing
     - csv files (or a dataframe) that start at the beginning of the climate history.
     - A dataframe (or the dictionary made from it by ``RCMpkl_to_spin.makeSpinFiles``).
   * - Who does the spin-up
     - ``firn_density_spin.py`` (class ``FirnDensitySpin``). It runs the firn column to a steady state at a constant climate; the result is the initial condition for the transient run in ``firn_density_nospin.py``.
     - ``firn_density_nospin.py`` (class ``FirnDensityNoSpin``). A representative stretch of the climate is repeated many times ahead of the real record, and the whole series is run as one continuous transient run. ``FirnDensitySpin`` only supplies a rough initial profile.
   * - Assumption
     - The firn is in steady state with the climate at the start of the record (a constant temperature and accumulation rate).
     - The firn is in dynamic equilibrium with the climate of the RCI: it has no net change over the RCI, but it does respond to the seasonal cycle, the weather and (if melt is on) the melt and refreezing of the RCI.

The RCI spin-up is the one that most current CFM work uses, and it is what ``run_CFM_example.py`` and the example notebook do (see :doc:`running`). The steady-state spin-up is described at the end of this page.

Spin-up with a reference climate interval
-----------------------------------------

The method follows the one used for the GSFC Firn Densification Model; see Section 2.1.2 of Medley et al. (2022, *The Cryosphere*, 16, 3971-4011, https://doi.org/10.5194/tc-16-3971-2022). A dry-snow steady-state spin-up builds a column that is in equilibrium with dry conditions only. If melt is then imposed, the meltwater processes produce large changes in the firn thickness and air content that are artifacts of the spin-up. Repeating a reference climate interval that includes the melt (and rain, sublimation, and so on) builds a column that is in equilibrium with the full climate of that interval, wet or dry.

The steps are done in ``RCMpkl_to_spin.makeSpinFiles`` (which prepares the forcing) and in ``FirnDensityNoSpin`` (which runs it):

1. **You choose the RCI.** It is a stretch of the climate record, given by the decimal dates ``spin_date_st`` (its start) and ``spin_date_end`` (the start of its last year). Both are inclusive of whole years: ``spin_date_st=1980.0, spin_date_end=1995.0`` is 1 January 1980 through 31 December 1995, which is 16 years. The RCI should be a period that you consider representative of the long-term climate; see "Choosing the RCI" below.

2. **The length of the spin-up is estimated.** The Herron and Langway (1980) analytic solution is evaluated with the mean temperature and the mean accumulation rate (in m ice equivalent per year) of the climate record, and it gives the depth at which the density reaches ``rho_bottom`` (default 916 kg m\ :sup:`-3`; ``run_CFM_example.py`` uses 900). That depth is divided by the mean accumulation rate to give the time that a parcel of snow takes to reach the bottom of the column. This is an overestimate of the age of the firn at that depth, which is intended: it makes sure that the whole column is refreshed. If you give the model domain depth as ``desired_depth``, that depth is used in place of the Herron and Langway depth (this is what ``main.py`` does).

3. **The RCI is repeated.** The number of repeats is the estimated spin-up time divided by the length of the RCI, rounded to the nearest integer (and at least one). You can give the number yourself with ``num_reps``. For example, if the estimate is 255 years and the RCI is 5 years long, the RCI is repeated 51 times.

4. **The repeats are put in front of the real record.** The repeated copies of the RCI are shifted back in time so that they end where the real record begins, and the real record follows them. The model therefore sees one continuous forcing series. Its time axis starts before the first date of your record, for example at 1740 for the DYE-2 example (40 repeats of the six-year RCI 1980-1985, then the record from 1980). The first date of the real record is stored as ``forcing_data_start`` in the forcing dictionary.

5. **The model is run through the whole series.** ``FirnDensitySpin`` first builds a rough initial profile (a Herron and Langway density profile for the mean climate, run for ``yearSpin`` years), and then ``FirnDensityNoSpin.time_evolve`` runs the transient model over the entire series, repeats and record, in one loop. The repeats are not a separate stage, and the spin-up repeats are the same physics and the same time step as the rest of the run, so they take as long to compute as the same number of years of the real record. In the DYE-2 example, 240 of the 285 model years are spin-up.

The timeline of the DYE-2 example (RCI 1980-1985, six years) looks like this::

   model time   1740                                    1980                              2025
                |<---- 40 repeats of the RCI (240 yr) ---->|<------- the real record ------->|
                | RCI | RCI | RCI |   ...   | RCI | RCI    | RCI (1980-85) | after the RCI ...
                ^                                          ^
                start of the run                           forcing_data_start
                (initial profile from FirnDensitySpin)

Because the last repeat ends where the real record starts, and the first years of the record are the RCI itself, the firn enters the record in equilibrium with the climate of the RCI. After the RCI, the column responds to the real climate.

Setting it up
~~~~~~~~~~~~~

The call that prepares the forcing is:

.. code-block:: python

    import RCMpkl_to_spin as RCM

    climateTS, stpsPerYear, depth_S1, depth_S2, grid_bottom, SEBfluxes = RCM.makeSpinFiles(
        df,                      # the climate dataframe (or the name of a pickle)
        timeres='1d',            # model time step
        spin_date_st=1980.0,     # start of the RCI
        spin_date_end=1985.0,    # start of the last year of the RCI (RCI = 1980-01-01 to 1985-12-31)
        melt=True,
        desired_depth=None,      # None: use the Herron and Langway depth to rho_bottom
        rho_bottom=900,
        num_reps=None)           # None: estimate the number of repeats

    climateTS['forcing_data_start'] = 1980.0   # first date of the real record

The values that the function returns are used to set the configuration (see :doc:`running` for a complete script): ``stpsPerYear`` is the number of time steps per year, ``depth_S1`` and ``depth_S2`` are suggestions for the depths at which the grid of the model is coarsened (``grid1bottom`` and ``grid2bottom``), and ``grid_bottom`` is the depth to which the domain should extend (``HbaseSpin`` = ``H`` - ``grid_bottom``).

The configuration keys that matter for this spin-up are:

* **yearSpin** (:ref:`json keys <json-page>`): the length of the initial dry spin-up in ``FirnDensitySpin``. In this mode it only has to give a reasonable first guess, since the repeats do the real work; a small value (for example 1) is enough, and a long one wastes computing time.
* **NewSpin**: whether ``FirnDensitySpin`` runs again to make the initial profile. It is cheap in this mode, so ``true`` is a good choice. (The repeated RCI is part of the forcing, so it is always run.)
* **TWriteStart**: the time at which outputs start being written. Set it to the first date of the real record to leave out the spin-up repeats from the output file.
* **stpsPerYear**, **DFresample**, **nodestocombine**, **multnodestocombine**: the time step, and the merging of the deep nodes, which depends on it.

Choosing the RCI
~~~~~~~~~~~~~~~~

The RCI is an assumption about what the long-term climate of the site is, and what the column is in equilibrium with. Some things to consider:

* It should be representative of the climate the firn has experienced for the last decades to centuries. Medley et al. (2022) use 1980-1995 for the Greenland ice sheet, which is before the large increase in temperature and melt of the following years, and 1980-2019 for Antarctica, where there was no appreciable change in climate.
* The firn is assumed to have no net change over the RCI. If the RCI includes a trend, the column will be in equilibrium with an average that no part of the interval resembles.
* It should be long enough to include the variability that matters (the year-to-year variability in melt is an example). The longer it is, the fewer times it is repeated for the same spin-up time.
* ``spin_date_st`` should be the first date of the climate record and be at the start of a year. The repeats are taken from the RCI and tiled back from the first date of the record. If ``spin_date_st`` is not the first date, ``makeSpinFiles`` prints a warning, and the repeats will overlap or leave gaps where they meet the record.

Things to know
~~~~~~~~~~~~~~

* The mean temperature and accumulation rate that are used to estimate the spin-up time are means over the **whole** climate record that you give to ``makeSpinFiles``, not only over the RCI. A record that is much wetter or drier than the RCI changes the number of repeats.
* The accumulation rate that sets the spin-up time is snowfall plus sublimation (``bdm_sublim=True``, the default; sublimation is negative). Melt and rain are not subtracted or added. In Medley et al. (2022) the burial rate also includes melt.
* The number of repeats is rounded to the nearest integer, so the spin-up can be shorter than the estimate. Medley et al. (2022) round up. Use ``num_reps`` to set it yourself.
* If the mean accumulation rate of the record is zero the estimate is not defined and ``makeSpinFiles`` fails. For a site with no net accumulation (snowfall less than sublimation and melt), set ``num_reps`` and ``desired_depth`` yourself, and consider starting from a column of solid ice with ``iceblock`` (see :doc:`json`) so that the model builds the column from the seasonal snow.
* The rate at which ice leaves the bottom of the column (``iceout``) is set from the mean mass flux into the column over the spin-up part of the series, and is held for the whole run. The model prints it at the start of the run, and prints the mean runoff and a mass-consistent value at the end. See **manual_iceout** and **iceout**.
* The totals that are printed at the end of the run (melt, rain, refreezing, runoff) are for the real record only, not the repeats.
* ``main.py`` with ``"input_type": "dataframe"`` uses the default RCI, 1980-1995 (``spin_date_st=1980.0, spin_date_end=1995.0``), and the depth of the domain (``H`` - ``HbaseSpin``) in place of the Herron and Langway depth. They are not read from the configuration file. To use another RCI, call ``makeSpinFiles`` in your own script (see :doc:`running`).

.. _spinupdate-section:

Saving the spin-up and restarting a run (``spinUpdate``)
--------------------------------------------------------

The reference-climate-interval spin-up can be long (in the DYE-2 example, 240 of the 285 model years). With ``"spinUpdate": true`` the model saves its state while it runs, so that a run does not have to start again from the beginning. It was added in version 3.0.0. The state (density, temperature, age, depth, the mean accumulation rate, and the liquid water content and the grain size if they are modeled) is saved as a **snapshot** in the spin file (**spinFileName**, usually ``CFMspin.hdf5``, in **resultsFolder**). The snapshots are written:

* at the first time step after each 100 years of model time (for example at 1900, 2000, and so on; the interval is fixed in the code), counting from the start of the run, and
* at the time given by **spinUpdateDate**, which should be the end of the spin-up: the first date of the real record, which is the start of the reference climate interval (for example ``1980.0``). The snapshot is made at the first time step at or after this date.

The first snapshot in the file is the initial profile that ``FirnDensitySpin`` made. All of the snapshots are kept, each as a new column of the arrays in the spin file, so the file grows with the length of the run. The first row of each array holds the model time of the snapshot, so you can list the snapshots that are available with

.. code-block:: python

    import h5py
    with h5py.File('CFMoutput_example/CFMresults_66.5_-46.25_GSFC2020/CFMspin.hdf5', 'r') as f:
        print(f['densitySpin'][0, :])   # model times of the snapshots; the first is the initial profile

To restart from a snapshot, run the model again with the same configuration and forcing, with **spinUpdate** true, **NewSpin** false (do not use ``-n`` with ``main.py``), and **spinUpdateDate** set to the time of the snapshot that you want to restart from. The model uses the snapshot closest to **spinUpdateDate** and continues from the next time step, with the forcing from that date on. A restarted run does not write new snapshots. Its results file (``CFMresults.hdf5`` in **resultsFolder**) holds the restarted run only and replaces the earlier one, so rename or copy the earlier results if you want to keep them; the first row of the output is repeated (see the changelog for version 4.0.0). In a script it is

.. code-block:: python

    # config.json has "spinUpdate": true and "spinUpdateDate": 1900.0
    firn = FirnDensityNoSpin('config.json', climateTS=climateTS, NewSpin=False, SEBfluxes=SEBfluxes)
    firn.time_evolve()

This gives two uses:

* **Recovering a run that stopped** before it finished, for example because it ran out of time on a cluster. Restart it from the last snapshot that it wrote before it stopped.
* **Extending a run with more forcing without repeating the spin-up.** Do a first run, with **spinUpdate** true, so that it writes the snapshot at the end of the spin-up (**spinUpdateDate** = the first date of the record). When more climate data are available, restart from that snapshot with the longer forcing. The model starts at the end of the spin-up and runs through the whole, longer record. Only the forcing after the snapshot is used, but ``makeSpinFiles`` estimates the number of repeats from the means of the whole record, so give ``num_reps`` explicitly to keep the time axis (and so the times of the snapshots) the same.

The snapshots belong to the configuration that made them (the physics, the grid, the time step), so a restart should use the same configuration. Spin files written by older versions are not interchangeable with those of 4.0.0 and later, because the snapshot is now written after the regridding step.

Steady-state spin-up
--------------------

This is the original spin-up of the CFM, for runs forced by a long climate history, such as one derived from an ice core, at a resolution of years to centuries. It is used when the forcing is not prepared with an RCI, for example when it is read from csv files (see :doc:`inputs`).

1. ``FirnDensitySpin`` takes a temperature and an accumulation rate from the forcing: the mean of the forcing, or its first value (**spinup_climate_type**: ``mean`` or ``initial``). You can give your own long-term values with **manual_climate**, **deepT** and **bdot_long**, which is useful if the forcing is short and not representative of the long-term climate.
2. The initial column is a Herron and Langway density profile for that temperature and accumulation rate, and the surface density **rhos0**.
3. The model is run at that constant climate until the column is in steady state. The length of the run is **yearSpin** years, or, if **AutoSpinUpTime** is true, the depth at which the initial profile reaches 850 kg m\ :sup:`-3` divided by the accumulation rate. The spin-up uses the time step set by **stpsPerYear**.
4. The state at the end is written to the spin file (**spinFileName**, usually ``CFMspin.hdf5``, in the results folder). It is the initial condition for the transient run, which ``FirnDensityNoSpin`` then starts at the beginning of the forcing.

If the results folder already holds a spin file, the spin-up is not run again, and the transient run starts from the existing spin file, unless **NewSpin** is true (or ``-n`` is given to ``main.py``). The spin file is made from the configuration at the time (climate, physics, grid), so it should be made again if you change them.

To run a model without any spin-up, you can start from a uniform column with **iceblock** (a solid-ice column, or any density given by **iceblock_rho**).
