CFM outputs
===========

A CFM run writes three kinds of files into the folder **resultsFolder**:

* the **results file** (**resultsFileName**, ``CFMresults.hdf5`` in the examples), with the model output described on this page;
* the **spin file** (**spinFileName**, ``CFMspin.hdf5``): the initial profile and the snapshots that the model restarts from (see :doc:`spinup`);
* the **forcing file** (**forcingFileName**, ``CFMforcing.hdf5``), written when the forcing is given as a dataframe or ``climateTS`` dictionary: a record of the forcing that the model used, so that a run can be repeated.

The configuration file is also copied into the folder by ``main.py``. All of the files are in the HDF5 format and are read with ``h5py`` (below), or with any other HDF5 reader.

Which outputs are saved
-----------------------

The key **outputs** in the configuration file is a list of the names of the profiles to save (for example ``["density", "depth", "temperature", "age", "DIP", "LWC", "meltoutputs", "climate"]``). The model also adds some outputs by itself:

* If **MELT** is true, ``meltoutputs`` is added to the list (a message is printed). It stands for ``refreeze``, ``runoff``, ``meltvol`` and ``LWC``.
* If **MELT** is false, ``LWC`` is removed, since there is no liquid water. ``isotopes`` is removed if **isoDiff** is false, and ``gasses`` if **FirnAir** is false.
* With **FirnAir** true, the outputs of the firn-air module that are listed in ``AirConfig.json`` are added (see :doc:`json`).
* With **isoDiff** true, ``isotopes_<isotope>`` and ``iso_sig2_<isotope>`` are written for each isotope in **iso** (for example ``isotopes_d18O`` and ``iso_sig2_d18O``).
* ``grainsize`` is written under the name ``r2``, and ``temp_Hx`` is the temperature history of the Morris (2014) model.

The outputs are written at the times set by **TWriteInt** and **TWriteStart** (see :doc:`json`). Only the times that are written are in the file.

Layout of the file
------------------

Each output is a dataset in the results file. There are two kinds.

**Profiles** (2-D arrays, one row per time that was written). Column 0 of a row is the time (decimal date, model time), and the other columns are the values in the nodes (the layers of the model grid) at that time, from the surface down. The values correspond to the depths in the same row of ``depth``.

* **Row 0 is the initial state**, written at the model start time. With a reference-climate-interval spin-up that time is before the first date of the record (see :doc:`spinup`), so it is not at **TWriteStart**; leave row 0 out when you analyze the run. After a restart the first two rows have the same time.
* **Depth.** With **grid_outputs** false, the model nodes are written, and ``depth`` holds their depths at each time. With **grid_outputs** true (the usual choice), the profiles are interpolated onto a regular grid with the spacing **grid_output_res** that goes down to **grid_output_max_depth**. The grid is rebuilt at every write, so ``depth`` is a 2-D array too, with one row per write. **Rows are padded with trailing NaN** where the column is shorter than the grid, so use ``np.nanmax``, ``np.nanmean`` and similar. If the domain gets deeper than **grid_output_max_depth**, the output is cut off at that depth and a warning is printed once.
* Amounts per layer (``LWC``, ``refrozen``, ``compaction``, ``PLWC_mem``) are regridded so that their total is conserved. Intensive quantities (density, temperature, age) are interpolated linearly. ``Dcon`` is interpolated with the nearest value.
* **truncate_outputs** true saves only every fifth row of ``density``, ``temperature``, ``LWC`` and ``age`` (and row 0), to reduce the file size, so their rows no longer line up with the other outputs.
* The numbers are written with the precision of **output_bits** (``float32`` by default). Decimal dates around 2000 are therefore only precise to about 1e-4 years in float32; this is why the times in the examples look like 1980.0819.

**Vectors** (2-D arrays with a few columns: ``DIP``, ``BCO``, ``Modelclimate``, ``meltvol``, ``refreeze``, ``runoff``): column 0 is the time, and the other columns are the values at that time, as described below.

The profile outputs
-------------------

``depth``
    (m) The depth of each node (or grid point) below the surface.

``density``
    (kg m\ :sup:`-3`) Density.

``temperature``
    (K) Temperature.

``age``
    (years) Age of the firn. Column 0 is the decimal date, as in the other outputs (before 4.0.0 it was divided by the number of seconds in a year by mistake).

``r2``
    (mm\ :sup:`2`) Grain size (the mean squared grain radius), from the ``grainsize`` output. Needs **physGrain** true.

``temp_Hx``
    The temperature history of the firn, used by the ``Morris2014`` model (Morris and Wingham, 2014).

``LWC``
    (m) Volume of liquid water per node (m of water). Only with **MELT** true.

``refrozen``
    (m w.e.) Liquid water that refroze in that node during the time step (not accumulated over time). For each time step the nodes add up to the total in ``refreeze``; with the ``darcy`` scheme, water that refroze in a surface node and was melted again later in the same step is not included. It includes the refreezing in the melt scheme and in the heat diffusion that follows it.

``PLWC_mem``
    (m) The water that was in the preferential-flow domain (``prefsnowpack`` scheme); zero for the other schemes.

``compaction``
    (m) Thinning of each node by densification since the previous time step, positive when the node got thinner. (Horizontal strain is not included.) To get a compaction rate, divide by the time step; to get the compaction over an interval, sum the values of the rows (following the layers, not the fixed grid points). For each time step the nodes add up to minus ``DIP`` column 4, except that the compaction of nodes that melt or sublimate away during the same step is not in the nodes.

``viscosity``
    (Pa s) The viscosity of a Newtonian viscous fluid with the densification rate of the model, :math:`\eta = \rho\sigma/(2\,d\rho/dt)`, where :math:`\sigma` is the overburden stress. See :doc:`../extras/densification`.

``bdot_mean``
    (m a\ :sup:`-1` ice equivalent) The mean accumulation rate over the lifetime of each parcel of firn.

``Dcon``
    (advanced use) A layer-tracking routine. To use it you need to edit the code and program it how you want; for example you could set each node that has liquid water to 1 and all others to 0. Gridded outputs use the nearest value.

``DIPc``
    (m) The depth-integrated porosity from the surface down to each node (the cumulative sum of porosity times layer thickness); its last value is the firn air content in ``DIP`` column 1.

``isotopes_<isotope>``, ``iso_sig2_<isotope>``
    The isotope value (per mil) and the diffusion length variance (:math:`\sigma^2`) of each isotope in **iso**. Only with **isoDiff** true.

Firn-air outputs (the concentration of each gas under its name, ``diffusivity``, ``gas_age``, ``w_air``, ``w_firn``) are written when **FirnAir** is true; see :doc:`../files/firn_air`.

The vector outputs
------------------

``Modelclimate``
    (written when ``climate`` is in **outputs**) The forcing at each written time step. Six columns: time, accumulation (m a\ :sup:`-1` ice equivalent), skin temperature (K), melt (m a\ :sup:`-1` ice equivalent), rain (m a\ :sup:`-1` ice equivalent) and sublimation (m a\ :sup:`-1` ice equivalent, negative = mass loss). With the surface energy balance the melt is the melt that the SEB module calculated. These are the values of the time step; for the whole forcing, use the forcing file.

``DIP``
    The depth-integrated porosity (firn air content) and the change of the surface elevation. The columns (with Python's 0-based indexing) are: (0) time, (1) depth-integrated porosity, a.k.a. firn air content (m), to the bottom of the model domain, (2) surface elevation change during the last time step (m), (3) cumulative elevation change since the start of the model run (m), (4) change in thickness (m) of the firn column by densification alone since the last time step, negative when the column got thinner (it is minus the sum of ``compaction`` over the nodes; zero in the first written row), (5) 'corrected' surface elevation change during the last time step (m), (6) 'corrected' cumulative elevation change since the start of the model run (m), (7) firn air content to the horizon given by **DIPhorizon** (m).

    The surface elevation change of a step is the sum of the thickness change by densification, the thickness change by horizontal strain (if included), the thickness of the snow added at the surface, and the thickness removed from the top by melt and sublimation, minus the thickness of ice that leaves at the bottom of the column (the mean mass flux into the column, **iceout**, over the reference period). The thickness added by the bottom of the grid when layers melt away is not part of it. The uncorrected fields (2 and 3) take the ice that leaves at the bottom as a thickness of ice (917 kg m\ :sup:`-3`); that is only right if the bottom of the domain is at ice density. The 'corrected' fields (5 and 6) use the density of the layer that actually leaves, i.e. they account for the compaction that occurs between the bottom of the model domain and the depth where the column reaches ice density, in the case that you are not modeling to the ice density. In a column in equilibrium (the same forcing repeated for centuries) the corrected cumulative elevation change levels off at zero (to within about half a percent of the accumulation in the tests with Summit and DYE-2 forcing), while the uncorrected one keeps a constant drift when the bottom is lighter than ice. Use the corrected fields unless your domain reaches ice density.

    Firn air to a horizon (7) is there because, due to the model's Lagrangian format, the bottom of the domain varies in time, and if the bottom is at a density less than ice density it can be useful to consider the FAC to a specified depth.

    Use caution with the modeled elevation change fields: deeper nodes are merged to save computing time, and that can make the modeled dH have blips in it. Analyses based on the FAC (column 1 or 7) are more robust.

``BCO``
    Bubble close-off properties. Ten columns: time, Martinerie close-off age, Martinerie close-off depth, age of the 830 kg m\ :sup:`-3` density horizon, depth of the 830 kg m\ :sup:`-3` density horizon, Martinerie lock-in age, Martinerie lock-in depth, age of the 815 kg m\ :sup:`-3` density horizon, depth of the 815 kg m\ :sup:`-3` density horizon, depth of zero closed porosity. The Martinerie fields refer to parameterizations for close-off properties published by Patricia Martinerie in the early 1990s (see the references below).

``meltvol``, ``refreeze``, ``runoff``
    Two columns each: time, and the total melt volume, the total liquid water refreezing, or the total runoff at that time step (m w.e.). All four liquid schemes (``bucket``, ``darcy``, ``prefsnowpack``, ``resingledomain``) write them.

The forcing file
----------------

``CFMforcing.hdf5`` holds the group ``main`` with the arrays the model used: ``time`` (decimal date), ``TSKIN`` (K), ``BDOT``, ``SUBLIM``, ``RAIN`` and ``SMELT`` (m ice equivalent per year; see :doc:`inputs`), and ``forcing_data_start``, the first date of the record. The reference-climate-interval repeats of the spin-up are not in it: it starts at the beginning of the record. With the surface energy balance there is also a group ``SEB`` with the sub-daily SEB forcing (``dtRATIO`` is the number of SEB steps per model step). It is not written when the forcing is read from csv files.

The spin file
-------------

``CFMspin.hdf5`` has the datasets ``densitySpin``, ``tempSpin``, ``ageSpin``, ``depthSpin``, ``bdot_meanSpin`` and, when they are used, ``LWCSpin``, ``r2Spin``, ``HxSpin``, ``gridSpin`` (the tracking of the merged grid) and the isotope profiles. Each column is a snapshot (the first is the end of the spin-up, the others are the ``spinUpdate`` snapshots), and the first row is the time of the snapshot. See :doc:`spinup` for how restarts use them.

Reading the results
-------------------

.. code-block:: python

    import h5py
    import numpy as np

    with h5py.File('CFMoutput_example/df/CFMresults.hdf5', 'r') as f:
        print(list(f.keys()))        # which outputs are in the file
        t     = f['density'][:, 0]   # time of each row (decimal date)
        rho   = f['density'][:, 1:]  # density, one row per time, one column per node
        z     = f['depth'][:, 1:]    # depth; 2-D, with trailing NaN if the grid is shorter
        fac   = f['DIP'][:, 1]       # firn air content (m)
        mclim = f['Modelclimate'][:] # time, accumulation, Ts, melt, rain, sublimation

    rec = t >= 1980.0                # leave out row 0 (the initial state) and any spin-up
    # depth of the 830 kg m-3 level at the last time step (the first grid point is the surface)
    z830 = np.interp(830, rho[-1][~np.isnan(rho[-1])], z[-1][~np.isnan(z[-1])])

With ``grid_outputs`` false and an older results file (before 4.0.0), ``depth`` is a single row (a 1-D axis); see the changelog for the format changes of each version.

**References**
*Goujon, C., Barnola, J.-M., and Ritz, C. (2003), Modeling the densification of polar firn including heat diffusion: Application to close-off characteristics and gas isotopic fractionation for Antarctica and Greenland sites, J. Geophys. Res., 108, 4792, doi:10.1029/2002JD003319, D24.*

*Martinerie, P., Lipenkov, V. Y., Raynaud, D., Chappellaz, J., Barkov, N. I., and Lorius, C. (1994), Air content paleo record in the Vostok ice core (Antarctica): A mixed record of climatic and glaciological parameters, J. Geophys. Res., 99( D5), 10565– 10576, doi:10.1029/93JD03223.*

*Martinerie, P., Raynaud, D., Etheridge, D. M., Barnola, J.-M., & Mazaudier, D. (1992, August). Physical and climatic parameters which influence the air content in polar ice. Earth and Planetary Science Letters. Elsevier BV. https://doi.org/10.1016/0012-821x(92)90002-d*

