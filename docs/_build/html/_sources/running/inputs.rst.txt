Inputs for the CFM
==================

The CFM is forced by surface-temperature and accumulation-rate boundary conditions. Additionally, the user can specify the surface-melt, surface-density and water-isotope values. These files are .csv formatted. The first row of these files is time (decimal date, i.e. 2015.3487) and the second row is the corresponding temperature/accumulation rate/boundary condition value at that time. Time must be going forward, i.e. the first column is a date some time ago and the last column is the most recent. (If the model is being forced with ice-core data, the user must be careful to ensure this is the case as ice-core data are often presented as years before present.) The times in the various input files do not need to be the same; they are interpolated onto a common axis. The units for temperature can be K or C. The CFM uses K, but it will change the temperature to K if you use C. The units for accumulation rate/surface mass balance are m ice equivalent per year (see note in section 5.3)

:ref:`Test Linking Pages <json-page>`

