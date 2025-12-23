CFM outputs
===========

The CFM writes its outputs to a single.hdf5-format file. In older CFM versions, all nodes (aka model grid layers) were written to file. This is still an option, but there is also now an option to interpolate the results onto a regular grid (using "grid_outputs": [true/false] in the .json file, setting the grid resolution using "grid_output_res").

The outputs are only saved at the time steps specified by the user with the variable TWrite. Most of the outputs should be self-explanatory. Many of them are big 2D matrices; the first column is time, and the values throughout are the particular values at the depth found in the corresponding cell in the depth output. Set the outputs you want in the .json file. The available outputs are:

depth: (m) The depth of each model node. If "grid_outputs" = True, this will be a vector of depths. (The first value is the decimal date of the model initialization time, which you should ignore.)

MATRIX OUTPUTS:

density: (kg m-3) The density at the depths in ‘depth'

temperature: (K) Temperature at the depths in ‘depth’

age: (years) Age of the firn at the depths in ‘depth’

grainsize: (mm2) the grain size of the firn, corresponds to ‘depth’
temp_Hx: the temperature history of the firn (See Morris and Wingham, 2014)

isotopes: (per mil) water isotope values, corresponds to ‘depth’

LWC: (m3) volume of liquid present in that node, corresponds to ‘depth’

compaction: (m) Total compaction of each node since the previous time step; corresponds to ‘depth’. To get compaction rate you need to divide by the time-step size. To get compaction over an interval you need to sum numerous boxes.

dcon: (advanced use) Dcon is a layer-tracking routine; to use it you need to dig into the code a bit and program it how you want, but for example you could set it up so that each model node that has liquid water gets a 1 and all others get a zero. Corresponds to depth vector/matrix.

bdot_mean: (m a-1 ice equivalent) the mean accumulation rate over the lifetime of each parcel of firn, corresponds with ‘depth’

FirnAir: only works if FirnAir is true in the config.json. Saves gas concentrations, diffusivity profile, gas age, and advection rates of air and firn, all corresponding to ‘depth’.

VECTOR OUTPUTS:

forcing: This is the climate forcing inputs, i.e. it is just keeping a record of the data that is fed into the CFM. I included it to enable reproducibility; i.e., if you need to, you should be able to use these to create new forcing files to repeat a model run if you need to. 5 columns (as of this documentation update, though you can check writer.py, in the 'if forcing_dict' block, for any updates there.): decimal date, skin temperature, accumulation, snowmelt, rain.

Modelclimate: 3 columns: timestep, temperature (K), and accumulation rate (m a-1 ice equivalent). This is useful if using interpolation to find determine the climate (e.g. long model runs for paleoclimate). This is probably not useful if you are doing runs with 'timesetup'='exact', as would be expected e.g. for runs investigating modern height and mass change.

DIP: the depth-integrated porosity and change in surface elevation. The columns in DIP are (with python's 0-based indexing):
(0) timestep, (1) Depth Integrated Porosity, aka Firn Air Content (m), to the bottom of the model domain, (2) surface elevation change since last time step (m), (3) cumulative elevation change since start of model run (m), (4) total compaction (m) of the firn column since last time step, (5) 'corrected' cumulative elevation change since start of model run (m), (6) 'corrected' cumulative elevation change since start of model run (m), (7) Firn Air Content above a specified horizon.

The 'corrected' fields (5 and 6) tries to account for any compaction that occurs between the bottom of the model domain and the depth where the column reaches ice density (917), in the case that you are not modeling to the ice density.

Firn air to a horizon (7) is there because due to the model's lagrangian format, the bottom of the domain varies in time, and if the bottom is at a density less than ice density it can be useful to consider the FAC to a specified depth,.

(I'd advise using caution if you use the modeled elevation change fields - we have made changes to the code to merge deeper nodes together to save computing time, but in doing so that can make the model-outputed dH have blips in it. If you want to work with this field let me know and I can do a bit of testing. I'd recommend structuring your analyses around the FAC variable.)

BCO: bubble close-off properties. 10 columns: time, Martinerie close-off age, Martinerie close-off depth, age of 830 kg m-3 density horizon, depth of 830 kg m-3 density horizon, Martinerie lock-in age, Martinerie lock-in depth, age of 815 kg m-3 density horizon, depth of 815 kg m-3 density horizon, depth of zero closed porosity. The Martinerie fields refer to paramaterizations for close-off properties published by Patricia Martinerie in the early 1990s. See references below. 

meltvol: two columnes: decimal time, total melt volume at that time step [m w.e.]

refreeze: two columns: decimal time, total liquid water refreezing at that time step [m w.e.]

runoff: two columns: decimal time, total liquid water runoff at that time step [m w.e.]

**References**
*Goujon, C., Barnola, J.-M., and Ritz, C. (2003), Modeling the densification of polar firn including heat diffusion: Application to close-off characteristics and gas isotopic fractionation for Antarctica and Greenland sites, J. Geophys. Res., 108, 4792, doi:10.1029/2002JD003319, D24.*

*Martinerie, P., Lipenkov, V. Y., Raynaud, D., Chappellaz, J., Barkov, N. I., and Lorius, C. (1994), Air content paleo record in the Vostok ice core (Antarctica): A mixed record of climatic and glaciological parameters, J. Geophys. Res., 99( D5), 10565– 10576, doi:10.1029/93JD03223.*

*Martinerie, P., Raynaud, D., Etheridge, D. M., Barnola, J.-M., & Mazaudier, D. (1992, August). Physical and climatic parameters which influence the air content in polar ice. Earth and Planetary Science Letters. Elsevier BV. https://doi.org/10.1016/0012-821x(92)90002-d*

