CFM outputs
===========

The CFM writes its outputs to a single.hdf5-format file. By default, all nodes are written to file. The output is only saved at the time steps specified by the user with the variable TWrite. Most of the outputs should be self-explanatory. Many of them are big 2D matrices; the first column is time, and the values throughout are the particular values at the depth found in the corresponding cell in the depth output. Set the outputs you want in the .json file. The available outputs are: 
depth: (m) The depth of each model node.
density: (kg m-3) The density at the depths in ‘depth 
temperature: (K) Temperature at the depths in ‘depth’
age: (years) Firn Age at the depths in ‘depth’
dcon: Dcon is a layer-tracking routine; to use it you need to dig into the code a bit and program it how you want, but for example you could set it up so that each model node that has liquid water gets a 1 and all others get a zero. Corresponds to depth.

bdot_mean: (m a-1 ice equivalent) the mean accumulation rate over the lifetime of each parcel of firn, corresponds with ‘depth’

climate: The temperature (K) and accumulation rate (m a-1 ice equivalent) at each time step – useful if using interpolation to find determine the climate.

compaction: (m) Total compaction of each node since the previous time step; corresponds to ‘depth’. To get compaction rate you need to divide by the time-step size. To get compaction over an interval you need to sum numerous boxes.

grainsize: (mm2) the grain size of the firn, corresponds to ‘depth’
temp_Hx: the temperature history of the firn (See Morris and Wingham, 2014)

isotopes: (per mil) water isotope values, corresponds to ‘depth’

LWC: (m3) volume of liquid present in that node, corresponds to ‘depth’

DIP: the depth-integrated porosity and change in surface elevation. 4 columns: The first is time, second is DIP to the bottom of the model domain (m), third is change in domain thickness since last time step (m), fourth is change in domain thickness since start of model run (m). 
DIP also saves a variable called DIPc, which is a matrix of the cumulative porosity to the depth in ‘depth’

BCO: bubble close-off properties. 10 columns: time, Martinerie close-off age, Marinerie close-off depth, age of 830 kg m-3 density horizon, depth of 830 kg m-3 density horizon, Martinerie lock-in age, Marinerie lock-in depth, age of 815 kg m-3 density horizon, depth of 815 kg m-3 density horizon, depth of zero porosity.

FirnAir: only works if FirnAir is true in the config.json. Saves gas concentrations, diffusivity profile, gas age, and advection rates of air and firn, all corresponding to ‘depth’.
