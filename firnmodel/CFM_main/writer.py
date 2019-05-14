#!/usr/bin/env python
import csv
import os
import numpy as np
import h5py

def write_nospin_hdf5(self):

    f4 = h5py.File(os.path.join(self.c['resultsFolder'], self.c['resultsFileName']),'w')
    
    if 'density' in self.output_list:
        f4.create_dataset('density',data = self.rho_out)
    if 'temperature' in self.output_list:
        f4.create_dataset('temperature',data = self.Tz_out)
    if 'age' in self.output_list:
        f4.create_dataset('age',data = self.age_out[-1,:])
        # f4.create_dataset('age',data = self.age_out) # use this one if you want a matrix of ages
    if 'depth' in self.output_list:    
        f4.create_dataset('depth',data = self.z_out)
    if 'dcon' in self.output_list:    
        f4.create_dataset('Dcon',data = self.D_out)
    if 'bdot_mean' in self.output_list:    
        f4.create_dataset('bdot',data = self.bdot_out)
    if 'climate' in self.output_list:    
        f4.create_dataset('Modelclimate',data = self.Clim_out)
    if 'compaction' in self.output_list:    
        f4.create_dataset('compaction', data = self.comp_out)
    if self.c['FirnAir']:
        if "gasses" in self.cg['outputs']:
            for gas in self.cg['gaschoice']:       
                f4.create_dataset(gas, data = self.gas_out[gas])
        if "diffusivity" in self.cg['outputs']:
            f4.create_dataset('diffusivity', data = self.diffu_out)
        if "gas_age" in self.cg['outputs']:
            f4.create_dataset('gas_age', data = self.gas_age_out)
        if "advection_rate" in self.cg['outputs']:
            f4.create_dataset('w_air', data = self.w_air_out)
            f4.create_dataset('w_firn', data = self.w_firn_out)
    if 'grainsize' in self.output_list:
        f4.create_dataset('r2', data = self.r2_out)
        f4.create_dataset('dr2_dt', data = self.dr2_dt_out)
    if 'temp_Hx' in self.output_list:    
        f4.create_dataset('Hx',data = self.Hx_out)
    if 'isotopes' in self.output_list:    
        f4.create_dataset('isotopes',data = self.iso_out)
    if 'LWC' in self.output_list:
        f4.create_dataset('LWC',data = self.LWC_out)
    if 'DIP' in self.output_list:
        f4.create_dataset('DIP',data = self.DIP_out)
        f4.create_dataset('DIPc', data = self.DIPc_out)
    if 'BCO' in self.output_list:
        f4.create_dataset('BCO',data = self.BCO_out) 
    if 'LIZ' in self.output_list:
        f4.create_dataset('LIZ',data = self.LIZ_out)
    if 'refrozen' in self.output_list:
        f4.create_dataset('refrozen',data = self.refrozen_out)
    if 'runoff' in self.output_list:
        f4.create_dataset('runoff',data = self.runoff_out)
    if 'viscosity' in self.output_list:
        f4.create_dataset('viscosity',data = self.viscosity_out)

    f4.close()

def write_spin_hdf5(self):

    f5 = h5py.File(os.path.join(self.c['resultsFolder'], self.c['spinFileName']), 'w')

    f5.create_dataset('densitySpin', data = self.rho_time)
    f5.create_dataset('tempSpin', data = self.Tz_time)
    f5.create_dataset('ageSpin', data = self.age_time)
    f5.create_dataset('depthSpin', data = self.z_time)
    if self.c['physGrain']:
        f5.create_dataset('r2Spin', data = self.r2_time)
    if self.THist:
        f5.create_dataset('HxSpin', data = self.Hx_time)
    if self.c['isoDiff']:
        f5.create_dataset('IsoSpin', data = self.iso_time)
    if self.doublegrid:
        f5.create_dataset('gridSpin', data = self.grid_time)
    if self.c['MELT']: #VV
        f5.create_dataset('LWCSpin', data = self.LWC_time)
    f5.close()
