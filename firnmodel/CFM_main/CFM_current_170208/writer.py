import csv
import os
import numpy as np
import h5py

def write_nospin_hdf5(self):

    f4 = h5py.File(os.path.join(self.c['resultsFolder'], self.c['resultsFileName']),'w')
    
    if 'density' in self.output_list:
        f4.create_dataset('density',data=self.rho_out)
    if 'temperature' in self.output_list:
        f4.create_dataset('temperature',data=self.Tz_out)
    if 'age' in self.output_list:
        f4.create_dataset('age',data=self.age_out)
    if 'depth' in self.output_list:    
        f4.create_dataset('depth',data=self.z_out)
    if 'dcon' in self.output_list:    
        f4.create_dataset('Dcon',data=self.D_out)
    if 'bdot_mean' in self.output_list:    
        f4.create_dataset('bdot',data=self.bdot_out)
    if 'climate' in self.output_list:    
        f4.create_dataset('Modelclimate',data=self.Clim_out)
    if 'compaction' in self.output_list:    
        f4.create_dataset('compaction_rate', data=self.crate_out)
    if 'gasses' in self.output_list:       
        f4.create_dataset('gasses', data=self.gas_out)
    if 'grainsize' in self.output_list:
        f4.create_dataset('r2',data=self.r2_out)
    if 'temp_Hx' in self.output_list:    
        f4.create_dataset('Hx',data=self.Hx_out)
    if 'isotopes' in self.output_list:    
        f4.create_dataset('isotopes',data=self.iso_out)
    if 'LWC' in self.output_list:
        f4.create_dataset('LWC',data=self.LWC_out)
    if 'DIP' in self.output_list:
        f4.create_dataset('DIP',data = self.DIP_out)  
    if 'BCO' in self.output_list:
        f4.create_dataset('BCO',data = self.BCO_out) 
    if 'LIZ' in self.output_list:
        f4.create_dataset('LIZ',data = self.LIZ_out)

    f4.close()

def write_spin_hdf5(folder, spinFileName, physGrain, THist, isoDiff, rho_time, Tz_time, age_time, z_time, r2_time, Hx_time, iso_time):

    f5 = h5py.File(os.path.join(folder, spinFileName), 'w')

    f5.create_dataset('densitySpin', data = rho_time)
    f5.create_dataset('tempSpin', data = Tz_time)
    f5.create_dataset('ageSpin', data = age_time)
    f5.create_dataset('depthSpin', data = z_time)
    if physGrain:
        f5.create_dataset('r2Spin', data = r2_time)
    if THist:
        f5.create_dataset('HxSpin', data = Hx_time)
    if isoDiff:
        f5.create_dataset('IsoSpin', data = iso_time)
    f5.close()
