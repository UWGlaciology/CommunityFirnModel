import csv
import os
from string import join
import numpy as np
import h5py

def write_nospin_init(folder, physGrain, THist, rho_time, Tz_time, age_time, z_time, D_time, Clim_time, bdot_time, r2_time, Hx_time):
    '''
    Writes the results of density, temperature, age, depth, Dcon, Clim, bdot, and r2 (if specified)
    Not used in spin mode

    :param folder: the name of the folder that the results will be written to
    :param physGrain: specifies whether the grain growth is turned on, and whether r2 exists
    :param rho_time:
    :param Tz_time:
    :param age_time:
    :param z_time:
    :param D_time:
    :param Clim_time:
    :param bdot_time:
    :param r2_time:
    '''

    densityPath     = os.path.join(folder, 'density.csv')
    tempPath        = os.path.join(folder, 'temp.csv')
    agePath         = os.path.join(folder, 'age.csv')
    depthPath       = os.path.join(folder, 'depth.csv')
    DconPath        = os.path.join(folder, 'Dcon.csv')
    ClimPath        = os.path.join(folder, 'Clim.csv')
    bdotPath        = os.path.join(folder, 'bdot_mean.csv')

    with open(densityPath, "w") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(rho_time)
    with open(tempPath, "w") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(Tz_time)
    with open(agePath, "w") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(age_time)
    with open(depthPath, "w") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(z_time)
    with open(DconPath, "w") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(D_time)
    with open(ClimPath, "w") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(Clim_time)
    with open(bdotPath, "w") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(bdot_time)
    if physGrain:
        r2Path = os.path.join(folder, 'r2.csv')
        with open(r2Path, "w") as f:
            csvwriter = csv.writer(f)
            csvwriter.writerow(r2_time)
    if THist:
        HxPath = os.path.join(folder, 'Hx.csv')
        with open(HxPath, "w") as f:
            csvwriter = csv.writer(f)
            csvwriter.writerow(Hx_time)

# def write_nospin_hdf5(folder, physGrain, THist, rho_out, Tz_out, age_out, z_out, D_out, Clim_out, bdot_out, r2_out, Hx_out):
def write_nospin_hdf5(self):

    f4 = h5py.File(os.path.join(self.c['resultsFolder'], self.c['resultsFileName']),'w')
    f4.create_dataset('density',data=self.rho_out)
    f4.create_dataset('temperature',data=self.Tz_out)
    f4.create_dataset('age',data=self.age_out)
    f4.create_dataset('depth',data=self.z_out)
    f4.create_dataset('Dcon',data=self.D_out)
    f4.create_dataset('bdot',data=self.bdot_out)
    f4.create_dataset('Modelclimate',data=self.Clim_out)
    if self.c['physGrain']:
        f4.create_dataset('r2',data=self.r2_out)
    if self.THist:
        f4.create_dataset('Hx',data=self.Hx_out)
    if self.c['isoDiff']:
        f4.create_dataset('isotopes',data=self.iso_out)

    # timewrite = np.append(self.modeltime[0],self.TWrite[:len(self.intPhiAll)])
    # timewrite = self.TWrite_out
    # timewrite = np.append(self.modeltime[0],self.TWrite_out)

    # DIPwrite=np.vstack((timewrite, self.intPhiAll, self.dHOut, self.dHOutC))
    f4.create_dataset('DIP',data = self.DIP_out)

    # BCOwrite=np.vstack((timewrite, self.bcoAgeMartAll, self.bcoDepMartAll, self.bcoAge815All, self.bcoDep815All))
    f4.create_dataset('BCO',data = self.BCO_out)

    # LIZwrite=np.vstack((timewrite, self.LIZAgeAll, self.LIZDepAll))
    f4.create_dataset('LIZ',data = self.LIZ_out)

    f4.close()

# def write_nospin(folder, physGrain, THist, rho_time, Tz_time, age_time, z_time, D_time, Clim_time, bdot_time, r2_time, Hx_time):
#     '''
#     Writes the results of density, temperature, age, depth, Dcon, Clim, bdot, and r2 (if specified)
#     Not used in spin mode

#     :param folder: the name of the folder that the results will be written to
#     :param physGrain: specifies whether the grain growth is turned on, and whether r2 exists
#     :param rho_time:
#     :param Tz_time:
#     :param age_time:
#     :param z_time:
#     :param D_time:
#     :param Clim_time:
#     :param bdot_time:
#     :param r2_time:
#     '''

#     densityPath     = os.path.join(folder, 'density.csv')
#     tempPath        = os.path.join(folder, 'temp.csv')
#     agePath         = os.path.join(folder, 'age.csv')
#     depthPath       = os.path.join(folder, 'depth.csv')
#     DconPath        = os.path.join(folder, 'Dcon.csv')
#     ClimPath        = os.path.join(folder, 'Clim.csv')
#     bdotPath        = os.path.join(folder, 'bdot_mean.csv')

#     with open(densityPath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(rho_time)
#     with open(tempPath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(Tz_time)
#     with open(agePath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(age_time)
#     with open(depthPath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(z_time)
#     with open(DconPath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(D_time)
#     with open(ClimPath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(Clim_time)
#     with open(bdotPath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(bdot_time)
#     if physGrain:
#         r2Path = os.path.join(folder, 'r2.csv')
#         with open(r2Path, "a") as f:
#             csvwriter = csv.writer(f)
#             csvwriter.writerow(r2_time)
#     if THist:
#         HxPath = os.path.join(folder, 'Hx.csv')
#         with open(HxPath, "a") as f:
#             csvwriter = csv.writer(f)
#             csvwriter.writerow(Hx_time)

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

# def write_spin(folder, physGrain, THist, rho_time, Tz_time, age_time, z_time, r2_time, Hx_time):
#     '''
#     Writes the results of density, temperature, age, depth, and r2 (if specified)
#     Used in spin mode

#     :param folder: the name of the folder that the results will be written to
#     :param physGrain: specifies whether the grain growth is turned on, and whether r2 should be written
#     :param rho_time:
#     :param Tz_time:
#     :param age_time:
#     :param z_time:
#     :param r2_time:
#     '''

#     densityPath = os.path.join(folder, 'densitySpin.csv')
#     tempPath    = os.path.join(folder, 'tempSpin.csv')
#     agePath     = os.path.join(folder, 'ageSpin.csv')
#     depthPath   = os.path.join(folder, 'depthSpin.csv')

#     if physGrain:
#         r2Path = os.path.join(folder, 'r2Spin.csv')

#     with open(densityPath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(rho_time)
#     with open(tempPath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(Tz_time)
#     with open(agePath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(age_time)
#     with open(depthPath, "a") as f:
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(z_time)
#     if physGrain:
#         r2Path = os.path.join(folder, 'r2Spin.csv')
#         with open(r2Path, "a") as f:
#             csvwriter = csv.writer(f)
#             csvwriter.writerow(r2_time)
#     if THist:
#         HxPath = os.path.join(folder, 'HxSpin.csv')
#         with open(HxPath, "a") as f:
#             csvwriter = csv.writer(f)
#             csvwriter.writerow(Hx_time)

# def write_nospin_BCO(folder, bcoAgeMartAll, bcoDepMartAll, bcoAge815All, bcoDep815All,modeltime,TWrite):
#     '''
#     Writes the results of bubble close-off depth and age
#     Not used in spin mode

#     :param folder: the name of the folder that the results will be written to
#     :param bcoAgeMartAll:
#     :param bcoDepMartAll:
#     :param bcoAge815All:
#     :param bcoDep815All:
#     '''

#     bcoPath = os.path.join(folder, 'BCO.csv')
#     with open(bcoPath, "w") as f: #write BCO.csv file. rows are: time, BCO age (mart), BCO depth (mart),BCO age (815), BCO depth (815)
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(np.append(modeltime[0],TWrite[:len(bcoAge815All)]))
#         csvwriter.writerow(bcoAgeMartAll)
#         csvwriter.writerow(bcoDepMartAll)
#         csvwriter.writerow(bcoAge815All)
#         csvwriter.writerow(bcoDep815All)

# def write_nospin_LIZ(folder, LIZAgeAll, LIZDepAll,modeltime,TWrite):
#     '''
#     Writes the results of lock-in zone depth and age
#     Not used in spin mode

#     :param folder: the name of the folder that the results will be written to
#     :param LIZAgeAll:
#     :param LIZDepAll:
#     '''

#     lidPath = os.path.join(folder, 'LID.csv')
#     with open(lidPath, "w") as f: #write LIZ information. Rows: time, age, dep
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(np.append(modeltime[0],TWrite[:len(LIZDepAll)]))
#         csvwriter.writerow(LIZAgeAll)
#         csvwriter.writerow(LIZDepAll)

# # def write_nospin_DIP(folder, intPhiAll, dsurfOut, dsurfOutC, modeltime,TWrite):
# def write_nospin_DIP(folder, intPhiAll, modeltime,TWrite):
#     '''
#     Writes the results of depth-integrated porosity
#     Not used in spin mode

#     :param folder: the name of the folder that the results will be written to
#     :param intPhiAll: the depth-integrated porosity at that time step
#     :param dsurfOut: the change in surface elevation since the last time step
#     :param dsurfOutC: the total change in surface elevation since start of model run
#     '''

#     intPhiPath = os.path.join(folder, 'porosity.csv')
#     with open(intPhiPath, "w") as f: #depth-integrated porosity
#         csvwriter = csv.writer(f)
#         csvwriter.writerow(np.append(modeltime[0],TWrite[:len(intPhiAll)]))
#         csvwriter.writerow(intPhiAll)
#         # csvwriter.writerow(dsurfOut)
#         # csvwriter.writerow(dsurfOutC)
