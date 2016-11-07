import csv
import os
from string import join

def write_nospin(folder, physGrain, rho_time, Tz_time, age_time, z_time, D_time, Clim_time, bdot_time, r2_time):
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

    densityPath = os.path.join(folder, 'density.csv')
    tempPath    = os.path.join(folder, 'temp.csv')
    agePath     = os.path.join(folder, 'age.csv')
    depthPath   = os.path.join(folder, 'depth.csv')
    DconPath    = os.path.join(folder, 'Dcon.csv')
    ClimPath    = os.path.join(folder, 'Clim.csv')
    bdotPath    = os.path.join(folder, 'bdot_mean.csv')

    if c["physGrain"]:
        r2Path = os.path.join(folder, 'r2.csv')

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
    if c["physGrain"]:
        with open(r2Path, "w") as f:
            csvwriter = csv.writer(f)
            csvwriter.writerow(r2_time)

def write_spin(folder, physGrain, rho_time, Tz_time, age_time, z_time, r2_time):
    '''
    Writes the results of density, temperature, age, depth, and r2 (if specified)
    Used in spin mode

    :param folder: the name of the folder that the results will be written to
    :param physGrain: specifies whether the grain growth is turned on, and whether r2 should be written
    :param rho_time:
    :param Tz_time:
    :param age_time:
    :param z_time:
    :param r2_time:
    '''

    densityPath = os.path.join(folder, 'densitySpin.csv')
    tempPath    = os.path.join(folder, 'tempSpin.csv')
    agePath     = os.path.join(folder, 'ageSpin.csv')
    depthPath   = os.path.join(folder, 'depthSpin.csv')

    if c['physGrain']:
        r2Path = os.path.join(folder, 'r2Spin.csv')

    with open(densityPath, "a") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(rho_time)
    with open(tempPath, "a") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(Tz_time)
    with open(agePath, "a") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(age_time)
    with open(depthPath, "a") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(z_time)
    if c['physGrain']:
        with open(r2Path, "a") as f:
            csvwriter = csv.writer(f)
            csvwriter.writerow(r2_time)

def write_nospin_BCO(folder, bcoAgeMartAll, bcoDepMartAll, bcoAge815All, bcoDep815All):
    '''
    Writes the results of bubble close-off depth and age
    Not used in spin mode

    :param folder: the name of the folder that the results will be written to
    :param bcoAgeMartAll:
    :param bcoDepMartAll:
    :param bcoAge815All:
    :param bcoDep815All:
    '''

    bcoPath = os.path.join(folder, 'BCO.csv')
    with open(bcoPath, "w") as f: #write BCO.csv file. rows are: time, BCO age (mart), BCO depth (mart),BCO age (815), BCO depth (815)
        csvwriter = csv.writer(f)
        csvwriter.writerow(np.append(modeltime[0],TWrite[:len(bcoAge815All)]))
        csvwriter.writerow(bcoAgeMartAll)
        csvwriter.writerow(bcoDepMartAll)
        csvwriter.writerow(bcoAge815All)
        csvwriter.writerow(bcoDep815All)

def write_nospin_LIZ(folder, LIZAgeAll, LIZDepAll):
    '''
    Writes the results of lock-in zone depth and age
    Not used in spin mode

    :param folder: the name of the folder that the results will be written to
    :param LIZAgeAll:
    :param LIZDepAll:
    '''

    lidPath = os.path.join(folder, 'LID.csv')
    with open(lidPath, "w") as f: #write LIZ information. Rows: time, age, dep
        csvwriter = csv.writer(f)
        csvwriter.writerow(np.append(modeltime[0],TWrite[:len(LIZDepAll)]))
        csvwriter.writerow(LIZAgeAll)
        csvwriter.writerow(LIZDepAll)

def write_nospin_DIP(folder, intPhiAll):
    '''
    Writes the results of depth-integrated porosity
    Not used in spin mode

    :param folder: the name of the folder that the results will be written to
    :param intPhiAll:
    '''

    intPhiPath = os.path.join(folder, 'porosity.csv')
    with open(intPhiPath, "w") as f: #depth-integrated porosity
        csvwriter = csv.writer(f)
        csvwriter.writerow(np.append(modeltime[0],TWrite[:len(intPhiAll)]))
        csvwriter.writerow(intPhiAll)

