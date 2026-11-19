#!usr/bin/env python
'''
Functions to read model inputs.
'''

import os
import numpy as np
# from string import join
from constants import *
import h5py

def read_input(filename,StartDate=None):
    '''
    Read in data from csv input files

    :param filename: name of the file which holds the accumulation rate data

    :return input_data: vector of field of interest (e.g. temperature, accumulation rate from a specified csv file
    :return input_year: corresponding time vector (in years)
    '''

    spot = os.getcwd()

    FID        = os.path.join(spot, filename)
    data       = np.loadtxt(FID, delimiter=',') #changed 3/6/17 to loadtxt from genfromtxt; much faster
    xx,yy = np.shape(data)
    n_data = np.minimum(xx, yy)
    if n_data==2:
        if xx>yy:
            input_year = data[:, 0]
            input_data = data[:, 1]
        else:
            input_year = data[0, :]
            input_data = data[1, :]
    else:  # Load input with multiple data rows.
        if xx>yy:
            input_year = data[:, 0]
            input_data = np.transpose(data[:, 1:n_data])
        else:
            input_year = data[0, :]
            input_data = data[1:n_data, :]

    input_year_full = input_year.copy()
    input_data_full = input_data.copy()

    if StartDate==None:
        pass
    else:
        StartInd = np.where(input_year>=StartDate)[0]
        input_year = input_year[StartInd]
        input_data = input_data[StartInd]

    return input_data, input_year, input_data_full, input_year_full

def read_init(folder, resultsFileName, varname, udate = None):

    '''
    Read in data for initial depth, age, density, and temperature to run the model without spinup

    :param folder: the folder containing the files holding depth, age, density, and temperature

    '''
    f5          = h5py.File(os.path.join(folder, resultsFileName),'r')
    
    init_value_full  = f5[varname][:]
    if init_value_full.ndim==1:
        init_value = init_value_full
    else:
        restart_times = init_value_full[0,:]
        if udate is not None:
            res_index = (np.abs(restart_times-udate)).argmin()
        else:
            res_index = -1
        init_value = init_value_full[:,res_index]
    f5.close()

    return init_value



### Forcing file (CFMforcing.hdf5) ###
######################################

FORCING_FILENAME = 'CFMforcing.hdf5'
RCI_KEYS = ('sds', 'sde', 'num_reps')


def expand_rci(record, sds, sde, num_reps):
    '''
    Put the repeated reference climate interval (RCI) in front of a forcing record.

    This is the spin-up series that RCMpkl_to_spin.makeSpinFiles builds: the part of the
    record from sds to the end of the year that starts at sde is repeated num_reps times, shifted
    back in time so that the last repeat ends where the record starts, and the record follows. It is
    used by the forcing writer and reader, so that a forcing file does not need to hold the repeats.

    :param record: dictionary of the record. 'time' (decimal date) and arrays of the same length as
        'time'. Scalar values (for example forcing_data_start, dtRATIO) are passed on unchanged.
    :param sds: decimal date of the start of the RCI (spin_date_st of makeSpinFiles)
    :param sde: decimal date of the start of the last year of the RCI (spin_date_end of makeSpinFiles)
    :param num_reps: number of repeats of the RCI

    :return: dictionary with the same keys, with the repeats in front of the record
    '''
    t = np.asarray(record['time'])
    RCI_length = sde - sds + 1
    years = num_reps * RCI_length
    sub = np.arange(-1 * years, 0, RCI_length)
    msk = (t >= sds) & (t < sde + 1)

    out = {}
    for k, v in record.items():
        if np.ndim(v) == 0:
            out[k] = v
        elif k == 'time':
            out[k] = np.concatenate(((sub[:, np.newaxis] + t[msk]).flatten(), t))
        else:
            out[k] = np.concatenate((np.tile(v[msk], len(sub)), v))
    return out


def _read_forcing_group(g):
    '''Read one group of the forcing file (main or SEB), putting back the repeats if they were removed.'''
    d = {k: g[k][()] for k in g.keys() if isinstance(g[k], h5py.Dataset) and k not in RCI_KEYS}
    if g.attrs.get('rci_trimmed', False):
        d = expand_rci(d, float(g['sds'][()]), float(g['sde'][()]), int(g['num_reps'][()]))
        if 'patch' in g:  # the few elements that the repeats do not reproduce (see writer._trim_rci)
            for k in g['patch'].keys():
                d[k] = d[k].copy()
                d[k][g['patch'][k]['index'][()]] = g['patch'][k]['value'][()]
    return d


def read_forcing(path, filename=FORCING_FILENAME):
    '''
    Read the forcing that a model run used, from its forcing file. With the configuration file of the
    run (and the same random_seed), this is what is needed to repeat the run::

        climateTS, SEBfluxes, strain_rates = read_forcing('CFMoutput_example/CFMresults_66.5_-46.25_GSFC2020')
        firn = FirnDensityNoSpin(config, climateTS=climateTS, SEBfluxes=SEBfluxes, strain_rates=strain_rates, NewSpin=True)

    The repeats of the reference climate interval are not stored in the file (if the writer could remove
    them without changing the series); they are put back here.

    :param path: the results folder, or the forcing file
    :param filename: name of the forcing file in the folder (the configuration key forcingFileName)

    :return climateTS: dictionary of the forcing (as returned by makeSpinFiles)
    :return SEBfluxes: dictionary of the surface energy balance forcing, or None if the run did not have one
    :return strain_rates: the strain_rates argument of the run, or None
    '''
    if os.path.isdir(path):
        path = os.path.join(path, filename)
    with h5py.File(path, 'r') as f:
        climateTS = _read_forcing_group(f['main'])
        SEBfluxes = _read_forcing_group(f['SEB']) if 'SEB' in f else None
        strain_rates = None
        if 'strain_rates' in f:
            if isinstance(f['strain_rates'], h5py.Group):
                strain_rates = {k: f['strain_rates'][k][()] for k in f['strain_rates'].keys()}
            else:
                strain_rates = f['strain_rates'][()]
    return climateTS, SEBfluxes, strain_rates
