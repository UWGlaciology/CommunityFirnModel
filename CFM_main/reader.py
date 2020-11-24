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
    if xx>yy:
        input_year = data[:, 0]
        input_data = data[:, 1]
    else:        
        input_year = data[0, :]
        input_data = data[1, :]

    if StartDate==None:
        pass
    else:
        StartInd = np.where(input_year>=StartDate)[0]
        input_year = input_year[StartInd]
        input_data = input_data[StartInd]

    return input_data, input_year

def read_init(folder, resultsFileName, varname):

    '''
    Read in data for initial depth, age, density, and temperature to run the model without spinup

    :param folder: the folder containing the files holding depth, age, density, and temperature

    '''
    f5          = h5py.File(os.path.join(folder, resultsFileName),'r')
    init_value  = f5[varname][:]
    f5.close()

    return init_value


# def read_snowmelt(file):
#     '''
#     Read in data for initial melt rates

#     :param file: name of the file which holds the accumulation rate data

#     :return input_bdot: accumulation rate vector from a specified csv file
#     :return input_year_bdot: corresponding time vector (in years)
#     '''

#     spot = os.getcwd()

#     FID_melt        = os.path.join(spot, file)
#     data_melt       = np.genfromtxt(FID_melt, delimiter=',')
#     input_year_melt = data_melt[0, :]
#     input_melt      = data_melt[1, :]

#     return input_snowmelt, input_year_snowmelt

# def read_snowmelt(file):
#     '''
#     Read in data for initial melt rates

#     :param file: name of the file which holds the accumulation rate data

#     :return input_bdot: accumulation rate vector from a specified csv file
#     :return input_year_bdot: corresponding time vector (in years)
#     '''

#     spot = os.getcwd()

#     FID_melt        = os.path.join(spot, file)
#     data_melt       = np.genfromtxt(FID_melt, delimiter=',')
#     input_year_melt = data_melt[0, :]
#     input_melt      = data_melt[1, :]

#     return input_snowmelt, input_year_snowmelt


