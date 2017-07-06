import os
import numpy as np
from string import join
from constants import *
import h5py

# def read_temp(file):
#     '''
#     Read in data for initial temperatures

#     :param file: name of the file which holds the temperature data

#     :return input_temp: temperature vector from a specified csv file
#     :return input_year_temp: corresponding time vector (in years)
#     '''
#     # spot = os.path.dirname(sys.argv[0])
#     spot = os.getcwd()

#     FID_temp        = os.path.join(spot, file)
#     data_temp       = np.genfromtxt(FID_temp, delimiter=',')
#     input_year_temp = data_temp[0, :]
#     input_temp      = data_temp[1, :]
#     if input_temp[0] < 0.0:
#         input_temp = input_temp + K_TO_C

#     return input_temp, input_year_temp

# def read_bdot(file):
#     '''
#     Read in data for initial accumulation rates

#     :param file: name of the file which holds the accumulation rate data

#     :return input_bdot: accumulation rate vector from a specified csv file
#     :return input_year_bdot: corresponding time vector (in years)
#     '''

#     spot = os.getcwd()

#     FID_bdot        = os.path.join(spot, file)
#     data_bdot       = np.genfromtxt(FID_bdot, delimiter=',')
#     input_year_bdot = data_bdot[0, :]
#     input_bdot      = data_bdot[1, :]

#     return input_bdot, input_year_bdot

def read_input(filename):
    '''
    Read in data from csv input files

    :param file: name of the file which holds the accumulation rate data

    :return input_bdot: accumulation rate vector from a specified csv file
    :return input_year_bdot: corresponding time vector (in years)
    '''

    spot = os.getcwd()

    FID        = os.path.join(spot, filename)
    data       = np.loadtxt(FID, delimiter=',') #changed 3/6/17 to loadtxt from genfromtxt; much faster
    input_year = data[0, :]
    input_data = data[1, :]

    return input_data, input_year


# def read_input2(file):
#     '''
#     Read in data for initial accumulation rates

#     :param file: name of the file which holds the accumulation rate data

#     :return input_bdot: accumulation rate vector from a specified csv file
#     :return input_year_bdot: corresponding time vector (in years)
#     '''

#     spot = os.getcwd()

#     FID        = os.path.join(spot, file)
#     data       = np.genfromtxt(FID, delimiter=',')
#     input_year = data[0, :]
#     input_data = data[1, :]

#     return input_data, input_year

def read_snowmelt(file):
    '''
    Read in data for initial melt rates

    :param file: name of the file which holds the accumulation rate data

    :return input_bdot: accumulation rate vector from a specified csv file
    :return input_year_bdot: corresponding time vector (in years)
    '''

    spot = os.getcwd()

    FID_melt        = os.path.join(spot, file)
    data_melt       = np.genfromtxt(FID_melt, delimiter=',')
    input_year_melt = data_melt[0, :]
    input_melt      = data_melt[1, :]

    return input_snowmelt, input_year_snowmelt

def read_snowmelt(file):
    '''
    Read in data for initial melt rates

    :param file: name of the file which holds the accumulation rate data

    :return input_bdot: accumulation rate vector from a specified csv file
    :return input_year_bdot: corresponding time vector (in years)
    '''

    spot = os.getcwd()

    FID_melt        = os.path.join(spot, file)
    data_melt       = np.genfromtxt(FID_melt, delimiter=',')
    input_year_melt = data_melt[0, :]
    input_melt      = data_melt[1, :]

    return input_snowmelt, input_year_snowmelt

# def read_init(folder):
#     '''
#     Read in data for initial depth, age, density, and temperature to run the model without spin

#     :param folder: the folder containing the files holding depth, age, density, and temperature

#     :return initDepth: initial depth vector from a specified csv file
#     :return initAge: initial age vector from a specified csv file
#     :return initDensity: initial density vector from a specified csv file
#     :return initTemp: initial temperature vector from a specified csv file
#     '''

#     densityPath = os.path.join(folder, 'densitySpin.csv')
#     tempPath    = os.path.join(folder, 'tempSpin.csv')
#     agePath     = os.path.join(folder, 'ageSpin.csv')
#     depthPath   = os.path.join(folder, 'depthSpin.csv')

#     initDepth   = np.genfromtxt(depthPath, delimiter = ',')
#     initAge     = np.genfromtxt(agePath, delimiter = ',' )
#     initDensity = np.genfromtxt(densityPath, delimiter = ',')
#     initTemp    = np.genfromtxt(tempPath, delimiter = ',')

#     return initDepth, initAge, initDensity, initTemp

# def read_init(folder, resultsFileName):
#     '''
#     Read in data for initial depth, age, density, and temperature to run the model without spin

#     :param folder: the folder containing the files holding depth, age, density, and temperature

#     :return initDepth: initial depth vector from a specified csv file
#     :return initAge: initial age vector from a specified csv file
#     :return initDensity: initial density vector from a specified csv file
#     :return initTemp: initial temperature vector from a specified csv file
#     '''

#     f5 = h5py.File(os.path.join(folder, resultsFileName),'r')

#     initDensity = f5['densitySpin'][:]
#     initAge = f5['ageSpin'][:]
#     initDepth = f5['depthSpin'][:]
#     initTemp = f5['tempSpin'][:]

#     # densityPath = os.path.join(folder, 'densitySpin.csv')
#     # tempPath    = os.path.join(folder, 'tempSpin.csv')
#     # agePath     = os.path.join(folder, 'ageSpin.csv')
#     # depthPath   = os.path.join(folder, 'depthSpin.csv')

#     # initDepth   = np.genfromtxt(depthPath, delimiter = ',')
#     # initAge     = np.genfromtxt(agePath, delimiter = ',' )
#     # initDensity = np.genfromtxt(densityPath, delimiter = ',')
#     # initTemp    = np.genfromtxt(tempPath, delimiter = ',')

#     return initDepth, initAge, initDensity, initTemp

def read_init(folder, resultsFileName, varname):
    '''
    Read in data for initial depth, age, density, and temperature to run the model without spin

    :param folder: the folder containing the files holding depth, age, density, and temperature

    :return initDepth: initial depth vector from a specified csv file
    :return initAge: initial age vector from a specified csv file
    :return initDensity: initial density vector from a specified csv file
    :return initTemp: initial temperature vector from a specified csv file
    '''

    f5 = h5py.File(os.path.join(folder, resultsFileName),'r')
    init_value = f5[varname][:]


    # initAge = f5['ageSpin'][:]
    # initDepth = f5['depthSpin'][:]
    # initTemp = f5['tempSpin'][:]

    # initDensity = f5['densitySpin'][:]
    # initAge = f5['ageSpin'][:]
    # initDepth = f5['depthSpin'][:]
    # initTemp = f5['tempSpin'][:]

    # densityPath = os.path.join(folder, 'densitySpin.csv')
    # tempPath    = os.path.join(folder, 'tempSpin.csv')
    # agePath     = os.path.join(folder, 'ageSpin.csv')
    # depthPath   = os.path.join(folder, 'depthSpin.csv')

    # initDepth   = np.genfromtxt(depthPath, delimiter = ',')
    # initAge     = np.genfromtxt(agePath, delimiter = ',' )
    # initDensity = np.genfromtxt(densityPath, delimiter = ',')
    # initTemp    = np.genfromtxt(tempPath, delimiter = ',')

    return init_value
