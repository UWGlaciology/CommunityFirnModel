import sys
import os
import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import math
import csv
import json


# Used to pull data from csv files and import in the data
# FOR NOW hard coded to just take in years and another set of data. in order by line: 1) years 2) other data
def file_impt(configName):
	temp_data = np.genfromtxt(configName, delimiter=',')
	imp_yrs = temp_data[0,:]
	imp_data = temp_data[1,:]

	return imp_yrs, imp_data


# 1-d Linear Interp, stps defined by user in JSON
def linear_interp(stps, years_t, temp, years_b, bDot):
	max_t = max(years_t)
	max_b = max(years_b)
	min_t = min(years_t)
	min_b = min(years_b)

	# Finds interp range
	interp_max = max_t if max_t <= max_b else max_b
	interp_min = min_t if min_t >= min_b else min_b


	temp_func = interpolate.interp1d(years_t, temp)
	bDot_func = interpolate.interp1d(years_b, bDot)

	temp_years = np.arange(interp_min, interp_max + stps, stps)

	if temp_years[len(temp_years) - 1] > interp_max:
		temp_years[len(temp_years) - 1] = interp_max

	interp_temp = temp_func.__call__(temp_years)
	interp_bDot = bDot_func.__call__(temp_years)

	name_t = 'intp_test_temp.csv'
	name_b = 'intp_test_bDot.csv'

	write_interp(name_t, name_b, temp_years, interp_temp, interp_bDot)

	return name_t, name_b


# Need to figure out what I want to add here in terms of the JSON file, and do we pass in the reader????
# for now this is pretty broken...
# ONLY DOES TEMPERATURE FOR NOW!!! AKA one at a time in terms of which data to synthetically create
# tt is for years/time the synth runs
#
def synth_data(typ, synth_info):
	# these may have to go inside a certai function... may not work well with other functions
	t_range = synth_info[0]
	t_pert = synth_info[1]
	T_pert = synth_info[2]
	T_0 = synth_info[3]
	b_pert = synth_info[4]
	b_0 = synth_info[5]
	#syn_step = synth_info[6]


	# hard coded for dev, normally in the list 'synth_info':
	syn_step = 15

	# time array in years:
	tt = np.arange(t_range[0], t_range[1] + t_pert, t_pert)

	if typ == 'step':
		# used to get number of data points
		T_1 = np.ones(len(tt)) * T_0
		b_1 = np.ones(len(tt)) * b_0
		ii = np.nonzero(tt >= t_pert + t_range[0])
		while (len(ii[0]) != 0):
			T_1[ii] = T_1[ii] + T_pert
			b_1[ii] = b_1[ii] + b_pert
			# Catches the index out of bounds... there aren't enough time steps to fill out the user wanted step length
			if len(ii[0]) >= syn_step:
				# this syn_step value.... should this be user generated? as in this is how long the number should be at the value...
				ii = np.nonzero(tt >= t_pert + tt[ii[0][syn_step - 1]])
			else:
				ii = np.nonzero(tt >= t_pert + tt[ii[0][len(ii[0]) - 1]])

		
		# for testing:
		print 'these are the years: ' + str(tt)
		print 'these are the temperature steps: ' + str(T_1)
		print 'these are the accu steps: ' + str(b_1)
		print
		print 'these are the lengths:'
		print 'length of tt: ' + str(len(tt))
		print 'length of Temp: ' + str(len(T_1))
		print 'length of Accu: ' + str(len(b_1))


# Helper function to write out the csv files
# Should i make it only one file write??????
# better naming of the file
def write_interp(name_t, name_b, years, intp_temp, intp_bDot):
	with open(name_t, 'wb') as csvfile:
		interp_writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
		interp_writer.writerow(years)
		interp_writer.writerow(intp_temp)

	with open(name_b, 'wb') as csvfile:
		interp_writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
		interp_writer.writerow(years)
		interp_writer.writerow(intp_bDot)
