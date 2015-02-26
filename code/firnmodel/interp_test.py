# THIS FILE IS FOR TESTING THE data_interp.py SCRIPT. DOESN'T ACTUALLY DO ANYTHING JUST SEES IF THE FUNCTION WORKS


import sys
import os
import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import math
import csv
import json
import data_interp as IntpData

years_t, data_temp = IntpData.file_impt('GISP2_temp.csv')
years_b, data_bDot = IntpData.file_impt("GISP2_accu.csv")

# Here, hard coded for ONLY linear interpolation. Future, add more interp choices ** 2/4/15 mikeoon
intpStp = 10.
intp_name_t, intp_name_b = IntpData.linear_interp(intpStp, years_t, data_temp, years_b, data_bDot)

input_year_temp, input_temp = IntpData.file_impt(intp_name_t)
input_year_accu, input_AR = IntpData.file_impt(intp_name_b)