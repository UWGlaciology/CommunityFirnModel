import sys
import os
import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import math
import csv
import json
import data_interp


# this file is used to test the Synthetic Step Data creator
testA = [[55., 1000.], 10, 2, 240, 0.02, 0.1]

data_interp.synth_data('step', testA)
