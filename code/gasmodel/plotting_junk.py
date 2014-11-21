import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

nodes=d['nodes']
d15N2=d['d15N2']

fig1=plt.figure(1)
plt.plot(nodes,d15N2[:,-1])
plt.show()