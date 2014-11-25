'''
Created on Aug 7, 2013

@author: Paul
'''
import csv

import json

import sys

import numpy as np

#import scipy as sp

from scipy import interpolate 
from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr, spsolve
from collections import deque
import matplotlib.pyplot as plt
import math
import os
import sys

spot = os.path.dirname(sys.argv[0]) #Add resultsFolder 
sys.path.insert(1,os.path.join(spot,'Results'))

def plotData(GG, resultsFolder, plotsFolder):
        
    densityPath = os.path.join(resultsFolder, 'density.csv')
    tempPath = os.path.join(resultsFolder, 'temp.csv')
    agePath = os.path.join(resultsFolder, 'age.csv')
    depthPath = os.path.join(resultsFolder, 'depth.csv')
    
    densityPlotPath = os.path.join(plotsFolder, 'density_depth.png')
    agePlotPath = os.path.join(plotsFolder, 'age_depth.png')
    tempPlotPath = os.path.join(plotsFolder, 'temp_depth.png')
    
    initDepth = np.genfromtxt(depthPath, delimiter = ',')[1]
    initAge = np.genfromtxt(agePath, delimiter = ',')[1]
    initDensity = np.genfromtxt(densityPath, delimiter = ',')[1]
    initTemp = np.genfromtxt(tempPath, delimiter = ',')[1]
    
    finalDepth = np.genfromtxt(depthPath, delimiter = ',')[-1]
    finalAge = np.genfromtxt(agePath, delimiter = ',')[-1]
    finalDensity = np.genfromtxt(densityPath, delimiter = ',')[-1]
    finalTemp = np.genfromtxt(tempPath, delimiter = ',')[-1]

    plt.figure(0)
    p1, = plt.plot(initDensity[1:], initDepth[1:], 'ro')
    p2, = plt.plot(finalDensity[1:], finalDepth[1:], 'g^')
    plt.xlabel('Density (kg m$^{-3}$)',fontsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.title('Density and Depth')
    plt.legend([p1,p2], ["Initial", "Final"])
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.savefig(densityPlotPath)
    
    plt.figure(1)
    p1, = plt.plot(initAge[1:], initDepth[1:],'ro')
    p2, = plt.plot(finalAge[1:], finalDepth[1:], 'g^')
    plt.xlabel('Age (year)',fontsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.title('Age and Depth')
    plt.legend([p1,p2], ["Initial", "Final"])
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.savefig(agePlotPath)
    
    plt.figure(2)
    p1, = plt.plot(initTemp[1:], initDepth[1:],'ro')
    p2, = plt.plot(finalTemp[1:], finalDepth[1:], 'g^')
    plt.xlabel('Temperature (K)',fontsize=14)
    plt.ylabel('Depth (m)',fontsize=14)
    plt.title('Temperature and Depth')
    plt.legend([p1,p2], ["Initial", "Final"])
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.savefig(tempPlotPath)
    
    if GG == 'on':
        r2Path = os.path.join(resultsFolder, 'r2.csv')
        r2PlotPath = os.path.join(plotsFolder, 'r2_depth.png')
        initR2 = np.genfromtxt(r2Path, delimiter = ',')[1]
        finalR2 = np.genfromtxt(r2Path, delimiter = ',')[-1]
        plt.figure(3)
        p1, = plt.plot(initR2[1:], initDepth[1:],'ro')
        p2, = plt.plot(finalR2[1:], finalDepth[1:], 'g^')
        plt.xlabel('r$^{2}$',fontsize=14)
        plt.ylabel('Depth (m)',fontsize=14)
        plt.title('r$^{2}$ and Depth')
        plt.legend([p1,p2], ["Initial", "Final"])
        ax = plt.gca()
        ax.set_ylim(ax.get_ylim()[::-1])
        plt.savefig(r2PlotPath)
       
    plt.show()
        
if __name__ == '__main__':
    plotData('on', 'Results', 'Plots')