#! /usr/bin/env python                                                                                                             
# -*- coding: utf-8 -*-                                                                                                            
from scipy.io import netcdf
from scipy.io import loadmat
import numpy as np
from pylab import clf, plot, show, floor, ceil, imshow
import matplotlib
import matplotlib.pyplot as plt
import os
import csv
import sys
import glob
sys.path.append('/noc/users/am8e13/Python/python_functions/')
from barotropic import *
from topostrophy import *
from rho import *

def grid_read(res):
    ''' use: grid = grid_read(res) '''
    path = "/scratch/general/am8e13/results"+str(res)+"km/grid.nc"
    file2read = netcdf.NetCDFFile(path,'r')
    # Bathy is 1 on land and 0 over sea                                                                                                
    grid = {}
    for var in ['HFacC','HFacW','HFacS','YC','XC','Z','Y','X','rA']:
        temp = file2read.variables[var]
        grid[var] = temp[:]*1
    file2read.close()        
    area = np.zeros_like(grid['HFacC'])
    for z in range(len(grid['Z'])):
        area[z,:,:] = grid['HFacC'][z,:,:]*grid['rA']*(grid['Z'][z]-grid['Z'][z-1])
    grid['Area'] = area
    return grid
