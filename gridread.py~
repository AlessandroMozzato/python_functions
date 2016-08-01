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
    for var in ['HFacC','HFacW','HFacS','YC','XC','Z','Y','X','rA','Depth']:
        temp = file2read.variables[var]
        grid[var] = temp[:]*1
    file2read.close()
    area = np.zeros_like(grid['HFacC'])
    Zp = [10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
             10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04 , 19.82, 24.85,
             31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
             93.96, 96.58, 98.25, 99.25,100.01,101.33,104.56,111.33,122.83,
             139.09,158.94,180.83,203.55,226.50,249.50,272.50,295.50,318.50,
             341.50,364.50,387.50,410.50,433.50,456.50 ]
    for z in range(len(grid['Z'])):
        area[z,:,:] = grid['HFacC'][z,:,:]*grid['rA']*Zp[z]
    grid['Area'] = area
    grid['Zp'] = Zp
    return grid
