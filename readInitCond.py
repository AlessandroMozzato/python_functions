#! /usr/bin/env python                                                                                                                 
# -*- coding: utf-8 -*-     

from scipy.io import netcdf
import numpy as np
from pylab import clf, plot, show, floor, ceil, imshow
import os
import csv
import sys
import glob
from netCDF4 import Dataset

from IPython.display import display, Math, Latex, Image
from mpl_toolkits.basemap import Basemap
import mpl_toolkits

sys.path.append('/noc/users/am8e13/PyNGL-1.4.0/lib/python2.7/site-packages/PyNGL/')
import Ngl
sys.path.append('/noc/users/am8e13/Python/')
import komod
sys.path.append('/noc/users/am8e13/Python/PyNIO-1.4.0/')
import Nio


class Woa:
    def __init__(self):
        self.lat = run36.lat
        self.lon = run36.lon
        file1 = '/hpcdata/scratch/am8e13/cs_36km_tutorial/run_input/WOA05_THETA_JAN_210x192x50_arctic'
        self.T = komod.mitbin(file1,xdim=192,ydim=210,zdim=50,datatype='float32')
        self.T = np.squeeze(self.T,axis=0)
        self.T[run36.bathy == 0] = np.nan
        file1 = '/hpcdata/scratch/am8e13/cs_36km_tutorial/run_input/WOA05_SALT_JAN_210x192x50_arctic'
        self.S = komod.mitbin(file1,xdim=192,ydim=210,zdim=50,datatype='float32')
        self.S = np.squeeze(self.S,axis=0)
        self.S[run36.bathy == 0] = np.nan
        self.depth = run36.depth
        self.title = 'Woa'

class Phc:
    def __init__(self):
        self.lat = run36.lat
        self.lon = run36.lon
        file1 = '/hpcdata/scratch/am8e13/cs_36km_tutorial/run_input/theta.phc_210_192_50_12'
        self.T = komod.mitbin(file1,xdim=192,ydim=210,zdim=50,datatype='float32')
        self.T = np.squeeze(self.T,axis=0)
        self.T[run36.bathy == 0] = np.nan
        file1 = '/hpcdata/scratch/am8e13/cs_36km_tutorial/run_input/salt.phc_210_192_50_12'
        self.S = komod.mitbin(file1,xdim=192,ydim=210,zdim=50,datatype='float32')
        self.S = np.squeeze(self.S,axis=0)
        self.S[run36.bathy == 0] = np.nan
        self.depth = run36.depth
        self.title = 'PHC'

