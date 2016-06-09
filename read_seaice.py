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

def read_seaicedata():
    class seaice():
        def __init__(self):
            self.seaice = {}
            self.title = 'NSIDC'
            self.lat = []
            self.lon = []
    nsidc = seaice()
    path = '/scratch/general/am8e13/seaice/'
    filename = 'seaice.nc'
    file2read = netcdf.NetCDFFile(path+filename,'r')
    temp = file2read.variables['SIarea']
    nsidc.seaice['SIarea'] = temp[:]*1
    temp = file2read.variables['lat']
    nsidc.lat = temp[:]*1
    temp = file2read.variables['lon']
    nsidc.lon = temp[:]*1
    return nsidc