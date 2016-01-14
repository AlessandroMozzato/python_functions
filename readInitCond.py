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

from rho import *
from jmd95 import *

class Woa:
    def __init__(self):
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/results36km/grid.nc",'r')
        bathy = file2read.variables['HFacC']
        self.bathy = bathy[:]*1
        Z = file2read.variables['Z']
        Z = Z[:]*1
        self.Z = Z
        lat = file2read.variables['YC']
        self.lat = lat[:]*1
        lon = file2read.variables['XC']
        self.lon = lon[:]*1
        file1 = '/hpcdata/scratch/am8e13/cs_36km_tutorial/run_input/WOA05_THETA_JAN_210x192x50_arctic'
        self.T = komod.mitbin(file1,xdim=192,ydim=210,zdim=50,datatype='float32')
        self.T = np.squeeze(self.T,axis=0)
        self.T[self.bathy == 0] = np.nan
        file1 = '/hpcdata/scratch/am8e13/cs_36km_tutorial/run_input/WOA05_SALT_JAN_210x192x50_arctic'
        self.S = komod.mitbin(file1,xdim=192,ydim=210,zdim=50,datatype='float32')
        self.S = np.squeeze(self.S,axis=0)
        self.S[self.bathy == 0] = np.nan
        self.rho = rho(self.S,self.T)
        self.rhop = rhop(self.S,self.T)
        self.depth = Z
        # calculate jmd 95 density                                                                                                                                       
        self.rho_jmd = np.zeros_like(self.rho)
        for z in range(len(self.Z)):
            self.rho_jmd[z,:,:] = densjmd95(self.S[z,:,:],self.T[z,:,:],-9.81*self.Z[z]*1025)

        self.title = 'Woa'
        self.dataDyn = {}
        temp_lv = np.nanmean(np.nanmean(self.T,axis=1),axis=1)
        self.dataDyn['theta_lv_mean'] = (np.ones((400,1))*temp_lv)
        temp_lv = np.nanmean(np.nanmean(self.S,axis=1),axis=1)
        self.dataDyn['salt_lv_mean'] = (np.ones((400,1))*temp_lv)
        self.dataDyn['rho_lv_mean'] = rhop(self.dataDyn['salt_lv_mean'],self.dataDyn['theta_lv_mean'])

class Phc:
    def __init__(self):
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/results36km/grid.nc",'r')
        bathy = file2read.variables['HFacC']
        self.bathy = bathy[:]*1
        Z = file2read.variables['Z']
        Z = Z[:]*1
        self.Z = Z
        lat = file2read.variables['XC']
        self.lat = lat[:]*1
        lon = file2read.variables['YC']
        self.lon = lon[:]*1
        file1 = '/hpcdata/scratch/am8e13/cs_36km_tutorial/run_input/theta.phc_210_192_50_12'
        self.T = komod.mitbin(file1,xdim=192,ydim=210,zdim=50,datatype='float32')
        self.T = np.squeeze(self.T,axis=0)
        self.T[self.bathy == 0] = np.nan
        file1 = '/hpcdata/scratch/am8e13/cs_36km_tutorial/run_input/salt.phc_210_192_50_12'
        self.S = komod.mitbin(file1,xdim=192,ydim=210,zdim=50,datatype='float32')
        self.S = np.squeeze(self.S,axis=0)
        self.S[self.bathy == 0] = np.nan
        self.rho = rho(self.S,self.T)
        self.rhop = rhop(self.S,self.T)
        self.depth = Z
        # calculate jmd 95 density                                                                                                                                       
        self.rho_jmd = np.zeros_like(self.rho)
        for z in range(len(self.Z)):
            self.rho_jmd[z,:,:] = densjmd95(self.S[z,:,:],self.T[z,:,:],-9.81*self.Z[z]*1025)

        self.title = 'PHC'
        self.dataDyn = {}
        temp_lv = np.nanmean(np.nanmean(self.T,axis=1),axis=1)
        self.dataDyn['theta_lv_mean'] = (np.ones((400,1))*temp_lv)
        temp_lv = np.nanmean(np.nanmean(self.S,axis=1),axis=1)
        self.dataDyn['salt_lv_mean'] = (np.ones((400,1))*temp_lv)
        self.dataDyn['rho_lv_mean'] = rhop(self.dataDyn['salt_lv_mean'],self.dataDyn['theta_lv_mean'])
