#! /usr/bin/env python                                                                                                                         
# -*- coding: utf-8 -*-                                                                                                                        
                                                                                                                                               
from scipy.io import netcdf
from scipy.io import loadmat
import numpy as np
from pylab import clf, plot, show, floor, ceil, imshow
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import csv
import sys
import glob
sys.path.append('/noc/users/am8e13/Python/python_functions/')
from barotropic import *
from topostrophy import *
from rho import *
from jmd95 import *
from regions_def import *

def read_mxld_data():
    file2read = netcdf.NetCDFFile('/scratch/general/am8e13/mxldepth/MIMOC_ML_v2.2_PT_S_MLP_month01.nc','r')
    temp = file2read.variables['LATITUDE']
    lat = temp[:]*1
    temp = file2read.variables['LONGITUDE']
    lon = temp[:]*1
    lon, lat = np.meshgrid(lon,lat)
    mxldepth = np.zeros((12,lon.shape[0],lon.shape[1]))
    for j in range(1,13,1):
        if j < 10:
            mon = '0'+str(j)
        else:
            mon = str(j)
        file2read = netcdf.NetCDFFile('/scratch/general/am8e13/mxldepth/MIMOC_ML_v2.2_PT_S_MLP_month'+mon+'.nc','r')
        temp = file2read.variables['DEPTH_MIXED_LAYER']
        mxldepth[j-1,:,:] = temp[:]*1
    mxldepth_mean = np.nanmean(mxldepth,axis=0)
    return mxldepth,mxldepth_mean,lon,lat
def mxld_dic():
    class mixl():
        def __init__(self):
            self.lat = []
            self.lon = []
    mxl = mixl()
    lat_ts=90.0 ; lat_0=90.0 ; lon_0=-45.0
    sgn=1 ; width=7000000. ;height=7000000.0
    fig, axes = plt.subplots(1,1)
    m = Basemap(ax=axes,width=width,height=height,resolution='l',\
                    projection='stere',lat_ts=lat_ts,lat_0=lat_0,lon_0=lon_0)
    mxldepth,mxldepth_mean,lon,lat = read_mxld_data()
    lon_or = np.array(lon)
    mxl.mxldepth = np.zeros_like(mxldepth)
    for j in range(12): 
        lons, mxl.mxldepth[j,:,:] = m.shiftdata(lon_or, datain = mxldepth[j,:,:], lon_0=0)
    mxl.mxldepth_mean = np.nanmean(mxl.mxldepth[j,:,:],axis=0)
    mxl.lat = lat
    mxl.lon = lons
    mxl.title = 'MIMOC'
    return mxl
