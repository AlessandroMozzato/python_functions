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

def baro_stream(vel):
    ### This function calculates the barotropic streamfunction for a nt nz ny nx array returning a nt ny nx array                  
    ### The function works for both 18 and 36km and 9 km                                                                                    
    if vel.shape[2] == 192:
        x="/scratch/general/am8e13/results36km"
    elif vel.shape[2] == 384:
        x="/scratch/general/am8e13/results18km"
    elif vel.shape[2] == 768:
        x="/scratch/general/am8e13/results9km"
    
    os.chdir(x)
    file2read = netcdf.NetCDFFile("grid.nc",'r')
    hfacw = file2read.variables['HFacW']
    hfacw = hfacw[:]*1
    dyg = file2read.variables['dyG']
    dyg = dyg[:]*1
    drf = file2read.variables['drF']
    drf = drf[:]*1
    dydz = np.zeros_like(hfacw)
    psi = np.zeros_like(vel[:,1,:,:])
    # Volume calculation   
    for k in range(drf.shape[0]):                                                                                                           
        dydz[k,:,:] = drf[k]*np.multiply(dyg,hfacw[k,:,:])
    for temp in range(vel.shape[0]):
        utemp = np.zeros_like(hfacw)
        utemp[:,:,:] = dydz[:,:,:]*vel[temp,:,:,:]

        # vertical integration                                                                                                         
        utemp = np.nansum( utemp, 0 );
        # Integration horizontally                                                                                                     
        psi[temp,:,:] = np.cumsum( -utemp, axis = 0 );
        psi[:,hfacw[0,:,:]==0]=np.nan
        
    return psi / 10**6

def baro_stream_old(vel):
    ### This function calculates the barotropic streamfunction for a nt nz ny nx array returning a nt ny nx array                   
    ### The function works for both 18 and 36km and 9km                                                                                   
    if vel.shape[2] == 192:
        x="/scratch/general/am8e13/results36km"
    elif vel.shape[2] == 384:
        x="/scratch/general/am8e13/results18km"
    elif vel.shape[2] == 768:
        x="/scratch/general/am8e13/results9km"

    os.chdir(x)
    file2read = netcdf.NetCDFFile("grid.nc",'r')
    hfacw = file2read.variables['HFacW']
    hfacw = hfacw[:]*1
    dyg = file2read.variables['dyG']
    dyg = dyg[:]*1
    drf = file2read.variables['drF']
    drf = drf[:]*1
    dydz = np.zeros_like(hfacw)
    psi = np.zeros_like(vel[:,1,:,:])
    # Volume calculation                                                                                                                             
    for i in range(dyg.shape[0]):
        for j in range(dyg.shape[1]):
            for k in range(drf.shape[0]):
                dydz[k,i,j] = drf[k]*dyg[i,j]*hfacw[k,i,j]

    for temp in range(vel.shape[0]):
        utemp = np.zeros_like(hfacw)
        for i in range(dyg.shape[0]):
            for j in range(dyg.shape[1]):
                for k in range(drf.shape[0]):
                    utemp[k,i,j] = dydz[k,i,j]* vel[temp,k,i,j]

        # vertical integration                                                                                                                       
        utemp = np.nansum( utemp, 0 );
        # Integration horizontally                                                                                                                   
        psi[temp,:,:] = np.nancumsum( -utemp, axis = 0 );

    return psi / 10**6
