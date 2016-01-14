#! /usr/bin/env python                                                                                                                     
# -*- coding: utf-8 -*-                                                                                                                   

from scipy.io import netcdf
import numpy as np

def rho(s,t):
    # This function calculates the density temperature and salinity                                                                        
    s0 = 35
    t0 = 5
    alpha = 0.0002
    beta = 0.0008
    rho0 = 1027.5
    return rho0*(1 - alpha*(t - t0) + beta*(s - s0))

def pressure(s,t):
    file2read = netcdf.NetCDFFile('/scratch/general/am8e13/results36km/grid.nc','r')
    Z=file2read.variables['Z']
    Z=Z[:]*1
    pres = np.zeros_like(s)
    g = 9.81
    for z in range(len(Z)):
        pres[z,:] = -g*Z[z]*1027.5
    return pres

def rhop(s,t):
    # This function calculates the density temperature and salinity                                                                                                      
    s0 = 35
    t0 = 5
    alpha = 0.0002
    beta = 0.0008
    k = 4.1*10**(-11)
    rho0 = 1027.5
    pres = pressure(s,t)
    return rho0*(1 - alpha*(t - t0) + beta*(s - s0) + k*pres )
