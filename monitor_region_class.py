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

def monitor(x,iter_list):
    theta_mean_tot =[]
    theta_max_tot =[]
    theta_min_tot =[]
    eta_mean_tot =[]
    eta_max_tot =[]
    eta_min_tot =[]
    salt_mean_tot =[]
    salt_max_tot =[]
    salt_min_tot =[]
    sss_mean_tot =[]
    sss_max_tot =[]
    sss_min_tot =[]
    sst_mean_tot =[]
    sst_max_tot =[]
    sst_min_tot =[]
    vvel_mean_tot =[]
    vvel_max_tot =[]
    vvel_min_tot =[]
    uvel_mean_tot =[]
    uvel_max_tot =[]
    uvel_min_tot =[]
    ke_mean_tot =[]
    ke_max_tot =[]
    ke_vol_tot =[]
    time_seconds_tot = []
    
    for iter in iter_list:
        zerosto_add = '0'
        for zeros in range(9 - len(str(iter))):
            zerosto_add = '0'+zerosto_add
        
        file2read = netcdf.NetCDFFile(x+"monitor."+zerosto_add+str(iter)+".t001.nc",'r')

        time_seconds = file2read.variables['time_secondsf']
        time_seconds=time_seconds[:]*1    
        theta_mean=file2read.variables['dynstat_theta_mean']
        theta_mean=theta_mean[:]*1    
        theta_max=file2read.variables['dynstat_theta_max']
        theta_max=theta_max[:]*1
        theta_min=file2read.variables['dynstat_theta_min']
        theta_min=theta_min[:]*1    
        salt_mean=file2read.variables['dynstat_salt_mean']
        salt_mean=salt_mean[:]*1
        salt_max=file2read.variables['dynstat_salt_max']
        salt_max=salt_max[:]*1
        salt_min=file2read.variables['dynstat_salt_min']
        salt_min=salt_min[:]*1    
        sst_mean=file2read.variables['dynstat_sst_mean']
        sst_mean=sst_mean[:]*1
        sst_max=file2read.variables['dynstat_sst_max']
        sst_max=sst_max[:]*1
        sst_min=file2read.variables['dynstat_sst_max']
        sst_min=sst_min[:]*1    
        sss_mean=file2read.variables['dynstat_sss_mean']
        sss_mean=sss_mean[:]*1
        sss_mim=file2read.variables['dynstat_sss_min']
        sss_min=sss_mean[:]*1
        sss_max=file2read.variables['dynstat_sss_max']
        sss_max=sss_max[:]*1    
        eta_mean=file2read.variables['dynstat_eta_mean']
        eta_mean=eta_mean[:]*1
        eta_min=file2read.variables['dynstat_eta_min']
        eta_min=eta_min[:]*1
        eta_max=file2read.variables['dynstat_eta_max']
        eta_max=eta_max[:]*1    
        uvel_mean=file2read.variables['dynstat_uvel_mean']
        uvel_mean=uvel_mean[:]*1
        uvel_max=file2read.variables['dynstat_uvel_max']
        uvel_max=uvel_max[:]*1
        uvel_min=file2read.variables['dynstat_uvel_min']
        uvel_min=uvel_min[:]*1    
        vvel_mean=file2read.variables['dynstat_vvel_mean']
        vvel_mean=vvel_mean[:]*1
        vvel_max=file2read.variables['dynstat_vvel_max']
        vvel_max=vvel_max[:]*1
        vvel_min=file2read.variables['dynstat_vvel_min']
        vvel_min=vvel_min[:]*1    
        ke_mean=file2read.variables['ke_mean']
        ke_mean=ke_mean[:]*1
        ke_max=file2read.variables['ke_max']
        ke_max=ke_max[:]*1
        ke_vol=file2read.variables['ke_vol']
        ke_vol=ke_vol[:]*1
        
        time_seconds_tot =np.concatenate([time_seconds_tot , time_seconds])
        theta_mean_tot =np.concatenate([theta_mean_tot , theta_mean])
        theta_min_tot = np.concatenate([theta_min_tot , theta_min])
        theta_max_tot =np.concatenate([theta_max_tot , theta_max])    
        salt_mean_tot =np.concatenate([salt_mean_tot , salt_mean])
        salt_min_tot =np.concatenate([salt_min_tot , salt_min])
        salt_max_tot =np.concatenate([salt_max_tot , salt_max])
        sst_mean_tot =np.concatenate([sst_mean_tot , sst_mean])
        sst_min_tot =np.concatenate([sst_min_tot , sst_min])
        sst_max_tot =np.concatenate([sst_max_tot , sst_max])    
        sss_mean_tot =np.concatenate([sss_mean_tot , sss_mean])
        sss_min_tot =np.concatenate([sss_min_tot , sss_min])
        sss_max_tot =np.concatenate([sss_max_tot , sss_max])    
        vvel_mean_tot =np.concatenate([vvel_mean_tot , vvel_mean])
        vvel_min_tot =np.concatenate([vvel_min_tot , vvel_min])
        vvel_max_tot =np.concatenate([vvel_max_tot , vvel_max])    
        uvel_mean_tot =np.concatenate([uvel_mean_tot , uvel_mean])
        uvel_min_tot =np.concatenate([uvel_min_tot , uvel_min])
        uvel_max_tot =np.concatenate([uvel_max_tot , uvel_max])
        eta_mean_tot =np.concatenate([eta_mean_tot , eta_mean])
        eta_min_tot = np.concatenate([eta_min_tot , eta_min])
        eta_max_tot =np.concatenate([eta_max_tot , eta_max])
        ke_mean_tot =np.concatenate([ke_mean_tot , ke_mean])
        ke_vol_tot = np.concatenate([ke_vol_tot , ke_vol])
        ke_max_tot =np.concatenate([ke_max_tot , ke_max])
    
    #averages
    ave_theta_mean = np.zeros(32)
    ave_sss_mean = np.zeros(32)
    ave_sst_mean = np.zeros(32)
    for i in range(32):
        ave_theta_mean[i] = np.mean(theta_mean_tot[i*36+0:i*36+36])
        ave_sss_mean[i] = np.mean(sss_mean_tot[i*36+0:i*36+36])
        ave_sst_mean[i] = np.mean(sst_mean_tot[i*36+0:i*36+36])
        
    return     theta_mean_tot, theta_max_tot, theta_min_tot, eta_mean_tot, eta_max_tot, eta_min_tot, \
                salt_mean_tot, salt_max_tot, salt_min_tot, sss_mean_tot, sss_max_tot, sss_min_tot, \
                sst_mean_tot, sst_max_tot, sst_min_tot, vvel_mean_tot, vvel_max_tot, vvel_min_tot, \
                uvel_mean_tot, uvel_max_tot, uvel_min_tot, ke_mean_tot, ke_max_tot, ke_vol_tot, \
                time_seconds_tot


def monitor_seaice(x,iter_list):

    seaice_area_max_tot = []
    seaice_area_min_tot = []
    seaice_area_mean_tot = []
    seaice_heff_max_tot = []
    seaice_heff_min_tot = []
    seaice_heff_mean_tot = []
    time_seconds_ice_tot = []
    
    for iter in iter_list:
        zerosto_add = '0'
        for zeros in range(9 - len(str(iter))):
            zerosto_add = '0'+zerosto_add
            
        file2read3 = netcdf.NetCDFFile(x+"monitor_seaice."+zerosto_add+str(iter)+".t001.nc",'r')
        # Reading Seaice
        time_seconds_ice = file2read3.variables['seaice_time_sec']
        time_seconds_ice =time_seconds_ice[:]*1            
        seaice_area_max = file2read3.variables['seaice_area_max']
        seaice_area_max = seaice_area_max[:]*1
        seaice_area_min = file2read3.variables['seaice_area_min']
        seaice_area_min = seaice_area_min[:]*1
        seaice_area_mean = file2read3.variables['seaice_area_mean']
        seaice_area_mean = seaice_area_mean[:]*1
        seaice_heff_max = file2read3.variables['seaice_heff_max']
        seaice_heff_max = seaice_heff_max[:]*1
        seaice_heff_min = file2read3.variables['seaice_heff_min']
        seaice_heff_min = seaice_heff_min[:]*1
        seaice_heff_mean = file2read3.variables['seaice_heff_mean']
        seaice_heff_mean = seaice_heff_mean[:]*1
        
        time_seconds_ice_tot = np.concatenate([time_seconds_ice_tot , time_seconds_ice])
        seaice_area_max_tot = np.concatenate([seaice_area_max_tot , seaice_area_max])
        seaice_area_min_tot = np.concatenate([seaice_area_min_tot , seaice_area_min])
        seaice_area_mean_tot = np.concatenate([seaice_area_mean_tot , seaice_area_mean])    
        seaice_heff_max_tot = np.concatenate([seaice_heff_max_tot , seaice_heff_max])
        seaice_heff_min_tot = np.concatenate([seaice_heff_min_tot , seaice_heff_min])
        seaice_heff_mean_tot =np.concatenate([seaice_heff_mean_tot , seaice_heff_mean])        

    return      seaice_area_max_tot, seaice_area_min_tot, seaice_area_mean_tot, seaice_heff_max_tot, \
                seaice_heff_min_tot, seaice_heff_mean_tot, time_seconds_ice_tot


def dynStDiag(x,iter_list):
    theta_mean_tot =np.zeros([1,7,1])
    theta_max_tot =np.zeros([1,7,1])
    theta_min_tot =np.zeros([1,7,1])
    eta_mean_tot =np.zeros([1,7,1])
    eta_max_tot =np.zeros([1,7,1])
    eta_min_tot =np.zeros([1,7,1])
    salt_mean_tot =np.zeros([1,7,1])
    salt_max_tot =np.zeros([1,7,1])
    salt_min_tot =np.zeros([1,7,1])
    sss_mean_tot = [] 
    sss_max_tot =[]
    sss_min_tot =[]
    sst_mean_tot =[]
    sst_max_tot =[]
    sst_min_tot =[]
    vvel_mean_tot =np.zeros([1,7,1])
    vvel_max_tot =np.zeros([1,7,1])
    vvel_min_tot =np.zeros([1,7,1])
    uvel_mean_tot =np.zeros([1,7,1])
    uvel_max_tot =np.zeros([1,7,1])
    uvel_min_tot =np.zeros([1,7,1])
    ke_mean_tot =np.zeros([1,7,1])
    ke_max_tot =np.zeros([1,7,1])
    ke_vol_tot =np.zeros([1,7,1])
    time_seconds_tot = []
    
    # This script is meant to read dynStDiag files
    theta_lv_mean_tot =np.zeros([1,7,50])
    theta_lv_max_tot =np.zeros([1,7,50])
    theta_lv_min_tot =np.zeros([1,7,50])
    salt_lv_mean_tot =np.zeros([1,7,50])
    salt_lv_max_tot =np.zeros([1,7,50])
    salt_lv_min_tot =np.zeros([1,7,50])
    vvel_lv_mean_tot =np.zeros([1,7,50])
    vvel_lv_max_tot =np.zeros([1,7,50])
    vvel_lv_min_tot =np.zeros([1,7,50])
    uvel_lv_mean_tot =np.zeros([1,7,50])
    uvel_lv_max_tot =np.zeros([1,7,50])
    uvel_lv_min_tot =np.zeros([1,7,50])
    ke_lv_mean_tot =np.zeros([1,7,50])
    ke_lv_max_tot =np.zeros([1,7,50])
    time_lv_tot = []
    
    for iter in iter_list:
        zerosto_add = '0'
        for zeros in range(9 - len(str(iter))):
            zerosto_add = '0'+zerosto_add
    
        file2read = netcdf.NetCDFFile(x+"dynStDiag."+zerosto_add+str(iter)+".t001.nc",'r')
        # 2D fields
        time_lv = file2read.variables['T']
        time_lv = time_lv[:]*1
        theta_mean_lv=file2read.variables['THETA_lv_ave']
        theta_mean_lv=theta_mean_lv[:]*1    
        theta_max_lv=file2read.variables['THETA_lv_max']
        theta_max_lv=theta_max_lv[:]*1
        theta_min_lv=file2read.variables['THETA_lv_min']
        theta_min_lv=theta_min_lv[:]*1 
        salt_mean_lv=file2read.variables['SALT_lv_ave']
        salt_mean_lv=salt_mean_lv[:]*1
        salt_max_lv=file2read.variables['SALT_lv_max']
        salt_max_lv=salt_max_lv[:]*1
        salt_min_lv=file2read.variables['SALT_lv_min']
        salt_min_lv=salt_min_lv[:]*1
        uvel_mean_lv=file2read.variables['UVEL_lv_ave']
        uvel_mean_lv=uvel_mean_lv[:]*1
        uvel_max_lv=file2read.variables['UVEL_lv_max']
        uvel_max_lv=uvel_max_lv[:]*1
        uvel_min_lv=file2read.variables['UVEL_lv_min']
        uvel_min_lv=uvel_min_lv[:]*1
        vvel_mean_lv=file2read.variables['VVEL_lv_ave']
        vvel_mean_lv=vvel_mean_lv[:]*1
        vvel_max_lv=file2read.variables['VVEL_lv_max']
        vvel_max_lv=vvel_max_lv[:]*1
        vvel_min_lv=file2read.variables['VVEL_lv_min']
        vvel_min_lv=vvel_min_lv[:]*1
        ke_mean_lv=file2read.variables['momKE_lv_ave']
        ke_mean_lv=ke_mean_lv[:]*1
        ke_max_lv=file2read.variables['momKE_lv_max']
        ke_max_lv=ke_max_lv[:]*1
        
        theta_lv_mean_tot =np.concatenate([theta_lv_mean_tot , theta_mean_lv],axis=0)
        theta_lv_max_tot =np.concatenate([theta_lv_max_tot , theta_max_lv])
        theta_lv_min_tot =np.concatenate([theta_lv_min_tot , theta_min_lv])
        salt_lv_mean_tot =np.concatenate([salt_lv_mean_tot , salt_mean_lv])
        salt_lv_max_tot =np.concatenate([salt_lv_max_tot , salt_max_lv])
        salt_lv_min_tot =np.concatenate([salt_lv_min_tot , salt_min_lv])
        vvel_lv_mean_tot =np.concatenate([vvel_lv_mean_tot , vvel_mean_lv])
        vvel_lv_max_tot =np.concatenate([vvel_lv_max_tot , vvel_max_lv])
        vvel_lv_min_tot =np.concatenate([vvel_lv_min_tot , vvel_min_lv])
        uvel_lv_mean_tot =np.concatenate([uvel_lv_mean_tot , uvel_mean_lv])
        uvel_lv_max_tot =np.concatenate([uvel_lv_max_tot , uvel_max_lv])
        uvel_lv_min_tot =np.concatenate([uvel_lv_min_tot , uvel_min_lv])
        ke_lv_mean_tot =np.concatenate([ke_lv_mean_tot , ke_mean_lv])
        ke_lv_max_tot =np.concatenate([ke_lv_max_tot , ke_max_lv])
        time_lv_tot = np.concatenate([time_lv_tot , time_lv])
        
        # 1D fields
        time_seconds = file2read.variables['T']
        time_seconds=time_seconds[:]*1    
        theta_mean=file2read.variables['THETA_ave']
        theta_mean=theta_mean[:]*1    
        theta_max=file2read.variables['THETA_max']
        theta_max=theta_max[:]*1
        theta_min=file2read.variables['THETA_min']
        theta_min=theta_min[:]*1    
        salt_mean=file2read.variables['SALT_ave']
        salt_mean=salt_mean[:]*1
        salt_max=file2read.variables['SALT_max']
        salt_max=salt_max[:]*1
        salt_min=file2read.variables['SALT_min']
        salt_min=salt_min[:]*1    
        #sst_mean=file2read.variables['dynstat_sst_mean']
        #sst_mean=sst_mean[:]*1
        #sst_max=file2read.variables['dynstat_sst_max']
        #sst_max=sst_max[:]*1
        #sst_min=file2read.variables['dynstat_sst_max']
        #sst_min=sst_min[:]*1    
        #sss_mean=file2read.variables['dynstat_sss_mean']
        #sss_mean=sss_mean[:]*1
        #sss_mim=file2read.variables['dynstat_sss_min']
        #sss_min=sss_mean[:]*1
        #sss_max=file2read.variables['dynstat_sss_max']
        #sss_max=sss_max[:]*1    
        eta_mean=file2read.variables['ETAN_ave']
        eta_mean=eta_mean[:]*1
        eta_min=file2read.variables['ETAN_min']
        eta_min=eta_min[:]*1
        eta_max=file2read.variables['ETAN_max']
        eta_max=eta_max[:]*1    
        uvel_mean=file2read.variables['UVEL_ave']
        uvel_mean=uvel_mean[:]*1
        uvel_max=file2read.variables['UVEL_max']
        uvel_max=uvel_max[:]*1
        uvel_min=file2read.variables['UVEL_min']
        uvel_min=uvel_min[:]*1    
        vvel_mean=file2read.variables['VVEL_ave']
        vvel_mean=vvel_mean[:]*1
        vvel_max=file2read.variables['VVEL_max']
        vvel_max=vvel_max[:]*1
        vvel_min=file2read.variables['VVEL_min']
        vvel_min=vvel_min[:]*1    
        ke_mean=file2read.variables['momKE_ave']
        ke_mean=ke_mean[:]*1
        ke_max=file2read.variables['momKE_max']
        ke_max=ke_max[:]*1
        ke_vol=file2read.variables['momKE_vol']
        ke_vol=ke_vol[:]*1
        
        time_seconds_tot =np.concatenate([time_seconds_tot , time_seconds])
        theta_mean_tot =np.concatenate([theta_mean_tot , theta_mean])
        theta_min_tot = np.concatenate([theta_min_tot , theta_min])
        theta_max_tot =np.concatenate([theta_max_tot , theta_max])    
        salt_mean_tot =np.concatenate([salt_mean_tot , salt_mean])
        salt_min_tot =np.concatenate([salt_min_tot , salt_min])
        salt_max_tot =np.concatenate([salt_max_tot , salt_max])
        #sst_mean_tot =np.concatenate([sst_mean_tot , sst_mean])
        #sst_min_tot =np.concatenate([sst_min_tot , sst_min])
        #sst_max_tot =np.concatenate([sst_max_tot , sst_max])    
        #sss_mean_tot =np.concatenate([sss_mean_tot , sss_mean])
        #sss_min_tot =np.concatenate([sss_min_tot , sss_min])
        #sss_max_tot =np.concatenate([sss_max_tot , sss_max])    
        vvel_mean_tot =np.concatenate([vvel_mean_tot , vvel_mean])
        vvel_min_tot =np.concatenate([vvel_min_tot , vvel_min])
        vvel_max_tot =np.concatenate([vvel_max_tot , vvel_max])    
        uvel_mean_tot =np.concatenate([uvel_mean_tot , uvel_mean])
        uvel_min_tot =np.concatenate([uvel_min_tot , uvel_min])
        uvel_max_tot =np.concatenate([uvel_max_tot , uvel_max])
        eta_mean_tot =np.concatenate([eta_mean_tot , eta_mean])
        eta_min_tot = np.concatenate([eta_min_tot , eta_min])
        eta_max_tot =np.concatenate([eta_max_tot , eta_max])
        ke_mean_tot =np.concatenate([ke_mean_tot , ke_mean])
        ke_vol_tot = np.concatenate([ke_vol_tot , ke_vol])
        ke_max_tot =np.concatenate([ke_max_tot , ke_max])

    
    return     theta_mean_tot, theta_max_tot, theta_min_tot, eta_mean_tot, eta_max_tot, eta_min_tot, \
                salt_mean_tot, salt_max_tot, salt_min_tot, sss_mean_tot, sss_max_tot, sss_min_tot, \
                sst_mean_tot, sst_max_tot, sst_min_tot, vvel_mean_tot, vvel_max_tot, vvel_min_tot, \
                uvel_mean_tot, uvel_max_tot, uvel_min_tot, ke_mean_tot, ke_max_tot, ke_vol_tot, \
                time_seconds_tot,\
                theta_lv_mean_tot, theta_lv_max_tot, theta_lv_min_tot, \
                salt_lv_mean_tot, salt_lv_max_tot, salt_lv_min_tot, \
                vvel_lv_mean_tot, vvel_lv_max_tot, vvel_lv_min_tot, \
                uvel_lv_mean_tot, uvel_lv_max_tot, uvel_lv_min_tot, \
                ke_lv_mean_tot, ke_lv_max_tot, time_lv_tot

def titles():
    titles = {'theta_mean' : 'Temperature Mean', 'theta_min' : 'Temperature Min', 'theta_max' :  'Temperature Max', \
            'eta_mean' : 'ETA mean', 'eta_max' : 'ETA max', 'eta_min' : 'ETA min', 'salt_mean' : 'Salinity Mean' , \
            'salt_max' : 'Salinity Max' , 'salt_min' : 'Salinity Min' , 'sss_mean': 'SSS mean' , \
            'sss_max' : 'SSS max', 'sss_min' : 'SSS min', 'sst_mean' : 'SST mean', 'sst_max' : 'SST max' , \
            'sst_min' : 'SST min', 'vvel_mean' : 'V mean', 'vvel_max' : 'V max' , 'vvel_min' : 'V min', \
            'uvel_mean' : 'U mean', 'uvel_max' : 'U max', 'uvel_min' : 'U min', 'ke_mean' : 'Kinetic mean', \
            'ke_max' : 'Kinetic max', 'ke_vol' : 'Kinetic volume', 'seaice_area_max' : 'Seaice area max', \
            'seaice_area_min' : 'Seaice area min', 'seaice_area_mean' : 'Seaice area mean', \
            'seaice_heff_max' : 'Seaice thickness max', 'seaice_heff_min' : 'Seaicea thickness max', \
            'seaice_heff_mean' : 'Seaice thickness mean', 'time_seconds' : 'Time seconds' , \
          'time_years' : 'Time years'} 
    unity = {'theta_mean' : 'C', 'theta_min' : 'C', 'theta_max' :  'C', \
            'eta_mean' : 'm', 'eta_max' : 'm', 'eta_min' : 'm', 'salt_mean' : 'psu' , \
            'salt_max' : 'psu' , 'salt_min' : 'psu' , 'sss_mean': 'psu' , \
            'sss_max' : 'psu', 'sss_min' : 'psu', 'sst_mean' : 'C', 'sst_max' : 'C' , \
            'sst_min' : 'C', 'vvel_mean' : 'm/s', 'vvel_max' : 'm/s' , 'vvel_min' : 'm/s', \
            'uvel_mean' : 'm/s', 'uvel_max' : 'm/s', 'uvel_min' : 'm/s', 'ke_mean' : 'm^2/s^2', \
            'ke_max' : 'm^2/s^2', 'ke_vol' : 'm^2/s^2', 'seaice_area_max' : '%', \
            'seaice_area_min' : '%', 'seaice_area_mean' : '%', \
            'seaice_heff_max' : 'm', 'seaice_heff_min' : 'm', \
            'seaice_heff_mean' : 'm', 'time_seconds' : 's' , 'time_seconds_ice' : 's', 'time_years' : 'Years',\
          'time_years_ice' : 'Years'} 
    time = {'theta_mean' : 'time_years', 'theta_min' : 'time_years', 'theta_max' :  'time_years', \
        'eta_mean' : 'time_years', 'eta_max' : 'time_years', 'eta_min' : 'time_years', \
        'salt_mean' : 'time_years' , 'salt_max' : 'time_years' , 'salt_min' : 'time_years' , \
        'sss_mean': 'time_years' , 'sss_max' : 'time_years', 'sss_min' : 'time_years', \
        'sst_mean' : 'time_years', 'sst_max' : 'time_years' ,'sst_min' : 'time_years', \
        'vvel_mean' : 'time_years', 'vvel_max' : 'time_years' , 'vvel_min' : 'time_years', \
        'uvel_mean' : 'time_years', 'uvel_max' : 'time_years', 'uvel_min' : 'time_years', \
        'ke_mean' : 'time_years', 'ke_max' : 'time_years', 'ke_vol' : 'time_years', \
        'seaice_area_max' : 'time_years_ice', 'seaice_area_min' : 'time_years_ice', \
        'seaice_area_mean' : 'time_years_ice', 'seaice_heff_max' : 'time_years_ice', \
        'seaice_heff_min' : 'time_years_ice', 'seaice_heff_mean' : 'time_years_ice'} 
    return titles,unity,time

class DynStRead():
    def __init__(self):
        self.data = {'theta_mean' : [], 'theta_min' : [], 'theta_max' : [],\
                     'eta_mean' : [], 'eta_max' : [], 'eta_min' : [],\
                     'salt_mean' : [] , 'salt_max' : [] , 'salt_min' : [],\
                     'sss_mean': [] , 'sss_max' : [], 'sss_min' : [],\
                     'sst_mean' : [], 'sst_max' : [], 'sst_min' : [],\
                     'vvel_mean' : [], 'vvel_max' : [], 'vvel_min' : [],\
                     'uvel_mean' : [], 'uvel_max' : [], 'uvel_min' : [],\
                     'ke_mean' : [], 'ke_max' : [], 'ke_vol' : [],\
                     'seaice_area_max' : [], 'seaice_area_min' : [], 'seaice_area_mean' : [],\
                     'seaice_heff_max' : [], 'seaice_heff_min' : [], 'seaice_heff_mean' : [],\
                     'time_seconds' : [] , 'time_years' : [] ,\
                     'time_seconds_ice' : [],'time_years_ice' : []}
        self.dataDyn =  {'theta_mean' : [], 'theta_min' : [], 'theta_max' : [] ,\
                         'eta_mean' : [], 'eta_max' : [], 'eta_min' : [], \
                         'salt_mean' : [] , 'salt_max' : [] , 'salt_min' : [] ,\
                         'sss_mean': [] , 'sss_max' : [], 'sss_min' : [],\
                         'sst_mean' : [], 'sst_max' : [], 'sst_min' : [], \
                         'vvel_mean' : [], 'vvel_max' : [], 'vvel_min' : [],\
                         'uvel_mean' : [], 'uvel_max' : [], 'uvel_min' : [],\
                         'ke_mean' : [], 'ke_max' : [], 'ke_vol' : [],\
                         'theta_lv_mean' : [], 'theta_lv_min' : [], 'theta_lv_max' : [] , \
                         'salt_lv_mean' : [] , 'salt_lv_max' : [] , 'salt_lv_min' : [] , \
                         'vvel_lv_mean' : [], 'vvel_lv_max' : [], 'vvel_lv_min' : [], \
                         'uvel_lv_mean' : [], 'uvel_lv_max' : [], 'uvel_lv_min' : [], \
                         'ke_lv_mean' : [], 'ke_lv_max' : [], \
                         'time_lv' : [] , 'time_lv_years' : [] }
        
    def readData(self,path,iters):
        self.data['theta_mean'], self.data['theta_max'], self.data['theta_min'], self.data['eta_mean'],\
        self.data['eta_max'], self.data['eta_min'], self.data['salt_mean'], self.data['salt_max'],\
        self.data['salt_min'], self.data['sss_mean'], self.data['sss_max'], self.data['sss_min'],\
        self.data['sst_mean'], self.data['sst_max'], self.data['sst_min'], self.data['vvel_mean'],\
        self.data['vvel_max'], self.data['vvel_min'], self.data['uvel_mean'], self.data['uvel_max'],\
        self.data['uvel_min'], self.data['ke_mean'], self.data['ke_max'], self.data['ke_vol'], \
        self.data['time_seconds'] = monitor(path,iters)
        self.data['time_years'] = (self.data['time_seconds']- self.data['time_seconds'][0])/(360*60*60*24)

    def readseaiceData(self,path,iters):
        self.data['seaice_area_max'], self.data['seaice_area_min'], self.data['seaice_area_mean'], \
        self.data['seaice_heff_max'], self.data['seaice_heff_min'], self.data['seaice_heff_mean'], \
        self.data['time_seconds_ice'] = monitor_seaice(path,iters)
        self.data['time_years_ice'] = (self.data['time_seconds_ice']- self.data['time_seconds_ice'][0])/(360*60*60*24)

    def readDynStDiag(self,path,iters):        
        self.dataDyn['theta_mean'], self.dataDyn['theta_max'], self.dataDyn['theta_min'],\
        self.dataDyn['eta_mean'], self.dataDyn['eta_max'], self.dataDyn['eta_min'],\
        self.dataDyn['salt_mean'], self.dataDyn['salt_max'], self.dataDyn['salt_min'],\
        self.dataDyn['sss_mean'], self.dataDyn['sss_max'], self.dataDyn['sss_min'],\
        self.dataDyn['sst_mean'], self.dataDyn['sst_max'], self.dataDyn['sst_min'],\
        self.dataDyn['vvel_mean'], self.dataDyn['vvel_max'], self.dataDyn['vvel_min'],\
        self.dataDyn['uvel_mean'], self.dataDyn['uvel_max'], self.dataDyn['uvel_min'],\
        self.dataDyn['ke_mean'], self.dataDyn['ke_max'], self.dataDyn['ke_vol'], self.dataDyn['time_seconds'], \
        self.dataDyn['theta_lv_mean'], self.dataDyn['theta_lv_max'], self.dataDyn['theta_lv_min'], \
        self.dataDyn['salt_lv_mean'], self.dataDyn['salt_lv_max'], self.dataDyn['salt_lv_min'], \
        self.dataDyn['vvel_lv_mean'], self.dataDyn['vvel_lv_max'], self.dataDyn['vvel_lv_min'], \
        self.dataDyn['uvel_lv_mean'], self.dataDyn['uvel_lv_max'], self.dataDyn['uvel_lv_min'], \
        self.dataDyn['ke_lv_mean'], self.dataDyn['ke_lv_max'],\
        self.dataDyn['time_lv'] = dynStDiag(path,iters)    
        self.dataDyn['rho_lv_mean'] = rhop(self.dataDyn['salt_lv_mean'],self.dataDyn['theta_lv_mean'])
        self.dataDyn['time_lv_years'] = (self.dataDyn['time_lv']- self.dataDyn['time_lv'][0])/(360*60*60*24)
        self.data['time_years'] = (self.data['time_seconds']- self.data['time_seconds'][0])/(360*60*60*24)
     
    def title(self,title,color):
        self.title = title
        self.color = color
