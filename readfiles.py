#! /usr/bin/env python                                                                                        # -*- coding: utf-8 -*-                                                                                                    

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

def monitor_extract(x,iter_list):
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
    seaice_area_max_tot = []
    seaice_area_min_tot = []
    seaice_area_mean_tot = []
    seaice_heff_max_tot = []
    seaice_heff_min_tot = []
    seaice_heff_mean_tot = []
    time_seconds_tot = []
    
    for iter in iter_list:
        zerosto_add = '0'
        for zeros in range(9 - len(str(iter))):
            zerosto_add = '0'+zerosto_add
        
        file2read = netcdf.NetCDFFile(x+"monitor."+zerosto_add+str(iter)+".t001.nc",'r')
        file2read3 = netcdf.NetCDFFile(x+"monitor_seaice."+zerosto_add+str(iter)+".t001.nc",'r')

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
    
        seaice_area_max_tot =np.concatenate([seaice_area_max_tot , seaice_area_max])
        seaice_area_min_tot = np.concatenate([seaice_area_min_tot , seaice_area_min])
        seaice_area_mean_tot =np.concatenate([seaice_area_mean_tot , seaice_area_mean])    
        seaice_heff_max_tot =np.concatenate([seaice_heff_max_tot , seaice_heff_max])
        seaice_heff_min_tot = np.concatenate([seaice_heff_min_tot , seaice_heff_min])
        seaice_heff_mean_tot =np.concatenate([seaice_heff_mean_tot , seaice_heff_mean])        
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
                seaice_area_max_tot, seaice_area_min_tot, seaice_area_mean_tot, seaice_heff_max_tot, \
                seaice_heff_min_tot, seaice_heff_mean_tot, time_seconds_tot

def dynStDiag_extract(x,iter_list):
    # This script is meant to read dynStDiag files
    os.chdir(x)

    theta_lv_mean_tot =np.zeros([1,1,50])
    theta_lv_max_tot =np.zeros([1,1,50])
    theta_lv_min_tot =np.zeros([1,1,50])
    salt_lv_mean_tot =np.zeros([1,1,50])
    salt_lv_max_tot =np.zeros([1,1,50])
    salt_lv_min_tot =np.zeros([1,1,50])
    vvel_lv_mean_tot =np.zeros([1,1,50])
    vvel_lv_max_tot =np.zeros([1,1,50])
    vvel_lv_min_tot =np.zeros([1,1,50])
    uvel_lv_mean_tot =np.zeros([1,1,50])
    uvel_lv_max_tot =np.zeros([1,1,50])
    uvel_lv_min_tot =np.zeros([1,1,50])
    ke_lv_mean_tot =np.zeros([1,1,50])
    ke_lv_max_tot =np.zeros([1,1,50])
    time_lv_tot = []

    for iter in iter_list:
        zerosto_add = '0'
        for zeros in range(9 - len(str(iter))):
            zerosto_add = '0'+zerosto_add
    
        file2read2 = netcdf.NetCDFFile("dynStDiag."+zerosto_add+str(iter)+".t001.nc",'r')
    
        # Dyn Stat
        time_lv = file2read2.variables['T']
        time_lv = time_lv[:]*1
    
        theta_mean_lv=file2read2.variables['THETA_lv_ave']
        theta_mean_lv=theta_mean_lv[:]*1    
        theta_max_lv=file2read2.variables['THETA_lv_max']
        theta_max_lv=theta_max_lv[:]*1
        theta_min_lv=file2read2.variables['THETA_lv_min']
        theta_min_lv=theta_min_lv[:]*1 
        salt_mean_lv=file2read2.variables['SALT_lv_ave']
        salt_mean_lv=salt_mean_lv[:]*1
        salt_max_lv=file2read2.variables['SALT_lv_max']
        salt_max_lv=salt_max_lv[:]*1
        salt_min_lv=file2read2.variables['SALT_lv_min']
        salt_min_lv=salt_min_lv[:]*1
        uvel_mean_lv=file2read2.variables['UVEL_lv_ave']
        uvel_mean_lv=uvel_mean_lv[:]*1
        uvel_max_lv=file2read2.variables['UVEL_lv_max']
        uvel_max_lv=uvel_max_lv[:]*1
        uvel_min_lv=file2read2.variables['UVEL_lv_min']
        uvel_min_lv=uvel_min_lv[:]*1
        vvel_mean_lv=file2read2.variables['VVEL_lv_ave']
        vvel_mean_lv=vvel_mean_lv[:]*1
        vvel_max_lv=file2read2.variables['VVEL_lv_max']
        vvel_max_lv=vvel_max_lv[:]*1
        vvel_min_lv=file2read2.variables['VVEL_lv_min']
        vvel_min_lv=vvel_min_lv[:]*1
        ke_mean_lv=file2read2.variables['momKE_lv_ave']
        ke_mean_lv=ke_mean_lv[:]*1
        ke_max_lv=file2read2.variables['momKE_lv_max']
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
    
    theta_lv_mean = theta_lv_mean_tot[1:,:,:].squeeze(axis=1)
    theta_lv_max = theta_lv_max_tot[1:,:,:].squeeze(axis=1)
    theta_lv_min = theta_lv_min_tot[1:,:,:].squeeze(axis=1)
    salt_lv_mean = salt_lv_mean_tot[1:,:,:].squeeze(axis=1)
    salt_lv_max = salt_lv_max_tot[1:,:,:].squeeze(axis=1)
    salt_lv_min = salt_lv_min_tot[1:,:,:].squeeze(axis=1)
    vvel_lv_mean = vvel_lv_mean_tot[1:,:,:].squeeze(axis=1)
    vvel_lv_max = vvel_lv_max_tot[1:,:,:].squeeze(axis=1)
    vvel_lv_min = vvel_lv_min_tot[1:,:,:].squeeze(axis=1)
    uvel_lv_mean = uvel_lv_mean_tot[1:,:,:].squeeze(axis=1)
    uvel_lv_max = uvel_lv_max_tot[1:,:,:].squeeze(axis=1)
    uvel_lv_min = uvel_lv_min_tot[1:,:,:].squeeze(axis=1)
    ke_lv_mean = ke_lv_mean_tot[1:,:,:].squeeze(axis=1)
    ke_lv_max = ke_lv_max_tot[1:,:,:].squeeze(axis=1)
    
    return     theta_lv_mean, theta_lv_max, theta_lv_min, \
                salt_lv_mean, salt_lv_max, salt_lv_min, \
                vvel_lv_mean, vvel_lv_max, vvel_lv_min, \
                uvel_lv_mean, uvel_lv_max, uvel_lv_min, \
                ke_lv_mean, ke_lv_max, time_lv_tot

class MonitorRead():
    def __init__(self):
        self.data = {'theta_mean' : [], 'theta_min' : [], 'theta_max' : [] , 'eta_mean' : [], 'eta_max' : [], \
                     'eta_min' : [], 'salt_mean' : [] , 'salt_max' : [] , 'salt_min' : [] , 'sss_mean': [] , \
                     'sss_max' : [], 'sss_min' : [], 'sst_mean' : [], 'sst_max' : [], 'sst_min' : [], \
                     'vvel_mean' : [], 'vvel_max' : [], 'vvel_min' : [], 'uvel_mean' : [], 'uvel_max' : [], \
                     'uvel_min' : [], 'ke_mean' : [], 'ke_max' : [], 'ke_vol' : [], 'seaice_area_max' : [], \
                     'seaice_area_min' : [], 'seaice_area_mean' : [], 'seaice_heff_max' : [], 'seaice_heff_min' : [], \
                     'seaice_heff_mean' : [], 'time_seconds' : [] , 'time_years' : [] }
        self.dataDyn =  {'theta_lv_mean' : [], 'theta_lv_min' : [], 'theta_lv_max' : [] , \
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
        self.data['seaice_area_max'], self.data['seaice_area_min'], self.data['seaice_area_mean'], \
        self.data['seaice_heff_max'], self.data['seaice_heff_min'], self.data['seaice_heff_mean'], \
        self.data['time_seconds'] = monitor_extract(path,iters)
        self.data['time_years'] = (self.data['time_seconds']- self.data['time_seconds'][0])/(360*60*60*24)
    
    def readDynStDiag(self,path,iters):
        self.dataDyn['theta_lv_mean'], self.dataDyn['theta_lv_max'], self.dataDyn['theta_lv_min'], \
        self.dataDyn['salt_lv_mean'], self.dataDyn['salt_lv_max'], self.dataDyn['salt_lv_min'], \
        self.dataDyn['vvel_lv_mean'], self.dataDyn['vvel_lv_max'], self.dataDyn['vvel_lv_min'], \
        self.dataDyn['uvel_lv_mean'], self.dataDyn['uvel_lv_max'], self.dataDyn['uvel_lv_min'], \
        self.dataDyn['ke_lv_mean'], self.dataDyn['ke_lv_max'],\
        self.dataDyn['time_lv'] = dynStDiag_extract(path,iters)
        self.dataDyn['time_lv_years'] = (self.dataDyn['time_lv']- self.dataDyn['time_lv'][0])/(360*60*60*24)
    
    def title(self,title,color):
        self.title = title
        self.color = color

class StateRead:
    def __init__(self):
        self.data = {'T' : [], 'V' : [], 'U' : [] , 'S' : [], 'days' : [], 'years' : []}
        self.psi = []
        self.psi_mean = []
        self.psi_max = []
        self.psi_min = []
        self.years = []
        self.res = []
        self.grid = []

    def readData(self,path,list_var):
        file2read = netcdf.NetCDFFile(path+'state.nc','r')
        Temp=file2read.variables['Temp']
        self.data['T']=Temp[list_var]*1
        V=file2read.variables['V']
        self.data['V']=V[list_var]*1
        U=file2read.variables['U']
        self.data['U']=U[list_var]*1
        S=file2read.variables['S']
        self.data['S']=S[list_var]*1
        days=file2read.variables['T']
        self.data['days']=days[list_var]*1
        self.years = (self.data['days'] - self.data['days'][0])/(60*60*24*360)
                                                                                    
        if self.data['T'].shape[3] == 210:
            self.grid = "/scratch/general/am8e13/results36km/grid.nc"
            self.res = 36
        elif self.data['T'].shape[3] == 420:
            self.grid = "/scratch/general/am8e13/results18km/grid.nc"
            self.res = 18
        elif self.data['T'].shape[3] == 840:
            self.grid = "/scratch/general/am8e13/results9km/grid.nc"
            self.res = 9
        file2read = netcdf.NetCDFFile(self.grid,'r')
        
        # Bathy is 1 on land and 0 over sea 
        hfacc = file2read.variables['HFacC']
        hfacc = hfacc[:]*1
        hfacw = file2read.variables['HFacW']
        hfacw = hfacw[:]*1
        hfacs = file2read.variables['HFacS']
        hfacs = hfacs[:]*1
        self.data['T'][:,hfacc==0] = np.nan
        self.data['U'][:,hfacw==0] = np.nan
        self.data['V'][:,hfacs==0] = np.nan
        self.data['S'][:,hfacc==0] = np.nan

    def fluxCalc(self):
        # This function calculates 
        #establish
        if self.res == 36:
            kk = 1
        elif self.res == 18:
            kk = 2
        elif self.res == 9:
            kk = 4
        self.grid
        file2read = netcdf.NetCDFFile(self.grid,'r')
        hfacc = file2read.variables['HFacC']
        hfacc = hfacc[:]*1
        drf = file2read.variables['drF']
        drf = drf[:]*1
        rA = file2read.variables['rA']
        rA = rA[:]*1
        dyF = file2read.variables['dyF']
        dyF = dyF[:]*1
        dxF = file2read.variables['dxF']
        dxF = dxF[:]*1
        dydx = np.zeros_like(hfacc)
        for k in range(len(drf)):
            dydx[k,:,:] = drf[k]*rA*hfacc[k,:,:]
        Area_x = dydx/dyF
        Area_y = dydx/dxF
        Area_x[hfacc==0]=np.nan
        Area_y[hfacc==0]=np.nan
        
        self.Fram = {'Flux' : np.zeros_like(self.data['V'][:,:,58*kk:80*kk,76*kk]) , \
                        'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                        'FluxS' : np.zeros_like(self.data['S'][:,:,58*kk:80*kk,76*kk]),\
                        'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                        'FluxT' : np.zeros_like(self.data['T'][:,:,58*kk:80*kk,76*kk]),\
                        'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                        'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                        'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                        'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
        
        self.Barents = {'Flux' : np.concatenate((self.data['V'][:,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                    self.data['U'][:,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=2),\
                        'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                        'FluxS' : np.concatenate((self.data['V'][:,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                    self.data['U'][:,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=2),\
                        'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                        'FluxT' : np.concatenate((self.data['V'][:,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                    self.data['U'][:,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=2),\
                        'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0])}
        
        for t in range(self.data['V'].shape[0]):
            self.Fram['Flux'][t,:,:] = self.data['V'][t,:,58*kk:80*kk,76*kk]*Area_x[:,58*kk:80*kk,76*kk]
            self.Fram['FluxSum'][t] = np.nansum(np.nansum(self.Fram['Flux'][t,:,:]))
            self.Fram['FluxT'][t,:,:] = self.Fram['Flux'][t,:,:]*self.data['T'][t,:,58*kk:80*kk,76*kk]*Area_x[:,58*kk:80*kk,76*kk]
            self.Fram['FluxSumT'][t] = np.nansum(np.nansum(self.Fram['FluxT'][t,:,:]))
            self.Fram['FluxS'][t,:,:] = self.Fram['Flux'][t,:,:]*self.data['S'][t,:,58*kk:80*kk,76*kk]*Area_x[:,58*kk:80*kk,76*kk]
            self.Fram['FluxSumS'][t] = np.nansum(np.nansum(self.Fram['FluxS'][t,:,:]))
            
            self.Barents['Flux'][t,:,:] = np.concatenate((self.data['V'][t,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                                      self.data['U'][t,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=1)
            self.Barents['FluxSum'][t] = np.nansum(np.nansum(self.Barents['Flux'][t,:,:]))
            self.Barents['FluxT'][t,:,:] = self.Barents['Flux'][t,:,:]*np.concatenate((self.data['T'][t,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                                      self.data['T'][t,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=1)
            self.Barents['FluxSumT'][t] = np.nansum(np.nansum(self.Barents['FluxT'][t,:,:]))
            self.Barents['FluxS'][t,:,:] = self.Barents['Flux'][t,:,:]*np.concatenate((self.data['T'][t,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                                      self.data['T'][t,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=1)
            self.Barents['FluxSumS'][t] = np.nansum(np.nansum(self.Barents['FluxS'][t,:,:]))

        n1 = 22
        n2 = 33
        self.Fram = {'FluxTop' : np.nansum(np.nansum(self.Fram['Flux'][:,0:n1,:],axis=2),axis=1),\
                    'FluxMid' : np.nansum(np.nansum(self.Fram['Flux'][:,n1:n2,:],axis=2),axis=1),\
                    'FluxBot' : np.nansum(np.nansum(self.Fram['Flux'][:,n2:50,:],axis=2),axis=1),\
                    'FluxTopS' : np.nansum(np.nansum(self.Fram['FluxS'][:,0:n1,:],axis=2),axis=1),\
                    'FluxMidS' : np.nansum(np.nansum(self.Fram['FluxS'][:,n1:n2,:],axis=2),axis=1),\
                    'FluxBotS' : np.nansum(np.nansum(self.Fram['FluxS'][:,n2:50,:],axis=2),axis=1),\
                    'FluxTopT' : np.nansum(np.nansum(self.Fram['FluxT'][:,0:n1,:],axis=2),axis=1),\
                    'FluxMidT' : np.nansum(np.nansum(self.Fram['FluxT'][:,n1:n2,:],axis=2),axis=1),\
                    'FluxBotT' : np.nansum(np.nansum(self.Fram['FluxT'][:,n2:50,:],axis=2),axis=1)}

    def baroCalc(self):
        self.psi = baro_stream(self.data['U'])
        self.psi_mean = np.nanmean(self.psi,axis = 1)
        self.psi_mean = np.nanmean(self.psi_mean,axis = 1)
        self.psi_max = np.nanmax(self.psi,axis=1)
        self.psi_max = np.nanmax(self.psi_max,axis=1)
        self.psi_min = np.nanmin(self.psi,axis=1)
        self.psi_min = np.nanmin(self.psi_min,axis=1)
        
    def topoCalc(self):
        topo,topo_opposit = Topostrophy(self.data['U'],self.data['V'])
        self.topo = topo
        self.topo_opposit = topo_opposit

    def title(self,title):
        self.title = title
    
    def print_title(self):
        print self.title
