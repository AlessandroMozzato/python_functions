#! /usr/bin/env python                                                                                                                     
# -*- coding: utf-8 -*-                                                                                                                   

from scipy.io import netcdf
import numpy as np
from regions_def import *
from gridread import *
from jmd95 import *
from rho import *

## Here we add a set of functions to calculate the overflow
# calculation of density
def density_calc(run):
    run.data['rho_jmd'] = np.zeros_like(run.data['S'])
    run.data['rhop'] = np.zeros_like(run.data['S'])
    for t in range(run.data['S'].shape[0]):
        for z in range(len(run.Z)):
            run.data['rho_jmd'][t,z,:,:] = densjmd95(run.data['S'][t,z,:,:],run.data['T'][t,z,:,:],-9.81*run.Z[z]*1025)
        run.data['rhop'][t,:,:,:] = rhop(run.data['S'][t,:,:,:],run.data['T'][t,:,:,:])

# calculation of fluxes including 
def fluxesCalculation(run):
    kdic = {36:1,18:2,9:4}
    kk = kdic[run.res]

    file2read = netcdf.NetCDFFile('/scratch/general/am8e13/results'+str(run.res)+'km/grid.nc','r')
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
    Area_x = dydx/dxF
    Area_y = dydx/dyF
    Area_x[hfacc==0]=np.nan
    Area_y[hfacc==0]=np.nan
    
    coords = {}
    coords = {'Fram': [58,80,76,76], 'Fram1' : [55,85,78,78] , 'Fram2' : [60,83,72,72] , 'Denmark': [100,100,37,48] , \
             'Bering' : [80,89,178,178] , 'Davis' : [113,135,75,75] , 'Davis1' : [135,135,52,73] , 'Davis2' : [113,135,75,75]}
    ax_d = {'Fram': 3, 'Fram1' : 3 , 'Fram2' : 3 , 'Denmark': 2 , 'Bering' : 3 , 'Davis' : 3 , 'Davis1' : 2 , 'Davis2' : 3}
    vel = {'Fram': 'U', 'Fram1' : 'V' , 'Fram2' : 'V' , 'Denmark': 'U' , \
             'Bering' : 'V' , 'Davis' : 'U' , 'Davis1' : 'U' , 'Davis2' : 'U' }
    
    def coord_calc(coord):
        if coord[0] == coord[1]:
            coord[0] = coord[0]*kk
            coord[1] = coord[0]+1
            coord[2] = coord[2]*kk
            coord[3] = coord[3]*kk
        elif coord[2] == coord[3]:
            coord[0] = coord[0]*kk
            coord[1] = coord[1]*kk
            coord[2] = coord[2]*kk
            coord[3] = coord[2]+1
        return coord

    run.fluxes2 = {}
    for var in ['Fram','Fram1','Fram2','Denmark','Bering','Davis','Davis1','Davis2']:
        coord = coords[var]
        coord = coord_calc(coord)
        run.fluxes2[var] = {'Flux' : np.zeros_like(run.data[vel[var]][:,:,coord[0]:coord[1],coord[2]:coord[3]]),\
                            'FluxSum' : np.zeros_like(run.data[vel[var]][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(run.data[vel[var]][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(run.data[vel[var]][:,0,0,0]),\
                            'FluxS' : np.zeros_like(run.data['S'][:,:,coord[0]:coord[1],coord[2]:coord[3]]),\
                            'FluxSumS' : np.zeros_like(run.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(run.data['T'][:,:,coord[0]:coord[1],coord[2]:coord[3]]),\
                            'FluxSumT' : np.zeros_like(run.data['T'][:,0,0,0]),\
                            'FluxFW' : np.zeros_like(run.data['T'][:,:,coord[0]:coord[1],coord[2]:coord[3]]),\
                            'FluxSumFW' : np.zeros_like(run.data['T'][:,0,0,0]),\
                            'FluxFW1' : np.zeros_like(run.data['T'][:,:,coord[0]:coord[1],coord[2]:coord[3]]),\
                            'FluxSumFW1' : np.zeros_like(run.data['T'][:,0,0,0]),\
                            'FluxOverFlow' : np.zeros_like(run.data['T'][:,:,coord[0]:coord[1],coord[2]:coord[3]]),\
                            'FluxSumOverFlow' : np.zeros_like(run.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
        
    S0 = 34.8 # reference salinity
    rho0 = 1027.8
    # this is to calculate the FW flux correctly
    tmp = np.ones_like(run.data['S'])
    tmp[:,31:,:,:] = 0
    tmpof = np.zeros_like(run.data['S'])
    tmpof[run.data['rhop']>rho0] = 1

    # this is to calculate the FW flux correctly
    tmp1 = np.ones_like(run.data['S'])
    tmp1[run.data['S']>S0] = 0
    for var in ['Fram','Fram1','Fram2','Denmark','Bering','Davis','Davis1','Davis2']:
        coord = coords[var]
        coord = coord_calc(coord)
        for t in range(run.data['V'].shape[0]):    
            # Fram fillign
            run.fluxes2[var]['Flux'][t,:,:] = run.data[vel[var]][t,:,coord[0]:coord[1],coord[2]:coord[3]]*\
                    Area_y[:,coord[0]:coord[1],coord[2]:coord[3]]/10**6
            run.fluxes2[var]['FluxSum'][t] = np.nansum(np.nansum(run.fluxes2[var]['Flux'][t,:,:]))
            run.fluxes2[var]['FluxInSum'][t] = np.nansum(np.nansum(run.fluxes2[var]['Flux'][t,run.fluxes2[var]['Flux'][t,:,:]>0]))
            run.fluxes2[var]['FluxOutSum'][t] = np.nansum(np.nansum(run.fluxes2[var]['Flux'][t,run.fluxes2[var]['Flux'][t,:,:]<0]))
            run.fluxes2[var]['FluxT'][t,:,:] = run.fluxes2[var]['Flux'][t,:,:]*\
                                                run.data['T'][t,:,coord[0]:coord[1],coord[2]:coord[3]]
            run.fluxes2[var]['FluxSumT'][t] = np.nansum(np.nansum(run.fluxes2[var]['FluxT'][t,:,:]))
            run.fluxes2[var]['FluxS'][t,:,:] = run.fluxes2[var]['Flux'][t,:,:]*\
                                                run.data['S'][t,:,coord[0]:coord[1],coord[2]:coord[3]]
            run.fluxes2[var]['FluxSumS'][t] = np.nansum(np.nansum(run.fluxes2[var]['FluxS'][t,:,:])) 
            run.fluxes2[var]['FluxFW'][t,:,:] = run.fluxes2[var]['Flux'][t,:,:]*\
                    (1 - run.data['S'][t,:,coord[0]:coord[1],coord[2]:coord[3]]/S0)*tmp[t,:,coord[0]:coord[1],coord[2]:coord[3]]
            run.fluxes2[var]['FluxSumFW'][t] = np.nansum(np.nansum(run.fluxes2[var]['FluxFW'][t,:,:]))                 
            run.fluxes2[var]['FluxFW1'][t,:,:] = run.fluxes2[var]['Flux'][t,:,:]*\
                    (1 - run.data['S'][t,:,coord[0]:coord[1],coord[2]:coord[3]]/S0)*tmp1[t,:,coord[0]:coord[1],coord[2]:coord[3]]
            run.fluxes2[var]['FluxSumFW1'][t] = np.nansum(np.nansum(run.fluxes2[var]['FluxFW1'][t,:,:]))   
            run.fluxes2[var]['FluxOverFlow'][t,:,:] = run.fluxes2[var]['Flux'][t,:,:]*tmpof[t,:,coord[0]:coord[1],coord[2]:coord[3]]
            run.fluxes2[var]['FluxSumOverFlow'][t] = np.nansum(np.nansum(run.fluxes2[var]['FluxOverFlow'][t,:,:]))  
            
        for flux in ['Flux','FluxT','FluxS','FluxFW','FluxFW1','FluxOverFlow']:
            run.fluxes2[var][flux] = np.squeeze(run.fluxes2[var][flux],axis=ax_d[var])

# function takes a runclass object and adds heat and freshwater content
def freshwater_content(run):
    mask1,mask2,mask3,mask4,mask5,mask6,mask7,mask8 = region_mask(run.res)
    mask0 = np.ones_like(mask1)
    grid = grid_read(run.res)
    areamasks = {0:mask0*grid['Area'], 1 : mask1*grid['Area'], 2 : mask2*grid['Area'] ,\
                 3 : mask3*grid['Area'], 4 : mask4*grid['Area'], 5 : mask5*grid['Area'],\
                 6 : mask6*grid['Area'], 7: mask7*grid['Area'], 8 : mask8*grid['Area']}
    freshwater_levels = np.zeros((len(areamasks),run.data['S'].shape[0],len(run.Z)))
    freshwater_total = np.zeros((len(areamasks),run.data['S'].shape[0]))
    freshwater_levels1 = np.zeros((len(areamasks),run.data['S'].shape[0],len(run.Z)))
    freshwater_total1 = np.zeros((len(areamasks),run.data['S'].shape[0]))
    heat_levels = np.zeros((len(areamasks),run.data['S'].shape[0],len(run.Z)))
    heat_total = np.zeros((len(areamasks),run.data['S'].shape[0]))
    for j in range(len(areamasks)):
        maskcalc = np.tile(areamasks[j],(run.data['S'].shape[0],1,1,1))     
        datam = (1 - run.data['S']/34.8)*maskcalc 
        freshwater_levels[j,:,:] = np.nansum(np.nansum(datam,axis=2),axis=2)
        freshwater_total[j,:] = np.nansum(np.nansum(np.nansum(datam[:,0:30,:,:],axis=2),axis=2),axis=1)
        run.freshwater_levels = freshwater_levels
        run.freshwater_total  = freshwater_total
        if j == 0:
            freshwater_content = np.nansum(datam[:,0:30,:,:],axis=1)
            run.freswater_content = freshwater_content
            
        maskcalc = np.tile(areamasks[j],(run.data['S'].shape[0],1,1,1))     
        datam = (1 - run.data['S']/34.8)*(run.data['S']<34.8)*maskcalc 
        freshwater_levels1[j,:,:] = np.nansum(np.nansum(datam,axis=2),axis=2)
        freshwater_total1[j,:] = np.nansum(np.nansum(np.nansum(datam[:,0:30,:,:],axis=2),axis=2),axis=1)
        run.freshwater_levels1 = freshwater_levels1
        run.freshwater_total1  = freshwater_total1
        if j == 0:
            freshwater_content1 = np.nansum(datam[:,0:30,:,:],axis=1)
            run.freswater_content1 = freshwater_content1
        
        datam = run.data['T']*maskcalc 
        heat_levels[j,:,:] = np.nansum(np.nansum(datam,axis=2),axis=2)
        heat_total[j,:] = np.nansum(np.nansum(np.nansum(datam,axis=2),axis=2),axis=1)
        if j == 0:
            heat_content = np.nansum(datam,axis=1)
            run.heat_content = heat_content
        run.heat_levels = heat_levels
        run.heat_total  = heat_total

