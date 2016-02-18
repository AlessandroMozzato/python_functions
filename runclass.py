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
from jmd95 import *

class RunRead:
    def __init__(self):
        self.data = {'T' : [], 'V' : [], 'U' : [] , 'S' : [], 'Eta' : [] , 'days' : [], 'years' : [],\
                    'theta_mean' : [], 'theta_min' : [], 'theta_max' : [],\
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
        
        self.mean = {}
        self.psi = []
        self.psi_mean = []
        self.psi_max = []
        self.psi_min = []
        self.years = []
        self.res = []
        self.grid = []
        self.hfacc = []
        self.hfacs = []
        self.hfacw = []
        self.lon = []
        self.lat = []
        self.depth = []
        self.X = []
        self.Y = []
        self.fluxes = {}
        self.totalFluxes = {}
        
    def readMonitorData(self,iters):
        self.data['theta_mean'], self.data['theta_max'], self.data['theta_min'], self.data['eta_mean'],\
        self.data['eta_max'], self.data['eta_min'], self.data['salt_mean'], self.data['salt_max'],\
        self.data['salt_min'], self.data['sss_mean'], self.data['sss_max'], self.data['sss_min'],\
        self.data['sst_mean'], self.data['sst_max'], self.data['sst_min'], self.data['vvel_mean'],\
        self.data['vvel_max'], self.data['vvel_min'], self.data['uvel_mean'], self.data['uvel_max'],\
        self.data['uvel_min'], self.data['ke_mean'], self.data['ke_max'], self.data['ke_vol'], \
        self.data['time_seconds'] = monitor(self.path,iters)
        self.data['time_years'] = (self.data['time_seconds']- self.data['time_seconds'][0])/(360*60*60*24)

    def readMonitorSeaiceData(self,iters):
        self.data['seaice_area_max'], self.data['seaice_area_min'], self.data['seaice_area_mean'], \
        self.data['seaice_heff_max'], self.data['seaice_heff_min'], self.data['seaice_heff_mean'], \
        self.data['time_seconds_ice'] = monitor_seaice(self.path,iters)
        self.data['time_years_ice'] = (self.data['time_seconds_ice']- self.data['time_seconds_ice'][0])/(360*60*60*24)

    def readDynStData(self,iters):        
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
        self.dataDyn['time_lv'] = dynStDiag(self.path,iters)    
        self.dataDyn['rho_lv_mean'] = rhop(self.dataDyn['salt_lv_mean'],self.dataDyn['theta_lv_mean'])
        self.dataDyn['time_lv_years'] = (self.dataDyn['time_lv']- self.dataDyn['time_lv'][0])/(360*60*60*24)
        self.data['time_years'] = (self.data['time_seconds']- self.data['time_seconds'][0])/(360*60*60*24)
     
    def title(self,title,color):
        self.title = title
        self.color = color
        
    def getPath(self,path):
        self.path = path

    def readStateData(self,list_var):
        file2read = netcdf.NetCDFFile(self.path+'state.nc','r')
        Temp=file2read.variables['Temp']
        self.data['T']=Temp[list_var]*1
        V=file2read.variables['V']
        self.data['V']=V[list_var]*1
        U=file2read.variables['U']
        self.data['U']=U[list_var]*1
        S=file2read.variables['S']
        self.data['S']=S[list_var]*1
        Eta=file2read.variables['Eta']
        self.data['Eta']=Eta[list_var]*1
        days=file2read.variables['T']
        self.data['days']=days[list_var]*1
        self.years = (self.data['days'] - self.data['days'][0])/(60*60*24*360)
        file2read.close()

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
        self.hfacc = hfacc[:]*1
        hfacw = file2read.variables['HFacW']
        self.hfacw = hfacw[:]*1
        hfacs = file2read.variables['HFacS']
        self.hfacs = hfacs[:]*1
        lat = file2read.variables['YC']
        self.lat = lat[:]*1        
        lon = file2read.variables['XC']
        self.lon = lon[:]*1
        depth = file2read.variables['Z']
        self.depth = depth[:]*1
        X = file2read.variables['X']
        self.X = X[:]*1
        Y = file2read.variables['Y']
        self.Y = Y[:]*1
        Z = file2read.variables['Z']
        self.Z = Z[:]*1
        file2read.close()
        
        self.data['T'][:,self.hfacc==0] = np.nan
        self.data['U'][:,self.hfacw==0] = np.nan
        self.data['V'][:,self.hfacs==0] = np.nan
        self.data['S'][:,self.hfacc==0] = np.nan
        self.data['Eta'][:,self.hfacc[0,:,:]==0] = np.nan
        self.T = np.nanmean(self.data['T'],axis=0)
        self.S = np.nanmean(self.data['S'],axis=0)
        self.rho = rho(self.S,self.T)
        self.rhop = rhop(self.S,self.T)

        # calculate jmd 95 density
        self.rho_jmd = np.zeros_like(self.rho)
        for z in range(len(self.Z)):
            self.rho_jmd[z,:,:] = densjmd95(self.S[z,:,:],self.T[z,:,:],-9.81*self.Z[z]*1027.5)
        
        self.V = np.nanmean(self.data['V'],axis=0)
        self.U = np.nanmean(self.data['U'],axis=0)
        self.Vda = np.nanmean(np.nanmean(self.data['V'],axis=1),axis=0)
        self.Uda = np.nanmean(np.nanmean(self.data['U'],axis=1),axis=0)

        print 'Data read from '+self.path

    def saveflux(self,fold):
        save_flux(self,fold)

    def read_flux(self,fold):
        self.fluxes['Fram'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }
        self.fluxes['Fram1'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }
        self.fluxes['Fram2'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  } 
        self.fluxes['Davis'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }
        self.fluxes['Davis1'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }
        self.fluxes['Davis2'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }
        self.fluxes['Barents'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }
        self.fluxes['Barents1'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }
        self.fluxes['Barents2'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }        
        self.fluxes['Bering'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }
        self.fluxes['Denmark'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }
        self.fluxes['Norwice'] = { 'FluxSum' : [] , 'FluxInSum' : [] , 'FluxOutSum' : [] ,
                            'FluxSumS'  : [] , 'FluxSumT' : []  }
        
        file2read = netcdf.NetCDFFile('/scratch/general/am8e13/results_saved/'+fold+'/fluxes.nc', 'r') 
        for var in ['Davis' , 'Davis1' , 'Davis2' , 'Fram' , 'Fram1', 'Fram2', 'Norwice' , \
                        'Denmark' , 'Barents' , 'Barents1', 'Barents2','Bering', ]:
            for flux in ['FluxSum' , 'FluxSumS' , 'FluxSumT' , 'FluxInSum' , 'FluxOutSum' ]:
                temp = file2read.variables[var+'_'+flux]
                self.fluxes[var][flux] = temp[:]*1
        temp = file2read.variables['years']
        self.fluxes['years'] = temp[:]*1

    def fluxCalc(self):
            def fluxTransport(data):
                ave_pos = np.zeros_like(data)
                ave_pos[data<0]=0
                ave_pos[data>0] = data[data>0] 
                ave_neg = np.zeros_like(data)
                ave_neg[data>0]=0
                ave_neg[data<0] = data[data<0]
                n_pos = np.mean(np.sum(np.sum(ave_pos,axis=2),axis=1),axis=0)
                n_neg = np.mean(np.sum(np.sum(ave_neg,axis=2),axis=1),axis=0)
                n_all = np.nanmean(np.nansum(np.nansum(data,axis=2),axis=1),axis=0)
                return {'Inflow' : round(n_pos/10**6,2), 'Outflow' : round(n_neg/10**6,2) , 'Total flow' : round(n_all/10**6,2)}

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
            Area_x = dydx/dxF
            Area_y = dydx/dyF
            Area_x[self.hfacc==0]=np.nan
            Area_y[self.hfacc==0]=np.nan
            
            barents = [40,41,53,68,40,58,68,69]
            fram = [58,80,76,77,80,81,76,77] 

            self.fluxes['Fram'] = {'Flux' : np.zeros_like(self.data['V'][:,:,58*kk:80*kk,76*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,58*kk:80*kk,76*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,58*kk:80*kk,76*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
            
            self.fluxes['Fram1'] = {'Flux' : np.zeros_like(self.data['V'][:,:,55*kk:85*kk,78*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,55*kk:85*kk,78*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,55*kk:85*kk,78*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
                    
            self.fluxes['Fram2'] = {'Flux' : np.zeros_like(self.data['V'][:,:,60*kk:83*kk,72*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,60*kk:83*kk,72*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,60*kk:83*kk,72*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
            
            self.fluxes['Bering'] = {'Flux' : np.zeros_like(self.data['U'][:,:,80*kk:89*kk,178*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['U'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,80*kk:89*kk,178*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,80*kk:89*kk,178*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
            
            self.fluxes['Davis'] = {'Flux' : np.zeros_like(self.data['U'][:,:,113*kk:135*kk,75*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['U'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,113*kk:135*kk,75*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,113*kk:135*kk,75*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}   
            
            self.fluxes['Davis1'] = {'Flux' : np.zeros_like(self.data['U'][:,:,135*kk,52*kk:73*kk]), \
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,135*kk,52*kk:73*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,135*kk,52*kk:73*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}   
            
            self.fluxes['Davis2'] = {'Flux' : np.zeros_like(self.data['U'][:,:,113*kk:135*kk,75*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['U'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,113*kk:135*kk,75*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,113*kk:135*kk,75*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []} 
            
            self.fluxes['Denmark'] = {'Flux' : np.zeros_like(self.data['U'][:,:,100*kk,37*kk:48*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['U'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,100*kk,37*kk:48*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,100*kk,37*kk:48*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
 
            self.fluxes['Barents'] = {'Flux' : np.concatenate((\
                                        self.data['V'][:,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                        self.data['U'][:,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=2),\
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.concatenate((\
                                        self.data['S'][:,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                        self.data['S'][:,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=2),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.concatenate((\
                                        self.data['T'][:,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                        self.data['T'][:,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=2),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0])}

            
            # Barents 1
            self.fluxes['Barents1'] = {'Flux' : np.concatenate((\
                                        self.data['V'][:,:,45*kk,49*kk:66*kk]*Area_x[:,45*kk,49*kk:66*kk],\
                                        self.data['U'][:,:,45*kk:58*kk,66*kk]*Area_y[:,45*kk:58*kk,66*kk]),axis=2),\
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.concatenate((\
                                        self.data['S'][:,:,45*kk,49*kk:66*kk]*Area_x[:,45*kk,49*kk:66*kk],\
                                        self.data['S'][:,:,45*kk:58*kk,66*kk]*Area_y[:,45*kk:58*kk,66*kk]),axis=2),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.concatenate((\
                                        self.data['T'][:,:,45*kk,49*kk:66*kk]*Area_x[:,45*kk,49*kk:66*kk],\
                                        self.data['T'][:,:,45*kk:58*kk,66*kk]*Area_y[:,45*kk:58*kk,66*kk]),axis=2),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0])}
        
            # Barents 2
            self.fluxes['Barents2'] = {'Flux' : np.concatenate((\
                                        self.data['V'][:,:,58*kk,49*kk:66*kk]*Area_x[:,58*kk,49*kk:66*kk],\
                                        self.data['U'][:,:,45*kk:58*kk,49*kk]*Area_y[:,45*kk:58*kk,49*kk]),axis=2),\
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.concatenate((\
                                        self.data['S'][:,:,58*kk,49*kk:66*kk]*Area_x[:,58*kk,49*kk:66*kk],\
                                        self.data['S'][:,:,45*kk:58*kk,49*kk]*Area_y[:,45*kk:58*kk,49*kk]),axis=2),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.concatenate((\
                                        self.data['T'][:,:,58*kk,49*kk:66*kk]*Area_x[:,58*kk,49*kk:66*kk],\
                                        self.data['T'][:,:,45*kk:58*kk,49*kk]*Area_y[:,45*kk:58*kk,49*kk]),axis=2),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0])}
        
            self.fluxes['Norwice'] = {'Flux' : np.concatenate((\
                                        self.data['U'][:,:,60*kk:95*kk,15*kk]*Area_y[:,60*kk:95*kk,15*kk],\
                                        self.data['V'][:,:,95*kk,15*kk:30*kk]*Area_x[:,95*kk,15*kk:30*kk]),axis=2),\
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxInSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxOutSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.concatenate((\
                                        self.data['S'][:,:,60*kk:95*kk,15*kk]*Area_y[:,60*kk:95*kk,15*kk],\
                                        self.data['S'][:,:,95*kk,15*kk:30*kk]*Area_x[:,95*kk,15*kk:30*kk]),axis=2),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.concatenate((\
                                        self.data['T'][:,:,60*kk:95*kk,15*kk]*Area_y[:,60*kk:95*kk,15*kk],\
                                        self.data['T'][:,:,95*kk,15*kk:30*kk]*Area_x[:,95*kk,15*kk:30*kk]),axis=2),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0])}

            for t in range(self.data['V'].shape[0]):
                # Fram fillign
                self.fluxes['Fram']['Flux'][t,:,:] = self.data['U'][t,:,58*kk:80*kk,76*kk]*Area_y[:,58*kk:80*kk,76*kk]
                self.fluxes['Fram']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Fram']['Flux'][t,:,:]))
                self.fluxes['Fram']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Fram']['Flux'][t,self.fluxes['Fram']['Flux'][t,:,:]>0]))
                self.fluxes['Fram']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Fram']['Flux'][t,self.fluxes['Fram']['Flux'][t,:,:]<0]))
                self.fluxes['Fram']['FluxT'][t,:,:] = self.fluxes['Fram']['Flux'][t,:,:]*\
                                            self.data['T'][t,:,58*kk:80*kk,76*kk]
                self.fluxes['Fram']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Fram']['FluxT'][t,:,:]))
                self.fluxes['Fram']['FluxS'][t,:,:] = self.fluxes['Fram']['Flux'][t,:,:]*\
                                            self.data['S'][t,:,58*kk:80*kk,76*kk]
                self.fluxes['Fram']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Fram']['FluxS'][t,:,:]))                    
                self.totalFluxes['Fram'] = fluxTransport(self.fluxes['Fram']['Flux'])

                # Fram 1 fillign
                self.fluxes['Fram1']['Flux'][t,:,:] = self.data['U'][t,:,55*kk:85*kk,78*kk]*Area_y[:,55*kk:85*kk,78*kk]
                self.fluxes['Fram1']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Fram1']['Flux'][t,:,:]))
                self.fluxes['Fram1']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Fram1']['Flux'][t,self.fluxes['Fram1']['Flux'][t,:,:]>0]))
                self.fluxes['Fram1']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Fram1']['Flux'][t,self.fluxes['Fram1']['Flux'][t,:,:]<0]))
                self.fluxes['Fram1']['FluxT'][t,:,:] = self.fluxes['Fram1']['Flux'][t,:,:]*\
                                            self.data['T'][t,:,55*kk:85*kk,78*kk]
                self.fluxes['Fram1']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Fram1']['FluxT'][t,:,:]))
                self.fluxes['Fram1']['FluxS'][t,:,:] = self.fluxes['Fram1']['Flux'][t,:,:]*\
                                            self.data['S'][t,:,55*kk:85*kk,78*kk]
                self.fluxes['Fram1']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Fram1']['FluxS'][t,:,:]))                    
                self.totalFluxes['Fram1'] = fluxTransport(self.fluxes['Fram1']['Flux'])
                
                # Fram 2 fillign
                self.fluxes['Fram2']['Flux'][t,:,:] = self.data['U'][t,:,60*kk:83*kk,72*kk]*Area_y[:,60*kk:83*kk,72*kk]
                self.fluxes['Fram2']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Fram2']['Flux'][t,:,:]))
                self.fluxes['Fram2']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Fram2']['Flux'][t,self.fluxes['Fram2']['Flux'][t,:,:]>0]))
                self.fluxes['Fram2']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Fram2']['Flux'][t,self.fluxes['Fram2']['Flux'][t,:,:]<0]))
                self.fluxes['Fram2']['FluxT'][t,:,:] = self.fluxes['Fram2']['Flux'][t,:,:]*\
                                            self.data['T'][t,:,60*kk:83*kk,72*kk]
                self.fluxes['Fram2']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Fram2']['FluxT'][t,:,:]))
                self.fluxes['Fram2']['FluxS'][t,:,:] =self.fluxes['Fram2']['Flux'][t,:,:]*\
                                            self.data['S'][t,:,60*kk:83*kk,72*kk]
                self.fluxes['Fram2']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Fram2']['FluxS'][t,:,:]))                    
                self.totalFluxes['Fram2'] = fluxTransport(self.fluxes['Fram2']['Flux'])                    
                
                # Bering fillign
                self.fluxes['Bering']['Flux'][t,:,:] = -self.data['U'][t,:,80*kk:89*kk,178*kk]*Area_y[:,80*kk:89*kk,178*kk]
                self.fluxes['Bering']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Bering']['Flux'][t,:,:]))
                self.fluxes['Bering']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Bering']['Flux'][t,self.fluxes['Bering']['Flux'][t,:,:]>0]))
                self.fluxes['Bering']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Bering']['Flux'][t,self.fluxes['Bering']['Flux'][t,:,:]<0]))
                self.fluxes['Bering']['FluxT'][t,:,:] = self.fluxes['Bering']['Flux'][t,:,:]*\
                                            self.data['T'][t,:,80*kk:89*kk,178*kk]
                self.fluxes['Bering']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Bering']['FluxT'][t,:,:]))
                self.fluxes['Bering']['FluxS'][t,:,:] = self.fluxes['Bering']['Flux'][t,:,:]*\
                                            self.data['S'][t,:,80*kk:89*kk,178*kk]
                self.fluxes['Bering']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Bering']['FluxS'][t,:,:]))
                self.totalFluxes['Bering'] = fluxTransport(self.fluxes['Bering']['Flux'])
                
                # Davis fillign
                self.fluxes['Davis']['Flux'][t,:,:] = self.data['U'][t,:,113*kk:135*kk,75*kk]*Area_y[:,113*kk:135*kk,75*kk]
                self.fluxes['Davis']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Davis']['Flux'][t,:,:]))
                self.fluxes['Davis']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Davis']['Flux'][t,self.fluxes['Davis']['Flux'][t,:,:]>0]))
                self.fluxes['Davis']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Davis']['Flux'][t,self.fluxes['Davis']['Flux'][t,:,:]<0]))
                self.fluxes['Davis']['FluxT'][t,:,:] = self.fluxes['Davis']['Flux'][t,:,:]*\
                                            self.data['T'][t,:,113*kk:135*kk,75*kk]
                self.fluxes['Davis']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Davis']['FluxT'][t,:,:]))
                self.fluxes['Davis']['FluxS'][t,:,:] = self.fluxes['Davis']['Flux'][t,:,:]*\
                                            self.data['S'][t,:,113*kk:135*kk,75*kk]
                self.fluxes['Davis']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Davis']['FluxS'][t,:,:]))
                self.totalFluxes['Davis'] = fluxTransport(self.fluxes['Davis']['Flux']) 
                
                # Davis1 fillign
                self.fluxes['Davis1']['Flux'][t,:,:] = self.data['V'][t,:,135*kk,52*kk:73*kk]*Area_x[:,135*kk,52*kk:73*kk]
                self.fluxes['Davis1']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Davis1']['Flux'][t,:,:]))
                self.fluxes['Davis1']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Davis1']['Flux'][t,self.fluxes['Davis1']['Flux'][t,:,:]>0]))
                self.fluxes['Davis1']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Davis1']['Flux'][t,self.fluxes['Davis1']['Flux'][t,:,:]<0]))                
                self.fluxes['Davis1']['FluxT'][t,:,:] = self.fluxes['Davis1']['Flux'][t,:,:]*\
                                            self.data['T'][t,:,135*kk,52*kk:73*kk]
                self.fluxes['Davis1']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Davis1']['FluxT'][t,:,:]))
                self.fluxes['Davis1']['FluxS'][t,:,:] = self.fluxes['Davis1']['Flux'][t,:,:]*\
                                            self.data['S'][t,:,135*kk,52*kk:73*kk]
                self.fluxes['Davis1']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Davis1']['FluxS'][t,:,:]))
                self.totalFluxes['Davis1'] = fluxTransport(self.fluxes['Davis1']['Flux']) 
                
                # Davis 2 fillign
                self.fluxes['Davis2']['Flux'][t,:,:] = self.data['U'][t,:,113*kk:135*kk,72*kk]*Area_y[:,113*kk:135*kk,72*kk]
                self.fluxes['Davis2']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Davis2']['Flux'][t,:,:]))
                self.fluxes['Davis2']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Davis2']['Flux'][t,self.fluxes['Davis2']['Flux'][t,:,:]>0]))
                self.fluxes['Davis2']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Davis2']['Flux'][t,self.fluxes['Davis2']['Flux'][t,:,:]<0]))                
                self.fluxes['Davis2']['FluxT'][t,:,:] = self.fluxes['Davis2']['Flux'][t,:,:]*\
                                            self.data['T'][t,:,113*kk:135*kk,72*kk]
                self.fluxes['Davis2']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Davis2']['FluxT'][t,:,:]))
                self.fluxes['Davis2']['FluxS'][t,:,:] = self.fluxes['Davis2']['Flux'][t,:,:]*\
                                            self.data['S'][t,:,113*kk:135*kk,72*kk]
                self.fluxes['Davis2']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Davis2']['FluxS'][t,:,:]))
                self.totalFluxes['Davis2'] = fluxTransport(self.fluxes['Davis2']['Flux']) 

                # Barents fillng
                self.fluxes['Barents']['Flux'][t,:,:] = np.concatenate((\
                                            self.data['V'][t,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                            self.data['U'][t,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=1)
                self.fluxes['Barents']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Barents']['Flux'][t,:,:]))
                self.fluxes['Barents']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Barents']['Flux'][t,self.fluxes['Barents']['Flux'][t,:,:]>0]))
                self.fluxes['Barents']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Barents']['Flux'][t,self.fluxes['Barents']['Flux'][t,:,:]<0]))                
                self.fluxes['Barents']['FluxT'][t,:,:] = self.fluxes['Barents']['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,40*kk,53*kk:68*kk],\
                                            self.data['T'][t,:,40*kk:58*kk,68*kk]),axis=1)
                self.fluxes['Barents']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Barents']['FluxT'][t,:,:]))
                self.fluxes['Barents']['FluxS'][t,:,:] = self.fluxes['Barents']['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,40*kk,53*kk:68*kk],\
                                            self.data['T'][t,:,40*kk:58*kk,68*kk]),axis=1)
                self.fluxes['Barents']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Barents']['FluxS'][t,:,:]))
                self.totalFluxes['Barents'] = fluxTransport(self.fluxes['Barents']['Flux'])

                # Barents 1 fillng
                y1 = 45; y2 = 58 ; x1 = 49 ; x2 = 66 # barents 1
                self.fluxes['Barents1']['Flux'][t,:,:] = np.concatenate((\
                                            self.data['V'][t,:,45*kk,49*kk:66*kk]*Area_x[:,45*kk,49*kk:66*kk],\
                                            self.data['U'][t,:,45*kk:58*kk,66*kk]*Area_y[:,45*kk:58*kk,66*kk]),axis=1)
                self.fluxes['Barents1']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Barents1']['Flux'][t,:,:]))
                self.fluxes['Barents1']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Barents1']['Flux'][t,self.fluxes['Barents1']['Flux'][t,:,:]>0]))
                self.fluxes['Barents1']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Barents1']['Flux'][t,self.fluxes['Barents1']['Flux'][t,:,:]<0]))                  
                self.fluxes['Barents1']['FluxT'][t,:,:] = self.fluxes['Barents1']['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,45*kk,49*kk:66*kk],\
                                            self.data['T'][t,:,45*kk:58*kk,66*kk]),axis=1)
                self.fluxes['Barents1']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Barents1']['FluxT'][t,:,:]))
                self.fluxes['Barents1']['FluxS'][t,:,:] = self.fluxes['Barents1']['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,45*kk,49*kk:66*kk],\
                                            self.data['T'][t,:,45*kk:58*kk,66*kk]),axis=1)
                self.fluxes['Barents1']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Barents1']['FluxS'][t,:,:]))
                self.totalFluxes['Barents1'] = fluxTransport(self.fluxes['Barents1']['Flux'])
                
                # Barents 2 fillng
                y1 = 45; y2 = 58 ; x1 = 49 ; x2 = 66 # barents 1
                self.fluxes['Barents2']['Flux'][t,:,:] = np.concatenate((\
                                            self.data['V'][t,:,58*kk,49*kk:66*kk]*Area_x[:,58*kk,49*kk:66*kk],\
                                            self.data['U'][t,:,45*kk:58*kk,49*kk]*Area_y[:,45*kk:58*kk,49*kk]),axis=1)
                self.fluxes['Barents2']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Barents2']['Flux'][t,:,:]))
                self.fluxes['Barents']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Barents']['Flux'][t,self.fluxes['Barents']['Flux'][t,:,:]>0]))
                self.fluxes['Barents']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Barents']['Flux'][t,self.fluxes['Barents']['Flux'][t,:,:]<0]))                  
                self.fluxes['Barents2']['FluxT'][t,:,:] = self.fluxes['Barents2']['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,58*kk,49*kk:66*kk],\
                                            self.data['T'][t,:,45*kk:58*kk,49*kk]),axis=1)
                self.fluxes['Barents2']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Barents2']['FluxT'][t,:,:]))
                self.fluxes['Barents2']['FluxS'][t,:,:] = self.fluxes['Barents2']['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,58*kk,49*kk:66*kk],\
                                            self.data['T'][t,:,45*kk:58*kk,49*kk]),axis=1)
                self.fluxes['Barents2']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Barents2']['FluxS'][t,:,:]))
                self.totalFluxes['Barents2'] = fluxTransport(self.fluxes['Barents2']['Flux'])
               
                # Denmark filling
                self.fluxes['Denmark']['Flux'][t,:,:] = self.data['V'][t,:,100*kk,37*kk:48*kk]*Area_x[:,100*kk,37*kk:48*kk]
                self.fluxes['Denmark']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Denmark']['Flux'][t,:,:]))
                self.fluxes['Denmark']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Denmark']['Flux'][t,self.fluxes['Denmark']['Flux'][t,:,:]>0]))
                self.fluxes['Denmark']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Denmark']['Flux'][t,self.fluxes['Denmark']['Flux'][t,:,:]<0]))  
                self.fluxes['Denmark']['FluxT'][t,:,:] = self.fluxes['Denmark']['Flux'][t,:,:]*\
                                            self.data['T'][t,:,100*kk,37*kk:48*kk]
                self.fluxes['Denmark']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Denmark']['FluxT'][t,:,:]))
                self.fluxes['Denmark']['FluxS'][t,:,:] = self.fluxes['Denmark']['Flux'][t,:,:]*\
                                            self.data['S'][t,:,100*kk,37*kk:48*kk]
                self.fluxes['Denmark']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Denmark']['FluxS'][t,:,:]))
                self.totalFluxes['Denmark'] = fluxTransport(self.fluxes['Denmark']['Flux'])

                # Norwey-Iceland strait filling
                self.fluxes['Norwice']['Flux'][t,:,:] = np.concatenate((\
                                        self.data['V'][t,:,60*kk:95*kk,15*kk]*Area_y[:,60*kk:95*kk,15*kk],\
                                        self.data['U'][t,:,95*kk,15*kk:30*kk]*Area_x[:,95*kk,15*kk:30*kk]),axis=1)
                self.fluxes['Norwice']['FluxSum'][t] = np.nansum(np.nansum(self.fluxes['Norwice']['Flux'][t,:,:]))
                self.fluxes['Norwice']['FluxInSum'][t] = np.nansum(np.nansum(self.fluxes['Norwice']['Flux'][t,self.fluxes['Norwice']['Flux'][t,:,:]>0]))
                self.fluxes['Norwice']['FluxOutSum'][t] = np.nansum(np.nansum(self.fluxes['Norwice']['Flux'][t,self.fluxes['Norwice']['Flux'][t,:,:]<0]))  
                self.fluxes['Norwice']['FluxT'][t,:,:] = self.fluxes['Norwice']['Flux'][t,:,:]*np.concatenate((\
                                        self.data['T'][t,:,60*kk:95*kk,15*kk],\
                                        self.data['T'][t,:,95*kk,15*kk:30*kk]),axis=1)
                self.fluxes['Norwice']['FluxSumT'][t] = np.nansum(np.nansum(self.fluxes['Norwice']['FluxT'][t,:,:]))
                self.fluxes['Norwice']['FluxS'][t,:,:] = self.fluxes['Norwice']['Flux'][t,:,:]*np.concatenate((\
                                        self.data['T'][t,:,60*kk:95*kk,15*kk],\
                                        self.data['T'][t,:,95*kk,15*kk:30*kk]),axis=1)
                self.fluxes['Norwice']['FluxSumS'][t] = np.nansum(np.nansum(self.fluxes['Norwice']['FluxS'][t,:,:]))
                self.totalFluxes['Norwice'] = fluxTransport(self.fluxes['Norwice']['Flux'])

    def baroCalc(self):
        self.psi = baro_stream(self.data['U'])
        self.psi_mean = np.nanmean(self.psi,axis = 1)
        self.psi_mean = np.nanmean(self.psi_mean,axis = 1)
        self.psi_max = np.nanmax(self.psi,axis=1)
        self.psi_max = np.nanmax(self.psi_max,axis=1)
        self.psi_min = np.nanmin(self.psi,axis=1)
        self.psi_min = np.nanmin(self.psi_min,axis=1)
        self.psi_ave = np.nanmean(self.psi,axis=0)

    # We save the barotropic streamfunction
    def savepsi(self,fold):
        save_psi(self,fold)

    # We read the saved barotropic streamfunction
    def read_psi(self,fold):
        file2read = netcdf.NetCDFFile('/scratch/general/am8e13/results_saved/'+fold+'/barostream.nc', 'r') 
        self.barostream = {}
        temp = file2read.variables['psi_max']
        self.barostream['psi_max']= temp[:]*1
        temp = file2read.variables['psi_min']
        self.barostream['psi_min']= temp[:]*1
        temp = file2read.variables['psi_mean']
        self.barostream['psi_mean']= temp[:]*1
        temp = file2read.variables['psi_years']
        self.barostream['years']= temp[:]*1

    def topoCalc(self):
        topo,topo_opposit = Topostrophy(self.data['U'],self.data['V'])
        self.topo = topo
        self.topo_opposit = topo_opposit
    
    def print_title(self):
        print self.title
        
    def vorticity(self,list_var):
        # this functions reads vorticity/divergence/kynetic energy from the files
        self.vort = {'vorticity' : [] , 'div' : [] , 'KE' : [] , 'seconds' : [] , 'years_vort' : []}
        file2read = netcdf.NetCDFFile(self.path+'vorticity.nc','r')
        momHDiv = file2read.variables['momHDiv']
        self.vort['div'] = momHDiv[list_var]*1
        momKE = file2read.variables['momKE']
        self.vort['KE'] = momKE[list_var]*1
        momVort3 = file2read.variables['momVort3']
        self.vort['vorticity'] = momVort3[list_var]*1
        T = file2read.variables['T']
        self.vort['seconds'] = T[list_var]*1
        self.vort['years_vort'] = self.vort['seconds']/(60*60*24*360) - self.vort['seconds'][0]/(60*60*24*360)
        
        # put nans
        self.vort['vorticity'][:,self.hfacc==0] = np.nan
        self.vort['div'][:,self.hfacc==0] = np.nan
        self.vort['KE'][:,self.hfacc==0] = np.nan

    def seaiceread(self,list_var):
        self.seaice = { 'SIarea' : [] , 'SIheff' : [] , 'SIuice' : [] , 'SIvice' : [] , 'seconds' : [] , 'years_seaice' : [] }
        file2read = netcdf.NetCDFFile(self.path+'seaice.nc','r')
        SIarea = file2read.variables['SIarea']
        self.seaice['SIarea'] = SIarea[list_var]*1
        SIheff = file2read.variables['SIheff']
        self.seaice['SIheff'] = SIheff[list_var]*1
        SIuice = file2read.variables['SIuice']
        self.seaice['SIuice'] = SIuice[list_var]*1
        SIvice = file2read.variables['SIvice']
        self.seaice['SIvice'] = SIvice[list_var]*1
        T = file2read.variables['T']
        self.seaice['seconds'] = T[list_var]*1
        self.seaice['years_seaice'] = self.seaice['seconds']/(60*60*24*360) - self.seaice['seconds'][0]/(60*60*24*360)
        
        # put nans
        self.seaice['SIarea'][:,0,self.hfacc[0,:,:]==0] = np.nan
        self.seaice['SIheff'][:,0,self.hfacc[0,:,:]==0] = np.nan
        self.seaice['SIuice'][:,0,self.hfacc[0,:,:]==0] = np.nan
        self.seaice['SIvice'][:,0,self.hfacc[0,:,:]==0] = np.nan

    def mxldepthread(self,list_var):
        file2read = netcdf.NetCDFFile(self.path+'MXLDEPTH.nc','r')
        MXLDEPTH = file2read.variables['MXLDEPTH']
        self.mxldepth = MXLDEPTH[list_var]*1
        self.mxldepth = self.mxldepth[:,0,:,:]
        self.mxldepth[:,self.hfacc[0,:,:]==0]=np.nan
        # Read time
        days=file2read.variables['T']
        days=days[list_var]*1
        self.mxldepth_years = (days - days[0])/(60*60*24*360)
        self.mxldepth_mean = np.nanmean(np.nanmean(self.mxldepth,axis=1),axis=1)
        self.mxldepth_max = np.nanmax(np.nanmax(self.mxldepth,axis=1),axis=1)
        self.mxldepth_min = np.nanmin(np.nanmin(self.mxldepth,axis=1),axis=1)

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
        
        # option added to make the old spinup readable
        if theta_mean_lv.shape[1] == 1:
            theta_mean_lv = np.tile(theta_mean_lv,(1,7,1))
            theta_max_lv = np.tile(theta_max_lv,(1,7,1))
            theta_min_lv = np.tile(theta_min_lv,(1,7,1))
            salt_mean_lv = np.tile(salt_mean_lv,(1,7,1))
            salt_max_lv = np.tile(salt_max_lv,(1,7,1))
            salt_min_lv = np.tile(salt_min_lv,(1,7,1))
            vvel_mean_lv = np.tile(vvel_mean_lv,(1,7,1))
            vvel_max_lv = np.tile(vvel_max_lv,(1,7,1))
            vvel_min_lv = np.tile(vvel_min_lv,(1,7,1))
            uvel_mean_lv = np.tile(uvel_mean_lv,(1,7,1))
            uvel_max_lv = np.tile(uvel_max_lv,(1,7,1))
            uvel_min_lv = np.tile(uvel_min_lv,(1,7,1))
            ke_mean_lv = np.tile(ke_mean_lv,(1,7,1))
            ke_max_lv = np.tile(ke_max_lv,(1,7,1))
            
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

        # This is to make old spinup readable
        if theta_mean.shape[1] == 1:                                                                                                   
            theta_mean=np.tile(theta_mean,(1,7,1))
            theta_max=np.tile(theta_max,(1,7,1))
            theta_min=np.tile(theta_min,(1,7,1))
            salt_mean=np.tile(salt_mean,(1,7,1))
            salt_max=np.tile(salt_max,(1,7,1))
            salt_min=np.tile(salt_min,(1,7,1))
            eta_mean=np.tile(eta_mean,(1,7,1))
            eta_max=np.tile(eta_max,(1,7,1))
            eta_min=np.tile(eta_min,(1,7,1))
            uvel_mean=np.tile(uvel_mean,(1,7,1))
            uvel_max=np.tile(uvel_max,(1,7,1))
            uvel_min=np.tile(uvel_min,(1,7,1))
            vvel_mean=np.tile(vvel_mean,(1,7,1))
            vvel_max=np.tile(vvel_max,(1,7,1))
            vvel_min=np.tile(vvel_min,(1,7,1))
            ke_mean=np.tile(ke_mean,(1,7,1))
            ke_max=np.tile(ke_max,(1,7,1))
            ke_vol=np.tile(ke_vol,(1,7,1))                                                                          
 
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

def mon_titles():
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

def save_flux(data,fold): 
    # this files saves fluxes
    # generate the folder
    if not(os.path.isdir('/scratch/general/am8e13/results_saved/'+fold)):
        os.mkdir('/scratch/general/am8e13/results_saved/'+fold)
    
    years = len(data.years)
    # make file
    f = netcdf.NetCDFFile('/scratch/general/am8e13/results_saved/'+fold+'/fluxes.nc', 'w')  
    # dimension
    f.createDimension('time', years)
    
    for var in data.fluxes:
        for flux in ['FluxSum' , 'FluxSumS' , 'FluxSumT' , 'FluxInSum' , 'FluxOutSum' ]:
            # variable
            name = var+'_'+flux
            temp = f.createVariable(name, 'float', ('time',))
            # data
            temp[:] = data.fluxes[var][flux]
    temp = f.createVariable('years', 'float', ('time',))
    temp[:] = data.years
    f.close()

def save_psi(data,fold): 
    # generate the folder
    if not(os.path.isdir('/scratch/general/am8e13/results_saved/'+fold)):
        os.mkdir('/scratch/general/am8e13/results_saved/'+fold)
    
    years = len(data.years)
    # make file
    f = netcdf.NetCDFFile('/scratch/general/am8e13/results_saved/'+fold+'/barostream.nc', 'w')  
    # dimension
    f.createDimension('time', years)
    # variable
    temp = f.createVariable('psi_max', 'float', ('time',))
    temp[:] = data.psi_max
    temp = f.createVariable('psi_min', 'float', ('time',))
    temp[:] = data.psi_min
    temp = f.createVariable('psi_mean', 'float', ('time',))
    temp[:] = data.psi_mean
    temp = f.createVariable('psi_years', 'float', ('time',))
    temp[:] = data.years
    f.close()
