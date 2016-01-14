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

class StateRead:
    def __init__(self):
        self.data = {'T' : [], 'V' : [], 'U' : [] , 'S' : [], 'Eta' : [] , 'days' : [], 'years' : []}
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
        self.totalFluxes = {}
        
    def getPath(self,path):
        self.path = path

    def readData(self,list_var):
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

    def getMeans(self,list_iter,list_var):
        file2read = netcdf.NetCDFFile(self.path+'state.nc','r')
        for var in list_var:
            self.mean[var] = np.ndarray(self.hfacc.shape)
        for i in list_iter:
            for var in list_var:
                tmp = file2read.variables[var]
                self.mean[var] = self.mean[var] + tmp[i,:,:,:]*float(1)/float(len(list_iter))
        file2read.close()

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

            self.Fram = {'Flux' : np.zeros_like(self.data['V'][:,:,58*kk:80*kk,76*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,58*kk:80*kk,76*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,58*kk:80*kk,76*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
            
            self.Fram1 = {'Flux' : np.zeros_like(self.data['V'][:,:,55*kk:85*kk,78*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,55*kk:85*kk,78*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,55*kk:85*kk,78*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
                    
            self.Fram2 = {'Flux' : np.zeros_like(self.data['V'][:,:,60*kk:83*kk,72*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,60*kk:83*kk,72*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,60*kk:83*kk,72*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
            
            self.Bering = {'Flux' : np.zeros_like(self.data['U'][:,:,80*kk:89*kk,178*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['U'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,80*kk:89*kk,178*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,80*kk:89*kk,178*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
            
            self.Davis = {'Flux' : np.zeros_like(self.data['U'][:,:,113*kk:135*kk,75*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['U'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,113*kk:135*kk,75*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,113*kk:135*kk,75*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}   
            
            self.Davis1 = {'Flux' : np.zeros_like(self.data['U'][:,:,135*kk,52*kk:73*kk]), \
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,135*kk,52*kk:73*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,135*kk,52*kk:73*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}   
            
            self.Davis2 = {'Flux' : np.zeros_like(self.data['U'][:,:,113*kk:135*kk,75*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['U'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,113*kk:135*kk,75*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,113*kk:135*kk,75*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []} 
            
            self.Denmark = {'Flux' : np.zeros_like(self.data['U'][:,:,100*kk,37*kk:48*kk]) , \
                            'FluxSum' : np.zeros_like(self.data['U'][:,0,0,0]),\
                            'FluxS' : np.zeros_like(self.data['S'][:,:,100*kk,37*kk:48*kk]),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.zeros_like(self.data['T'][:,:,100*kk,37*kk:48*kk]),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0]),\
                            'FluxTop' : [], 'FluxMid' : [] , 'FluxBot' : [], \
                            'FluxTopS' : [], 'FluxMidS' : [] , 'FluxBotS' : [], \
                            'FluxTopT' : [], 'FluxMidT' : [] , 'FluxBotT' : []}
 
            self.Barents = {'Flux' : np.concatenate((\
                                        self.data['V'][:,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                        self.data['U'][:,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=2),\
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.concatenate((\
                                        self.data['S'][:,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                        self.data['S'][:,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=2),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.concatenate((\
                                        self.data['T'][:,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                        self.data['T'][:,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=2),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0])}

            
            # Barents 1
            self.Barents1 = {'Flux' : np.concatenate((\
                                        self.data['V'][:,:,45*kk,49*kk:66*kk]*Area_x[:,45*kk,49*kk:66*kk],\
                                        self.data['U'][:,:,45*kk:58*kk,66*kk]*Area_y[:,45*kk:58*kk,66*kk]),axis=2),\
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.concatenate((\
                                        self.data['S'][:,:,45*kk,49*kk:66*kk]*Area_x[:,45*kk,49*kk:66*kk],\
                                        self.data['S'][:,:,45*kk:58*kk,66*kk]*Area_y[:,45*kk:58*kk,66*kk]),axis=2),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.concatenate((\
                                        self.data['T'][:,:,45*kk,49*kk:66*kk]*Area_x[:,45*kk,49*kk:66*kk],\
                                        self.data['T'][:,:,45*kk:58*kk,66*kk]*Area_y[:,45*kk:58*kk,66*kk]),axis=2),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0])}
        
            # Barents 2
            self.Barents2 = {'Flux' : np.concatenate((\
                                        self.data['V'][:,:,58*kk,49*kk:66*kk]*Area_x[:,58*kk,49*kk:66*kk],\
                                        self.data['U'][:,:,45*kk:58*kk,49*kk]*Area_y[:,45*kk:58*kk,49*kk]),axis=2),\
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
                            'FluxS' : np.concatenate((\
                                        self.data['S'][:,:,58*kk,49*kk:66*kk]*Area_x[:,58*kk,49*kk:66*kk],\
                                        self.data['S'][:,:,45*kk:58*kk,49*kk]*Area_y[:,45*kk:58*kk,49*kk]),axis=2),\
                            'FluxSumS' : np.zeros_like(self.data['S'][:,0,0,0]),\
                            'FluxT' : np.concatenate((\
                                        self.data['T'][:,:,58*kk,49*kk:66*kk]*Area_x[:,58*kk,49*kk:66*kk],\
                                        self.data['T'][:,:,45*kk:58*kk,49*kk]*Area_y[:,45*kk:58*kk,49*kk]),axis=2),\
                            'FluxSumT' : np.zeros_like(self.data['T'][:,0,0,0])}
        
            self.Norwice = {'Flux' : np.concatenate((\
                                        self.data['U'][:,:,60*kk:95*kk,15*kk]*Area_y[:,60*kk:95*kk,15*kk],\
                                        self.data['V'][:,:,95*kk,15*kk:30*kk]*Area_x[:,95*kk,15*kk:30*kk]),axis=2),\
                            'FluxSum' : np.zeros_like(self.data['V'][:,0,0,0]),\
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
                self.Fram['Flux'][t,:,:] = self.data['U'][t,:,58*kk:80*kk,76*kk]*Area_y[:,58*kk:80*kk,76*kk]
                self.Fram['FluxSum'][t] = np.nansum(np.nansum(self.Fram['Flux'][t,:,:]))
                self.Fram['FluxT'][t,:,:] = self.Fram['Flux'][t,:,:]*\
                                            self.data['T'][t,:,58*kk:80*kk,76*kk]*Area_y[:,58*kk:80*kk,76*kk]
                self.Fram['FluxSumT'][t] = np.nansum(np.nansum(self.Fram['FluxT'][t,:,:]))
                self.Fram['FluxS'][t,:,:] = self.Fram['Flux'][t,:,:]*\
                                            self.data['S'][t,:,58*kk:80*kk,76*kk]*Area_y[:,58*kk:80*kk,76*kk]
                self.Fram['FluxSumS'][t] = np.nansum(np.nansum(self.Fram['FluxS'][t,:,:]))                    
                self.totalFluxes['Fram'] = fluxTransport(self.Fram['Flux'])

                # Fram 1 fillign
                self.Fram1['Flux'][t,:,:] = self.data['U'][t,:,55*kk:85*kk,78*kk]*Area_y[:,55*kk:85*kk,78*kk]
                self.Fram1['FluxSum'][t] = np.nansum(np.nansum(self.Fram1['Flux'][t,:,:]))
                self.Fram1['FluxT'][t,:,:] = self.Fram1['Flux'][t,:,:]*\
                                            self.data['T'][t,:,55*kk:85*kk,78*kk]*Area_y[:,55*kk:85*kk,78*kk]
                self.Fram1['FluxSumT'][t] = np.nansum(np.nansum(self.Fram1['FluxT'][t,:,:]))
                self.Fram1['FluxS'][t,:,:] = self.Fram1['Flux'][t,:,:]*\
                                            self.data['S'][t,:,55*kk:85*kk,78*kk]*Area_y[:,55*kk:85*kk,78*kk]
                self.Fram1['FluxSumS'][t] = np.nansum(np.nansum(self.Fram1['FluxS'][t,:,:]))                    
                self.totalFluxes['Fram1'] = fluxTransport(self.Fram1['Flux'])
                
                # Fram 2 fillign
                self.Fram2['Flux'][t,:,:] = self.data['U'][t,:,60*kk:83*kk,72*kk]*Area_y[:,60*kk:83*kk,72*kk]
                self.Fram2['FluxSum'][t] = np.nansum(np.nansum(self.Fram['Flux'][t,:,:]))
                self.Fram2['FluxT'][t,:,:] = self.Fram2['Flux'][t,:,:]*\
                                            self.data['T'][t,:,60*kk:83*kk,72*kk]*Area_y[:,60*kk:83*kk,72*kk]
                self.Fram2['FluxSumT'][t] = np.nansum(np.nansum(self.Fram2['FluxT'][t,:,:]))
                self.Fram2['FluxS'][t,:,:] = self.Fram2['Flux'][t,:,:]*\
                                            self.data['S'][t,:,60*kk:83*kk,72*kk]*Area_y[:,60*kk:83*kk,72*kk]
                self.Fram2['FluxSumS'][t] = np.nansum(np.nansum(self.Fram2['FluxS'][t,:,:]))                    
                self.totalFluxes['Fram2'] = fluxTransport(self.Fram2['Flux'])                    
                
                # Bering fillign
                self.Bering['Flux'][t,:,:] = -self.data['U'][t,:,80*kk:89*kk,178*kk]*Area_y[:,80*kk:89*kk,178*kk]
                self.Bering['FluxSum'][t] = np.nansum(np.nansum(self.Bering['Flux'][t,:,:]))
                self.Bering['FluxT'][t,:,:] = self.Bering['Flux'][t,:,:]*\
                                            self.data['T'][t,:,80*kk:89*kk,178*kk]*Area_y[:,80*kk:89*kk,178*kk]
                self.Bering['FluxSumT'][t] = np.nansum(np.nansum(self.Bering['FluxT'][t,:,:]))
                self.Bering['FluxS'][t,:,:] = self.Bering['Flux'][t,:,:]*\
                                            self.data['S'][t,:,80*kk:89*kk,178*kk]*Area_y[:,80*kk:89*kk,178*kk]
                self.Bering['FluxSumS'][t] = np.nansum(np.nansum(self.Bering['FluxS'][t,:,:]))
                self.totalFluxes['Bering'] = fluxTransport(self.Bering['Flux'])
                
                # Davis fillign
                self.Davis['Flux'][t,:,:] = self.data['U'][t,:,113*kk:135*kk,75*kk]*Area_y[:,113*kk:135*kk,75*kk]
                self.Davis['FluxSum'][t] = np.nansum(np.nansum(self.Davis['Flux'][t,:,:]))
                self.Davis['FluxT'][t,:,:] = self.Davis['Flux'][t,:,:]*\
                                            self.data['T'][t,:,113*kk:135*kk,75*kk]*Area_y[:,113*kk:135*kk,75*kk]
                self.Davis['FluxSumT'][t] = np.nansum(np.nansum(self.Davis['FluxT'][t,:,:]))
                self.Davis['FluxS'][t,:,:] = self.Davis['Flux'][t,:,:]*\
                                            self.data['S'][t,:,113*kk:135*kk,75*kk]*Area_y[:,113*kk:135*kk,75*kk]
                self.Davis['FluxSumS'][t] = np.nansum(np.nansum(self.Davis['FluxS'][t,:,:]))
                self.totalFluxes['Davis'] = fluxTransport(self.Davis['Flux']) 
                
                # Davis1 fillign
                self.Davis1['Flux'][t,:,:] = self.data['V'][t,:,135*kk,52*kk:73*kk]*Area_x[:,135*kk,52*kk:73*kk]
                self.Davis1['FluxSum'][t] = np.nansum(np.nansum(self.Davis1['Flux'][t,:,:]))
                self.Davis1['FluxT'][t,:,:] = self.Davis1['Flux'][t,:,:]*\
                                            self.data['T'][t,:,135*kk,52*kk:73*kk]*Area_x[:,135*kk,52*kk:73*kk]
                self.Davis1['FluxSumT'][t] = np.nansum(np.nansum(self.Davis1['FluxT'][t,:,:]))
                self.Davis1['FluxS'][t,:,:] = self.Davis1['Flux'][t,:,:]*\
                                            self.data['S'][t,:,135*kk,52*kk:73*kk]*Area_x[:,135*kk,52*kk:73*kk]
                self.Davis1['FluxSumS'][t] = np.nansum(np.nansum(self.Davis1['FluxS'][t,:,:]))
                self.totalFluxes['Davis1'] = fluxTransport(self.Davis1['Flux']) 
                
                # Davis 2 fillign
                self.Davis2['Flux'][t,:,:] = self.data['U'][t,:,113*kk:135*kk,72*kk]*Area_y[:,113*kk:135*kk,72*kk]
                self.Davis2['FluxSum'][t] = np.nansum(np.nansum(self.Davis2['Flux'][t,:,:]))
                self.Davis2['FluxT'][t,:,:] = self.Davis2['Flux'][t,:,:]*\
                                            self.data['T'][t,:,113*kk:135*kk,72*kk]*Area_y[:,113*kk:135*kk,72*kk]
                self.Davis2['FluxSumT'][t] = np.nansum(np.nansum(self.Davis2['FluxT'][t,:,:]))
                self.Davis2['FluxS'][t,:,:] = self.Davis2['Flux'][t,:,:]*\
                                            self.data['S'][t,:,113*kk:135*kk,72*kk]*Area_y[:,113*kk:135*kk,72*kk]
                self.Davis2['FluxSumS'][t] = np.nansum(np.nansum(self.Davis2['FluxS'][t,:,:]))
                self.totalFluxes['Davis2'] = fluxTransport(self.Davis2['Flux']) 

                # Barents fillng
                self.Barents['Flux'][t,:,:] = np.concatenate((\
                                            self.data['V'][t,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                            self.data['U'][t,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=1)
                self.Barents['FluxSum'][t] = np.nansum(np.nansum(self.Barents['Flux'][t,:,:]))
                self.Barents['FluxT'][t,:,:] = self.Barents['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                            self.data['T'][t,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=1)
                self.Barents['FluxSumT'][t] = np.nansum(np.nansum(self.Barents['FluxT'][t,:,:]))
                self.Barents['FluxS'][t,:,:] = self.Barents['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,40*kk,53*kk:68*kk]*Area_x[:,40*kk,53*kk:68*kk],\
                                            self.data['T'][t,:,40*kk:58*kk,68*kk]*Area_y[:,40*kk:58*kk,68*kk]),axis=1)
                self.Barents['FluxSumS'][t] = np.nansum(np.nansum(self.Barents['FluxS'][t,:,:]))
                self.totalFluxes['Barents'] = fluxTransport(self.Barents['Flux'])

                # Barents 1 fillng
                y1 = 45; y2 = 58 ; x1 = 49 ; x2 = 66 # barents 1
                self.Barents1['Flux'][t,:,:] = np.concatenate((\
                                            self.data['V'][t,:,45*kk,49*kk:66*kk]*Area_x[:,45*kk,49*kk:66*kk],\
                                            self.data['U'][t,:,45*kk:58*kk,66*kk]*Area_y[:,45*kk:58*kk,66*kk]),axis=1)
                self.Barents1['FluxSum'][t] = np.nansum(np.nansum(self.Barents1['Flux'][t,:,:]))
                self.Barents1['FluxT'][t,:,:] = self.Barents1['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,45*kk,49*kk:66*kk]*Area_x[:,45*kk,49*kk:66*kk],\
                                            self.data['T'][t,:,45*kk:58*kk,66*kk]*Area_y[:,45*kk:58*kk,66*kk]),axis=1)
                self.Barents1['FluxSumT'][t] = np.nansum(np.nansum(self.Barents1['FluxT'][t,:,:]))
                self.Barents1['FluxS'][t,:,:] = self.Barents1['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,45*kk,49*kk:66*kk]*Area_x[:,45*kk,49*kk:66*kk],\
                                            self.data['T'][t,:,45*kk:58*kk,66*kk]*Area_y[:,45*kk:58*kk,66*kk]),axis=1)
                self.Barents1['FluxSumS'][t] = np.nansum(np.nansum(self.Barents1['FluxS'][t,:,:]))
                self.totalFluxes['Barents1'] = fluxTransport(self.Barents1['Flux'])
                
                # Barents 2 fillng
                y1 = 45; y2 = 58 ; x1 = 49 ; x2 = 66 # barents 1
                self.Barents2['Flux'][t,:,:] = np.concatenate((\
                                            self.data['V'][t,:,58*kk,49*kk:66*kk]*Area_x[:,58*kk,49*kk:66*kk],\
                                            self.data['U'][t,:,45*kk:58*kk,49*kk]*Area_y[:,45*kk:58*kk,49*kk]),axis=1)
                self.Barents2['FluxSum'][t] = np.nansum(np.nansum(self.Barents2['Flux'][t,:,:]))
                self.Barents2['FluxT'][t,:,:] = self.Barents1['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,58*kk,49*kk:66*kk]*Area_x[:,58*kk,49*kk:66*kk],\
                                            self.data['T'][t,:,45*kk:58*kk,49*kk]*Area_y[:,45*kk:58*kk,49*kk]),axis=1)
                self.Barents2['FluxSumT'][t] = np.nansum(np.nansum(self.Barents2['FluxT'][t,:,:]))
                self.Barents2['FluxS'][t,:,:] = self.Barents1['Flux'][t,:,:]*np.concatenate((\
                                            self.data['T'][t,:,58*kk,49*kk:66*kk]*Area_x[:,58*kk,49*kk:66*kk],\
                                            self.data['T'][t,:,45*kk:58*kk,49*kk]*Area_y[:,45*kk:58*kk,49*kk]),axis=1)
                self.Barents2['FluxSumS'][t] = np.nansum(np.nansum(self.Barents2['FluxS'][t,:,:]))
                self.totalFluxes['Barents2'] = fluxTransport(self.Barents2['Flux'])
                
                # Denmark filling
                self.Denmark['Flux'][t,:,:] = self.data['V'][t,:,100*kk,37*kk:48*kk]*Area_x[:,100*kk,37*kk:48*kk]
                self.Denmark['FluxSum'][t] = np.nansum(np.nansum(self.Denmark['Flux'][t,:,:]))
                self.Denmark['FluxT'][t,:,:] = self.Denmark['Flux'][t,:,:]*\
                                            self.data['T'][t,:,100*kk,37*kk:48*kk]*Area_x[:,100*kk,37*kk:48*kk]
                self.Denmark['FluxSumT'][t] = np.nansum(np.nansum(self.Denmark['FluxT'][t,:,:]))
                self.Denmark['FluxS'][t,:,:] = self.Denmark['Flux'][t,:,:]*\
                                            self.data['S'][t,:,100*kk,37*kk:48*kk]*Area_x[:,100*kk,37*kk:48*kk]
                self.Denmark['FluxSumS'][t] = np.nansum(np.nansum(self.Denmark['FluxS'][t,:,:]))
                self.totalFluxes['Denmark'] = fluxTransport(self.Denmark['Flux'])

                # Norwey-Iceland strait filling
                self.Norwice['Flux'][t,:,:] = np.concatenate((\
                                        self.data['V'][t,:,60*kk:95*kk,15*kk]*Area_y[:,60*kk:95*kk,15*kk],\
                                        self.data['U'][t,:,95*kk,15*kk:30*kk]*Area_x[:,95*kk,15*kk:30*kk]),axis=1)
                self.Norwice['FluxSum'][t] = np.nansum(np.nansum(self.Barents['Flux'][t,:,:]))
                self.Norwice['FluxT'][t,:,:] = self.Norwice['Flux'][t,:,:]*np.concatenate((\
                                        self.data['T'][t,:,60*kk:95*kk,15*kk]*Area_y[:,60*kk:95*kk,15*kk],\
                                        self.data['T'][t,:,95*kk,15*kk:30*kk]*Area_x[:,95*kk,15*kk:30*kk]),axis=1)
                self.Norwice['FluxSumT'][t] = np.nansum(np.nansum(self.Barents['FluxT'][t,:,:]))
                self.Norwice['FluxS'][t,:,:] = self.Norwice['Flux'][t,:,:]*np.concatenate((\
                                        self.data['T'][t,:,60*kk:95*kk,15*kk]*Area_x[:,60*kk:95*kk,15*kk],\
                                        self.data['T'][t,:,95*kk,15*kk:30*kk]*Area_y[:,95*kk,15*kk:30*kk]),axis=1)
                self.Norwice['FluxSumS'][t] = np.nansum(np.nansum(self.Barents['FluxS'][t,:,:]))
                self.totalFluxes['Norwice'] = fluxTransport(self.Norwice['Flux'])


            n1 = 10
            n2 = 24
            self.FramPart = {'FluxTop' : np.nansum(np.nansum(self.Fram['Flux'][:,0:n1,:],axis=2),axis=1),\
                        'FluxMid' : np.nansum(np.nansum(self.Fram['Flux'][:,n1:n2,:],axis=2),axis=1),\
                        'FluxBot' : np.nansum(np.nansum(self.Fram['Flux'][:,n2:50,:],axis=2),axis=1),\
                        'FluxTopS' : np.nansum(np.nansum(self.Fram['FluxS'][:,0:n1,:],axis=2),axis=1),\
                        'FluxMidS' : np.nansum(np.nansum(self.Fram['FluxS'][:,n1:n2,:],axis=2),axis=1),\
                        'FluxBotS' : np.nansum(np.nansum(self.Fram['FluxS'][:,n2:50,:],axis=2),axis=1),\
                        'FluxTopT' : np.nansum(np.nansum(self.Fram['FluxT'][:,0:n1,:],axis=2),axis=1),\
                        'FluxMidT' : np.nansum(np.nansum(self.Fram['FluxT'][:,n1:n2,:],axis=2),axis=1),\
                        'FluxBotT' : np.nansum(np.nansum(self.Fram['FluxT'][:,n2:50,:],axis=2),axis=1)}            

            n1 = 10
            n2 = 20
            self.DenmarkPart = {'FluxTop' : np.nansum(np.nansum(self.Denmark['Flux'][:,0:n1,:],axis=2),axis=1),\
                        'FluxMid' : np.nansum(np.nansum(self.Denmark['Flux'][:,n1:n2,:],axis=2),axis=1),\
                        'FluxBot' : np.nansum(np.nansum(self.Denmark['Flux'][:,n2:50,:],axis=2),axis=1),\
                        'FluxTopS' : np.nansum(np.nansum(self.Denmark['FluxS'][:,0:n1,:],axis=2),axis=1),\
                        'FluxMidS' : np.nansum(np.nansum(self.Denmark['FluxS'][:,n1:n2,:],axis=2),axis=1),\
                        'FluxBotS' : np.nansum(np.nansum(self.Denmark['FluxS'][:,n2:50,:],axis=2),axis=1),\
                        'FluxTopT' : np.nansum(np.nansum(self.Denmark['FluxT'][:,0:n1,:],axis=2),axis=1),\
                        'FluxMidT' : np.nansum(np.nansum(self.Denmark['FluxT'][:,n1:n2,:],axis=2),axis=1),\
                        'FluxBotT' : np.nansum(np.nansum(self.Denmark['FluxT'][:,n2:50,:],axis=2),axis=1)}
            
            


    def baroCalc(self):
        self.psi = baro_stream(self.data['U'])
        self.psi_mean = np.nanmean(self.psi,axis = 1)
        self.psi_mean = np.nanmean(self.psi_mean,axis = 1)
        self.psi_max = np.nanmax(self.psi,axis=1)
        self.psi_max = np.nanmax(self.psi_max,axis=1)
        self.psi_min = np.nanmin(self.psi,axis=1)
        self.psi_min = np.nanmin(self.psi_min,axis=1)
        self.psi_ave = np.nanmean(self.psi,axis=0)

    def topoCalc(self):
        topo,topo_opposit = Topostrophy(self.data['U'],self.data['V'])
        self.topo = topo
        self.topo_opposit = topo_opposit

    def title(self,title):
        self.title = title
    
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
        self.mxldepth[:,0,self.hfacc[0,:,:]==0]=np.nan
