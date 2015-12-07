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

class Obcs:
    def __init__(self):

        self.data = {'Ev' : [] , 'Wv' : [] , 'Nv' : [] , \
                     'Eu' : [] , 'Wu' : [] , 'Nu' : [] , \
                     'Et' : [] , 'Wt' : [] , 'Nt' : [] , \
                     'Es' : [] , 'Ws' : [] , 'Ns' : [] }
        self.data_av = {'Ev' : [] , 'Wv' : [] , 'Nv' : [] , \
                     'Eu' : [] , 'Wu' : [] , 'Nu' : [] , \
                     'Et' : [] , 'Wt' : [] , 'Nt' : [] , \
                     'Es' : [] , 'Ws' : [] , 'Ns' : [] }
        self.mean = {'Ev' : [] , 'Wv' : [] , 'Nv' : [] , \
                     'Eu' : [] , 'Wu' : [] , 'Nu' : [] , \
                     'Et' : [] , 'Wt' : [] , 'Nt' : [] , \
                     'Es' : [] , 'Ws' : [] , 'Ns' : [] }
        self.title = 'OBCS'

    def ReadData(self,path,tdim,endingstr,res):
        obcs_data = ['Ev' , 'Wv' , 'Nv' , 'Eu' , 'Wu' , 'Nu' , 'Et' , 'Wt' , 'Nt' , 'Es' , 'Ws' , 'Ns' ]
        
        if res == 18:
            xdim = [420,1]
            ydim = [384,1]
            file_res = "_arctic_420x384."
            grid = '/scratch/general/am8e13/results18km/'
            kk = 2
        elif res == 36:
            xdim = [210,1]
            ydim = [192,1]
            file_res = "_arctic_210x192."
            grid = '/scratch/general/am8e13/results36km/'
            kk = 1
        elif res == 9:
            xdim = [840,1]
            ydim = [768,1]
            file_res = "_arctic_840x768."
            grid = '/scratch/general/am8e13/results9km/'
            kk = 4
        zdim = 50
        
        for var in obcs_data:
            if var == 'Nu' or var == 'Nv' or var == 'Ns' or var == 'Nt':
                xdim_i = 0
                ydim_i = 1
            elif var == 'Wu' or var == 'Wv' or var == 'Ws' or var == 'Wt' or \
                    var == 'Eu' or var == 'Ev' or var == 'Es' or var == 'Et':
                xdim_i = 1
                ydim_i = 0
            
            if var == 'Es' or var == 'Ns' or var == 'Ws':
                stable = 'stable'
                ending = ''
            elif var == 'Et' or var == 'Nt' or var == 'Wt':
                stable = 'stable'
                ending = endingstr
            else:
                stable = 'bin'
                ending = ''
                
            name = path+'/OB'+str(var)+file_res+str(stable)+'_mean'+str(ending)
            data = komod.mitbin(name,xdim=xdim[xdim_i],ydim=ydim[ydim_i],zdim=50,tdim=tdim,datatype='float32')
            data_av_temp = np.mean(np.mean(np.mean(data,axis=3),axis=2),1)    
            self.data[var] = data #_av_temp
            #self.data_av[var] = data_av_temp
            data_climy_33_years = []
            for year in range(10):
                data_climy_33_years = np.concatenate([data_climy_33_years, data_av_temp])                
            self.mean[var] = np.mean(self.data[var],axis = 0)
            self.data_av[var] = data_climy_33_years 
        
        file2read = netcdf.NetCDFFile(grid+'grid.nc','r')
        bathy=file2read.variables['HFacC']
        bathy=bathy[:]*1
        
        self.T = np.ndarray([50,192*kk,210*kk])
        self.T[:,:,209*kk] = self.mean['Et'][:,0,:]
        self.T[:,:,0] = self.mean['Wt'][:,0,:]
        self.T[:,191*kk,:]= self.mean['Nt'][:,:,0]
        self.T[bathy==0] = np.nan
        self.S = np.ndarray([50,192*kk,210*kk])
        self.S[:,:,209*kk] = self.mean['Es'][:,0,:]
        self.S[:,:,0] = self.mean['Ws'][:,0,:]
        self.S[:,191*kk,:]= self.mean['Ns'][:,:,0]
        self.S[bathy==0] = np.nan
