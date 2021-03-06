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

import random
import calendar

sys.path.append('/noc/users/am8e13/PyNGL-1.4.0/lib/python2.7/site-packages/PyNGL/')
import Ngl
sys.path.append('/noc/users/am8e13/Python/')
import komod
sys.path.append('/noc/users/am8e13/Python/PyNIO-1.4.0/')
import Nio

def read_forcing(path,dataset,init_year,tot_years):
    # This function is meant to read and store the forcing fields
    # Load adata files
    xdim = 320
    ydim = 160
    tdim = [365*4,366*4]

    data = {'v10m': np.zeros(360*tot_years) , 'u10m' : np.zeros(360*tot_years), \
            'rain': np.zeros(360*tot_years) , 'dlw' : np.zeros(360*tot_years) , \
            'dsw' : np.zeros(360*tot_years) , 'tmp2m_degC' : np.zeros(360*tot_years) , \
            'spfh2m' : np.zeros(360*tot_years), 'rain': np.zeros(360*tot_years)}

    data_all = {'v10m': np.zeros(360*tot_years*4) , 'u10m' : np.zeros(360*tot_years*4), \
                'rain': np.zeros(360*tot_years*4) , 'dlw' : np.zeros(360*tot_years*4) , \
                'dsw' : np.zeros(360*tot_years*4) , 'tmp2m_degC' : np.zeros(360*tot_years*4) , \
                'spfh2m' : np.zeros(360*tot_years*4), 'rain': np.zeros(360*tot_years*4)}


    data_ave = {'v10m': np.zeros(360) , 'u10m' : np.zeros(360), 'rain': np.zeros(360), \
                'dlw' : np.zeros(360) , 'dsw' : np.zeros(360) , 'tmp2m_degC' : np.zeros(360) , \
                'spfh2m' : np.zeros(360), 'rain': np.zeros(360*33),}

    data_ave_all = {'v10m': np.zeros(360*4) , 'u10m' : np.zeros(360*4), 'rain': np.zeros(360*4), \
                'dlw' : np.zeros(360*4) , 'dsw' : np.zeros(360*4) , 'tmp2m_degC' : np.zeros(360*4) , \
                'spfh2m' : np.zeros(360*4), 'rain': np.zeros(360*4),}

    for var in data:
        data_av = []
        #print "Now reading:"+str(var)
        for year in range(tot_years):
            if calendar.isleap(init_year+year):
                tdim_ly = 1
                n_pop = 6
            else:
                tdim_ly = 0
                n_pop = 5

            name = path+dataset+str(var)+'_'+str(init_year+year)
            data_read = komod.mitbin(name,xdim=xdim,ydim=ydim,zdim=1,tdim=tdim[tdim_ly],datatype='float32')
            data_av_temp = np.mean(np.mean(data_read,axis = 3),axis = 2).squeeze(axis = 1)        
            to_pop = random.sample(range(len(data_av_temp)), 4*n_pop)
            data_av_temp = np.delete(data_av_temp,to_pop,0) 
            data_all[var]=data_av_temp
            # Daily average for plotting purpose
            data_av_temp_daily = np.zeros(360)
            for day in range(360):
                data_av_temp_daily[day] = np.mean(data_av_temp[day*4 : day*4 +4])
                
            data_av = np.concatenate([data_av,data_av_temp_daily])
            data[var] = data_av

        name = path+dataset+str(var)+'_average'
        data_read = komod.mitbin(name,xdim=xdim,ydim=ydim,zdim=1,tdim=360*4,datatype='float32')
        data_climy_av = np.mean(np.mean(data_read,axis = 3),axis = 2).squeeze(axis = 1)
        data_ave_all[var]=data_climy_av  
        data_climy_av_daily = np.zeros(360)
        for day in range(360):
            data_climy_av_daily[day] = np.mean(data_climy_av[day*4 : day*4 +4])

        data_climy_33_years = []
        for year in range(tot_years):
            data_climy_33_years = np.concatenate([data_climy_33_years, data_climy_av_daily])

        data_ave[var] = data_climy_33_years
        
    clim = {'data' : data , 'data_all' : data_all , 'data_ave' : data_ave , 'data_ave_all' : data_ave_all}
    
    print 'read '+dataset
        
    return clim

def read_core(path):
    # This function is meant to read and store the forcing fields
    # Load adata files
    xdim = 320
    ydim = 160
    
    data_ave = {'v10m': np.zeros(366*4) , 'u10m' : np.zeros(366*4), 'rain': np.zeros(366), \
                'dlw' : np.zeros(366) , 'dsw' : np.zeros(366) , 'tmp2m_degC' : np.zeros(366*4) , \
                'spfh2m' : np.zeros(366*4), 'rain': np.zeros(12),}
    tdim_v = {'v10m': 366*4 , 'u10m' : 366*4, 'rain': 366, 'dlw' : 366 , 'dsw' : 366 , 'tmp2m_degC' : 366*4 , \
                'spfh2m' : 366*4, 'rain': 12,}
    
    time_core = {'v10m': 366*4 , 'u10m' : 366*4, 'rain': 366, 'dlw' : 366 , 'dsw' : 366 , 'tmp2m_degC' : 366*4 , \
                'spfh2m' : 366*4, 'rain': 12,}
    for var in data_ave: 
        name = path+'CORE2_'+str(var)
        data_read = komod.mitbin(name,xdim=192,ydim=94,tdim=tdim_v[var],datatype='float32')
        data_av = np.mean(np.mean(data_read,axis = 3),axis = 2).squeeze(axis = 1)
        
        if var == 'dlw' or var == 'dsw':
            data_ave[var] = -data_av
        else: 
            data_ave[var] = data_av
        
        time_core[var] = np.array(range(len(data_av)))/float(tdim_v[var])

    clim = {'data_ave' : data_ave}  
    print 'read CORE2'
    return clim,time_core

def dataset_unity():
    titles = {'v10m': 'Meridional wind' , 'u10m' : 'Zonal wind', 'rain_new33': 'Precipitation new33', \
            'dlw' : 'Long wave radiation' , 'dsw' : 'Short wave radiation' , 'tmp2m_degC' : 'Temperature at 2m' , \
            'spfh2m' : 'Specific humidity at 2m', 'rain': 'Precipitation'}
    unity = {'v10m': r"$m/s$" , 'u10m' : r"$m/s$", 'rain_new33': r"$mm/day$", \
            'dlw' : r"$W/m^2$" , 'dsw' : r"$W/m^2$" , 'tmp2m_degC' : r"$^{\circ}C$" , \
            'spfh2m' : r"$kg/kg", 'rain': r"$mm/day$"}
    return titles,unity
titles, unity = dataset_unity()
