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
import matplotlib.pyplot as plt 

sys.path.append('/noc/users/am8e13/PyNGL-1.4.0/lib/python2.7/site-packages/PyNGL/')
import Ngl
sys.path.append('/noc/users/am8e13/Python/')
import komod
sys.path.append('/noc/users/am8e13/Python/PyNIO-1.4.0/')
import Nio

def plot_obcs(data1,data2,data3,data4,data5,loc,var,vmin,vmax,row=1,col=5):
    vmin = vmin
    vmax = vmax
    tempbounds = range(vmin,vmax,2)
    var = var
    if var == 'S':
        if loc == 'N':
            field1 = np.array(data1.S[:,191,:])
            field2 = np.array(data2.S[:,191,:])
            field3 = np.array(data3.S[:,191,:])
            field4 = np.array(data4.S[:,191,:])
            field5 = np.array(data5.S[:,191,:])
        elif loc == 'E':
            field1 = np.array(data1.S[:,:,209])
            field2 = np.array(data2.S[:,:,209])
            field3 = np.array(data3.S[:,:,209])
            field4 = np.array(data4.S[:,:,209])
            field5 = np.array(data5.S[:,:,209])
        elif loc == 'W':
            field1 = np.array(data1.S[:,:,0])
            field2 = np.array(data2.S[:,:,0])
            field3 = np.array(data3.S[:,:,0])
            field4 = np.array(data4.S[:,:,0])
            field5 = np.array(data5.S[:,:,0])
    elif var == 'T':
        if loc == 'N':
            field1 = np.array(data1.T[:,191,:])
            field2 = np.array(data2.T[:,191,:])
            field3 = np.array(data3.T[:,191,:])
            field4 = np.array(data4.T[:,191,:])
            field5 = np.array(data5.T[:,191,:])
        elif loc == 'E':
            field1 = np.array(data1.T[:,:,209])
            field2 = np.array(data2.T[:,:,209])
            field3 = np.array(data3.T[:,:,209])
            field4 = np.array(data4.T[:,:,209])
            field5 = np.array(data5.T[:,:,209])
        elif loc == 'W':
            field1 = np.array(data1.T[:,:,0])
            field2 = np.array(data2.T[:,:,0])
            field3 = np.array(data3.T[:,:,0])
            field4 = np.array(data4.T[:,:,0])
            field5 = np.array(data5.T[:,:,0])        
        
    cs = plt.contourf(field4,vmin=vmin,vmax=vmax,level=tempbounds,extend='both')
    plt.colorbar(cs)
    plt.close()

    fig, ((ax1, ax2 , ax3 ,ax4 , ax5)) = plt.subplots(1,5)
    ax1 = plt.subplot(1,5,1)
    im1 = plt.contourf(field1,vmin=vmin,vmax=vmax,level=tempbounds,extend='both')
    ax1.set_title(data1.title+' '+var)
    ax1.title.set_fontsize('14')

    ax2 = plt.subplot(1,5,2)
    im2 = plt.contourf(field2,vmin=vmin,vmax=vmax,level=tempbounds,extend='both')
    ax2.set_title(data2.title+' '+var)
    ax2.title.set_fontsize('14')

    ax3 = plt.subplot(1,5,3)
    im3 = plt.contourf(field3,vmin=vmin,vmax=vmax,level=tempbounds,extend='both')
    ax3.set_title(data3.title+' '+var)
    ax3.title.set_fontsize('14')
    
    ax4 = plt.subplot(1,5,4)
    im4 = plt.contourf(field4,vmin=vmin,vmax=vmax,level=tempbounds,extend='both')
    ax4.set_title(data4.title+' '+var)
    ax4.title.set_fontsize('14')
    
    ax5 = plt.subplot(1,5,5)
    im5 = plt.contourf(field5,vmin=vmin,vmax=vmax,level=tempbounds,extend='both')
    ax5.set_title(data5.title+' '+var)
    ax5.title.set_fontsize('14')


    cbar_ax = fig.add_axes([1.9, 0.4, 0.045, .6])
    cbar = plt.colorbar(cs, cax=cbar_ax,)
    cbar.ax.set_ylabel('C')

    fig.subplots_adjust(right=1.7,top=1.)
