#! /usr/bin/env python                                                                                                                 
# -*- coding: utf-8 -*-  

from scipy.io import netcdf
import numpy as np
from pylab import clf, plot, show, floor, ceil, imshow
import matplotlib
import matplotlib.pyplot as plt 
import os
import csv
import sys
import glob
from netCDF4 import Dataset

from IPython.html.widgets import interact, interactive
from IPython.display import clear_output, display, HTML

def interact(field,cmap1,vimin1,vimin2,vimax1,vimax2,time1,time2,z1,z2):
    def plot_field2(vimin=10,vimax=-10,time=0,Z=0):
        # This function plots a 2D field, the field is meant to have NaNs on the land place                                                
        # vimin is the minimun, vimax is the maximum, setbad is the NaN color, unity is the unity in the colorbar                          
        # cmap is a colomap                                                                                                               
        fig,axes = plt.subplots(1,1)
        if len(field.shape) == 3:
            masked_array = np.ma.array(field[time,:,:], mask=np.isnan(field[time,:,:]))
        elif len(field.shape) == 4:
            masked_array = np.ma.array(field[time,Z,:,:], mask=np.isnan(field[time,Z,:,:]))
        elif len(field.shape) == 2:
            masked_array = np.ma.array(field[:,:], mask=np.isnan(field[:,:]))
        cmap= cmap1
        cmap.set_bad('grey',1.)
        ca = imshow(masked_array,vmin = vimin, vmax = vimax, interpolation='nearest',cmap = cmap ,origin="upper")
        cbar = fig.colorbar(ca )
        fig.subplots_adjust(right=2.4,top=2)
        return plot_field2
        
    w = interactive(plot_field2,vimin=(vimin1,vimin2),vimax=(vimax1,vimax2),time=(time1,time2),Z=(z1,z2))
    display(w)
