#! /usr/bin/env python                                                                                                              
# -*- coding: utf-8 -*-  

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def plot_field(field,years,vimin,vimax,cmap,setbad,unity):
    # This function plots a 2D field, the field is meant to have NaNs on the land place
    # vimin is the minimun, vimax is the maximum, setbad is the NaN color, unity is the unity in the colorbar
    # cmap is a colomap
    fig, axes = plt.subplots(nrows=2, ncols=3)
    i=0
    vimin = vimin
    vimax = vimax
    
    for ax in axes.flat:        
        ax.set_title("T = "+str(np.round(years[i],2))+" years")
        masked_array = np.ma.array(field[i,:,:], mask=np.isnan(field[1,:,:]))
        cmap= cmap
        cmap.set_bad(setbad,1.)
        ca = ax.imshow(masked_array,vmin = vimin, vmax = vimax, interpolation='nearest',cmap = cmap ,origin="lower")
        cbar = fig.colorbar(ca , ax=ax, )
        cbar.ax.set_ylabel(unity)       
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
            item.set_fontsize(14)
        i=i+1
    fig.subplots_adjust(right=2.4,top=2)

def plot_field_gif(field,vimin,vimax,cmap,setbad,unity,title,years,origin='upper'):
    # This function produces a gif animation of the field                                                                                                      
    # All the parameters are the same as the previous                                                                                                          
    vimin = vimin
    vimax = vimax
    fig, ax = plt.subplots(nrows=1, ncols=1)
    for i in range(len(years)):
        ax.set_title("T = "+str(np.round(years[i],2))+" years")
        masked_array = np.ma.array(field[i,:,:], mask=np.isnan(field[1,:,:]))
        cmap= cmap
        cmap.set_bad(setbad,1.)
        ca = ax.imshow(masked_array,vmin = vimin, vmax = vimax, interpolation='nearest',cmap = cmap,aspect='auto',origin=origin)
        if i == 0:
            cbar = fig.colorbar(ca , ax=ax, )
            cbar.ax.set_ylabel(unity)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
            item.set_fontsize(14)

        if i < 10:
            zeros = '000'
        elif i < 100:
            zeros = '00'
        else:
            zeros = '0'
        sout = '/scratch/general/am8e13/gifs/'
        fig.savefig(sout+str(title)+zeros+str(i)+'.png',dpi=300)
        clf
    os.system("convert -delay 1 -dispose Background +page " + str(sout) \
          + "/*.png -loop 0 " + str(sout) + "/animation_"+title+".gif")
    os.system("rm "+str(sout)+"*png")
    sout = '/scratch/general/am8e13/gifs/'

def plot_dynmon(var_list,monitor_list,row=3,col=3):
    # this function is meant to plot 2D monitor stats
    fig, axes = plt.subplots(row,col)
    ax_ind = 0
    for var in var_list:  
        for mon in monitor_list:
            ca = axes.flat[ax_ind].contourf(monitor_list[mon].dataDyn['time_lv_years'],Z,monitor_list[mon].dataDyn[var].T,40,vmin=vimins[var],vmax=vimaxs[var])
            axes.flat[ax_ind].set_title(monitor_list[mon].title)
            cbar = fig.colorbar(ca , ax=axes.flat[ax_ind],boundaries=[1,1])
            #cbar.set_clim(vimins[var],vimaxs[var])
            ax_ind += 1
        if ((ax_ind//col) == 0):
            axes.flat[ax_ind].set_ylabel(var)
            
            
    fig.subplots_adjust(right=2.4,top=1.7)

def plot_monitor(var_list,monitor_list,row=2,col=2):
    # this function is meant to plot monitor variables
    fig, axes = plt.subplots(row,col)
    ax_ind = 0
    for var in var_list:  
        for mon in monitor_list:
            axes.flat[ax_ind].plot(monitor_list[mon].data['time_years'],monitor_list[mon].data[var],monitor_list[mon].color)
            if ax_ind == 0:
                print monitor_list[mon].title, monitor_list[mon].color
        axes.flat[ax_ind].set_ylabel(unity[var])
        axes.flat[ax_ind].set_title(titles[var])
        ax_ind += 1
        
    fig.subplots_adjust(right=2.4,top=2.4)

def plot_flux(var_list,monitor_list,location):
    fig, axes = plt.subplots(2,3)
    ax_ind = 0
    for var in var_list:  
        for mon in monitor_list:
            axes.flat[ax_ind].plot(monitor_list[mon].years,monitor_list[mon].Fram[var])
            axes.flat[ax_ind+3].plot(monitor_list[mon].years,monitor_list[mon].Barents[var])
            #if ax_ind == 0:
            #    print monitor_list[mon].title, monitor_list[mon].color
        #axes.flat[ax_ind].set_ylabel(unity[var])
        axes.flat[ax_ind].set_title(var)
        axes.flat[ax_ind+3].set_title(var)
        ax_ind += 1
        
    fig.subplots_adjust(right=2.4,top=2.4)
