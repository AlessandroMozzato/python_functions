#! /usr/bin/env python                                                                                                              
# -*- coding: utf-8 -*-  

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def plot_field(field,years,vimin,vimax,cmap,setbad,unity,row=2,col=3,origin="lower"):
    # This function plots a 2D field, the field is meant to have NaNs on the land place
    # vimin is the minimun, vimax is the maximum, setbad is the NaN color, unity is the unity in the colorbar
    # cmap is a colomap
    fig, axes = plt.subplots(nrows=row, ncols=col)
    i=0
    vimin = vimin
    vimax = vimax
    
    for ax in axes.flat:        
        ax.set_title("T = "+str(np.round(years[i],2))+" years")
        masked_array = np.ma.array(field[i,:,:], mask=np.isnan(field[1,:,:]))
        cmap= cmap
        cmap.set_bad(setbad,1.)
        ca = ax.imshow(masked_array,vmin = vimin, vmax = vimax, interpolation='nearest',cmap = cmap ,origin=origin)
        cbar = fig.colorbar(ca , ax=ax, )
        cbar.ax.set_ylabel(unity)       
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
            item.set_fontsize(14)
        i=i+1
    fig.subplots_adjust(right=2.4,top=2)

def plot_field_gif(field,years,vimin,vimax,cmap,setbad,unity,title,row=1,col=1,origin='upper'):
    # This function produces a gif animation of the field                                                                         
    # All the parameters are the same as the previous                                                                                                          
    vimin = vimin
    vimax = vimax
    fig, ax = plt.subplots(nrows=row, ncols=col)
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

def plot_dynmon(monitor_list,var_list,row=3,col=3):
    # this function is meant to plot 2D monitor stats
    file2read1 = netcdf.NetCDFFile("/scratch/general/am8e13/results36km/grid.nc",'r')
    Z = file2read1.variables['Zp1']
    Z = Z[0:50]*1
    fig, axes = plt.subplots(row,col)
    ax_ind = 0
    for var in var_list:  
        for mon in monitor_list:
            ca = axes.flat[ax_ind].contourf(monitor_list[mon].dataDyn['time_lv_years'],Z,\
                                            monitor_list[mon].dataDyn[var].T,40,vmin=vimins[var],vmax=vimaxs[var])
            axes.flat[ax_ind].set_title(monitor_list[mon].title)
            cbar = fig.colorbar(ca , ax=axes.flat[ax_ind],boundaries=[1,1])
            #cbar.set_clim(vimins[var],vimaxs[var])
            ax_ind += 1
        if ((ax_ind//col) == 0):
            axes.flat[ax_ind].set_ylabel(var)
            
    fig.subplots_adjust(right=2.4,top=2.0)

def plot_monitor(monitor_list,var_list,row=2,col=2):
    # this function is meant to plot monitor variables                                                                                                          
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
                 'seaice_heff_mean' : 'm', 'time_seconds' : 's' , \
                 'time_years' : 'Years'}
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

    fig, axes = plt.subplots(row,col)
    ax_ind = 0
    for var in var_list:
        for mon in monitor_list:
            axes.flat[ax_ind].plot(monitor_list[mon].data[time[var]],monitor_list[mon].data[var],monitor_list[mon].color)
            if ax_ind == 0:
                print monitor_list[mon].title, monitor_list[mon].color
        axes.flat[ax_ind].set_ylabel(unity[var])
        axes.flat[ax_ind].set_title(titles[var])
        ax_ind += 1

    fig.subplots_adjust(right=2.4,top=2.4)

def plot_flux(var_list,monitor_list,location):
    fig, axes = plt.subplots(4,3)
    ax_ind = 0
    for var in var_list:  
        for mon in monitor_list:
            axes.flat[ax_ind].plot(monitor_list[mon].years,monitor_list[mon].Fram[var])
            axes.flat[ax_ind+3].plot(monitor_list[mon].years,monitor_list[mon].Barents[var])
            axes.flat[ax_ind+6].plot(monitor_list[mon].years,monitor_list[mon].Denmark[var])
            axes.flat[ax_ind+9].plot(monitor_list[mon].years,monitor_list[mon].Norwice[var])

        axes.flat[ax_ind].set_title(var)
        axes.flat[ax_ind+3].set_title(var)
        ax_ind += 1
        
    fig.subplots_adjust(right=2.4,top=2.4)
