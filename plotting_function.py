#! /usr/bin/env python                                                                                                              
# -*- coding: utf-8 -*-  

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def runningMeanFast(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

def reg_titles():
    region_titles = {0 : 'Global' , 1 : 'Arctic' , 2 : 'Nord Seas' , 3 : 'North Atl'}
    return region_titles

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

def plot_lv(monitor_list,var,reg,vimin,vimax,cmap,row=4):
    unit_titles = { 'theta_lv_mean' : 'T' , 'theta_lv_min' : 'T' , 'theta_lv_max' : 'T' , 
                'salt_lv_mean' : 'S' , 'salt_lv_min' : 'S' , 'salt_lv_max' : 'S' ,
                'rho_lv_mean' : 'Rho' , 'rho_lv_min' : 'Rho' , 'rho_lv_max' : 'Rho' , 
                'uvel_lv_mean' : 'U' ,'uvel_lv_min' : 'U' ,'uvel_lv_max' : 'U' , 
                'vvel_lv_mean' : 'V' ,'vvel_lv_min' : 'V' ,'vvel_lv_max' : 'V' , 
                 'ke_lv_mean' : 'KE' ,}
    reg_title = { 0 : 'Global' , 1 : 'Arctic' , 2 : 'Nordic Seas' , 3 : 'North Atl'}
    
    levels = { 'theta_lv_mean' : [-1,0,1,2,3,5,7] ,\
              'theta_lv_min' :  [-1,0,1,2,3,5,7], \
              'theta_lv_max' : [-1,0,1,2,3,5,7] , \
              'salt_lv_mean' : [33,34,35,37,40,45,50,60,70] ,\
              'salt_lv_min' : [33,34,35,37,40,45,50,60,70],\
              'salt_lv_max' : [33,34,35,37,40,45,50,60,70],\
              'rho_lv_mean' : [1026,1028,1030,1032,1034,1036,1038],\
              'rho_lv_min' : [1026,1028,1030,1032,1034,1036,1038],\
              'rho_lv_max' : [1026,1028,1030,1032,1034,1036,1038],\
              'ke_lv_mean' : [0, 0.01, 0.02, 0.03,0.04, 0.05,],\
              'uvel_lv_mean' : [-0.01, 0, 0.01, 0.01],\
              'uvel_lv_min' : [-0.01, 0, 0.01, 0.01],\
              'uvel_lv_max' : [-0.01, 0, 0.01, 0.01],\
              'vvel_lv_mean' : [-0.01, 0, 0.01, 0.01],\
              'vvel_lv_min' : [-0.01, 0, 0.01, 0.01],\
              'vvel_lv_max' : [-0.01, 0, 0.01, 0.01]
             }
    
    col = len(monitor_list)
     
    fig, axes = plt.subplots(nrows=row, ncols=col, sharex=True, sharey=True)
    ax_ind = 0
    ticks = np.linspace(vimin,vimax,5)
    for mon in monitor_list:
        reg = 0
        regz = 0
        ca = axes.flat[ax_ind].pcolor(monitor_list[mon].dataDyn['time_lv_years'],monitor_list[mon].Z,\
                    monitor_list[mon].dataDyn[var][:monitor_list[mon].dataDyn['time_lv_years'].shape[0],reg,:].T,\
                    vmin = vimin, vmax = vimax, cmap = cmap,)        
        cbar = fig.colorbar(ca , ax=axes.flat[ax_ind], ticks=ticks )
        axes.flat[ax_ind].set_title(unit_titles[var]+' '+reg_title[reg]+' '+monitor_list[mon].title)
        
        reg = 1
        regz = 46
        ca = axes.flat[ax_ind + col].pcolor(monitor_list[mon].dataDyn['time_lv_years'],monitor_list[mon].Z[0:regz],\
                    monitor_list[mon].dataDyn[var][:monitor_list[mon].dataDyn['time_lv_years'].shape[0],reg,0:regz].T,\
                    vmin = vimin, vmax = vimax, cmap = cmap,)
        cbar = fig.colorbar(ca , ax=axes.flat[ax_ind + col], ticks=ticks)
        axes.flat[ax_ind + col].set_title(unit_titles[var]+' '+reg_title[reg]+' '+monitor_list[mon].title)
        
        reg = 2
        regz = 45
        ca = axes.flat[ax_ind + col*2].pcolor(monitor_list[mon].dataDyn['time_lv_years'],monitor_list[mon].Z[0:regz],\
                    monitor_list[mon].dataDyn[var][:monitor_list[mon].dataDyn['time_lv_years'].shape[0],reg,0:regz].T,\
                    vmin = vimin, vmax = vimax, cmap = cmap,)
        cbar = fig.colorbar(ca , ax=axes.flat[ax_ind + col*2], ticks=ticks)
        axes.flat[ax_ind + col*2].set_title(unit_titles[var]+' '+reg_title[reg]+' '+monitor_list[mon].title)
        
        reg = 3
        regz = 50
        ca = axes.flat[ax_ind + col*3].pcolor(monitor_list[mon].dataDyn['time_lv_years'],monitor_list[mon].Z[0:regz],\
                    monitor_list[mon].dataDyn[var][:monitor_list[mon].dataDyn['time_lv_years'].shape[0],reg,0:regz].T,\
                    vmin = vimin, vmax = vimax, cmap = cmap,)
        cbar = fig.colorbar(ca , ax=axes.flat[ax_ind + col*3], ticks=ticks)
        axes.flat[ax_ind + col*3].set_title(unit_titles[var]+' '+reg_title[reg]+' '+monitor_list[mon].title)
        
        axes.flat[ax_ind + col*3].set_xlabel('Years')
        axes.flat[0].set_ylabel('m')
        axes.flat[col].set_ylabel('m')
        axes.flat[col*2].set_ylabel('m')
        axes.flat[col*3].set_ylabel('m')
        ax_ind += 1
        
    fig.subplots_adjust(right=2.4,top=2.8)

def plot_dynSt(monitor_list,var_list,reg,row=2,col=3):
    # this function is meant to plot monitor variables   
    
    titles = {'theta_mean' : 'Temperature Mean', 'theta_min' : 'Temperature Min',\
              'theta_max' :  'Temperature Max', 'eta_mean' : 'ETA mean', 'eta_max' : 'ETA max',\
              'eta_min' : 'ETA min', 'salt_mean' : 'Salinity Mean' , 'salt_max' : 'Salinity Max',\
              'salt_min' : 'Salinity Min' , 'sss_mean': 'SSS mean' ,'sss_max' : 'SSS max',\
              'sss_min' : 'SSS min', 'sst_mean' : 'SST mean', 'sst_max' : 'SST max' ,'sst_min' : 'SST min',\
              'vvel_mean' : 'V mean', 'vvel_max' : 'V max' , 'vvel_min' : 'V min', 'uvel_mean' : 'U mean',\
              'uvel_max' : 'U max', 'uvel_min' : 'U min', 'ke_mean' : 'Kinetic mean', 'ke_max' : 'Kinetic max',\
              'ke_vol' : 'Kinetic volume', 'seaice_area_max' : 'Seaice area max',\
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
    region = { 0 : 'Global' , 1 : 'Arctic' , 2 : 'Nordic Seas' , 3 : 'North Atl' }
    
    fig, axes = plt.subplots(row,col)
    ax_ind = 0
    for var in var_list:
        for mon in monitor_list:
            if reg == 'all':
                axes.flat[ax_ind].plot(monitor_list[mon].dataDyn['time_lv_years'],\
                                       monitor_list[mon].dataDyn[var][1:,0,0],monitor_list[mon].color,alpha=0.2)
                axes.flat[ax_ind].plot(monitor_list[mon].dataDyn['time_lv_years'],\
                                       monitor_list[mon].dataDyn[var][1:,1,0],monitor_list[mon].color,alpha=0.2)
                axes.flat[ax_ind].plot(monitor_list[mon].dataDyn['time_lv_years'],\
                                       monitor_list[mon].dataDyn[var][1:,2,0],monitor_list[mon].color,alpha=0.2)
                axes.flat[ax_ind].plot(monitor_list[mon].dataDyn['time_lv_years'],\
                                       monitor_list[mon].dataDyn[var][1:,3,0],monitor_list[mon].color,alpha=0.2)
                axes.flat[ax_ind].plot(monitor_list[mon].dataDyn['time_lv_years'][0:-12*3],\
                        runningMeanFast(monitor_list[mon].dataDyn[var][1:,0,0],12*3)[0:-12*3],monitor_list[mon].color)
                axes.flat[ax_ind].plot(monitor_list[mon].dataDyn['time_lv_years'][0:-12*3],\
                        runningMeanFast(monitor_list[mon].dataDyn[var][1:,1,0],12*3)[0:-12*3],monitor_list[mon].color)
                axes.flat[ax_ind].plot(monitor_list[mon].dataDyn['time_lv_years'][0:-12*3],\
                        runningMeanFast(monitor_list[mon].dataDyn[var][1:,2,0],12*3)[0:-12*3],monitor_list[mon].color)
                axes.flat[ax_ind].plot(monitor_list[mon].dataDyn['time_lv_years'][0:-12*3],\
                        runningMeanFast(monitor_list[mon].dataDyn[var][1:,3,0],12*3)[0:-12*3],monitor_list[mon].color,\
                                       label=monitor_list[mon].title)
                axes.flat[ax_ind].set_title(titles[var])
            else:
                axes.flat[ax_ind].plot(monitor_list[mon].dataDyn['time_lv_years'],\
                                       monitor_list[mon].dataDyn[var][1:,reg,0],monitor_list[mon].color,alpha=0.2)
                axes.flat[ax_ind].plot(monitor_list[mon].dataDyn['time_lv_years'][0:-12*3],\
                                       runningMeanFast(monitor_list[mon].dataDyn[var][1:,reg,0],12*3)[0:-12*3]\
                                       ,monitor_list[mon].color,label=monitor_list[mon].title)
                axes.flat[ax_ind].set_title(titles[var]+' '+region[reg])

            plt.legend(bbox_to_anchor=(-2.5, 1.02, 2., 1.5),
                   ncol=4, mode="expand", borderaxespad=0.)
            axes.flat[ax_ind].set_ylabel(unity[var])
            monitor_list[mon]
        ax_ind += 1

    fig.subplots_adjust(right=2.4,top=2.4)

def plot_psi(monitor_list,vimin,vimax,cmap,row=1,col=3):
    fig, axes = plt.subplots(nrows=row, ncols=col)
    for mon_v in monitor_list:        #print mon
        axes.flat[0].plot(monitor_list[mon_v].barostream['years'],monitor_list[mon_v].barostream['psi_mean'],\
                          monitor_list[mon_v].color,alpha=0.2)
        axes.flat[0].plot(monitor_list[mon_v].barostream['years'][0:-12],\
                          runningMeanFast(monitor_list[mon_v].barostream['psi_mean'],12)[0:-12],\
                          monitor_list[mon_v].color)
        axes.flat[0].set_xlabel('Month')
        axes.flat[0].set_ylabel('SV')
        axes.flat[0].set_title('Average barotropic streamfunction')
        
        axes.flat[1].plot(monitor_list[mon_v].barostream['years'],monitor_list[mon_v].barostream['psi_min'],\
                          monitor_list[mon_v].color,alpha=0.2)
        axes.flat[1].plot(monitor_list[mon_v].barostream['years'][0:-12],\
                          runningMeanFast(monitor_list[mon_v].barostream['psi_min'],12)[0:-12],\
                          monitor_list[mon_v].color)
        axes.flat[1].set_xlabel('Month')
        axes.flat[1].set_ylabel('SV')
        axes.flat[1].set_title('Minimum barotropic streamfunction')
        
        axes.flat[2].plot(monitor_list[mon_v].barostream['years'],monitor_list[mon_v].barostream['psi_max'],\
                          monitor_list[mon_v].color,alpha=0.2)
        axes.flat[2].plot(monitor_list[mon_v].barostream['years'][0:-12],\
                          runningMeanFast(monitor_list[mon_v].barostream['psi_max'],12)[0:-12],\
                          monitor_list[mon_v].color,label = monitor_list[mon_v].title)
        axes.flat[2].set_xlabel('Month')
        axes.flat[2].set_ylabel('SV')
        axes.flat[2].set_title('Maximum barotropic streamfunction')
        
        plt.legend(bbox_to_anchor=(-2.5, 1.02, 2., .3),
           ncol=4, mode="expand", borderaxespad=0.)
        
    fig.subplots_adjust(right=2.5,top=1.2)

def plot_flux_in_out(monitor_list,var_list,flux,row=2,col=3):
    # this function is meant to plot monitor variables   
    fig, axes = plt.subplots(row,col)
    ax_ind = 0
    for var in var_list:
        for mon in monitor_list:
            axes.flat[ax_ind].plot(monitor_list[mon].fluxes['years'],monitor_list[mon].fluxes[var]['FluxSum'],\
                                   monitor_list[mon].color,alpha=0.2)
            axes.flat[ax_ind].plot(monitor_list[mon].fluxes['years'][0:-12],\
                runningMeanFast(monitor_list[mon].fluxes[var]['FluxSum'],12)[0:-12],monitor_list[mon].color)
            axes.flat[ax_ind].plot(monitor_list[mon].fluxes['years'],\
                            monitor_list[mon].fluxes[var]['FluxInSum'],monitor_list[mon].color,alpha=0.2)
            axes.flat[ax_ind].plot(monitor_list[mon].fluxes['years'][0:-12],\
                runningMeanFast(monitor_list[mon].fluxes[var]['FluxInSum'],12)[0:-12],monitor_list[mon].color)
            axes.flat[ax_ind].plot(monitor_list[mon].fluxes['years'],\
                            monitor_list[mon].fluxes[var]['FluxOutSum'],monitor_list[mon].color,alpha=0.2)
            axes.flat[ax_ind].plot(monitor_list[mon].fluxes['years'][0:-12],\
                runningMeanFast(monitor_list[mon].fluxes[var]['FluxOutSum'],12)[0:-12],\
                                   monitor_list[mon].color,label = monitor_list[mon].title)
            axes.flat[ax_ind].set_title(var)
            axes.flat[ax_ind].set_ylabel('Sv')
            axes.flat[ax_ind].set_xlabel('Yers')
        ax_ind += 1
        plt.legend(bbox_to_anchor=(-2.5, 1.02, 2., 1.5),
                   ncol=4, mode="expand", borderaxespad=0.)
            
    fig.subplots_adjust(right=2.4,top=2.4)

def plot_flux(monitor_list,var_list,flux,row=2,col=3):
    # this function is meant to plot monitor variables   
    fig, axes = plt.subplots(row,col)
    ax_ind = 0
    for var in var_list:
        for mon in monitor_list:
            axes.flat[ax_ind].plot(monitor_list[mon].fluxes['years'],monitor_list[mon].fluxes[var][flux],\
                                   monitor_list[mon].color,alpha=0.2)
            axes.flat[ax_ind].plot(monitor_list[mon].fluxes['years'][0:-12],\
                        runningMeanFast(monitor_list[mon].fluxes[var][flux],12)[0:-12],\
                                   monitor_list[mon].color,label = monitor_list[mon].title)
            axes.flat[ax_ind].set_title(var)
            axes.flat[ax_ind].set_ylabel('Sv')
            axes.flat[ax_ind].set_xlabel('Yers')
        ax_ind += 1
        plt.legend(bbox_to_anchor=(-2.5, 1.02, 2., 1.5),
                   ncol=4, mode="expand", borderaxespad=0.)
    fig.subplots_adjust(right=2.4,top=2.4)

def plot_flux_total(monitor_list,var_list,row=2,col=3):
    # this function is meant to plot monitor variables                                                                                 
    ax_ind = 0
    row = len(var_list)
    flux_list = ['FluxSumFW','FluxSumS','FluxInSum','FluxOutSum','FluxSum']
    col = len(flux_list)
    fig, axes = plt.subplots(row,col, sharex=True,)
    for var in var_list:
        flux_ind = 0
        for flux in flux_list:
            for mon in monitor_list:
                #axes.flat[ax_ind + flux_ind].plot(monitor_list[mon].fluxes['years'],monitor_list[mon].fluxes[var][flux],\
                #                  monitor_list[mon].color,alpha=0.2)
                axes.flat[ax_ind + flux_ind].plot(monitor_list[mon].fluxes['years'][0:-12],\
                        runningMeanFast(monitor_list[mon].fluxes[var][flux],12)[0:-12],\
                                  monitor_list[mon].color)
            axes.flat[ax_ind + flux_ind].set_title(var+' '+flux[4:])
            axes.flat[ax_ind + flux_ind].set_ylabel('Sv')  
            axes.flat[(row-1)*col + flux_ind].set_xlabel('Years')  
            flux_ind += 1
        ax_ind += 1*len(flux_list)
        plt.legend(bbox_to_anchor=(-5, 6, 3.5, 1.5),
                ncol=4, mode="expand", borderaxespad=0.)
    fig.subplots_adjust(right=2.4,top=4.4)
    
def plot_mxldepth(monitor_list,reg,row=1,col=3):
    fig, axes = plt.subplots(nrows=row, ncols=col)
    regions = reg_titles()
    for mon_v in monitor_list:
        axes.flat[0].plot(monitor_list[mon_v].mxldepth_years, -monitor_list[mon_v].mxldepth_mean[reg],\
                          monitor_list[mon_v].color,alpha=0.2)
        axes.flat[0].plot(monitor_list[mon_v].mxldepth_years[0:-12],\
                runningMeanFast(-monitor_list[mon_v].mxldepth_mean[reg],12)[0:-12],
                monitor_list[mon_v].color)
        axes.flat[0].set_xlabel('Month')
        axes.flat[0].set_ylabel('m')
        axes.flat[0].set_title('Mean MXLDEPTH '+regions[reg])
        
        axes.flat[1].plot(monitor_list[mon_v].mxldepth_years, -monitor_list[mon_v].mxldepth_min[reg],
                monitor_list[mon_v].color,alpha=0.2)
        axes.flat[1].plot(monitor_list[mon_v].mxldepth_years[0:-12],\
                runningMeanFast(-monitor_list[mon_v].mxldepth_min[reg],12)[0:-12],
                monitor_list[mon_v].color)
        axes.flat[1].set_xlabel('Month')
        axes.flat[1].set_ylabel('m')
        axes.flat[1].set_title('Min MXLDEPTH '+regions[reg])
        
        axes.flat[2].plot(monitor_list[mon_v].mxldepth_years,-monitor_list[mon_v].mxldepth_min[reg],
                monitor_list[mon_v].color,alpha=0.2)
        axes.flat[2].plot(monitor_list[mon_v].mxldepth_years[0:-12],\
                runningMeanFast(-monitor_list[mon_v].mxldepth_max[reg],12)[0:-12],
                monitor_list[mon_v].color,label = monitor_list[mon_v].title)
        axes.flat[2].set_xlabel('Month')
        axes.flat[2].set_ylabel('m')
        axes.flat[2].set_title('Max MXLDEPTH '+regions[reg])
        
        plt.legend(bbox_to_anchor=(-2.5, 1.02, 2., .3),
           ncol=4, mode="expand", borderaxespad=0.)
        
    fig.subplots_adjust(right=2.5,top=1.2)


# this function plots different angles of ptracers
def plot_ptracer(data,ptracer,timings):
    ax_ind = 0
    col = 3
    row = len(timings)
    fig, axes = plt.subplots(row,col)
    vimin = 0.1
    vimax = 100
    
    cmap= matplotlib.cm.hot
    cmap.set_bad('grey',1.)
    
    for t in timings:
        masked_array = np.ma.array(np.nanmean(data.ptracers[ptracer][t,:,:,:],axis=0),\
                                   mask=np.isnan(np.nanmean(data.ptracers[ptracer][t,:,:,:],axis=0)))
        ca = axes.flat[ax_ind*col ].imshow(masked_array,vmin = vimin, vmax = vimax, interpolation='nearest',\
                                       cmap = cmap ,origin="left",aspect='auto',norm = matplotlib.colors.LogNorm())

        masked_array = np.ma.array(np.nanmean(data.ptracers[ptracer][t,:,:,0:200],axis=1),\
                                   mask=np.isnan(np.nanmean(data.ptracers[ptracer][t,:,:,0:200],axis=1)))
        ca = axes.flat[ax_ind*col +1].pcolormesh(data.X[0:200],data.Z,masked_array,vmin = vimin, vmax = vimax,cmap=cmap,
                                                 norm = matplotlib.colors.LogNorm())
        axes.flat[ax_ind*col+1].set_title("Perturbation "+ptracer+" after "+str(round(data.ptracers['years'][t],2))+" years")
        
        masked_array = np.ma.array(np.nanmean(data.ptracers[ptracer][t,:,:,:],axis=2),\
                                   mask=np.isnan(np.nanmean(data.ptracers[ptracer][t,:,:,:],axis=2)))
        ca = axes.flat[ax_ind*col + 2].pcolormesh(data.Y,data.Z,masked_array,vmin = vimin, vmax = vimax, cmap=cmap,
                                              norm = matplotlib.colors.LogNorm())
        cbar = fig.colorbar(ca , ax=axes.flat[ax_ind*col + 2])
        
        ax_ind += 1 
    fig.subplots_adjust(right=2.4,top=5.8)

