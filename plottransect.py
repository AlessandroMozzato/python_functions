# This script contains the functions required to plot a transect
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
from mpl_toolkits.basemap import Basemap
import sys
sys.path.append('/noc/users/am8e13/PyNGL-1.4.0/lib/python2.7/site-packages/PyNGL/')
import Ngl
sys.path.append('/noc/users/am8e13/Python/')
import komod
sys.path.append('/noc/users/am8e13/Python/PyNIO-1.4.0/')
import Nio

import cmocean

def lastnan(data,Z):
    # This function finds the last layer with some non-NaN value
    z = 0
    last_nonnan = len(Z)-1
    while z < len(Z):
        count = np.count_nonzero(~np.isnan(data[z,:]))
        if count == 0 or count == 1:
            last_nonnan = z
            z = len(Z)                    
        z = z + 1
    z2 = last_nonnan
    Z = np.array(Z[0:z2])
    return Z, z2

def pointfind2(plat, plon, lat, lon, pdif=1):
    """ return indeces and values of the grid point closest to lat/lon point
    the same as pointfind but could be faster
    Usage: 
        pointfind(plat, plon, lat, lon, pdif = 0.5)
    Input:
        plat - latitude of the point
        plon - longitude of the point
        lat  - 2d array of latutudes
        lon  - 2d array of longitudes
        pdif  - we don't need it but leave it to have the same input as pointfind
    Output:
        indeces and values of the points that fulfil conditions  
    """

    dist_min = 1000000.
    
    
    for i in range(lon.shape[0]):
        for j in range(lon.shape[1]):
            dist = Ngl.gc_dist(plat,plon,lat[i,j],lon[i,j])
            if dist_min > dist:
                dist_min = dist
                i_min = i
                j_min = j
                lat_min = lat[i,j]
                lon_min = lon[i,j]
    
#    print(i_min,j_min,lat_min,lon_min)
    gg1 = i_min, j_min
    
    return(gg1, lat_min, lon_min)

def get_transect1(lat, lon, data, lat1, lon1, lat2, lon2, npoints = 10, pdif = 1, norep=False):
    import numpy as numpy

    plat, plon = Ngl.gc_interp(lat1,lon1,lat2,lon2,npoints)
    
    m_grid = numpy.zeros((plon.shape[0]))
    n_grid = numpy.zeros((plat.shape[0]))
    grid_lats = numpy.zeros((plat.shape[0]))
    grid_lons = numpy.zeros((plat.shape[0]))
    
    for k in range(plon.shape[0]):
        coord, trash1, trash1 = pointfind2(plat[k], plon[k], lat, lon, pdif)
        m_grid[k] = coord[0]
        n_grid[k] = coord[1]
        grid_lats[k] = lat[m_grid[k],n_grid[k]]
        grid_lons[k] = lon[m_grid[k],n_grid[k]]
    
    
    if len(data.shape) == 4:
        data = data[0,:,:,:]
    elif len(data.shape) == 5:
        data = data[0,0,:,:,:]
    
        
    data_prof = numpy.zeros(([data.shape[0],npoints]))   
    
    for j in range(data_prof.shape[1]):
        data_prof[:,j] = data[:,m_grid[j],n_grid[j]]
    
    
    x_kilometers = numpy.zeros((grid_lats.shape[0]))
    
    for j in range(grid_lats.shape[0] - 1):
        angle = Ngl.gc_dist(grid_lats[j], grid_lons[j], grid_lats[j+1], grid_lons[j+1] )
        x_kilometers[j+1] = Ngl.gc_convert(angle, 2) + x_kilometers[j]
    
    
    #numpy.savez("outfile.npz", data_prof, x_kilometers)
    
    if norep == True:
        data_prof_norep = numpy.zeros(([data_prof.shape[0],1]))
        x_kilometers_norep = numpy.array(([]))
        m_grid_norep       = numpy.array(([]))
        n_grid_norep       = numpy.array(([]))
        present_location = -999999.999

        for k in range(x_kilometers.shape[0]):
            if int(x_kilometers[k]) != int(present_location):
                present_location = x_kilometers[k]
                data_prof_norep = numpy.hstack((data_prof_norep, data_prof[:,k:k+1])) 
                x_kilometers_norep = numpy.append(x_kilometers_norep, x_kilometers[k:k+1])
                m_grid_norep = numpy.append(m_grid_norep, m_grid[k:k+1])
                n_grid_norep = numpy.append(n_grid_norep, n_grid[k:k+1])
            

        data_prof_norep = data_prof_norep[:,1:] 
        data_prof = data_prof_norep
        x_kilometers = x_kilometers_norep
        m_grid = m_grid_norep
        n_grid = n_grid_norep
        
    return(data_prof, x_kilometers, m_grid, n_grid )

def comp_plot(data, region,npoints=15):
    import matplotlib as mpl
    # This function plots time evolution of temperature/salinity on a transect.                                                                                           
    # Transect are: Fram Strait, Bering Strait ...                                                                                                                        
    npl = len(data)

    lat1 = region[0]
    lon1 = region[1]
    lat2 = region[2]
    lon2 = region[3]
    fig, axes = plt.subplots(3,npl,sharex='col', sharey='row')
    ind = 0
    
    t_min = -2 ; t_max = 8 ; tempbounds = np.linspace(t_min,t_max,15) ; tempbounds1 = np.linspace(t_min,t_max,5)
    s_min = 30 ; s_max = 35.5 ; saltbounds = np.linspace(s_min,s_max,15) ; saltbounds1 = np.linspace(s_min,s_max,5)
    r_min = 26 ; r_max = 29 ; rhobounds = np.linspace(r_min,r_max,15) ; rhobounds1 = np.linspace(r_min,r_max,5)
    
    for run in data:
        # plot temperature
        data_prof, x_kilometers, m_grid, n_grid  = \
        get_transect1(data[run].lat, data[run].lon, data[run].T , lat1, lon1, lat2, lon2,\
                                                                 npoints = npoints, pdif = 1, norep=False)
        Z,z2 = lastnan(data_prof,data[run].depth)
        if ind == 0:
            imT = axes.flat[ind].contourf(x_kilometers,Z,data_prof[0:z2,:],vmin=t_min,vmax=t_max,levels = tempbounds,\
                    extend = 'both', cmap = cmocean.cm.temperature)
        else:
            axes.flat[ind].contourf(x_kilometers,Z,data_prof[0:z2,:],vmin=t_min,vmax=t_max,levels = tempbounds,\
                    extend = 'both', cmap = cmocean.cm.temperature)
        #axes.flat[ind].contour(x_kilometers,Z,data_prof[0:z2,:],colors='k',levels = tempbounds1,\
        #            extend = 'both')
        axes.flat[ind].set_title("T "+data[run].title)
        axes.flat[ind].title.set_fontsize('14')
        if ind == 0:
            axes.flat[ind].set_ylabel('m')

        # plot salinity
        data_prof, x_kilometers, m_grid, n_grid  = \
        get_transect1(data[run].lat, data[run].lon, data[run].S , lat1, lon1, lat2, lon2,\
                                                                 npoints = npoints, pdif = 1, norep=False)
        if ind == 0:
            imS = axes.flat[ind+npl].contourf(x_kilometers,Z,data_prof[0:z2,:],vmin=s_min,vmax=s_max,levels = saltbounds,\
                    extend = 'both' , cmap = cmocean.cm.salt)
        else:
            axes.flat[ind+npl].contourf(x_kilometers,Z,data_prof[0:z2,:],vmin=s_min,vmax=s_max,levels = saltbounds,\
                    extend = 'both' , cmap = cmocean.cm.salt)
        #axes.flat[ind+npl].contour(x_kilometers,Z,data_prof[0:z2,:],colors='k',levels = saltbounds1,\
        #            extend = 'both')
        axes.flat[ind+npl].set_title("S "+data[run].title)
        axes.flat[ind+npl].title.set_fontsize('14')
        if ind == 0:
            axes.flat[ind+npl].set_ylabel('m')
        
        # plot density
        data_prof, x_kilometers, m_grid, n_grid  = \
        get_transect1(data[run].lat, data[run].lon, data[run].rhop - 1000 , lat1, lon1, lat2, lon2,\
                                                                 npoints = npoints, pdif = 1, norep=False)
        if ind == 0:
            imrho = axes.flat[ind+npl*2].contourf(x_kilometers,Z,data_prof[0:z2,:],\
                    vmin=r_min,vmax=r_max,levels = rhobounds,\
                    extend = 'both',cmap = cmocean.cm.rho)
        else:
            axes.flat[ind+npl*2].contourf(x_kilometers,Z,data_prof[0:z2,:],\
                    vmin=r_min,vmax=r_max,levels = rhobounds,\
                    extend = 'both',cmap = cmocean.cm.rho)
        #axes.flat[ind+npl*2].contour(x_kilometers,Z,data_prof[0:z2,:],colors='k',levels = rhobounds,\
        #            extend = 'both')
        axes.flat[ind+npl*2].set_title("rho "+data[run].title)
        axes.flat[ind+npl*2].title.set_fontsize('14')
        if ind == 0:
            axes.flat[ind+npl*2].set_ylabel('m')
        axes.flat[ind+npl*2].set_xlabel('km')

        #for item in ([axes.flat[ind].title, axes.flat[ind].xaxis.label, axes.flat[ind].yaxis.label]):                                                                    
         #   item.set_fontsize(14)                                                                                                                                        
        ind = ind + 1

    cbar_ax = fig.add_axes([2.15, 2.2, 0.045, 0.7])
    cbar = plt.colorbar(imT, cax=cbar_ax,ticks=[-1,0,1,2,3,4,5,6,7])
    cbar.ax.set_ylabel('C')

    cbar_ax2 = fig.add_axes([2.15, 1.2, 0.045, 0.7])
    cbar2 = plt.colorbar(imS, cax=cbar_ax2,ticks=[30,31,32,33,34,35])
    cbar2.ax.set_ylabel('psu')

    cbar_ax3 = fig.add_axes([2.15, 0.2, 0.045, 0.7])
    cbar3 = plt.colorbar(imrho, cax=cbar_ax3,ticks=[26,26.5,27,27.5,28,28.5,29])
    cbar3.ax.set_ylabel('kg/m^3 - 1000')
    
    fig.subplots_adjust(right=2.1,top=3.)

# This function plots surface data

import numpy as np
from mpl_toolkits.basemap import Basemap

def regbase(region):
    '''Takes name of the region and returns dictionary with
    information necessary for creation of the Basemap instance
    '''

    mapDict = {}

    if region == 'Arctic':
        mapDict['projection'] = 'npstere'
        mapDict['boundinglat'] = 60
        mapDict['lon_0'] = 0
        mapDict['resolution'] = 'l'

    return mapDict

def bp(lon, lat, data, yescbar, region = 'Arctic', ptype = 'contourf', **kwargs):
    
    '''Basic Basemap plot function. Use coordinates (1d or 2d), data and name of the region
     as an input and plot data. Region defines in the "regbase" function.

     You can also provide any argument for matplotlib plotting functions.

     Usage:
         bp(lon, lat, data, region = 'Arctic', ptype = 'contourf', **kwargs)
     
     Input:
        lon         - 2D or 1D array of longitudes
        lat         - 2D or 1D array of latitudes
        data        - 2D array of scalar data.
        region      - one of the predefined regions (for list of regions see the "regbase" function)
        ptype       - plot type (contour, contourf, pcolor, pcolormesh)
        **kwargs    - arguments for plotting functions

     Output:
        Basemap instance.
    '''
    
    mapDict = regbase(region)

    # Create Basemap instance
    if mapDict['projection'] == 'npstere':
        m = Basemap(projection=mapDict['projection'],boundinglat=mapDict['boundinglat'],\
                    lon_0=mapDict['lon_0'],resolution=mapDict['resolution'])
    
    # Check if we have proper number of dimensions for lon (and hopefully lat as well)
    if lon.shape.__len__() == 1:
        lon, lat = np.meshgrid(lon, lat)
    elif lon.shape.__len__() > 2:
        raise Exception("Coordinate variables (lon) has too many dimensions")
    
    # Convert lat/lon to map coordinates
    x, y = m(lon, lat)

    # Make the map look better
    m.fillcontinents(color='gray',lake_color='gray')
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='white')
    
    # Draw values on the map
    if ptype == 'contourf':
        cs = m.contourf(x,y,data,**kwargs)
        if yescbar == True:
            cbar3 = plt.colorbar(cs)
        return m
    elif ptype == 'pcolormesh':
        cs = m.pcolormesh(x,y,data,**kwargs)
    elif ptype == 'contour':
        cs = m.contour(x,y,data,**kwargs)
    elif ptype == 'pcolor':
        cs = m.pcolor(x,y,data,**kwargs)
    else:
        raise Exception("Plot type not supported. Valid plot types are: contour, contourf, pcolor, pcolormesh ")
    
    return m

def plot_seaice(data1,data2,data3):
    tempbounds = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

    cs = plt.contourf(data3.ice,vmin=0,vmax=1,level=tempbounds,extend='both')
    plt.colorbar(cs)
    plt.close()

    fig, (ax1, ax2 , ax3) = plt.subplots(1,3)
    ax1 = plt.subplot(1,3,1)
    im1 = bp(data1.lon, data1.lat, data1.ice, yescbar = False,  vmin=0, vmax=1, level=tempbounds, extend='both')
    ax1.set_title('Sea-ice '+data1.title)
    ax1.title.set_fontsize('14')

    ax2 = plt.subplot(1,3,2)
    im2 = bp(data2.lon, data2.lat, data2.ice, yescbar = False, vmin=0,vmax=1,level=tempbounds,extend='both')
    ax2.set_title('Sea-ice '+data2.title)
    ax2.title.set_fontsize('14')

    ax3 = plt.subplot(1,3,3)
    im3 = bp(data3.lon, data3.lat, data3.ice, yescbar = False, vmin=0,vmax=1,level=tempbounds,extend='both')
    ax3.set_title('Sea-ice '+data3.title)
    ax3.title.set_fontsize('14')

    cbar_ax = fig.add_axes([2.2, 0.4, 0.045, 0.9])
    cbar = plt.colorbar(cs, cax=cbar_ax,)
    cbar.ax.set_ylabel('Fraction')  

    fig.subplots_adjust(right=2.1,top=1.6)

def plot_sst(data1,data2,data3,data4):
    vmin = -2
    vmax = 12
    tempbounds = range(vmin,vmax,2)

    cs = plt.contourf(data4.T[0,:,:],vmin=vmin,vmax=vmax,level=tempbounds,extend='both')
    plt.colorbar(cs)
    plt.close()

    fig, ((ax1, ax2) , (ax3 ,ax4)) = plt.subplots(2,2)
    ax1 = plt.subplot(2,2,1)
    im1 = bp(data1.lon, data1.lat, data1.T[0,:,:], yescbar = False,  vmin=vmin, vmax=vmax, level=tempbounds, extend='both')
    ax1.set_title('SST '+data1.title)
    ax1.title.set_fontsize('14')

    ax2 = plt.subplot(2,2,2)
    im2 = bp(data2.lon, data2.lat, data2.T[0,:,:], yescbar = False, vmin=vmin,vmax=vmax,level=tempbounds,extend='both')
    ax2.set_title('SST '+data2.title)
    ax2.title.set_fontsize('14')

    ax3 = plt.subplot(2,2,3)
    im3 = bp(data3.lon, data3.lat, data3.T[0,:,:], yescbar = False, vmin=vmin,vmax=vmax,level=tempbounds,extend='both')
    ax3.set_title('SST '+data3.title)
    ax3.title.set_fontsize('14')

    ax4 = plt.subplot(2,2,4)
    im4 = bp(data4.lon, data4.lat, data4.T[0,:,:], yescbar = False, vmin=vmin,vmax=vmax,level=tempbounds,extend='both')
    ax4.set_title('SST '+data4.title)
    ax4.title.set_fontsize('14')

    cbar_ax = fig.add_axes([1.9, 0.4, 0.045, 1.4])
    cbar = plt.colorbar(cs, cax=cbar_ax,)
    cbar.ax.set_ylabel('C')  

    fig.subplots_adjust(right=1.7,top=2.)
