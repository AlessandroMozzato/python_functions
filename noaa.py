# Getting NOA data

from scipy.io import netcdf
from scipy.io import loadmat
import numpy as np
import mpl_toolkits

class noaa:
    def __init__(self):
        self.lat = []
        self.lon = []
        self.T = []
        self.S =[]
        self.depth = []
    
    def read_noaa(self):
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/noaa_climatology/s00_04.nc",'r')
        s_an=file2read.variables['s_an']
        s_an=s_an[:]*1
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/noaa_climatology/t00_04.nc",'r')
        t_an=file2read.variables['t_an']
        t_an=t_an[:]*1
        lon=file2read.variables['lon']
        lon=lon[:]*1
        lat=file2read.variables['lat']
        lat=lat[:]*1
        depth=file2read.variables['depth']
        depth=depth[:]*1
        
        # We interpolate the data for the new resolution
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/results18km/grid.nc",'r')
        XC_18km = file2read.variables['XC']
        XC_18km = XC_18km[:]*1
        YC_18km = file2read.variables['YC']
        YC_18km = YC_18km[:]*1
        
        T_noaa = np.zeros((102, 384, 420))
        S_noaa = np.zeros((102, 384, 420))
        # looping in the vertical coordinate
        for k in range(0,102):
            T_noaa[k,:,:] = mpl_toolkits.basemap.interp(t_an[0,k,:,:], lon, \
                                            lat,XC_18km, YC_18km, checkbounds=False, masked=False, order=1)
            S_noaa[k,:,:] = mpl_toolkits.basemap.interp(s_an[0,k,:,:], lon, \
                                            lat,XC_18km, YC_18km, checkbounds=False, masked=False, order=1)
    
        T_noaa[T_noaa > 100] = np.nan
        S_noaa[S_noaa > 100] = np.nan
        
        self.lat = lat
        self.lon = lon
        self.T = T_noaa
        self.S = S_noaa
        self.depth = depth
