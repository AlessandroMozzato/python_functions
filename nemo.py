# Getting NEMO data
from scipy.io import netcdf
from scipy.io import loadmat
import numpy as np
import mpl_toolkits
import sys
sys.path.append('/noc/users/am8e13/Python/')
import komod

class nemo:
    def __init__(self):
        self.lat = []
        self.lon = []
        self.T = []
        self.S = []
        self.depth = []
        self.mask = []
        
    def read_nemo(self,res):
        if res == '83':
            xdim = 3059
            ydim = 4322
            zdim = 75
            ycut = 2000
        elif res == '25':
            xdim = 1021
            ydim = 1442
            zdim = 64
            ycut = 2000

        path = '/scratch/general/am8e13/NEMO_data/'
        mask_nemo = komod.mitbin(path+'NEMO'+res+'_mask',xdim=xdim,ydim=ydim,zdim=zdim,datatype='float32')
        lat_nemo = komod.mitbin(path+'NEMO'+res+'_lat',xdim=xdim,ydim=ydim,datatype='float32')
        lon_nemo = komod.mitbin(path+'NEMO'+res+'_lon',xdim=xdim,ydim=ydim,datatype='float32')
        depth_nemo = komod.mitbin(path+'NEMO'+res+'_depth',xdim=1,ydim=1,zdim=zdim,datatype='float32')
        T_nemo = komod.mitbin(path+'NEMO'+res+'_temp',xdim=xdim,ydim=ydim,zdim=zdim,datatype='float32')
        S_nemo = komod.mitbin(path+'NEMO'+res+'_salt',xdim=xdim,ydim=ydim,zdim=zdim,datatype='float32')

        T_nemo[mask_nemo == 0] = np.nan
        T_nemo = T_nemo.squeeze(axis=0)
        T_nemo = T_nemo[:,ycut:,:]

        S_nemo[mask_nemo == 0] = np.nan
        S_nemo = S_nemo.squeeze(axis=0)
        S_nemo = S_nemo[:,ycut:,:]

        lat_nemo = lat_nemo.squeeze(axis=0)
        lat_nemo = lat_nemo.squeeze(axis=0)
        lat_nemo = np.array(lat_nemo[ycut:,:])
        lon_nemo = lon_nemo.squeeze(axis=0)
        lon_nemo = lon_nemo.squeeze(axis=0)
        lon_nemo = np.array(lon_nemo[ycut:,:])

        self.lat = lat_nemo
        self.lon = lon_nemo
        self.T = T_nemo
        self.S = S_nemo
        self.depth = depth_nemo
        self.mask = mask_nemo
