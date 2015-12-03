from scipy.io import netcdf
import numpy as np

class spinup:
    def __init__(self):
        self.T = []
        self.S = []
        self.lat = []
        self.lon = []
        self.depth = []
        self.bathy = []
        self.title = []
        self.ice = []

    def readSpin(self,path):
        file2read = netcdf.NetCDFFile(path,'r')
        T=file2read.variables['T']
        self.T=T[:]*1
        S=file2read.variables['S']
        self.S=S[:]*1
        lat=file2read.variables['lat']
        self.lat=lat[:]*1        
        lon=file2read.variables['lon']
        self.lon=lon[:]*1        
        depth=file2read.variables['depth']
        self.depth=depth[:]*1 
        bathy = file2read.variables['bathy']
        self.bathy = bathy[:]*1
        ice = file2read.variables['ice']
        self.ice = ice[:]*1
        file2read.close()
        
        
        self.S[self.bathy==0]=np.nan
        self.T[self.bathy==0]=np.nan
        self.ice[self.ice>2]=0
        self.ice[self.ice<0]=0
        self.ice[self.bathy[0,:,:]==0]=np.nan
        
        if T.shape[1] == 192:
            self.title = '36km run'
        elif T.shape[1] == 384:
            self.title = '18km run'
        elif T.shape[1] == 768:
            self.title = '9km run'

            
