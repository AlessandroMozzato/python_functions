#! /usr/bin/env python                                                                                                         # -*- coding: utf-8 -*-                                                                                                         
"""####This script is meant to calculate topostrophy. 
Topostrophy is defined in Holloway 2008 as follows:
\begin{equation}
\tau = \hat{\mathrm{u}} \cdot ( \hat{\mathrm{z}}\ \times\ \hat{\mathrm{s}})
\end{equation}
where $\hat{\mathrm{u}}$ is the velocity, $\hat{\mathrm{z}}$ is the vertical unit vector $(0,0,1)$ and $ \hat{\mathrm{s}}=\nabla H$ is the gradient of the Bathymetry.\\
\\
A similar definition is given in Merryfield and Scott 2007 as follows:
\begin{equation}
T(H_j,d_k)=\frac{\sum_i(\hat{\mathrm{u}}_i\cdot \hat{\mathrm{u}}_{n,i})W_i \delta_{jk,i}dV_i}{\sum_i W_i \delta_{jk,i}dV_i}
\end{equation}
where $\hat{\mathrm{u}}_{n} = (\hat{z}\times\hat{s})f/|f|$ with $\hat{\mathrm{z}}$ and $\hat{\mathrm{s}}$ defined above and $f$ the Coriolis parameter, $W_i = |s_i||f_i|$ are the weights, $\delta_{jk,i}$ is a delta function and $dV_i$ is the cell volume. $k$ is the level we are considering whereas $j$ is the depth level for the considered cell. $i$ represent the considered cell, sums for $i$ can be considered to consider various areas.

####Caluclation for the gradient of bathymetry
Gradient of bathymetry is calculated using the following the following formula:
\begin{equation}
(s_x)_n = -\frac{(H_{n+1}-H_{n-1})}{(x_{n+1}-x_{n-1})}
\end{equation}
the same formulat is applied for the $y$ coordinate."""

def CellVolume(res):
    if res == 36:
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/results36km/grid.nc",'r')
    elif res == 18:
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/results18km/grid.nc",'r')
    hfacc = file2read.variables['HFacC']
    hfacc = hfacc[:]*1  
    drf = file2read.variables['drF']
    drf = drf[:]*1
    rA = file2read.variables['rA']
    rA = rA[:]*1
    dydx = np.zeros_like(hfacc)
    for k in range(len(drf)):
        dydx[k,:,:] = drf[k]*rA*hfacc[k,:,:]    
    dydx[hfacc==0]=np.nan
    return dydx

def balanceCalc(data):
    if data['ADVr_TH'].shape[3] == 210:
        res = 36
    elif data['ADVr_TH'].shape[3] == 420:
        res = 18
    else:
        print "Dimension Error"
        
    CellVol = CellVolume(res)
    
    file2read = netcdf.NetCDFFile("/scratch/general/am8e13/results"+str(res)+"km/grid.nc",'r')
    RF = file2read.variables['RF']
    RF = RF[:]*1
    drf = file2read.variables['drF']
    drf = drf[:]*1
    hfacc = file2read.variables['HFacC']
    hfacc = hfacc[:]*1
    RAC = file2read.variables['RC']
    hfacc[hfacc==0]=np.nan
    
    nz = CellVol.shape[0]
    ny = CellVol.shape[1]
    nx = CellVol.shape[2]
    nt = data['ADVr_TH'].shape[0]
    swfrac = 0.62*np.exp(RF[0:nz-1]/0.6) + (1.0 - 0.62)*np.exp(RF[0:nz-1]/20.0)
    swfrac1 = 0.62*np.exp(RF[1:nz]/0.6) + (1.0 - 0.62)*np.exp(RF[1:nz]/20.0)
    swdiff = swfrac - swfrac1
    rhoConst = 999.8
    Cp = 4.210
    rhoCp = rhoConst*Cp
    
    Adv_tend = np.zeros_like(data['ADVr_TH'][:,:,:,:])
    Dif_tend = np.zeros_like(data['ADVr_TH'][:,:,:,:])
    Kpp_tend = np.zeros_like(data['ADVr_TH'][:,:,:,:])    
    Qsw_tend = np.zeros_like(data['ADVr_TH'][:,:,:,:])
    Tflx_tend = np.zeros_like(data['ADVr_TH'][:,:,:,:])    
    Heat_tend = np.zeros_like(data['ADVr_TH'][:,:,:,:]) 
    
    for it in range(nt):
        Adv_tend[it,0:49,:,:] = -(data['ADVr_TH'][it,0:nz-1,:,:] - data['ADVr_TH'][it,1:nz,:,:])/CellVol[0:49,:,:]
        Adv_tend[it,:,:,:] += -((data['ADVx_TH'][it,:,:,1:nx+1] - data['ADVx_TH'][it,:,:,0:nx]) / CellVol + \
                                (data['ADVy_TH'][it,:,1:ny+1,:] - data['ADVy_TH'][it,:,0:ny,:]) / CellVol )
        Dif_tend[it,0:49,:,:] = -(data['DFrE_TH'][it,0:nz-1,:,:] - data['DFrE_TH'][it,1:nz,:,:]) /CellVol[0:49,:,:]                                                                                                
        Dif_tend[it,0:49,:,:] = -(data['DFrI_TH'][it,0:nz-1,:,:] - data['DFrI_TH'][it,1:nz,:,:]) /CellVol[0:49,:,:]                                                                                                
        Dif_tend[it,:,:,:] += -((data['DFxE_TH'][it,:,:,1:nx+1] - data['DFxE_TH'][it,:,:,0:nx]) / CellVol + \
                                (data['DFyE_TH'][it,:,1:ny+1,:] - data['DFyE_TH'][it,:,0:ny,:]) / CellVol ) 
        Kpp_tend[it,0:49,:,:] = -(data['KPPg_TH'][it,0:nz-1,:,:] - data['KPPg_TH'][it,1:nz,:,:])/CellVol[0:49,:,:]
                                                                                                
        for iz in range(nz-1):
            Qsw_tend[it,iz,:,:] = data['oceQsw'][it,0,:,:] / rhoCp /((drf[iz])*hfacc[iz,:,:])*swdiff[iz]                                                            
            Tflx_tend[it, iz,:,:] = (data['TFLUX'][it,0,:,:] - data['oceQsw'][it,0,:,:])/(rhoCp*drf[0]*hfacc[iz,:,:])
    
    for it in range(nt):                                                                
        for iz in range(nz-1):
            if iz == 0:
                Heat_tend[it,iz,:,:] = Adv_tend[it,iz,:,:] + Dif_tend[it,iz,:,:] + Kpp_tend[it,iz,:,:]  \
                                        + Qsw_tend[it,iz,:,:] + Tflx_tend[it,iz,:,:]
            else:
                Heat_tend[it,iz,:,:] = Adv_tend[it,iz,:,:] + Dif_tend[it,iz,:,:] + Kpp_tend[it,iz,:,:]  \
                                        + Qsw_tend[it,iz,:,:]
    return Adv_tend, Dif_tend, Kpp_tend, Qsw_tend, Tflx_tend, Heat_tend

class FieldForBalance():
    def __init__(self):
        self.data = {'ADVr_TH' : [] , 'ADVx_TH' : [] , 'ADVy_TH' : [] , 'TOTTTEND' : [] , \
                    'DFrE_TH' : [] , 'DFxE_TH' : [] , 'DFyE_TH' : [] , 'DFrI_TH ' : [] , \
                    'KPPg_TH' : [] , 'TFLUX' : [] , 'oceQsw' : [] , 'WTHMASS' : [] , 'oceQnet' : [] }
        
    def ReadData(self,path,list_var):
        file2read = netcdf.NetCDFFile(path+'heatbal1.nc','r')
        ADVr_TH = file2read.variables['ADVr_TH']
        self.data['ADVr_TH'] = ADVr_TH[list_var]*1
        ADVx_TH = file2read.variables['ADVx_TH']
        self.data['ADVx_TH'] = ADVx_TH[list_var]*1
        ADVy_TH = file2read.variables['ADVy_TH']
        self.data['ADVy_TH'] = ADVy_TH[list_var]*1
        TOTTTEND = file2read.variables['TOTTTEND']
        self.data['TOTTTEND'] = TOTTTEND[list_var]*1
        DFrE_TH = file2read.variables['DFrE_TH']
        self.data['DFrE_TH'] = DFrE_TH[list_var]*1
        DFxE_TH = file2read.variables['DFxE_TH']
        self.data['DFxE_TH'] = DFxE_TH[list_var]*1
        DFyE_TH = file2read.variables['DFyE_TH']
        self.data['DFyE_TH'] = DFyE_TH[list_var]*1
        DFrI_TH = file2read.variables['DFrI_TH']
        self.data['DFrI_TH'] = DFrI_TH[list_var]*1
        KPPg_TH = file2read.variables['KPPg_TH']
        self.data['KPPg_TH'] = KPPg_TH[list_var]*1
        days = file2read.variables['T']
        self.data['days']=days[list_var]*1
        self.years = (self.data['days'] - self.data['days'][0])/(60*60*24*360)  
       
        file2read = netcdf.NetCDFFile(path+'heatbal2.nc','r')
        TFLUX = file2read.variables['TFLUX']
        self.data['TFLUX'] = TFLUX[list_var]*1        
        oceQsw = file2read.variables['oceQsw']
        self.data['oceQsw'] = oceQsw[list_var]*1
        WTHMASS = file2read.variables['WTHMASS']
        self.data['WTHMASS'] = WTHMASS[list_var]*1
        oceQnet = file2read.variables['oceQnet']
        self.data['oceQnet'] = oceQsw[list_var]*1             
        
    def title(self,title):
        self.title = title
    
    def print_title(self):
        print self.title
        
    def fluxCalc(self):
        Adv_tend, Dif_tend, Kpp_tend, Qsw_tend, Tflx_tend, Heat_tend = balanceCalc(self.data)
        self.heatbal = {'Adv_tend' : Adv_tend, 'Dif_tend' : Dif_tend, 'Kpp_tend' : Kpp_tend, \
                        'Qsw_tend' : Qsw_tend, 'Tflx_tend' : Tflx_tend, 'Heat_tend' : Heat_tend}
    
    def fluxAverage(self):
        self.heat_ave = {'Adv_tend' : [], 'Dif_tend' : [], 'Kpp_tend' : [], \
                        'Qsw_tend' : [], 'Tflx_tend' : [], 'Heat_tend' : [], 'TOTTTEND': [] }
        for var in self.heat_ave:
            if var == 'TOTTTEND':
                self.heat_ave[var] = np.nansum(np.nansum(np.nansum(self.data[var],axis=3),axis=2),axis=1)
            else:
                self.heat_ave[var] = np.nansum(np.nansum(np.nansum(self.heatbal[var],axis=3),axis=2),axis=1)
