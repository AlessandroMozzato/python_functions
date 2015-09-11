#! /usr/bin/env python    
# -*- coding: utf-8 -*-  

"""####This script is meant to calculate topostrophy.                                                                           
Topostrophy is defined in Holloway 2008 as follows:                                                                             
\begin{equation}                                                                                                                
\tau = \hat{\mathrm{u}} \cdot ( \hat{\mathrm{z}}\ \times\ \hat{\mathrm{s}})                                                     
\end{equation}                                                                                                                  
where $\hat{\mathrm{u}}$ is the velocity, $\hat{\mathrm{z}}$ is the vertical unit vector $(0,0,1)$ and $ \hat{\mathrm{s}}=\nabl\
a H$ is the gradient of the Bathymetry.\\                                                                                       
\\                                                                                                                              
A similar definition is given in Merryfield and Scott 2007 as follows:                                                          
\begin{equation}                                                                                                                
T(H_j,d_k)=\frac{\sum_i(\hat{\mathrm{u}}_i\cdot \hat{\mathrm{u}}_{n,i})W_i \delta_{jk,i}dV_i}{\sum_i W_i \delta_{jk,i}dV_i}     
\end{equation}                                                                                                                  
where $\hat{\mathrm{u}}_{n} = (\hat{z}\times\hat{s})f/|f|$ with $\hat{\mathrm{z}}$ and $\hat{\mathrm{s}}$ defined above and $f$\
 the Coriolis parameter, $W_i = |s_i||f_i|$ are the weights, $\delta_{jk,i}$ is a delta function and $dV_i$ is the cell volume.\
 $k$ is the level we are considering whereas $j$ is the depth level for the considered cell. $i$ represent the considered cell,\
 sums for $i$ can be considered to consider various areas.                                                                      
                                                                                                                                
####Caluclation for the gradient of bathymetry                                                                                  
Gradient of bathymetry is calculated using the following the following formula:                                                 
\begin{equation}                                                                                                                
(s_x)_n = -\frac{(H_{n+1}-H_{n-1})}{(x_{n+1}-x_{n-1})}                                                                          
\end{equation}                                                                                                                  
the same formulat is applied for the $y$ coordinate."""


from scipy.io import netcdf
from scipy.io import loadmat
import numpy as np
from pylab import clf, plot, show, floor, ceil, imshow
import matplotlib
import matplotlib.pyplot as plt
import os
import csv
import sys
import glob

def BathyGrad(res,depth_partial):
    ### This function calculates the gradient of the bathymetry
    if res == 18:
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/results18km/grid.nc",'r')
    elif res == 36:
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/results36km/grid.nc",'r')
    Depth = file2read.variables['Depth']
    Depth = Depth[:]*1
    hfacc=file2read.variables['HFacC']
    hfacc=hfacc[:]*1
    dxF = file2read.variables['dxF']
    dxF = dxF[:]*1
    dyF = file2read.variables['dyF']
    dyF = dyF[:]*1
    X = file2read.variables['X']
    X = X[:]*1
    Y = file2read.variables['Y']
    Y = Y[:]*1
    Z = file2read.variables['Z']
    Z = Z[:]*1
    
    if (depth_partial==True):
        depth_recalc = np.zeros_like(hfacc)
        Zp = np.zeros_like(Z)
        for z in range(len(Z))[1:]:
            Zp[z] = Z[z] - Z[z-1]
        Zp[0]=-10 
        Zp = abs(Zp)
        for x in range(len(X)):
            for y in range(len(Y)):
                depth_recalc[:,y,x] = hfacc[:,y,x]*Zp
        Depth = np.sum(depth_recalc,axis=0)

    Sx = np.zeros_like(hfacc[0,:,:])
    Sy = np.zeros_like(hfacc[0,:,:])
    for x in range(len(X)-1)[1:]:
        for y in range(len(Y)-1)[1:]:
            if x == 0 and y == 0:
                Sx[y,x] = -(Depth[y,1]-Depth[y,0])/(dxF[y,1]-dxF[y,0])
                Sy[y,x] = -(Depth[1,x]-Depth[0,x])/(dyF[1,x]-dyF[1,x])
            elif x == 0 and y == (len(Y)-1):
                Sx[y,x] = -(Depth[y,x+1]-Depth[y,0])/(dxF[y,1]-dxF[y,0])
                Sy[y,x] = -(Depth[(len(Y)-1),x]-Depth[y-1,x])/(dyF[(len(Y)-1),x]-dyF[y-1,x])
            elif x == (len(X)-1) and y == 0:
                Sx[y,x] = -(Depth[y,(len(X)-1)]-Depth[y,x-1])/(dxF[y,(len(X)-1)]-dxF[y,x-1])
                Sy[y,x] = -(Depth[y+1,x]-Depth[0,x])/(dyF[y+1,x]-dyF[0,x])
            elif x == (len(X)-1) and y == (len(Y)-1):
                Sx[y,x] = -(Depth[y,(len(X)-1)]-Depth[y,x-1])/(dxF[y,(len(X)-1)]-dxF[y,x-1])
                Sy[y,x] = -(Depth[y+1,x]-Depth[y-1,x])/(dyF[(len(Y)-1),x]-dyF[y-1,x])
            elif x == 0:
                Sx[y,x] = -(Depth[y,1]-Depth[y,0])/(dxF[y,1]-dxF[y,0])
                Sy[y,x] = -(Depth[y+1,x]-Depth[y-1,x])/(dyF[y+1,x]-dyF[y-1,x])
            elif x == (len(X)-1):
                Sx[y,x] = -(Depth[y,(len(X)-1)]-Depth[y,x-1])/(dxF[y,(len(X)-1)]-dxF[y,x-1])
                Sy[y,x] = -(Depth[y+1,x]-Depth[y-1,x])/(dyF[y+1,x]-dyF[y-1,x])
            elif y == 0:
                Sx[y,x] = -(Depth[y,x+1]-Depth[y,x-1])/(dxF[y,x+1]-dxF[y,x-1])
                Sy[y,x] = -(Depth[1,x]-Depth[0,x])/(dxF[1,x]-dxF[0,x])
            elif y == (len(Y)-1):
                Sx[y,x] = -(Depth[y,x+1]-Depth[y,x-1])/(dxF[y,x+1]-dxF[y,x-1])
                Sy[y,x] = -(Depth[(len(Y)-1),x]-Depth[y-1,x])/(dyF[(len(Y)-1),x]-dyF[y-1,x])
            elif res == 36 and y == 69:
                Sx[y,x] = -(Depth[y,x+1]-Depth[y,x-1])/(dxF[y,x+1]-dxF[y,x-1])
                Sy[y,x] = -(Depth[y+1,x]-Depth[y-1,x])/(dyF[y+2,x]-dyF[y,x])
            else:
                Sx[y,x] = -(Depth[y,x+1]-Depth[y,x-1])/(dxF[y,x+1]-dxF[y,x-1])
                Sy[y,x] = -(Depth[y+1,x]-Depth[y-1,x])/(dyF[y+1,x]-dyF[y-1,x])
                
    return Sx,Sy

def Topostrophy(Uvel,Vvel,norm=False):
    ### Calculation of topostrophy
    nt = Uvel.shape[0]
    if Uvel.shape[3]==421:
        res = 18
        nx = 420
        ny = 384
        nz = 50
    elif Uvel.shape[3]==211:
        res = 36
        nx = 210
        ny = 192
        nz = 50
    else:
        print "Dimension Error"
        
    if res == 18:
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/results18km/grid.nc",'r')
    elif res == 36:
        file2read = netcdf.NetCDFFile("/scratch/general/am8e13/results36km/grid.nc",'r')
    hfacc=file2read.variables['HFacC']
    hfacc=hfacc[:]*1
    Uvel = Uvel[:,:,0:ny,0:nx]
    Vvel = Vvel[:,:,0:ny,0:nx]
    Uvel[:,hfacc==0] = np.nan
    Vvel[:,hfacc==0] = np.nan
    
    Sx,Sy = BathyGrad(res,True)
    Snorm = np.sqrt(Sx**2 + Sy**2)
    velnorm = np.sqrt(Uvel**2 + Vvel**2)
    tau = np.zeros_like(velnorm)
    tau_opposit = np.zeros_like(velnorm)

    for t in range(nt):
        for z in range(nz):
            if norm == True:
                tau[t,z,:,:] = (-Uvel[t,z,:,:]*Sy + Vvel[t,z,:,:]*Sx) / (velnorm[t,z,:,:]*Snorm)
                tau_opposit[t,z,:,:] = (Uvel[t,z,:,:]*Sx + Vvel[t,z,:,:]*Sy) / (velnorm[t,z,:,:]*Snorm) 
            else:
                tau[t,z,:,:] = (-Uvel[t,z,:,:]*Sy + Vvel[t,z,:,:]*Sx)
                tau_opposit[t,z,:,:] = (Uvel[t,z,:,:]*Sx + Vvel[t,z,:,:]*Sy)
    return tau, tau_opposit
