import numpy as np

def region_mask(res):
    if res == 36:
        k = 1 ; nx = 210 ; ny = 192 ; nz = 50 ;
    elif res == 18:
        k = 2 ; nx = 420 ; ny = 384 ; nz = 50 ;
    elif res == 9:
        k = 4 ; nx = 840 ; ny = 768 ; nz = 50 ;

    # this function gives the region masks                                                                                                   
    # mask 1: Arctic
    mask1 = np.zeros([nz,ny,nx]) ;
    mask1[:,:,:] = np.nan
    mask1[:,10*k:40*k,53*k:75*k] = 1 ;
    mask1[:,40*k:58*k,68*k:75*k] = 1 ;
    mask1[:,:,75*k:178*k] = 1 ;
    # mask 2:  Nordic Seas                                                                                                                   
    mask2 = np.zeros([nz,ny,nx]) ;
    mask2[:,:,:] = np.nan
    mask2[:,40*k:100*k,36*k:68*k] = 1 ;
    mask2[:,50*k:90*k,25*k:40*k] = 1 ;
    mask2[:,58*k:84*k,68*k:75*k] = 1 ;
    # mask 3 : North Atlantic                                                                                                                
    mask3 = np.zeros([nz,ny,nx]) ;
    mask3[:,:,:] = np.nan
    mask3[:,99*k:191*k,0:75*k]= 1 ;
    mask3[:,50*k:90*k,0:25*k] = 1 ;
    mask3[:,90*k:100*k,0:30*k] =1 ;
    
    # mask 4 : Lofoten
    mask4 = np.zeros([nz,ny,nx]) ;
    mask4[:,:,:] = np.nan
    mask4[:,58*k:79*k,23*k:38*k] = 1
    
    # mask 5 : Greenland
    mask5 = np.zeros([nz,ny,nx]) ;
    mask5[:,:,:] = np.nan
    mask5[:,60*k:80*k,52*k:65*k] = 1
    
    # mask 6 : Norwegian
    mask6 = np.zeros([nz,ny,nx]) ;
    mask6[:,:,:] = np.nan
    mask6[:,52*k:67*k,40*k:52*k] = 1    
    
    # mask 7 : Iceland
    mask7 = np.zeros([nz,ny,nx])
    mask7[:,:,:] = np.nan
    mask7[:,82*k:98*k,41*k:51*k] = 1

    # mask 8 : Labrador
    mask8 = np.zeros([nz,ny,nx])
    mask8[:,:,:] = np.nan
    mask8[:,135*k:158*k,38*k:61*k] = 1

    return mask1,mask2,mask3,mask4,mask5,mask6,mask7,mask8

def reg_dic(res):
    mask1,mask2,mask3,mask4,mask5,mask6,mask7,mask8 = region_mask(res)
    regions = {0 : np.ones_like(mask1) , 1 : mask1 , 2 : mask2 , 3 : mask3 , 4 : mask4,\
               5 : mask5, 6 : mask6, 7 : mask7 , 8 : mask8 }
    return regions

def reg_titles():
    region_titles = {0 : 'Global' , 1 : 'Arctic' , 2 : 'Nord Seas' , 3 : 'North Atl', 4: 'Lofoten' ,\
                     5: 'Greenland' , 6 : 'Norwegian' , 7 : 'Iceland Sea' , 8 : 'Labrador'}
    return region_titles

def mask(reg,res):
    mask1,mask2,mask3,mask4,mask5,mask6 = region_mask(res)
    if reg == 1:
        return mask1
    elif reg == 2:
        return mask2
    elif reg == 3:
        return mask3
    elif reg == 4:
        return mask4
    elif reg == 5:
        return mask5   
    elif reg == 6:
        return mask6
    elif reg == 7:
        return mask7
    elif reg == 8:
        return mask8
    else:
        print "Error: invalid mask number"
