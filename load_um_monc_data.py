""" Module to  load up um and monc data
    into combined (stacked) arrays """
import iris
import numpy as np
from netCDF4 import Dataset
    
def load_data(resolutions, expt, um_dir, monc_dir, monc_times, nx):
    
      
    # Load UM data
    # stash_names dictionary
    stash_names = {'w': 'm01s00i150', 'th': 'm01s00i004'}
    stream='r'
    
    
    for counter, res in enumerate(resolutions):
        #Load UM data
        filename=um_dir+expt+'_'+res+'a_p'+stream+'000.nc' 
        names = {'w', 'th'}
        
        variablelist = {}
        for name in names:
            if name in stash_names:
                variablelist[name] = iris.load_cube(filename,
                                                    iris.AttributeConstraint(STASH=stash_names[name]))
            else:
                variablelist[name] = iris.load_cube(um_dir+ filename)
        
        
        w=variablelist['w']
        th=variablelist['th']
        w_dat=w.data
        th_dat=th.data
        
        if counter==0:
            # Expand array to include resolution as first dimension
            w_um=np.zeros((np.size(resolutions),)+np.shape(w_dat))
            th_um=np.zeros((np.size(resolutions),)+np.shape(th_dat))
        
            
        w_um[counter,:,:,:,:]=w_dat[:,:,:,:]
        th_um[counter,:,:,:,:]=th_dat[:,:,:,:]
            
        # Load MONC data 
        # At 4 h= 14400s
        for nm, monc_t in enumerate(monc_times):
            nc_file=monc_dir+res+'/cbl_'+monc_t+'.nc'  
            nc=Dataset(nc_file,'r')
            # u_mean=nc.variables['u_wind_mean'][:,:]
            th_mean=nc.variables['theta_mean'][:,:]
            # th_mean_2=nc_2.variables['theta_mean'][:,:]
            zn_monc=nc.variables['zn'][:]
            z_monc=nc.variables['z'][:]
    
            w_m=nc.variables['w'][:,:,:,:]
            w_mo=np.swapaxes(w_m, 1, 3)
            # print(counter,nm,np.shape(w_monc),np.shape(w_mo),
            #       np.shape(th_mean_monc), np.shape(th_mean))
            
            if counter==0 and nm==0:
                # Expand array to include resolution as first dimension
                w_monc=np.zeros((np.size(resolutions),np.size(monc_times), 
                               np.size(zn_monc), nx, nx))
                th_mean_monc=np.zeros((np.size(resolutions),np.size(monc_times),
                                np.size(zn_monc)))
            
            w_monc[counter, nm,:,:,:]=w_mo[0,:,:,:]
            th_mean_monc[counter, nm,:]=th_mean[0,:]
        
 
        
    
    
        
    
    return w_um, w_monc, th_um, th_mean_monc, z_monc, zn_monc

#comment out theta
