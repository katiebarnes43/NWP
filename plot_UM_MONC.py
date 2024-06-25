#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 15:37:43 2023

@author: robertbeare
"""
import iris
# import iris.plot as iplt
# import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
# import scipy.stats as stats

from netCDF4 import Dataset

import w_power_spectrum_UM_ummonc as sp
import w_power_spectrum_UM_2_ummonc as sp2

# Figure plotting
import figure_plotting_module as fp 
# Data loading
import load_um_monc_data as ld



# UM data location
#um_dir='/Users/robertbeare/Documents/python/UM/um_data_hires/'
um_dir='/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/um_data_hires/'
# Changing SGS coefficients to standard ones (Brown et al 1994)
# um_dir='/Users/robertbeare/Documents/python/UM/um_data_hires/diff_sgs/'
# UM expt
expt='da249'


# MONC data location
# monc_dir='/Users/robertbeare/Documents/python/MONC/monc_data/sul_pat/'
# Changing SGS coefficients to standard ones (Brown et al 1994)
#monc_dir='/Users/robertbeare/Documents/python/MONC/monc_data/sul_pat/diff_sgs/'
monc_dir='/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/monc_data/sul_pat/'
dz=30.0 # Vertical grid length
nz=100 # No. of points in vertical
z_max=3000 # Domain top
z_um_th= np.arange(dz, z_max+dz, dz)

# Resolutions and grid lengths   
resolutions=['1km','500m','250m','200m','100m','50m']
dx_vals=[1000,500,250,200,100,50]

#Times of MONC data to be used
monc_times=['3600','7200','10800','14400']  

n=1 # Plot counter

domain_max=3200 #Half domain width for highest resolution

#nx=w_um.shape[3] # No. of points in x
nx=128
# print("nx ",nx)

#load up data
w_um, w_monc, th_um, th_mean_monc, z_monc, zn_monc= \
ld.load_data(resolutions, expt, um_dir, monc_dir, monc_times, nx)

# Figure counter
n=1 
vert_lev=16 #Vertical level for plotting mid-boundary layer w
  
hours={'1h':4, '2h':8, '3h':12, '4h':16}

Spec_um_arr=np.zeros((np.size(resolutions),len(hours),nx//2))
k_d_um_arr=np.zeros((np.size(resolutions),len(hours)))
kvals_um_arr=np.zeros((np.size(resolutions),nx//2))

th_mean_um_arr=np.zeros((np.size(resolutions),len(hours),np.size(z_um_th)))

Spec_monc_arr=np.zeros((np.size(resolutions),len(hours),nx//2))
k_d_monc_arr=np.zeros((np.size(resolutions),len(hours)))
kvals_monc_arr=np.zeros((np.size(resolutions),nx//2))
for res_index, res in enumerate(resolutions):
    print("Resolution ",res)
    
    dx=dx_vals[res_index] # Horizontal grid length
     
    
    # Directory for plots
    # plot_dir='UM_MONC_plots/'
    #plot_dir='UM_MONC_plots/diff_sgs/' commented out
    plot_dir='/Users/klb236/OneDrive - University of Exeter/4th Year/NWP Internship/UM_MONK/UM_MONC_plots/'
    # plot_dir='UM_MONC_plots/diff_sgs_both/'
    
           
    domain=dx*nx
    x_um = np.arange(-0.5*domain, 0.5*domain, dx)
    y_um=x_um
    

    # Change to 2D coordinate arays
    # Horizontal slices grid
    xx_um, yy_um=np.meshgrid(x_um,y_um)
    # Vertical slice grid
    xx_um_vert, zz_um_th=np.meshgrid(x_um, z_um_th)
    
    
   
    
    # Update plot destination for res folder
    plot_dir +=res+'/'
    
    for hr_index, hr in enumerate(hours):
        
        time_index=hours[hr]
        file_string=hr
        
        title_end_string=file_string+', '+res
        file_string_monc= file_string+'_MONC_Conventional'+res
        file_string_two= file_string+'_UM_MONC_'+res
        file_string += '_UM_Conventional'+res
        
        # 2D Power spectrum
        # UM
        Spec_um, kvals_um, k_d_um, l_d, ww_lev=    \
        sp.ww_spectrum(w_um[res_index,:,:,:,:], time_index, vert_lev, dx)
        
        Spec_um_arr[res_index, hr_index, :]=Spec_um[:]
        k_d_um_arr[res_index, hr_index]=k_d_um
        
        # MONC
        Spec_monc, kvals_monc, k_d_monc, l_d, ww_lev=    \
        sp.ww_spectrum(w_monc[res_index,:,:,:,:], 1, vert_lev, dx)
        
        Spec_monc_arr[res_index, hr_index, :]=Spec_monc[:]
        k_d_monc_arr[res_index, hr_index]=k_d_monc
        
        
        if hr_index==0:
            
            kvals_um_arr[res_index,:]=kvals_um[:]
            kvals_monc_arr[res_index,:]=kvals_monc[:]
        # # 1D Power spectrum
        # Spec_1d, kvals_1d, k_d_1d, l_d_1d =    \
        # sp.ww_spectrum_1d(w_um[res_index,:,:,:,:], time_index, vert_lev, dx)
        
        # #Normalising factor for inertial sub-range
        # ind=12
        # fact=Spec[ind]/ (kvals[ind]**(-5/3))
        # inertial=fact*kvals[:]**(-5/3)
        
         #<ww> 
         #for k in range(nz):  
        #     w_section=w_um[res_index,time_index,k,:,:].data
         #    ww[k]=np.mean(w_section**2)
        
        # # <th>
        th_data=th_um[res_index,time_index,:,:,:].data
        th_mean_um_arr[res_index,hr_index,:]=np.mean(th_data, axis=(1,2))
        
        
        
        # Vertical cross section of w
  
        
        # Horizontal cross section of w
        n=fp.plot_w_horizontal_lim(w_um[res_index,:,:,:,:], xx_um, yy_um, vert_lev, time_index,
                              n, plot_dir, file_string, title_end_string+' UM', 1)
        n=fp.plot_w_horizontal_lim(w_monc[res_index,:,:,:,:], xx_um, yy_um, vert_lev, hr_index,
                              n, plot_dir, file_string_monc, title_end_string+' MONC', 1)

        # Vertical velocity spectra
        n=fp.plot_w_spectra(kvals_um_arr[res_index,:], Spec_um_arr[res_index, hr_index,:],
                            k_d_um_arr[res_index, hr_index], n, plot_dir, file_string, 
                            title_end_string+' UM', 1)

        n=fp.plot_w_spectra(kvals_monc_arr[res_index,:], Spec_monc_arr[res_index, hr_index,:],
                            k_d_monc_arr[res_index, hr_index], n, plot_dir, file_string_monc, 
                            title_end_string+' MONC', 1)
        
        n=fp.plot_w_spectra_two(kvals_monc_arr[res_index,:], Spec_um_arr[res_index, hr_index,:],
                            Spec_monc_arr[res_index, hr_index,:],
                            k_d_monc_arr[res_index, hr_index], n, plot_dir, file_string_two, 
                            title_end_string+'', 1)
        # Mean potential temperature profile
        
        n=fp.plot_th_prof( th_mean_um_arr[res_index, hr_index,:], z_um_th,
                          n, plot_dir, file_string, title_end_string, 0, 'k',
                          0, [''])
        
        z_monc_th=np.zeros(101)
        z_monc_th[1:101]=z_um_th[0:100]
        
        n=fp.plot_th_prof( th_mean_monc[res_index,hr_index,1:100], z_monc_th[1:100], n, plot_dir, 
                          file_string_monc, title_end_string, 1, 'r.',
                          1,['UM','MONC'])
        
        # <ww> profiles
        n=fp.plot_ww_prof(w_um[res_index,:,:,:,:], z_um_th, nz, time_index, n, 
                         plot_dir, file_string, title_end_string, 0, 
                         'k', 0, [''])
        # Sort out vertical levels
        n=fp.plot_ww_prof(w_monc[res_index,:,:,:,:], z_monc_th[0:100], nz, hr_index, n, 
                         plot_dir, file_string, title_end_string, 1, 
                         'r.',1, ['UM','MONC'])
        
        
        
#         # Additional plot with Coarse-graining
#         if res == '1km' and hr=='6h':
#             Spec_1km_CG_250=np.load('ww_spec_CG_to_1km_from_250m.npy')
#             kvals_CG=np.load('ww_spec_CG_to_1km_from_250m_kvals.npy')
#             plt.figure(n)
#             plt.loglog(kvals, Spec)
#             plt.loglog(kvals_CG, Spec_1km_CG_250,'b--')
#             # plt.loglog(kvals_1d, Spec_1d,'.')
#             # plt.loglog(kvals, inertial, '--')
#             # plt.loglog([k_d, k_d],[1e-1,1e5], 'k--')
#             plt.xlabel("$k$")
#             plt.ylabel("$S(k)$")
#             plt.legend(['1 km Power spectrum', '250 m coarse-grained to 1 km'],frameon=False)
#             plt.title("Vertical velocity power spectrum, including CG at "+title_end_string)
#             plt.tight_layout()
#             plt.savefig(plot_dir+'w_power_spectrum_plus_CG_from_250m_'+file_string+'.pdf')
#             n+=1
      
#         plt.figure(n)
#         plt.plot(ww, z_um_th)
#         plt.ylabel("z (m)")
#         plt.xlabel("<ww> (m$^{2}$s$^{-2}$)")
#         plt.title("Vertical velocity variance at "+title_end_string)
#         plt.savefig(plot_dir+'ww_vert_profile_'+file_string+'.pdf')
#         n+=1

#         plt.figure(n)
#         plt.plot(th_mean, z_um_th)
#         plt.ylim((0, 3000))
#         plt.xlim((303, 315))
#         plt.ylabel("z (m)")
#         plt.xlabel("Potential temperature (K)")
#         plt.title("Potential temperature at "+title_end_string)
#         plt.savefig(plot_dir+'th_mean_vert_profile_'+file_string+'.pdf')
#         n+=1

#         if res == '50m':
#             ww_50m=ww
#         if res == '1km':
#             ww_1km=ww
            
#         if res == '250m':
#             ww_250m=ww
#             """ Coarse graining/binning of a 2D array
#                 For a factor of 2
#                 w_2d_coarse_grained=w_2d.reshape(nx//2,2,nx//2,2).mean(axis=3).mean(1)
#             """ 
            
#             w_2d_250= w[time_index,16,:,:].data
#             w_2d_500_cg=w_2d_250.reshape(64,2,64,2).mean(axis=3).mean(1)
#             w_2d_1000_cg=w_2d_500_cg.reshape(32,2,32,2).mean(axis=3).mean(1)
            

#             xx_um_500= xx_um.reshape(64,2,64,2).mean(3).mean(1)
#             yy_um_500= yy_um.reshape(64,2,64,2).mean(3).mean(1)
            
#             xx_um_1000= xx_um_500.reshape(32,2,32,2).mean(3).mean(1)
#             yy_um_1000= yy_um_500.reshape(32,2,32,2).mean(3).mean(1)
            
#             Spec_1km_CG, kvals, k_d, l_d, ww_lev=    \
#             sp2.ww_spectrum(w_2d_1000_cg, 1000)
            
#             # np.save('ww_spec_CG_to_1km_from_250m.npy',Spec_1km_CG)
#             # np.save('ww_spec_CG_to_1km_from_250m_kvals.npy',kvals)
            
#             plt.figure(n)
#             plt.gcf().set_size_inches(5, 6.5)
#             plt.contourf(xx_um_500, yy_um_500, w_2d_500_cg, cmap=contour_map)
#             plt.ylabel("y (m)")
#             plt.xlabel("x (m)")
#             plt.title("CG Vertical velocity (ms$^{-1}$), height 510 m at 500 m")
#             plt.colorbar(location='bottom')
#             plt.xlim((-domain_max, domain_max))
#             plt.ylim((-domain_max, domain_max))
#             plt.savefig(plot_dir+'w_hor_x_section_2_COARG_500m.pdf')
#             n+=1
            
#             plt.figure(n)
#             plt.gcf().set_size_inches(5, 6.5)
#             plt.contourf(xx_um_1000, yy_um_1000, w_2d_1000_cg, cmap=contour_map)
#             plt.ylabel("y (m)")
#             plt.xlabel("x (m)")
#             plt.title("CG Vertical velocity (ms$^{-1}$), height 510 m at  1 km")
#             plt.colorbar(location='bottom')
#             plt.xlim((-domain_max, domain_max))
#             plt.ylim((-domain_max, domain_max))
#             plt.savefig(plot_dir+'w_hor_x_section_2_COARG_1km.pdf')
#             n+=1
            
#             plt.figure(n)
#             plt.loglog(kvals, Spec_1km_CG)
            
#             # plt.loglog(kvals, inertial, '--')
#             # plt.loglog([k_d, k_d],[1e-1,1e5], 'k--')
#             plt.xlabel("$k$")
#             plt.ylabel("$S(k)$")
#             plt.title("Vertical velocity power spectrum, height 510 m at CG to 1km")
#             plt.tight_layout()
#             plt.savefig(plot_dir+'w_power_spectrum_COARG_1km.pdf')
#             n+=1
            

# plt.figure(n)
# plt.plot(ww_50m, z_um_th)
# plt.plot(ww_250m, z_um_th,'r')
# plt.plot(ww_1km, z_um_th,'k')
# plt.legend(['50m', '250m', '1km'],frameon=False)
# plt.ylabel("z (m)")
# plt.xlabel("<ww> (m$^{2}$s$^{-2}$)")
# plt.title("Vertical velocity variance at 4h")
# plt.savefig(plot_dir+'ww_vert_profile_different_resolutions'+file_string+'.pdf')
# n+=1            