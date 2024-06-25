"""
Module to calculate power spectra on horizontal
cross sections of the 4D vertical velocity, w, a 
time index time_int and vertical level vert_lev 
Based on:

https://bertvandenbroucke.netlify.app/2019/05/24/computing-a-power-spectrum-in-python/

1. Take the Fourier transform of the 2D vertical
velocity slice (using np.fft.fft2). 
Then take it complex modulus squared to form
fourier_amplitudes.
3. Construct a corresponding array of wave 
vectors k with the same layout as the Fourier 
amplitudes (using np.fft.fftfreq and
            np.meshgrid) . Flaten k into a 1D
array.
4. Bin the amplitudes of the Fourier signal into
 kbins that run up to the Nyquist limit and compute the total variance within each bin.
 (using stats.binned_statistic)

"""

import numpy as np
import scipy.stats as stats

def ww_spectrum(w, time_ind, vert_lev, dx):
    
    # horizontal vertical velocity section
    w_section=w[time_ind,vert_lev,:,:].data

    # No. of points in x
    nx = w_section.shape[0]
    
    # a 2D discrete fft of w_section
    # fourier_image is a complex array with same
    # dimensions as w_section
    fourier_image = np.fft.fft2(w_section)
    # Fourier amplitudes is the complex modulus squared
    # of the fourier image
    fourier_amplitudes = np.abs(fourier_image)**2
    
    # The discrete Fourier transform wavenumbers
    # A 1D array that is linear from 0 to nx/2 up to nx/2
    kfreq = np.fft.fftfreq(nx) * nx
    # Form 2D array from kfreq
    kfreq2D = np.meshgrid(kfreq, kfreq)
    
    # Form a positive nx by nx  array from the ffreq2D(2,nx,nx)
    knrm = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)

    # Flatten nx by nx 2D arrays into a nx*nx 1D array
    knrm = knrm.flatten()
    fourier_amplitudes = fourier_amplitudes.flatten()

    #Wave number bits up to the Nyquist limit
    # kbins runs from 0.5 up nx//2 +1
    kbins = np.arange(0.5, nx//2+1, 1.)
    
    # Average neibouring kbins to 
    # make k values from 1 to nx//2
    kvals = 0.5 * (kbins[1:] + kbins[:-1])

    # Normalise k so that max value is pi/dx 
    kvals *= np.pi/(dx*0.5*nx) 

    #Wavenumber increment
    dk=kvals[5]-kvals[4]

    #Bin the fourier ampltiudes into the nx//2 kbins
    binned_fourier, _, _ = stats.binned_statistic(knrm, fourier_amplitudes,
                                         statistic = "mean",                                        
                                         bins = kbins)
    # Normalise binned fourier amplitude
    binned_fourier *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)

    
    # Vertial velocity variance
    ww_lev=np.mean(w_section**2)

    # Integral of power spectrum
    int_ps=np.sum(binned_fourier)*dk
    
    # Power spectrum normalised to givevertical
    # variance when integrated
    Spec = binned_fourier * ww_lev/int_ps
    int_ps=np.sum(Spec)*dk
    
    #Dissipation wavenumber and length scale
    k_d= np.sqrt( np.sum((kvals**2)*Spec)*dk/int_ps)
    l_d=2*np.pi/k_d


    return Spec, kvals, k_d, l_d, ww_lev

""" Construct the spectrum from
    1d power spectra """
    
def ww_spectrum_1d(w, time_ind, vert_lev, dx):
 
    # horizontal vertical velocity section
    w_section=w[time_ind,vert_lev,:,:].data

    # No. of points in x
    nx = w_section.shape[0]
    
    mean_fourier_amplitudes=np.zeros(nx//2-1)
    for j in range(nx):
        # a 1D discrete fft of w_section
        fourier_image = np.fft.fft(w_section[j,:])
        
        # Fourier amplitudes is the complex modulus squared
        # of the fourier image
        fourier_amplitudes = np.abs(fourier_image[1:nx//2+1])**2
        
        mean_fourier_amplitudes += fourier_amplitudes[1:nx//2+1]/nx


    #Wave number bits up to the Nyquist limit
    kvals = np.arange(2, nx//2+1, 1.)
    
    # Normalise k so that max value is pi/dx 
    kvals *= np.pi/(dx*0.5*nx) 

    #Wavenumber increment
    dk=kvals[5]-kvals[4]

    # Vertial velocity variance
    ww_lev=np.mean(w_section**2)
    
    # Set to match the 2D spectrum
    # mean_fourier_amplitudes *= kvals* mean_fourier_amplitudes

    # Integral of power spectrum
    int_ps=np.sum(mean_fourier_amplitudes)*dk
    
    # Power spectrum normalised to give vertical
    # variance when integrated
    
    Spec = mean_fourier_amplitudes * ww_lev/int_ps
    int_ps=np.sum(Spec)*dk
    
    #Dissipation wavenumber and length scale
    k_d= np.sqrt( np.sum((kvals**2)*Spec)*dk/int_ps)
    l_d=2*np.pi/k_d


    return Spec, kvals, k_d, l_d