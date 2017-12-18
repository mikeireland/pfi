""" 24.1 to 24.8 microns is Q 

Hi-5:

snr.maglim([3.5,4.1], [290], diam=8, eta_c=0.4, eta_w=0.25, t_int=3600, snr=5, n_tel=4)
"""

#bp=[3.3,4.2]
#"{0:5.1f} & {1:5.1f} & {2:5.1f}".format(maglim(bp,280),maglim(bp,270),maglim(bp,220))
from __future__ import division, print_function
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
#import astropy.units as ap_u
import astropy.constants as ap_c
import pdb

#Choose one of the following lines for 10 or 1mm PWV
#trans = np.concatenate((np.loadtxt('cptrans_zm_100_10.dat'),np.loadtxt('cptrans_nq_100_10.dat')))
trans = np.concatenate((np.loadtxt('mktrans_zm_10_10.dat'),np.loadtxt('mktrans_nq_10_10.dat')))


lmin = 10.0 #In microns
lmax = 11.0 #In microns

def vega_photonrate(wave_in):
    """Returns the photon rate based on vega magnitude, calibrated from 
    Tokunga and Vacca (2005, PASP, 117, 421), Hewett et al. 2006, 
    Cohen, Walker, Barlow, and Deacon, 1992, Astronomical Journal, vol. 104, 1650-1657
    and a long-wavelength extrapolation:
    m_AB(lambda1)-m_AB(lambda0) = 5 log_10 (lambda1/lambda0), from the
    Rayleigh-Jeans tail. """
    dwave = 1.0
    try:
        if len(wave_in)> 2:
            print("Error: min and max wavelengths or central wavelength only")
        elif len(wave_in)==2:
            dwave = np.abs(0.5*(wave_in[0] - wave_in[1]))
            wave = 0.5*(wave_in[0] + wave_in[1])
        else:
            wave = wave_in[0]
    except:
        wave = wave_in
    
    filter_waves = [0.355,0.467,0.5466,0.616,0.747,0.892,1.031,1.248,1.644,2.121,2.149,2.198,3.754,4.702,10.2,30]
    filter_ab = [0.927,-0.103,0.026,0.146,0.366,0.533,0.634,0.940,1.38,1.84,1.86,1.90,2.94,3.40,4.98,7.32] 
    
    #First, convert the vega magnitude into AB magnitude.
    vega_ab = np.interp(np.log10(wave),np.log10(filter_waves), filter_ab)
    
    #Now, convert to Jansky.
    flux_jy = 3631.0*10**(-0.4*vega_ab)
    
    #Now convert to photons/s/m^2/micron
    wave_m = wave/1e6
    photon_rate = flux_jy*1e-26*ap_c.c.value/wave_m**2/(ap_c.h.value*ap_c.c.value/wave_m)*1e-6
    
    return photon_rate
    
def atm_throughput(wave_in):
    """Find the mean atmospheric throughput. NB Doesn't work between 5.6 and 7 microns"""
    ww = np.where( (trans[:,0] > wave_in[0]) * (trans[:,0] < wave_in[1]) )[0]
    return np.mean(trans[ww,1])
    
def snr_heterodyne_constant(wave_in, printit=True):
    wave = np.mean(wave_in)*1e-6
    #NB the following constant is per polarisation.
    snr_const = atm_throughput(wave_in)*vega_photonrate(wave_in)*np.mean(wave_in)/(ap_c.c.value/wave)/2
    return snr_const
    
def nulling_requirement(wave_in, twarm, diam):
    """Figure out the magnitude equivalent of the background"""
    wave = np.mean(wave_in)*1e-6
    
    try:
        n_temps=len(twarm)
    except:
        twarm=[twarm]
        
    mags=[]
    for t in twarm:
        #Background in photons/mode/Hz of bandwidth/s
        background = 2.0/(np.exp(ap_c.h.value*ap_c.c.value/ap_c.k_B.value/t/wave) - 1)

        #Background in photons/mode/unit fractional bandwidth
        background *= ap_c.c.value/wave
    
        #Signal per unit fractional bandwidth
        signal = atm_throughput(wave_in)*vega_photonrate(wave_in)*np.mean(wave_in)*np.pi*(diam/2.0)**2
        
        mags.append(2.5*np.log10(signal/background))
    return mags
    
def snr_direct_constant(wave_in, twarm, printit=True):
    """Find the constant for the direct detection SNR equation as a function of
    wavelength and warm optics temperature. 
    
    Note that the warm optics throughput has to include loss in coupling to the beam 
    combiner 'modes', which could be e.g. just a limited square field of view. A
    complete calculation would have to include the effects of cold pupil stops etc, 
    which may make such a coupling loss not as extreme as a warm optics efficiency loss.
    
    Similarly, atmospheric throughput loss is a 'warm' optics loss.
    
    snr_direct_constant([10,11],270)*np.pi*2**2*np.sqrt(0.44/(1-0.48))*0.48*np.sqrt(0.1/10.5)*10**(-0.4*13.3)*sqrt(100*3600.)*sqrt(20)
    """
    
    wave = np.mean(wave_in)*1e-6
    
    try:
        n_temps=len(twarm)
    except:
        twarm=[twarm]
        
    snrs=[]
    for t in twarm:
        #Background in photons/mode/Hz of bandwidth/s
        background = 2.0/(np.exp(ap_c.h.value*ap_c.c.value/ap_c.k_B.value/t/wave) - 1)
    
        #Background in photons/mode/unit fractional bandwidth
        background *= ap_c.c.value/wave
    
        #Signal per unit fractional bandwidth
        signal = atm_throughput(wave_in)*vega_photonrate(wave_in)*np.mean(wave_in)
    
        #Background limited SNR
        snrs.append(signal/np.sqrt(background))
        
    if printit:
        fmt = "{" + ":8.2e} & {".join([str(r) for r in range(len(snrs))]) + ":8.2e} \\\\"
        print(fmt)
        print(fmt.format(*snrs))
        
    if len(snrs)==1:
        return snrs[0]
    else:
        return snrs
    
def eta_func(eta_w,eta_c):
    """Return the component of the SNR equation due to efficiency"""
    return eta_w*np.sqrt(eta_c/(1-eta_w))
    
def maglim(wave_in, twarm, diam=2.5, A=None, eta_c=0.44, eta_w=0.39, t_int=1e4, snr=5.0, n_tel=12,printit=True):
    """Find the point-source detection magnitude limit for PFI based on a series of assumptions"""
    if not A:
        A = np.pi*(diam/2)**2
    frac_wave = np.abs(wave_in[1]-wave_in[0])/0.5/(wave_in[0] + wave_in[1])
    maglims = []
    for t in twarm:
        magfunc = snr/(snr_direct_constant(wave_in, t,printit=False)*A*eta_func(eta_w,eta_c)*\
            np.sqrt(t_int)*np.sqrt(frac_wave)*np.sqrt(n_tel-1))
        maglims.append(-2.5*np.log10(magfunc))
    
    if printit:
        fmt = "{" + ":5.1f} & {".join([str(r) for r in range(len(maglims))]) + ":5.1f} \\\\"
        print(fmt)
        print(fmt.format(*maglims))
        
    if len(maglims)==1:
        return maglims[0]

def maglim_heterodyne(wave_in, diam=2.5, A=None, eta_w=0.35, t_int=1e4, snr=5.0, n_tel=12,
        printit=True, polz='dual',const='new'):
    """Find the point-source detection magnitude limit for PFI based on a series of assumptions"""
    if not A:
        A = np.pi*(diam/2)**2
#    if (wave_in[0] < 8) | (wave_in[1]>13):
#        print("ERROR: N band only for now...")
#        raise UserWarning
    frac_wave = np.abs(wave_in[1]-wave_in[0])/0.5/(wave_in[0] + wave_in[1])
    mean_wave = np.mean(wave_in)*1e-6
    delta_nu = frac_wave * ap_c.c.value/mean_wave
    if const=='new':
        magfunc = snr/(snr_heterodyne_constant(wave_in)*A*eta_w*np.sqrt(2*t_int*delta_nu)*(n_tel-1)) 
    else:
        magfunc = snr/(9.8e-6*A*eta_w*np.sqrt(2*t_int*delta_nu)*(n_tel-1))
    if polz == 'dual':
        magfunc /= np.sqrt(2)
    
    return -2.5*np.log10(magfunc)

        
def planck_photons_beam(temperature,wave=10.5e-6, omega=1e-16, target=0):
    return 2.0*ap_c.c.value/wave**4/(np.exp(ap_c.h.value*ap_c.c.value/ap_c.k_B.value/temperature/wave) - 1)*1e-6*omega - target
        
def surface_brightness(wave_in, baselines=[1e3,2e3,4e3],twarm=280, diam=2.5, A=None, \
    eta_c=0.44, eta_w=0.39, t_int=1e4, snr=3.0, n_tel=12,printit=True, type='direct'):
    """Find the point-source detection magnitude limit for PFI based on a series of assumptions"""
    if not A:
        A = np.pi*(diam/2)**2
    frac_wave = np.abs(wave_in[1]-wave_in[0])/0.5/(wave_in[0] + wave_in[1])
    mean_wave = np.mean(wave_in)*1e-6
    delta_nu = frac_wave * ap_c.c.value/mean_wave
    t_surf = []
    for bl in baselines:
        if type=='direct':
            magfunc = snr/(snr_direct_constant(wave_in, twarm,printit=False)*A*eta_w*np.sqrt(eta_c/(1-eta_w))*\
                np.sqrt(t_int)*np.sqrt(frac_wave)*np.sqrt(n_tel-1))
        else:
            #Just use the formula from Ireland (2014)
            magfunc = snr/(9.8e-6*A*eta_w*np.sqrt(2*t_int*delta_nu)*(n_tel-1))
        photons_per_beam = vega_photonrate(wave_in)*magfunc
        omega = np.pi*(mean_wave/bl/2)**2
        t_surf.append(op.bisect(planck_photons_beam, 20, 5000, args=(mean_wave,omega,photons_per_beam)))
    
    if printit:
        fmt = "{" + ":5.0f} & {".join([str(r) for r in range(len(t_surf))]) + ":5.0f} \\\\"
        print(fmt)
        print(fmt.format(*t_surf))
        
    if len(t_surf)==1:
        return t_surf[0]
    else:
        return t_surf 

if __name__=="__main__":
    maglim_1hr = maglim([3.4,4.0], [285], diam=8.0, eta_c=0.4, eta_w=0.25, t_int=3600*2, snr=5, n_tel=4)-2.5*np.log10(2)
    print("Magnitude limit for Hi-5 for 1 hour integration: {:5.2f}".format(maglim_1hr))