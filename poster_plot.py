from astropy.io import fits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
matplotlib.rc('text', usetex=True)
matplotlib.rc('font',size=24)
import sys

def line_load():
    wave_tab_raw = np.loadtxt('line_table.dat',dtype=bytes,skiprows=1,delimiter=',').astype(str)
    num_lines = len(wave_tab_raw)
    wave_table = np.zeros(num_lines,dtype=[('NAME_SHORT','U4'),('DISP_NAME','U5'),
                                           ('REST_WAVE','f4'),('OBS_WAVE','f4')])
    wave_table['NAME_SHORT'] = wave_tab_raw[:,0]
    wave_table['DISP_NAME'] = wave_tab_raw[:,1]
    wave_table['REST_WAVE'] = wave_tab_raw[:,2].astype(np.float32)

    return wave_table

def file_plot(infile):
    spec_hdr = fits.open(infile)[0].header
    spec_array = fits.open(infile)[1].data
    logwave = spec_array['loglam'] #Wavelength as log10(lambda)
    wave_arr = spec_array['lam'] #Recast the wavelengths in Angstroms
    flux_raw = spec_array['flux_raw'] #Flux at a given wavelength. Units of 10^(-17) erg/s/cm^(2)/Ang
    flux_err = spec_array['ivar']
    zobj = float(spec_hdr['Z_VI'])
    obj_class = int(spec_hdr['OBJ_CLAS'])

    #Fit the continuum to a polynomial using Legendre polynomials
    box_continuum = spec_array['cont_flux']
    param_lgd = np.array([float(spec_hdr['LGD_PRM0']),float(spec_hdr['LGD_PRM1']),
                          float(spec_hdr['LGD_PRM2']),float(spec_hdr['LGD_PRM3'])])
    continuum_fit = spec_array['cont_fit'] #continuum fit
    flux_reduced = spec_array['flux_reduced'] #Reduce the flux
    err_reduced = spec_array['err_reduced'] #Reduce the ivar

    #Pick out the emission lines
    line_flux,line_err = spec_array['line_flux'],spec_array['line_err']

    line_arr = line_load()
    for i in range(len(line_arr)):
        line_arr['OBS_WAVE'][i] = line_arr['REST_WAVE'][i] * (zobj + 1)
    wlines = np.where((line_arr['OBS_WAVE']>=3700)&(line_arr['OBS_WAVE']<=9000))[0]
    vis_lines = line_arr[wlines]


    x_lower = np.amin(wave_arr)
    x_upper = np.amax(wave_arr)
    vis_range = np.where((wave_arr>=3700)&(wave_arr<=9000))[0]
    y_range = np.amax(flux_raw[vis_range]) - np.amin(flux_raw[vis_range])
    y_pad = float(y_range) / 10
    y_lower = np.amin(flux_raw[vis_range]) - y_pad
    y_upper = np.amax(flux_raw[vis_range]) + y_pad

    #Plot some things for analysis
    fig2,ax2 = plt.subplots(nrows=2,ncols=1,figsize=(10,8),sharex='col')
    ax2[0].plot(wave_arr,flux_raw,color='0.75',linewidth=1.3)
    ax2[0].plot(wave_arr,box_continuum,color='blue',linewidth=1.8)
    ax2[0].plot(wave_arr,continuum_fit,color='orange',linewidth=1.3)
    ax2[1].plot(wave_arr,flux_reduced,color='0.75',linewidth=1.3)
    ax2[1].plot(wave_arr,line_flux,linewidth=1.0,color='black')
    ax2[0].set_xlim((x_lower,x_upper))
    ax2[1].set_xlabel(r'Wavelength (\AA)')
    for i in range(2):
        if i < 1:
            ax2[i].set_ylim((y_lower,y_upper))
            ax2[i].set_ylabel(r'F$_{\lambda}$ ($10^{-17}$ ergs s$^{-1}$ cm$^{-2}$\,\AA$^{-1}$)')
        else:
            ax2[i].set_ylim((0,5))
            ax2[i].set_ylabel(r'Rel. Flux')
        ax2[i].tick_params(axis='both',direction='in')
        ax2[i].tick_params(axis='both',which='minor',direction='in')
        ax2[i].xaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax2[i].yaxis.set_minor_locator(ticker.MultipleLocator(1))

    fig2.tight_layout()
    fig2.subplots_adjust(wspace=0,hspace=0)
    fig2.savefig('poster_star_smooth2.png')

file_in = sys.argv[1]
file_plot(file_in)
