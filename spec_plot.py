from astropy.io import fits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('text', usetex=True)
import sys

def p_spec(infile,p_in,m_in,f_in):
    loglam = infile['loglam'] #Extract log values for the x-axis values.
    lam = 10**loglam #Define the x-axis values for the plot.
    flux = infile['flux']

    #Define the axes limits
    x_lower = np.amin(lam)
    x_upper = np.amax(lam)
    vis_range = np.where((lam>=3800)&(lam<=9000))[0]
    y_range = np.amax(flux[vis_range]) - np.amin(flux[vis_range])
    y_pad = float(y_range) / 10
    y_lower = np.amin(flux[vis_range]) - y_pad
    y_upper = np.amax(flux[vis_range]) + y_pad

    fig1, ax1 = plt.subplots(figsize=(10,8))
    ax1.set_xlim((x_lower,x_upper))
    ax1.set_ylim((y_lower,y_upper))
    ax1.tick_params(axis='both',direction='in')
    ax1.plot(lam,flux,color='black',linewidth=0.5)
    ax1.set_xticks(np.arange(x_lower,x_upper,6))
    ax1.set_xticklabels(ax1.get_xticks(),fontsize=15)
    ax1.set_yticklabels(ax1.get_yticks(),fontsize=15)
    ax1.set_xlabel('Wavelength (\AA)',fontsize=15)
    ax1.set_ylabel('$f_{\lambda}$ ($10^{-17}$ ergs s$^{-1}$ cm$^{-2}$\,\AA$^{-1}$)',fontsize=15)
    plt.tight_layout()
    plt.show()

def file_open(ifile,plate,mjd,fiberid):
    spec_file = fits.open(ifile)[1].data
    p_spec(spec_file,plate,mjd,fiberid)

ptest, mtest, ftest = 7294, 56739, 21
test_file = '../ica_spec/7294/spec-7294-56739-0021.fits'
file_open(test_file,ptest,mtest,ftest)
