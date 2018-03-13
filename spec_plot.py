from astropy.io import fits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('text', usetex=True)
import sys

def p_spec(infile,p_in,m_in,f_in,smoothing_pct=0,ogs=1):
    loglam = infile['loglam'] #Extract log values for the x-axis values.
    lam = 10**loglam #Define the x-axis values for the plot.
    flux = infile['flux']

    #Define the axes limits
    x_lower = np.amin(lam)
    x_upper = np.amax(lam)
    vis_range = np.where((lam>=3700)&(lam<=9000))[0]
    y_range = np.amax(flux[vis_range]) - np.amin(flux[vis_range])
    y_pad = float(y_range) / 10
    y_lower = np.amin(flux[vis_range]) - y_pad
    y_upper = np.amax(flux[vis_range]) + y_pad

    if smoothing_pct != 0:
        num_x = len(lam)
        step_smooth = int(num_x / smoothing_pct)
        num_smooth = int(num_x/step_smooth)
        xs_adr = np.arange(0,num_x,num_smooth)
        x_smooth = lam[xs_adr]
        y_smooth = flux[xs_adr]


    fig1, ax1 = plt.subplots(figsize=(10,8))
    ax1.set_xlim((x_lower,x_upper))
    ax1.set_ylim((y_lower,y_upper))
    ax1.tick_params(axis='both',direction='in')
    ax1.tick_params(axis='both',which='minor',direction='in')
    if ogs == 1:
        ax1.plot(lam,flux,color='grey',linewidth=0.6,label='Original Spectrum')
    if smoothing_pct != 0:
        ax1.plot(x_smooth,y_smooth,color='black',linewidth=0.8,label='Smoothed Spectrum ({})'.format(smoothing_pct))
    ax1.set_xticks([4000,5000,6000,7000,8000,9000,10000])
    ax1.set_xticklabels(ax1.get_xticks(),fontsize=15)
    ax1.set_yticklabels(ax1.get_yticks(),fontsize=15)
    ax1.set_xlabel('Wavelength (\AA)',fontsize=15)
    ax1.set_ylabel('$f_{\lambda}$ ($10^{-17}$ ergs s$^{-1}$ cm$^{-2}$\,\AA$^{-1}$)',fontsize=15)
    plt.legend()
    plt.minorticks_on()
    plt.tight_layout()
    plt.show()

def file_open(ifile,plate,mjd,fiberid,spct=0,ogspec=1):
    spec_file = fits.open(ifile)[1].data
    p_spec(spec_file,plate,mjd,fiberid,spct,ogspec)

check_smooth = int(sys.argv[1])
og_check = int(sys.argv[2])
#ptest, mtest, ftest = 7294, 56739, 21
#test_file = '../ica_spec/7294/spec-7294-56739-0021.fits'
spec_file_in = sys.argv[3]
pin = int(sys.argv[4])
mmin = int(sys.argv[5])
fin = int(sys.argv[6])
#file_open(test_file,ptest,mtest,ftest,check_smooth,og_check)
file_open(spec_file_in,pin,mmin,fin,check_smooth,og_check)
