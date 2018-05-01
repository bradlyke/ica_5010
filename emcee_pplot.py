from astropy.io import fits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
matplotlib.rc('text', usetex=True)
matplotlib.rc('font',size=32)
import sys
import emcee
import corner
import scipy.optimize as op

def line_func(x,g,x0,A,m,b):
    return (m*x) + b + (A*np.exp(-((x-x0)**(2))/(2*g**(2))))

def lnlike(theta,x,y,yerr):
    g,x0,A,m,b = theta
    model = line_func(x,g,x0,A,m,b)
    inv_sigma2 = yerr
    inv_model = 1/ model
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprior(theta):
    g,x0,A,m,b = theta
    #g_lower,g_upper = 10.0,100.0
    #x0_lower,x0_upper = 6800.0,7200.0
    #A_lower,A_upper = -1.0,4.0
    #m_lower,m_upper = -0.5,10.0
    #b_lower,b_upper = -4.0,5.0
    if ((g_lower<g<g_upper)&(x0_lower<x0<x0_upper)&(A_lower<A<A_upper)&(m_lower<m<m_upper)&(b_lower<b<b_upper)):
        return 0.0
    return -np.inf

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    lcheck = lp + lnlike(theta, x, y, yerr)
    if np.isnan(lcheck):
        return -np.inf
    return lcheck

def line_load():
    wave_tab_raw = np.loadtxt('line_table.dat',dtype=bytes,skiprows=1,delimiter=',').astype(str)
    num_lines = len(wave_tab_raw)
    wave_table = np.zeros(num_lines,dtype=[('NAME_SHORT','U4'),('DISP_NAME','U5'),
                                           ('REST_WAVE','f4'),('OBS_WAVE','f4')])
    wave_table['NAME_SHORT'] = wave_tab_raw[:,0]
    wave_table['DISP_NAME'] = wave_tab_raw[:,1]
    wave_table['REST_WAVE'] = wave_tab_raw[:,2].astype(np.float32)

    return wave_table

def file_load(infile):
    match_flag = 0
    #g_lower,g_upper = 10.0,100.0
    #x0_lower,x0_upper = 6800.0,7200.0
    #A_lower,A_upper = -1.0,4.0
    #m_lower,m_upper = -0.5,10.0
    #b_lower,b_upper = -4.0,5.0
    global g_lower,g_upper,x0_lower,x0_upper,A_lower,A_upper,m_lower,m_upper,b_lower,b_upper
    g_lower,g_upper = 1.0,50.0
    A_lower,A_upper = -2.0,3.0
    m_lower,m_upper = -1.0,1.0
    b_lower,b_upper = -1.0,1.0

    spec_hdr = fits.open(infile)[0].header
    spec_array = fits.open(infile)[1].data
    logwave = spec_array['loglam'] #Wavelength as log10(lambda)
    wave_arr = spec_array['lam'] #Recast the wavelengths in Angstroms
    flux_raw = spec_array['flux_raw'] #Flux at a given wavelength. Units of 10^(-17) erg/s/cm^(2)/Ang
    flux_err = spec_array['ivar']
    zobj = float(spec_hdr['Z_VI'])
    obj_class = int(spec_hdr['OBJ_CLAS'])
    if zobj < 0.0:
        zobj = 0.0

    wave_shift = wave_arr / (1+zobj)

    #Fit the continuum to a polynomial using Legendre polynomials
    box_continuum = spec_array['cont_flux']
    param_lgd = np.array([float(spec_hdr['LGD_PRM0']),float(spec_hdr['LGD_PRM1']),
                          float(spec_hdr['LGD_PRM2']),float(spec_hdr['LGD_PRM3'])])
    continuum_fit = spec_array['cont_fit'] #continuum fit
    flux_reduced = spec_array['flux_reduced'] #Reduce the flux
    err_reduced = spec_array['err_reduced'] #Reduce the ivar

    flux_reducedS = flux_reduced - 1.0

    #Pick out the emission lines
    line_flux,line_err = spec_array['line_flux'],spec_array['line_err']
    line_fluxS = line_flux - 1.0
    line_errS = line_err - 1.0

    line_arr = line_load()
    lower_shift = 3700 / (1+zobj)
    upper_shift = 9500 / (1+zobj)
    for i in range(len(line_arr)):
        line_arr['OBS_WAVE'][i] = line_arr['REST_WAVE'][i] * (zobj + 1)
    wlines = np.where((line_arr['REST_WAVE']>=lower_shift)&(line_arr['REST_WAVE']<=upper_shift))[0]
    vis_lines = line_arr[wlines]


    x_lower = np.amin(wave_arr)
    x_upper = np.amax(wave_arr)
    vis_range = np.where((wave_arr>=3700)&(wave_arr<=9000))[0]
    y_range = np.amax(flux_raw[vis_range]) - np.amin(flux_raw[vis_range])
    y_pad = float(y_range) / 10
    y_lower = np.amin(flux_raw[vis_range]) - y_pad
    y_upper = np.amax(flux_raw[vis_range]) + y_pad

    for i in range(len(wlines)):
        lam_center = vis_lines['REST_WAVE'][i]
        x0_lower = lam_center - 30.0
        x0_upper = lam_center + 30.0
        lam_lower = lam_center - 100.0
        lam_upper = lam_center + 100.0
        g_true,x0_true,A_true,m_true,b_true = 20.0,lam_center,1.5,0.0,0.0
        result = [g_true,x0_true,A_true,m_true,b_true]
        ndim, nwalkers = 5, 50
        p0 = [result + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

        wem = np.where((wave_shift>=lam_lower)&(wave_shift<=lam_upper))[0]
        x_em = wave_shift[wem]
        y_em = line_fluxS[wem]
        err_em = line_errS[wem]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x_em,y_em,err_em),a=3.5)

        #Burn-in
        pos,prob,state = sampler.run_mcmc(p0, 100) #----Ran well with 500
        sampler.reset()
        pos, prob, state = sampler.run_mcmc(pos, 1000) #----Ran well with 1000
        samples = sampler.flatchain
        g_mcmc, x0_mcmc, A_mcmc, m_mcmc,b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                 zip(*np.percentile(samples, [16, 50, 84],
                                                    axis=0)))
        acmean = np.mean(sampler.acceptance_fraction)

        #print('\n')
        #print('Central Wavelength: {0:6.1f}'.format(lam_center))
        #print('----------------PARAMETERS------------------')
        #print('Parameter |   50th   |   84th   |   16th   |')
        #print('--------------------------------------------')
        #print('    g     |  {0:5.2f}   |  {1:5.2f}   | {2:5.2f}    |'.format(g_mcmc[0],g_mcmc[1],g_mcmc[2]))
        #print('    x0    |  {0:4d}    |  {1:4d}    | {2:4d}     |'.format(int(x0_mcmc[0]),int(x0_mcmc[1]),int(x0_mcmc[2])))
        #print('    A     |  {0:5.2f}   |  {1:5.2f}   | {2:5.2f}    |'.format(A_mcmc[0],A_mcmc[1],A_mcmc[2]))
        #print('    m     |  {0:6.2f}  |  {1:6.2f}  | {2:6.2f}   |'.format(m_mcmc[0],m_mcmc[1],m_mcmc[2]))
        #print('    b     |  {0:.4f} | {1:0.6f} | {2:0.6f} |'.format(b_mcmc[0],b_mcmc[1],b_mcmc[2]))
        #print('--------------------------------------------')
        #print('Mean Acceptance Fraction: {0:6.4f}'.format(acmean))

        y_lor = line_func(x_em,g_mcmc[0],x0_mcmc[0],A_mcmc[0],m_mcmc[0],b_mcmc[0])
        max_amp = np.amax(y_lor)
        '''
        try:
            y_lor = line_func(x_em,g_mcmc[0],x0_mcmc[0],A_mcmc[0],m_mcmc[0],b_mcmc[0])
            max_amp = np.amax(y_lor)
        except Exception:
            continue
        '''
        #print(max_amp)
        if ((g_mcmc[0]>=5.0)&(g_mcmc[0]<50.0)&(acmean>0.22)&(acmean<0.5)&(max_amp>0.35)):
            match_flag = 1
            #print('\n')
            #print('TARGET IS QSO')
            #print('\n')
            #print('Max. Relative Flux: {0:06.4f}'.format(max_amp))

            matplotlib.rc('font',size=32)
            x_lower = np.amin(wave_shift)
            x_upper = np.amax(wave_shift)
            fig1,ax1 = plt.subplots(figsize=(10,8))
            ax1.plot(wave_shift,flux_reducedS,color='0.75',linewidth=1.3)
            ax1.plot(wave_shift,line_fluxS,linewidth=1.0,color='black')
            ax1.plot(x_em,y_lor,linewidth=1.5,color='orange')
            ax1.set_xlim((x_lower,x_upper))
            ax1.set_xlabel(r'Wavelength (\AA)')
            ax1.set_ylim((-2.0,5))
            ax1.set_ylabel(r'Rel. Flux')
            ax1.tick_params(axis='both',direction='in')
            ax1.tick_params(axis='both',which='minor',direction='in')
            ax1.xaxis.set_minor_locator(ticker.MultipleLocator(100))
            ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))

            matplotlib.rc('font',size=15)
            fig2 = corner.corner(samples[:,0:3], labels=["$\sigma_1$", "$\mu_1$", "$A_1$"],
                        label_kwargs={"fontsize": 15},quantiles=[0.16, 0.5, 0.84],
                        show_titles=True, levels=(1-np.exp(-0.5),), title_kwargs={"fontsize": 15})
            '''
            fig2 = corner.corner(samples[:,:], labels=["$\sigma_1$", "$\mu_1$", "$A_1$","$m$","$b$"],
                        label_kwargs={"fontsize": 15},quantiles=[0.16, 0.5, 0.84],
                        show_titles=True, title_kwargs={"fontsize": 15})
            '''
            fig1.tight_layout()
            spec_name = 'poster_spec_fit{}.png'.format(i)
            corner_name = 'poster_corner{}.png'.format(i)
            fig1.savefig(spec_name)
            fig2.savefig(corner_name)

            #break

ifile = sys.argv[1]
file_load(ifile)
