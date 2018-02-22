from astropy.io import fits
import numpy as np
import sys
import progressBar as pb
import os
import tmark
import glob

#This gets the spectra png, one file at a time. DR12 is public data, so wget does
#not need a user/password.
def get_png(c_in,p_in,m_in,f_in):

    loc = 'https://data.sdss.org/sas/dr12/boss/spectro/redux/images/v5_7_0/v5_7_0'
    loc2 = 'https://data.sdss.org/sas/dr12/boss/spectro/redux/images/v5_7_2/v5_7_2'
    png_file = '{0}/{1:04d}-{2:05d}/spec-image-{3:04d}-{4:05d}-{5:04d}.png'.format(loc,p_in,m_in,p_in,m_in,f_in)
    png_file2 = '{0}/{1:04d}-{2:05d}/spec-image-{3:04d}-{4:05d}-{5:04d}.png'.format(loc2,p_in,m_in,p_in,m_in,f_in)
    get_cmd = 'wget -qP ../ica_png/{} {}'.format(c_in,png_file)
    get_cmd2 = 'wget -qP ../ica_png/{} {}'.format(c_in,png_file2)
    os.system(get_cmd)
    os.system(get_cmd2)


#Main program for getting the spectra by reading the superset file.
def png_get(infile,num_ran):
    sd = fits.open(infile)[1].data

    for j in range(4):
        if j == 0:
            clp = 'QSO'
            wobj = np.where((sd['CLASS_PERSON']==3)&(sd['Z_CONF_PERSON']==3))[0]
            tstr = 'Retrieving QSO PNG files'
        elif j == 1:
            clp = 'BAL'
            wobj = np.where((sd['CLASS_PERSON']==30)&(sd['Z_CONF_PERSON']==3))[0]
            tstr = 'Retrieving BAL PNG files'
        elif j == 2:
            clp = 'GAL'
            wobj = np.where((sd['CLASS_PERSON']==4)&(sd['Z_CONF_PERSON']==3))[0]
            tstr = 'Retrieving GAL PNG files'
        elif j == 3:
            clp = 'STR'
            wobj = np.where((sd['CLASS_PERSON']==1)&(sd['Z_CONF_PERSON']==3))[0]
            tstr = 'Retrieving STR PNG files'

        ran_w = np.random.randint(0,len(wobj),num_ran)

        tmark.tm(tstr)
        for i in range(num_ran):
            pt = sd['PLATE'][wobj[ran_w[i]]]
            mt = sd['MJD'][wobj[ran_w[i]]]
            ft = sd['FIBERID'][wobj[ran_w[i]]]

            get_png(clp,pt,mt,ft)
            pb.pbar(i,num_ran)

    tmark.tm('Files Retrieved')

#These lines let me run from the command line.
set_file = sys.argv[1]
num_get = int(sys.argv[2])
png_get(set_file,num_get)
