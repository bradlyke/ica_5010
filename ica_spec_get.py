from astropy.io import fits
import numpy as np
import sys
import progressBar as pb
import os
import tmark
import glob

def make_dirs(tplate):
	check_loc = '../ica_spec/*'
	folder_list = glob.glob(check_loc)
	check_plate = '../ica_spec/{}'.format(tplate)
	if check_plate not in folder_list:
		plate_dir = 'mkdir ../ica_spec/{}'.format(tplate)
		os.system(plate_dir)

def get_files(p_in,m_in,f_in):

	loc = 'https://data.sdss.org/sas/dr12/boss/spectro/redux/v5_7_0/spectra'
	spec_file = '{0}/{1}/spec-{2:04d}-{3:05d}-{4:04d}.fits'.format(loc,p_in,p_in,m_in,f_in)
	get_cmd = 'wget -qP ../ica_spec/{} {}'.format(p_in,spec_file)
	os.system(get_cmd)

def spec_get(infile,tflag):
	sd = fits.open(infile)[1].data

	if tflag == 0:
		numrec = len(sd)
	elif tflag == 1:
		numrec = 10

	tmark.tm('Retrieving SEQUELS Spectra')
	for i in range(numrec):
		pt = sd['PLATE'][i]
		mt = sd['MJD'][i]
		ft = sd['FIBERID'][i]

		make_dirs(pt)
		get_files(pt,mt,ft)
		pb.pbar(i,numrec)
	tmark.tm('Files Retrieved')


set_file = sys.argv[1]
tflag = int(sys.argv[2])
spec_get(set_file)
