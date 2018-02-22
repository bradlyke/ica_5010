"""
	You will need one of my other repositories: utilities. That includes
	the code for tmark and progressBar.

	NOTE: Running this file with tflag = 0 on the DR12Q superset will download
	~56 GB to your computer (~36,000 files at ~1.5 MB a piece).
"""

from astropy.io import fits
import numpy as np
import sys
import progressBar as pb
import os
import tmark
import glob

#This will make the directories if they aren't present.
def make_dirs(tplate):
	check_loc = '../ica_spec/*'
	folder_list = glob.glob(check_loc)
	check_plate = '../ica_spec/{}'.format(tplate)
	if check_plate not in folder_list:
		plate_dir = 'mkdir ../ica_spec/{}'.format(tplate)
		os.system(plate_dir)

#This gets the spectra, one file at a time. DR12 is public data, so wget does
#not need a user/password.
def get_files(p_in,m_in,f_in):

	loc = 'https://data.sdss.org/sas/dr12/boss/spectro/redux/v5_7_0/spectra'
	loc2 = 'https://data.sdss.org/sas/dr12/boss/spectro/redux/v5_7_2/spectra'
	spec_file = '{0}/{1}/spec-{2:04d}-{3:05d}-{4:04d}.fits'.format(loc,p_in,p_in,m_in,f_in)
	spec_file2 = '{0}/{1}/spec-{2:04d}-{3:05d}-{4:04d}.fits'.format(loc2,p_in,p_in,m_in,f_in)
	get_cmd = 'wget -qP ../ica_spec/{} {}'.format(p_in,spec_file)
	get_cmd2 = 'wget -qP ../ica_spec/{} {}'.format(p_in,spec_file2)
	os.system(get_cmd)
	os.system(get_cmd2)

#Main program for getting the spectra by reading the superset file. Tflag is a
#testing flag. If you want to get all the files, go with tflag = 0.
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

#These lines let me run from the command line.
set_file = sys.argv[1]
test_flag = int(sys.argv[2])
spec_get(set_file,test_flag)
