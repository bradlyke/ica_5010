from astropy.io import fits
import numpy as np
import sys
import progressBar as pb
import os
import tmark

def spec_get(infile):
	sd = fits.open(infile)[1].data
	numsd = len(sd)

	tmark.tm('Retrieving SEQUELS Spectra')
	for i in range(numsd):
		pt = sd['PLATE'][i]
		mt = sd['MJD'][i]
		ft = sd['FIBERID'][i]

		loc = 'https://data.sdss.org/sas/dr12/boss/spectro/redux/v5_7_0/spectra'
		spec_file = '{0}/{1}/spec-{2:04d}-{3:05d}-{4:04d}.fits'.format(loc,pt,pt,mt,ft)
		get_cmd = 'wget -qP ../ica_spec/ {}'.format(spec_file)
		os.system(get_cmd)
		pb.pbar(i,numsd)

	print('Spectra retrieved')

set_file = sys.argv[1]
spec_get(set_file)
