# -*- coding: utf-8 -*-
from astropy.io import fits
from glob import glob

filelist = glob('*.fit')

print filelist

for sci in filelist:
  hdu_sci = fits.open(sci)
  hdu_sci[0].header['REFSPEC1'] = ref1
#  hdu_sci.writeto(sci,clobber=True)