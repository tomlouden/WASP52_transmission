# -*- coding: utf-8 -*-
from astropy.io import fits
from glob import glob
from numpy import *

filelist = glob('apps_r*')

print filelist

filelist = sort(filelist)

for sci in filelist:
#  sci = 'spec_'+arc.split('arc_')[1]
#  arc = 'arc3_'+sci.split('apps_')[1]
  arc = 'arc_'+sci.split('apps_')[1]
  print sci
  ref1 = 'arc_'+filelist[0].split('apps_')[1].strip('.fits') + ' 1'
  ref2 = 'arc_'+filelist[0].split('apps_')[1].strip('.fits') + ' 2'
  print ref2
#  hdu_arc = fits.open(arc)
  hdu_sci = fits.open(sci)

#  print hdu_arc[0].header['REFSPEC1']
#  print hdu_arc[0].header['REFSPEC2']
#  print hdu_sci[0].header['REFSPEC1']
#  print hdu_sci[0].header['REFSPEC2']

#  hdu_arc[0].header['REFSPEC1'] = ref1
#  hdu_arc[0].header['REFSPEC2'] = ref2
  hdu_sci[0].header['REFSPEC1'] = ref1
  hdu_sci[0].header['REFSPEC2'] = ref2
#  hdu_arc.writeto(arc,clobber=True)
  hdu_sci.writeto(sci,clobber=True)
