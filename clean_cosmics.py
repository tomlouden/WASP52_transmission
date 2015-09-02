# -*- coding: utf-8 -*-
from glob import glob
import numpy as np
import os
from astropy.io import fits 

def main():
  flist = glob('*_cal.fit')

  flist.sort()

  for i in range(2,len(flist)-2.0):
    min2 = flist[i-2]
    min1 = flist[i-1]
    cen = flist[i]
    plus1 = flist[i+1]
    plus2 = flist[i+2]

    min2_hdu = fits.open(min2)
    min1_hdu = fits.open(min1)
    cen_hdu = fits.open(cen)
    plus1_hdu = fits.open(plus1)
    plus2_hdu = fits.open(plus2)

    newname = cen.strip('.fit') + '_dif.fit'

    diff1 = (cen_hdu[0].data - min2_hdu[0].data)
    diff2 = (cen_hdu[0].data - min1_hdu[0].data)
    diff3 = (cen_hdu[0].data - plus1_hdu[0].data)
    diff4 = (cen_hdu[0].data - plus2_hdu[0].data)

    stack = np.array([diff1] + [diff2]+ [diff3]+ [diff4])

    difference = np.median(stack,axis=0)

    hdu = fits.PrimaryHDU(difference)

    hdu.header = cen_hdu[0].header

    hdu.writeto(newname,clobber=True)

main()