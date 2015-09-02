# -*- coding: utf-8 -*-
from glob import glob
from astropy.io import fits 
import numpy as np

def main():
  data_dir = '../good_wavelength'
  imagelist = glob(data_dir+'/science/*.fit')
  biaslist = glob(data_dir+'/bias/*.fit')
  flatlist = glob(data_dir+'/sky_image_flats/*.fit')

#  make_master_bias(biaslist)
  master_bias = fits.open('master_bias.fits')[1].data
#  make_master_flat(flatlist,master_bias)
  master_flat = fits.open('master_flat.fits')[1].data
#  make_bad_mask(master_flat)
  bad_mask = fits.open('bad_mask.fits')[0].data

  callibrate_images(imagelist,master_bias,master_flat,bad_mask)

def make_bad_mask(master_flat):
  bad_mask = np.round(master_flat,0)

  hdu = fits.PrimaryHDU(bad_mask)
  hdu.writeto('bad_mask.fits', clobber=True)

def make_master_bias(filelist):

  bias_cube = []
  for filename in filelist:
    bias_cube += [fits.open(filename)[1].data]
  bias_cube = np.array(bias_cube)
  master_bias = np.median(bias_cube, axis=0)

  offset = float(fits.open(filelist[0])[1].header['BZERO'])
  master_bias = master_bias - offset

  hdulist = fits.open(filelist[0])
  hdulist[1].data = master_bias
  
  hdulist.writeto('master_bias.fits', clobber=True)

def make_master_flat(filelist,masterbias):

  flat_cube = []
  for filename in filelist:
    flat_cube += [fits.open(filename)[1].data - masterbias]

  flat_cube = np.array(flat_cube)
  master_flat = np.median(flat_cube, axis=0)

  norm_flat = master_flat / np.median(master_flat)

  offset = float(fits.open(filelist[0])[1].header['BZERO'])
  norm_flat = norm_flat - offset

  hdulist = fits.open(filelist[0])

  hdulist[1].data = norm_flat
  hdulist.writeto('master_flat.fits', clobber=True)

def callibrate_images(filelist,master_bias,master_flat,bad_mask):

  for filename in filelist:
    newname = filename.split('/')[-1].strip('.fit')+('_cal.fit')
    hdulist = fits.open(filename)
    raw = hdulist[1].data
    calibrated = ((raw - master_bias)/master_flat)*bad_mask
    hdulist[1].data = calibrated[:,50:-50]

    hdu = fits.PrimaryHDU(calibrated[:,50:-50])

    hdu.header = hdulist[0].header

    hdu.writeto(newname,clobber=True)


main()