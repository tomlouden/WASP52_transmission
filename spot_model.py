# -*- coding: utf-8 -*-
from pylab import *
from ATLAS import *

def niceplot():

  spot_coverage = 0.1
  transit_depth = 0.164
  R = 16
  P = 'P00'
  v = 2
  t = 5000
  g = 4.5

  wvl = array([4250, 4350, 4450, 4550, 4650, 4750, 4850, 4950, 5050, 5150, 5250, 5350, 5450, 5550, 5650, 5750, 5850, 5950, 6050, 6150, 6250, 6350, 6450, 6550, 6650, 6750, 6850, 6950, 7050, 7150, 7250, 7350, 7450, 7550, 7650, 7750, 7850, 7950, 8050, 8150, 8250, 8350, 8450])

  for i in range(1,7):
    contrast = 250*i
#    wvl, depth = spotty_transit_model(spot_coverage,transit_depth,R,P,v,t,g,contrast)
    depth = spotty_transit(spot_coverage,transit_depth,R,P,v,t,g,contrast,wvl)
    spot_temp = t - contrast
    plot(wvl,depth/min(depth),label=str(spot_temp)+' K')

  ylabel('Transit depth modifaction')
  xlabel('wavelength (A)')
  legend()
  show()

def spotty_transit(spot_coverage,transit_depth,R,P,v,t,g,contrast,wvl_bins,bin_width=100):

  wvl, depth = spotty_transit_model(spot_coverage,transit_depth,R,P,v,t,g,contrast)

  bin_depth = []

  for i in range(0,len(wvl_bins)):
    min_bin = wvl_bins[i] - bin_width/2.0
    max_bin = wvl_bins[i] + bin_width/2.0
    index = [(wvl > min_bin) & (wvl < max_bin)]
    bin_depth += [mean(depth[index])]

  bin_depth = array(bin_depth)

  return bin_depth


def spotty_transit_model(spot_coverage,transit_depth,R,P,v,t,g,contrast):

  star_d = load_spec(P,v,t,g)
  spot_d = load_spec(P,v,t-contrast,g)

  wvl = star_d['wvl']
  star = star_d['flux']
  spot = spot_d['flux']

  index = [(wvl > 3000) & (wvl < 12000)]

  wvl = wvl[index]
  star = star[index]
  spot = spot[index]

  out_transit = star*(1.0- spot_coverage) + spot*(spot_coverage)

  old_out_transit = out_transit.copy()

  in_transit = star*(1.0- spot_coverage - transit_depth**2.0 ) + spot*(spot_coverage)

  old_in_transit = in_transit.copy()

  out_transit = convolve(wvl,out_transit,R)
  in_transit = convolve(wvl,in_transit,R)

  index = [(wvl > 4000) & (wvl < 8500)]

  wvl = wvl[index]
  out_transit = out_transit[index]
  old_out_transit = old_out_transit[index]
  in_transit = in_transit[index]
  old_in_transit = old_in_transit[index]

  out_transit = out_transit*(sum(old_out_transit)/sum(out_transit))
  in_transit = in_transit*(sum(old_in_transit)/sum(in_transit))

  depth = sqrt(1.0 - (in_transit/out_transit))

  old_depth = sqrt(1.0 - (old_in_transit/old_out_transit))

  return wvl, depth

def convolve(x,y,R):

  new_y = []
  for i in range(0,len(x)):
    sigma = x[i]/R
    gauss = gaussian(x,x[i],sigma)
    new = gauss*y
    new_y += [sum(new)]
  new_y = array(new_y)
  new_y = new_y*(sum(y)/sum(new_y))
  return new_y

def gaussian(x,x0,sigma):
  return (1.0/(sigma*sqrt(2.0*pi)))*exp(-((x-x0)**2/(2.0*sigma**2.0)))