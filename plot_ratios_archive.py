# -*- coding: utf-8 -*-
from pylab import *
from glob import glob
from astropy.io import fits 
from astropy.wcs import wcs 
from scipy.optimize import leastsq
from scipy.optimize import minimize
from PyAstronomy.pyasl.asl.astroTimeLegacy import *
import transit
from numpy import *

def main():


  filelist = glob('t_app1*')

  hjd, airmass, app1_table, app2_table = load_data(filelist,cache='all_data.fits')

  index = argsort(hjd)

  airmass = airmass[index]
  hjd = hjd[index]
  app1_table = app1_table[index]
  app2_table = app2_table[index]

  divider = 2456895
  #flux_table = flux_table[hjd < divider]
  #wvl_table = wvl_table[hjd < divider]
  #airmass = airmass[hjd < divider]
  #flux_frac_errors = flux_frac_errors[hjd < divider]
  #app1_table = app1_table[hjd < divider]
  #app2_table = app2_table[hjd < divider]
  #hjd = hjd[hjd < divider]

  reg_min = 5800
  reg_max = 6200

  #correct_app1, center_1 = align_lightcurves(app1_table[:,1],app1_table[:,2],app1_table[:,0],reg_min,reg_max)
  #correct_app2, center_2 = align_lightcurves(app2_table[:,1],app2_table[:,2],app2_table[:,0],reg_min,reg_max)

  #app1_table[:,0] = correct_app1
  #app2_table[:,0] = correct_app2 + (center_1[0] - center_2[0])

  #for i in range(1,len(app1_table)):
    #dif = cross_correlate(app1_table[0],app1_table[i])
    #app1_table[i][0] += dif

  #for i in range(1,len(app2_table)):
    #dif = cross_correlate(app2_table[0],app2_table[i])
    #app2_table[i][0] += dif

  #for i in range(0,len(app1_table[:,0])):
    #plot(app1_table[i,0],app1_table[i,1]/median(app1_table[i,1][(app1_table[i,0] > 6000) & (app1_table[i,0] < 7000)]))
    #plot(app2_table[i,0],app2_table[i,1]/median(app2_table[i,1][(app2_table[i,0] > 6000) & (app2_table[i,0] < 7000)]))
  #show()
  #quit()

  reg_min = 7900
  reg_max = 8000

  regs_min = []
  regs_max = []
  regs_cen = []
  transits = []
  lightcurves = []
  frac_errors = []
  results = []
  depth = []
  errors = []
  depth_error = []


  width = 100
  zero = 4000

  for i in range(1,45):
#  for i in range(30,32):
    regs_min += [zero + i*width]
    regs_max += [zero + width + i*width]
    regs_cen += [zero + (width/2.0)+ i*width]
    transit, lightcurve, frac_error, result, result_errors = fit_a_lightcurve(app1_table,app2_table,hjd,airmass,regs_min[-1],regs_max[-1])
    transits += [transit]
    results += [result]
    lightcurves += [lightcurve]
    frac_errors += [frac_error]
    depth += [result[0]]
    errors += [result_errors]
    depth_error +=[result_errors[0]]
    plot_lightcurve(lightcurves[-1],frac_errors[-1],results[-1],airmass,hjd,name=str(int(regs_cen[-1]))+'.pdf')

    #plot_lightcurve(lightcurves[-1],frac_errors[-1],results[-1],airmass,hjd)
    #show()
  
  print depth
  print depth_error

  errorbar(regs_cen,depth,yerr=depth_error,xerr=50,fmt='ko')
  ylabel('Transit depth')
  xlabel('Wavelength (A)')

  ylim([0.15,0.18])

#  show()
  savefig('spectrum.pdf')

def plot_lightcurve(lightcurve,frac_error,x,airmass,hjd,name=False):
  # first remove the polynomial from both transit and lightcurve.

  index = argsort(hjd)

  lightcurve = lightcurve[index]
  airmass = airmass[index]
  hjd = hjd[index]

  poly_coefs = x.copy()
  poly_coefs[0] = 0

  poly_m, z = transit_model(poly_coefs,hjd,airmass)

  transit_m, z = transit_model(x,hjd,airmass)

  norm_m = transit_m/poly_m

  norm_lightcurve = lightcurve / poly_m

  f, axarr = plt.subplots(2, sharex=True)

  residuals = norm_lightcurve - norm_m

  zero = norm_m - norm_m

  dummy = arange(0,len(norm_lightcurve))

  first_min = min(hjd)
  second_min = min(hjd[abs(hjd - 2456892.54327) > 3.0])

  t = hjd.copy()
  for i in range(0,len(t)):
    if (abs(t[i] - 2456892.54327) > 3.0):
      t[i] = ((t[i] - second_min)*(24*60)) + 500
    else:
      t[i] = (t[i] - first_min)*(24*60)

  dummy = t

  axarr[0].errorbar(dummy,norm_lightcurve,yerr=frac_error*norm_lightcurve,fmt='r,')
  axarr[0].plot(dummy,norm_m,'k-')

  axarr[0].set_ylabel('Normalised flux ratio')

  axarr[1].errorbar(dummy,residuals,yerr=frac_error*norm_lightcurve,fmt='r,')
  axarr[1].plot(dummy,zero,'k-')

  axarr[1].set_ylabel('Residuals')

  axarr[1].set_xlabel('Frame number')

  if name:
    savefig(name)
    close()
  else:
    show()

def fit_a_lightcurve(app1_table,app2_table,hjd,airmass,reg_min,reg_max):

  divider = 2456895

  b_range = str(reg_min) + ' - ' + str(reg_max)

  app1_b_lightcurve, app1_b_errors = create_lightcurve(app1_table[:,1],app1_table[:,2],app1_table[:,0],reg_min,reg_max,hjd)
  app2_b_lightcurve, app2_b_errors = create_lightcurve(app2_table[:,1],app2_table[:,2],app2_table[:,0],reg_min,reg_max,hjd)

  b_lightcurve = app2_b_lightcurve / app1_b_lightcurve
  b_errors = sqrt((app1_b_errors/app1_b_lightcurve)**2.0 + (app2_b_errors/app2_b_lightcurve)**2.0)


  b_frac_errors = b_errors

  #b_lightcurve = b_lightcurve / k_lightcurve
  #b_frac_errors = ((b_errors)**2 + (k_errors)**2)**0.5

  poly_order = 2

  filt_time = -1

  start_transit_hjd = 2456892.50423
  mid_transit_hjd = 2456892.54327
  end_transit_hjd = 2456892.58231

  start_transit_hjd_2 = 2456899.50335
  mid_transit_hjd_2 = 2456899.54239
  end_transit_hjd_2 = 2456899.58143

  start_transit = (start_transit_hjd - min(hjd))*24*60
  mid_transit = (mid_transit_hjd - min(hjd))*24*60
  end_transit = (end_transit_hjd - min(hjd))*24*60

  print start_transit, mid_transit, end_transit

  first_min = min(hjd)
  second_min = min(hjd[abs(hjd - 2456892.54327) > 3.0])

  start_transit_2 = (start_transit_hjd_2 - second_min)*24*60 + 500
  mid_transit_2 = (mid_transit_hjd_2 - second_min)*24*60 +500 
  end_transit_2 = (end_transit_hjd_2 - second_min)*24*60 +500

  t = hjd.copy()
  for i in range(0,len(t)):
    if (abs(t[i] - 2456892.54327) > 3.0):
      t[i] = ((t[i] - second_min)*(24*60)) + 500
    else:
      t[i] = (t[i] - first_min)*(24*60)

  b_lightcurve_raw = b_lightcurve.copy()
  b_frac_errors_raw = b_frac_errors.copy()

  hjd_raw = hjd.copy()

  index = argsort(t)
  b_lightcurve = b_lightcurve[index]
  b_frac_errors = b_frac_errors[index]
  airmass = airmass[index]
  hjd = hjd[index]
  t = t[index]

  filt_time = logical_or(((hjd<divider) & (t<430)),(hjd>divider))

  b_lightcurve = b_lightcurve[filt_time]
  b_frac_errors = b_frac_errors[filt_time]
  airmass = airmass[filt_time]
  hjd = hjd[filt_time]
  t = t[filt_time]

  start_spot = 146
  end_spot = 187

#  start_spot = 250
#  end_spot = 95

  not_spot = logical_or((t<start_spot),(t>end_spot))
  b_lightcurve = b_lightcurve[not_spot]
  b_frac_errors = b_frac_errors[not_spot]
  airmass = airmass[not_spot]
  hjd = hjd[not_spot]
  t = t[not_spot]

  start_spot = 676
  end_spot = 710
  not_spot = logical_or((t<start_spot),(t>end_spot))
  b_lightcurve = b_lightcurve[not_spot]
  b_frac_errors = b_frac_errors[not_spot]
  airmass = airmass[not_spot]
  hjd = hjd[not_spot]
  t = t[not_spot]

  start_spot = 0
  end_spot = 30
  not_spot = logical_or((t<start_spot),(t>end_spot))
  b_lightcurve = b_lightcurve[not_spot]
  b_frac_errors = b_frac_errors[not_spot]
  airmass = airmass[not_spot]
  hjd = hjd[not_spot]
  t = t[not_spot]

  start_spot = 490
  end_spot = 530
  not_spot = logical_or((t<start_spot),(t>end_spot))
  b_lightcurve = b_lightcurve[not_spot]
  b_frac_errors = b_frac_errors[not_spot]
  airmass = airmass[not_spot]
  hjd = hjd[not_spot]
  t = t[not_spot]

  not_transit_1 = logical_or((t<start_transit),(t>end_transit))
  not_transit_2 = logical_or((t<start_transit_2),(t>end_transit_2))

  not_transit = logical_and(not_transit_1,not_transit_2)

  transit_window = [start_transit,end_transit]

  b_lightcurve_old = b_lightcurve.copy()[not_transit]

  hjd_c = hjd.copy()[not_transit]

  for i in range(0,len(hjd)):
    if hjd[i] < divider:
      b_lightcurve[i] = b_lightcurve[i]/median(b_lightcurve_old[hjd_c<divider])
    else:
      b_lightcurve[i] = b_lightcurve[i]/median(b_lightcurve_old[hjd_c>divider])

  for i in range(0,len(hjd_raw)):
    if hjd_raw[i] < divider:
      b_lightcurve_raw[i] = b_lightcurve_raw[i]/median(b_lightcurve_old[hjd_c<divider])
    else:
      b_lightcurve_raw[i] = b_lightcurve_raw[i]/median(b_lightcurve_old[hjd_c>divider])

  p0 = 0.164
  tt = mid_transit_hjd - 2.4568925e6
  per = 1.7497798
  inc = 85.35
  ars = 7.38
  gamma_0 = 0.5
  gamma_1 = 0.5

  x = [p0,tt,per,inc,ars,gamma_0,gamma_1,-0.3,1.00,1.00,1.0,0.002,1.0,0.002]

  poly_prior = [-10.0,10.0]

  x += poly_prior*2

  x = [1.62858863e-01,-1.25135820e-04,1.74977980e+00,8.53500000e+01,7.38000000e+00,3.17007330e-01,3.17007330e-01,-1.97555202e-01,9.72396652e-01,9.60859487e-01,1.03377789e+00,-8.28656790e-03,1.04746041e+00,-1.96077433e-02,-5.61442377e+00,-5.27598514e+00,-6.11726671e+00,-4.77924359e+00]

  b_transit, b_final, b_final_errors = fit_transit(hjd,airmass,b_lightcurve,b_frac_errors,x,log=True)

  b_transit, b_final, b_final_errors = fit_transit(hjd,airmass,b_lightcurve,b_frac_errors,b_final,log=False,fix=b_final)

  print 'finale!'

  print b_final

  return b_transit, b_lightcurve_raw, b_frac_errors_raw, b_final, b_final_errors

def some_plotting_stuff():

  binup = 1

  #k_lightcurve = binning(k_lightcurve,binup)
  #b_lightcurve = binning(b_lightcurve,binup)
  #g_lightcurve = binning(g_lightcurve,binup)
  #r_lightcurve = binning(r_lightcurve,binup)
  #t_bin = binning(t,binup)

  t_bin = t

  #k_frac_errors = (binning(k_frac_errors,binup))/sqrt(binup)
  #b_frac_errors = (binning(b_frac_errors,binup))/sqrt(binup)
  #g_frac_errors = (binning(g_frac_errors,binup))/sqrt(binup)
  #r_frac_errors = (binning(r_frac_errors,binup))/sqrt(binup)

  divider = 2456895

  guide = array([0]*len(t))

  guide_time = array([-5,max(t)+5])

  f, axarr = plt.subplots(2, sharex=True)

  axarr[0].plot([start_transit,start_transit],[-1.8,1.8],'k--')
  axarr[0].plot([mid_transit,mid_transit],[-1.8,1.8],'k--')
  axarr[0].plot([end_transit,end_transit],[-1.8,1.8],'k--')

  axarr[0].plot([start_transit_2,start_transit_2],[-1.8,1.8],'k--')
  axarr[0].plot([mid_transit_2,mid_transit_2],[-1.8,1.8],'k--')
  axarr[0].plot([end_transit_2,end_transit_2],[-1.8,1.8],'k--')

  axarr[0].errorbar(t_bin,(k_lightcurve+0.05)-1.0,yerr=k_frac_errors,fmt='ko',label=k_range)
  axarr[0].errorbar(t_bin,b_lightcurve-1.0,yerr=b_frac_errors,fmt='bo',label=b_range)
  axarr[0].errorbar(t_bin,(g_lightcurve-0.05)-1.0,yerr=g_frac_errors,fmt='go',label=g_range)
  axarr[0].errorbar(t_bin,(r_lightcurve-0.1)-1.0,yerr=r_frac_errors,fmt='ro',label=r_range)

  axarr[0].plot(t[hjd < divider],(k_transit[hjd < divider]+0.05)-1.0,'k-')
  axarr[0].plot(t[hjd < divider],b_transit[hjd < divider]-1.0,'k-')
  axarr[0].plot(t[hjd < divider],(g_transit[hjd < divider]-0.05)-1.0,'k-')
  axarr[0].plot(t[hjd < divider],(r_transit[hjd < divider]-0.1)-1.0,'k-')

  axarr[0].plot(t[hjd > divider],(k_transit[hjd > divider]+0.05)-1.0,'k-')
  axarr[0].plot(t[hjd > divider],b_transit[hjd > divider]-1.0,'k-')
  axarr[0].plot(t[hjd > divider],(g_transit[hjd > divider]-0.05)-1.0,'k-')
  axarr[0].plot(t[hjd > divider],(r_transit[hjd > divider]-0.1)-1.0,'k-')

  axarr[1].plot([start_transit,start_transit],[-1.8,1.8],'k--')
  axarr[1].plot([mid_transit,mid_transit],[-1.8,1.8],'k--')
  axarr[1].plot([end_transit,end_transit],[-1.8,1.8],'k--')

  axarr[1].plot([start_transit_2,start_transit_2],[-1.8,1.8],'k--')
  axarr[1].plot([mid_transit_2,mid_transit_2],[-1.8,1.8],'k--')
  axarr[1].plot([end_transit_2,end_transit_2],[-1.8,1.8],'k--')

  axarr[1].errorbar(t_bin,(k_lightcurve-k_transit)+0.05,yerr=k_frac_errors,fmt='ko',label=k_range)
  axarr[1].errorbar(t_bin,(b_lightcurve-b_transit),yerr=b_frac_errors,fmt='bo',label=b_range)
  axarr[1].errorbar(t_bin,(g_lightcurve-g_transit)-0.05,yerr=g_frac_errors,fmt='go',label=g_range)
  axarr[1].errorbar(t_bin,(r_lightcurve-r_transit)-0.1,yerr=r_frac_errors,fmt='ro',label=r_range)

  axarr[1].plot(t[hjd < divider],guide[hjd < divider]+0.05,'k-')
  axarr[1].plot(t[hjd < divider],guide[hjd < divider],'k-')
  axarr[1].plot(t[hjd < divider],guide[hjd < divider]-0.05,'k-')
  axarr[1].plot(t[hjd < divider],guide[hjd < divider]-0.1,'k-')

  axarr[1].plot(t[hjd > divider],guide[hjd > divider]+0.05,'k-')
  axarr[1].plot(t[hjd > divider],guide[hjd > divider],'k-')
  axarr[1].plot(t[hjd > divider],guide[hjd > divider]-0.05,'k-')
  axarr[1].plot(t[hjd > divider],guide[hjd > divider]-0.1,'k-')

  axarr[0].legend(loc='upper right')

  b_chi2 = sum(((b_lightcurve - 1)/b_frac_errors)**2)/(len(b_lightcurve) - 1 - poly_order)
  g_chi2 = sum(((g_lightcurve - 1)/g_frac_errors)**2)/(len(g_lightcurve) - 1 - poly_order)
  r_chi2 = sum(((r_lightcurve - 1)/r_frac_errors)**2)/(len(r_lightcurve) - 1 - poly_order)

  print "Chi2 for hypothesis of flat (b,g,r):"

  print b_chi2, g_chi2, r_chi2

  b_rms = std(b_lightcurve - b_transit)
  g_rms = std(g_lightcurve - g_transit)
  r_rms = std(r_lightcurve - r_transit)

  print "fractional rms (b,g,r):"

  print b_rms, g_rms, r_rms

  cadence = max(t)/len(t)

  transit_dur = 60*3.0

  sn_increase = sqrt(transit_dur)

  print "over the timescale of a transit, that's(b,g,r):"

  print b_rms/sn_increase, g_rms/sn_increase, r_rms/sn_increase

  print "giving OTOO this significance per scale height (assume scale height ~ 1e-4) (b,g,r)"

  scale_height = 1e-4

  print scale_height/(b_rms/sn_increase), scale_height/(g_rms/sn_increase), scale_height/(r_rms/sn_increase)

  axarr[0].set_ylabel('Normalized flux - 1')
  axarr[1].set_ylabel('Residuals')
  xlabel('Time from start (minutes)')

  axarr[0].set_xlim([-5,max(t)+5])
  axarr[0].set_ylim([-0.16,0.16])
  axarr[1].set_ylim([-0.14,0.09])

  show()

def old_main():

  filelist = glob('txt_*')

  hjd, airmass, app1_table, app2_table = load_data(filelist,cache='all_data.fits')

  divider = 2456895
  #flux_table = flux_table[hjd < divider]
  #wvl_table = wvl_table[hjd < divider]
  #airmass = airmass[hjd < divider]
  #flux_frac_errors = flux_frac_errors[hjd < divider]
  #app1_table = app1_table[hjd < divider]
  #app2_table = app2_table[hjd < divider]
  #hjd = hjd[hjd < divider]

  reg_min = 5800
  reg_max = 6200

  #correct_app1, center_1 = align_lightcurves(app1_table[:,1],app1_table[:,2],app1_table[:,0],reg_min,reg_max)
  #correct_app2, center_2 = align_lightcurves(app2_table[:,1],app2_table[:,2],app2_table[:,0],reg_min,reg_max)

  #app1_table[:,0] = correct_app1
  #app2_table[:,0] = correct_app2 - (center_1[0] - center_2[0])

  #reg_min = 4850 
  #reg_max = 5150

  #correct_app1_l2, center_1_l2 = align_lightcurves(app1_table[:,1].copy(),app1_table[:,2].copy(),app1_table[:,0].copy(),reg_min,reg_max)

  #reg_min = 5100 
  #reg_max = 5500

  #correct_app2_l2, center_2_l2 = align_lightcurves(app2_table[:,1].copy(),app2_table[:,2].copy(),app2_table[:,0].copy(),reg_min,reg_max)

  #delta = (center_1[0] - center_1_l2[0])
  #line2_errs = center_1_l2[0] - center_1_l2
  #gradients = 1.0 - ((delta + line2_errs)/delta)
  #for i in range(0, len(app1_table[:,0])):
    #dist = (center_1[0] - app1_table[i,0])
    #correction = dist*gradients[i]
    #print correction
    #app1_table[i,0] -= correction


  #delta = (center_2[0] - center_2_l2[0])
  #line2_errs = center_2_l2[0] - center_2_l2
  #gradients = 1.0 - ((delta + line2_errs)/delta)
  #for i in range(0, len(app2_table[:,0])):
    #dist = (center_2[0] - app2_table[i,0])
    #correction = dist*gradients[i]
    #print correction
    #app2_table[i,0] -= correction


  reg_min = 4000
  reg_max = 4900

  reg_min = 4850 
  reg_max = 5150

  reg_min = 4000 
  reg_max = 4100

  b_range = str(reg_min) + ' - ' + str(reg_max)

  app1_b_lightcurve, app1_b_errors = create_lightcurve(app1_table[:,1],app1_table[:,2],app1_table[:,0],reg_min,reg_max)
  app2_b_lightcurve, app2_b_errors = create_lightcurve(app2_table[:,1],app2_table[:,2],app2_table[:,0],reg_min,reg_max)

  b_lightcurve = app2_b_lightcurve / app1_b_lightcurve
  b_errors = sqrt((app1_b_errors/app1_b_lightcurve)**2.0 + (app2_b_errors/app2_b_lightcurve)**2.0)

  reg_min = 5450
  reg_max = 5950

  reg_min = 4100
  reg_max = 4200

  g_range = str(reg_min) + ' - ' + str(reg_max) 

  app1_g_lightcurve, app1_g_errors = create_lightcurve(app1_table[:,1],app1_table[:,2],app1_table[:,0],reg_min,reg_max)
  app2_g_lightcurve, app2_g_errors = create_lightcurve(app2_table[:,1],app2_table[:,2],app2_table[:,0],reg_min,reg_max)
  g_lightcurve = app2_g_lightcurve / app1_g_lightcurve
  g_errors = sqrt((app1_g_errors/app1_g_lightcurve)**2.0 + (app2_g_errors/app2_g_lightcurve)**2.0)

  reg_min = 7800
  reg_max = 8200

  reg_min = 4200
  reg_max = 4300

  r_range = str(reg_min) + ' - ' + str(reg_max)   

  app1_r_lightcurve, app1_r_errors = create_lightcurve(app1_table[:,1],app1_table[:,2],app1_table[:,0],reg_min,reg_max)
  app2_r_lightcurve, app2_r_errors = create_lightcurve(app2_table[:,1],app2_table[:,2],app2_table[:,0],reg_min,reg_max)
  r_lightcurve = app2_r_lightcurve / app1_r_lightcurve
  r_errors = sqrt((app1_r_errors/app1_r_lightcurve)**2.0 + (app2_r_errors/app2_r_lightcurve)**2.0)

  reg_min = 4000
  reg_max = 8200

  k_range = str(reg_min) + ' - ' + str(reg_max)   

  app1_k_lightcurve, app1_k_errors = create_lightcurve(app1_table[:,1],app1_table[:,2],app1_table[:,0],reg_min,reg_max)
  app2_k_lightcurve, app2_k_errors = create_lightcurve(app2_table[:,1],app2_table[:,2],app2_table[:,0],reg_min,reg_max)

  k_lightcurve = app2_k_lightcurve / app1_k_lightcurve
  k_errors = sqrt((app1_k_errors/app1_k_lightcurve)**2.0 + (app2_k_errors/app2_k_lightcurve)**2.0)

  b_frac_errors = b_errors
  g_frac_errors = g_errors
  r_frac_errors = r_errors
  k_frac_errors = k_errors

  #b_lightcurve = b_lightcurve / k_lightcurve
  #g_lightcurve = g_lightcurve / k_lightcurve
  #r_lightcurve = r_lightcurve / k_lightcurve
  #b_frac_errors = ((b_errors)**2 + (k_errors)**2)**0.5
  #g_frac_errors = ((g_errors)**2 + (k_errors)**2)**0.5
  #r_frac_errors = ((r_errors)**2 + (k_errors)**2)**0.5

  poly_order = 2

  filt_time = -1

  start_transit_hjd = 2456892.50423
  mid_transit_hjd = 2456892.54327
  end_transit_hjd = 2456892.58231

  start_transit_hjd_2 = 2456899.50335
  mid_transit_hjd_2 = 2456899.54239
  end_transit_hjd_2 = 2456899.58143

  start_transit = (start_transit_hjd - min(hjd))*24*60
  mid_transit = (mid_transit_hjd - min(hjd))*24*60
  end_transit = (end_transit_hjd - min(hjd))*24*60

  print start_transit, mid_transit, end_transit

  t = hjd.copy()

  first_min = min(hjd)
  second_min = min(hjd[abs(hjd - 2456892.54327) > 3.0])

  start_transit_2 = (start_transit_hjd_2 - second_min)*24*60 + 500
  mid_transit_2 = (mid_transit_hjd_2 - second_min)*24*60 +500 
  end_transit_2 = (end_transit_hjd_2 - second_min)*24*60 +500

  for i in range(0,len(t)):
    if (abs(t[i] - 2456892.54327) > 3.0):
      t[i] = ((t[i] - second_min)*(24*60)) + 500
    else:
      t[i] = (t[i] - first_min)*(24*60)

  b_lightcurve = b_lightcurve[t > filt_time]
  g_lightcurve = g_lightcurve[t > filt_time]
  r_lightcurve = r_lightcurve[t > filt_time]
  k_lightcurve = k_lightcurve[t > filt_time]
  b_frac_errors = b_frac_errors[t > filt_time]
  g_frac_errors = g_frac_errors[t > filt_time]
  r_frac_errors = r_frac_errors[t > filt_time]
  k_frac_errors = k_frac_errors[t > filt_time]

  airmass = airmass[t > filt_time]
  hjd = hjd[t > filt_time]
  t = t[t > filt_time]

  index = argsort(t)

  b_lightcurve = b_lightcurve[index]
  g_lightcurve = g_lightcurve[index]
  r_lightcurve = r_lightcurve[index]
  k_lightcurve = k_lightcurve[index]
  b_frac_errors = b_frac_errors[index]
  g_frac_errors = g_frac_errors[index]
  r_frac_errors = r_frac_errors[index]
  k_frac_errors = k_frac_errors[index]

  airmass = airmass[index]
  hjd = hjd[index]
  t = t[index]

  start_spot = 95
  end_spot = 250

  start_spot = 146
  end_spot = 187

#  start_spot = 250
#  end_spot = 95

  not_spot = logical_or((t<start_spot),(t>end_spot))

  b_lightcurve = b_lightcurve[not_spot]
  g_lightcurve = g_lightcurve[not_spot]
  r_lightcurve = r_lightcurve[not_spot]
  k_lightcurve = k_lightcurve[not_spot]
  b_frac_errors = b_frac_errors[not_spot]
  g_frac_errors = g_frac_errors[not_spot]
  r_frac_errors = r_frac_errors[not_spot]
  k_frac_errors = k_frac_errors[not_spot]

  airmass = airmass[not_spot]
  hjd = hjd[not_spot]
  t = t[not_spot]

  not_transit_1 = logical_or((t<start_transit),(t>end_transit))
  not_transit_2 = logical_or((t<start_transit_2),(t>end_transit_2))

  not_transit = logical_and(not_transit_1,not_transit_2)

  transit_window = [start_transit,end_transit]

  #k_lightcurve = detrend_poly(k_lightcurve,airmass,k_frac_errors,poly_order,transit_window,t,istime=True)
  #quit()

  plot_white = False

  if plot_white == True:
    plot(t,(app1_k_lightcurve/median(app1_k_lightcurve))-1.0,'ko')
    plot(t,(app2_k_lightcurve/median(app2_k_lightcurve))-1.0,'ro')
    #plot([start_transit,start_transit],[-0.1,0.1],'k--')
    #plot([mid_transit,mid_transit],[-0.1,0.1],'k--')
    #plot([end_transit,end_transit],[-0.1,0.1],'k--')
    show()
    quit()

  #b_lightcurve = detrend_poly(b_lightcurve,k_lightcurve,poly_order,transit_window,t)
  #g_lightcurve = detrend_poly(g_lightcurve,k_lightcurve,poly_order,transit_window,t)
  #r_lightcurve = detrend_poly(r_lightcurve,k_lightcurve,poly_order,transit_window,t)

  
  #b_lightcurve = detrend_poly(b_lightcurve,airmass,poly_order,transit_window,t)
  #g_lightcurve = detrend_poly(g_lightcurve,airmass,poly_order,transit_window,t)
  #r_lightcurve = detrend_poly(r_lightcurve,airmass,poly_order,transit_window,t)
  #k_lightcurve = detrend_poly(k_lightcurve,airmass,poly_order,transit_window,t)

  #b_lightcurve = b_lightcurve/median(b_lightcurve[not_transit])
  #g_lightcurve = g_lightcurve/median(g_lightcurve[not_transit])
  #r_lightcurve = r_lightcurve/median(r_lightcurve[not_transit])
  #k_lightcurve = k_lightcurve/median(k_lightcurve[not_transit])

  b_lightcurve_old = b_lightcurve.copy()[not_transit]
  g_lightcurve_old = g_lightcurve.copy()[not_transit]
  r_lightcurve_old = r_lightcurve.copy()[not_transit]
  k_lightcurve_old = k_lightcurve.copy()[not_transit]

  hjd_c = hjd.copy()[not_transit]

  for i in range(0,len(hjd)):
    if hjd[i] < divider:
      b_lightcurve[i] = b_lightcurve[i]/median(b_lightcurve_old[hjd_c<divider])
      g_lightcurve[i] = g_lightcurve[i]/median(g_lightcurve_old[hjd_c<divider])
      r_lightcurve[i] = r_lightcurve[i]/median(r_lightcurve_old[hjd_c<divider])
      k_lightcurve[i] = k_lightcurve[i]/median(k_lightcurve_old[hjd_c<divider])
    else:
      b_lightcurve[i] = b_lightcurve[i]/median(b_lightcurve_old[hjd_c>divider])
      g_lightcurve[i] = g_lightcurve[i]/median(g_lightcurve_old[hjd_c>divider])
      r_lightcurve[i] = r_lightcurve[i]/median(r_lightcurve_old[hjd_c>divider])
      k_lightcurve[i] = k_lightcurve[i]/median(k_lightcurve_old[hjd_c>divider])

  p0 = 0.164
  tt = mid_transit_hjd
  per = 1.7497798
  inc = 85.35
  ars = 7.38
  gamma_0 = -0.2
  gamma_1 = 0.7

  #prior = [p0,tt,per,inc,ars,gamma_0,gamma_1]
  #k_transit = fit_transit(hjd,2.0-k_lightcurve,k_frac_errors,prior)
  #r_transit = fit_transit(hjd,2.0-r_lightcurve,r_frac_errors,prior)
  #g_transit = fit_transit(hjd,2.0-g_lightcurve,g_frac_errors,prior)
  #b_transit = fit_transit(hjd,2.0-b_lightcurve,b_frac_errors,prior)


  x = [p0,tt,per,inc,ars,gamma_0,gamma_1,-0.3,0.998,0.1]

  poly_prior = [-10.0,10.0]

  x += poly_prior*2

  k_transit, k_final = fit_transit(hjd,airmass,k_lightcurve,k_frac_errors,x)
  r_transit, r_final = fit_transit(hjd,airmass,r_lightcurve,r_frac_errors,x)
  g_transit, g_final = fit_transit(hjd,airmass,g_lightcurve,g_frac_errors,x)
  b_transit, b_final = fit_transit(hjd,airmass,b_lightcurve,b_frac_errors,x)

  print "done fitting"

  print k_final
  print r_final
  print g_final
  print b_final

  print 'white depth =',k_final[0]
  print 'blue depth =',b_final[0]
  print 'green depth =',g_final[0]
  print 'red depth =',r_final[0]

#  k_transit, r_transit, g_transit, b_transit = fit_transit_colours(hjd,airmass,k_lightcurve,k_frac_errors,r_lightcurve,r_frac_errors,g_lightcurve,g_frac_errors,b_lightcurve,b_frac_errors,prior)

  binup = 1

  #k_lightcurve = binning(k_lightcurve,binup)
  #b_lightcurve = binning(b_lightcurve,binup)
  #g_lightcurve = binning(g_lightcurve,binup)
  #r_lightcurve = binning(r_lightcurve,binup)
  #t_bin = binning(t,binup)

  t_bin = t

  #k_frac_errors = (binning(k_frac_errors,binup))/sqrt(binup)
  #b_frac_errors = (binning(b_frac_errors,binup))/sqrt(binup)
  #g_frac_errors = (binning(g_frac_errors,binup))/sqrt(binup)
  #r_frac_errors = (binning(r_frac_errors,binup))/sqrt(binup)

  divider = 2456895

  guide = array([0]*len(t))

  guide_time = array([-5,max(t)+5])

  f, axarr = plt.subplots(2, sharex=True)

  axarr[0].plot([start_transit,start_transit],[-1.8,1.8],'k--')
  axarr[0].plot([mid_transit,mid_transit],[-1.8,1.8],'k--')
  axarr[0].plot([end_transit,end_transit],[-1.8,1.8],'k--')

  axarr[0].plot([start_transit_2,start_transit_2],[-1.8,1.8],'k--')
  axarr[0].plot([mid_transit_2,mid_transit_2],[-1.8,1.8],'k--')
  axarr[0].plot([end_transit_2,end_transit_2],[-1.8,1.8],'k--')

  axarr[0].errorbar(t_bin,(k_lightcurve+0.05)-1.0,yerr=k_frac_errors,fmt='ko',label=k_range)
  axarr[0].errorbar(t_bin,b_lightcurve-1.0,yerr=b_frac_errors,fmt='bo',label=b_range)
  axarr[0].errorbar(t_bin,(g_lightcurve-0.05)-1.0,yerr=g_frac_errors,fmt='go',label=g_range)
  axarr[0].errorbar(t_bin,(r_lightcurve-0.1)-1.0,yerr=r_frac_errors,fmt='ro',label=r_range)

  axarr[0].plot(t[hjd < divider],(k_transit[hjd < divider]+0.05)-1.0,'k-')
  axarr[0].plot(t[hjd < divider],b_transit[hjd < divider]-1.0,'k-')
  axarr[0].plot(t[hjd < divider],(g_transit[hjd < divider]-0.05)-1.0,'k-')
  axarr[0].plot(t[hjd < divider],(r_transit[hjd < divider]-0.1)-1.0,'k-')

  axarr[0].plot(t[hjd > divider],(k_transit[hjd > divider]+0.05)-1.0,'k-')
  axarr[0].plot(t[hjd > divider],b_transit[hjd > divider]-1.0,'k-')
  axarr[0].plot(t[hjd > divider],(g_transit[hjd > divider]-0.05)-1.0,'k-')
  axarr[0].plot(t[hjd > divider],(r_transit[hjd > divider]-0.1)-1.0,'k-')

  axarr[1].plot([start_transit,start_transit],[-1.8,1.8],'k--')
  axarr[1].plot([mid_transit,mid_transit],[-1.8,1.8],'k--')
  axarr[1].plot([end_transit,end_transit],[-1.8,1.8],'k--')

  axarr[1].plot([start_transit_2,start_transit_2],[-1.8,1.8],'k--')
  axarr[1].plot([mid_transit_2,mid_transit_2],[-1.8,1.8],'k--')
  axarr[1].plot([end_transit_2,end_transit_2],[-1.8,1.8],'k--')

  axarr[1].errorbar(t_bin,(k_lightcurve-k_transit)+0.05,yerr=k_frac_errors,fmt='ko',label=k_range)
  axarr[1].errorbar(t_bin,(b_lightcurve-b_transit),yerr=b_frac_errors,fmt='bo',label=b_range)
  axarr[1].errorbar(t_bin,(g_lightcurve-g_transit)-0.05,yerr=g_frac_errors,fmt='go',label=g_range)
  axarr[1].errorbar(t_bin,(r_lightcurve-r_transit)-0.1,yerr=r_frac_errors,fmt='ro',label=r_range)

  axarr[1].plot(t[hjd < divider],guide[hjd < divider]+0.05,'k-')
  axarr[1].plot(t[hjd < divider],guide[hjd < divider],'k-')
  axarr[1].plot(t[hjd < divider],guide[hjd < divider]-0.05,'k-')
  axarr[1].plot(t[hjd < divider],guide[hjd < divider]-0.1,'k-')

  axarr[1].plot(t[hjd > divider],guide[hjd > divider]+0.05,'k-')
  axarr[1].plot(t[hjd > divider],guide[hjd > divider],'k-')
  axarr[1].plot(t[hjd > divider],guide[hjd > divider]-0.05,'k-')
  axarr[1].plot(t[hjd > divider],guide[hjd > divider]-0.1,'k-')

  axarr[0].legend(loc='upper right')

  b_chi2 = sum(((b_lightcurve - 1)/b_frac_errors)**2)/(len(b_lightcurve) - 1 - poly_order)
  g_chi2 = sum(((g_lightcurve - 1)/g_frac_errors)**2)/(len(g_lightcurve) - 1 - poly_order)
  r_chi2 = sum(((r_lightcurve - 1)/r_frac_errors)**2)/(len(r_lightcurve) - 1 - poly_order)

  print "Chi2 for hypothesis of flat (b,g,r):"

  print b_chi2, g_chi2, r_chi2

  b_rms = std(b_lightcurve - b_transit)
  g_rms = std(g_lightcurve - g_transit)
  r_rms = std(r_lightcurve - r_transit)

  print "fractional rms (b,g,r):"

  print b_rms, g_rms, r_rms

  cadence = max(t)/len(t)

  transit_dur = 60*3.0

  sn_increase = sqrt(transit_dur)

  print "over the timescale of a transit, that's(b,g,r):"

  print b_rms/sn_increase, g_rms/sn_increase, r_rms/sn_increase

  print "giving OTOO this significance per scale height (assume scale height ~ 1e-4) (b,g,r)"

  scale_height = 1e-4

  print scale_height/(b_rms/sn_increase), scale_height/(g_rms/sn_increase), scale_height/(r_rms/sn_increase)

  axarr[0].set_ylabel('Normalized flux - 1')
  axarr[1].set_ylabel('Residuals')
  xlabel('Time from start (minutes)')

  axarr[0].set_xlim([-5,max(t)+5])
  axarr[0].set_ylim([-0.16,0.16])
  axarr[1].set_ylim([-0.14,0.09])

  show()

def fit_transit_colours(t,airmass,flux_k,error_k,flux_r,error_r,flux_g,error_g,flux_b,error_b,prior):

  f = minimize(transit_func_colour, prior, method='Powell', args=(flux_k,error_k,flux_r,error_r,flux_g,error_g,flux_b,error_b,t,airmass))
  x = f['x']

  #f, jacob = leastsq(transit_func_colour, prior, args=(flux_k,error_k,flux_r,error_r,flux_g,error_g,flux_b,error_b,t,airmass))
  #x = f

  x_k = [x[0],x[4],x[5],x[6],x[7],x[8],x[9],x[16],x[17],x[18],x[19],x[20],x[21]]
  x_r = [x[1],x[4],x[5],x[6],x[7],x[10],x[11],x[22],x[23],x[24],x[25],x[26],x[27]]
  x_g = [x[2],x[4],x[5],x[6],x[7],x[12],x[13],x[28],x[29],x[30],x[31],x[32],x[33]]
  x_b = [x[3],x[4],x[5],x[6],x[7],x[14],x[15],x[34],x[35],x[36],x[37],x[38],x[39]]

  model_k, z = transit_model(x_k,t,airmass)
  model_r, z = transit_model(x_r,t,airmass)
  model_g, z = transit_model(x_g,t,airmass)
  model_b, z = transit_model(x_b,t,airmass)

  print 'white depth =',x[0]
  print 'blue depth =',x[1]
  print 'green depth =',x[2]
  print 'red depth =',x[3]

  return model_k,model_r,model_g,model_b

def fit_transit(t,airmass,flux,error,prior,log=False,fix=[]):

  f = minimize(transit_func, prior, method='BFGS', args=(flux,error,t,airmass,log,fix))
  x = f['x']
  jacob = f['jac']
  hess = f['hess_inv']

  errors = []

  for i in range(0,len(hess)):
    errors += [sqrt(hess[i][i])]

  errors = array(errors)

  #f, jacob  = leastsq(transit_func, prior, args=(flux,error,t,airmass))
  #x = f

  model, z = transit_model(x,t,airmass)

  print 'white depth =',x[0]

  return model, x, errors

def transit_func(x,data,error,variable,airmass,log,fix,plotting=False):

  model, z = transit_model(x,variable,airmass,fix)

  dif = (data - model)/error

  chi2 = sum(dif**2.0)

  if plotting == True:
    print x
    plot(variable,model,'b-')
    plot(variable,data,'ro')
    show()
    plot(variable,dif,'ro')
    show()

  rchi2 = chi2/(len(data) - len(x) - 1)

  print rchi2

  if log == True:
    return log10(chi2)
  else:
    return chi2
  #return [chi2]*50

def transit_func_colour(x,data_k,error_k,data_r,error_r,data_g,error_g,data_b,error_b,variable,airmass,plotting=True):

  x_k = [x[0],x[4],x[5],x[6],x[7],x[8],x[9],x[16],x[17],x[18],x[19],x[20],x[21]]
  x_r = [x[1],x[4],x[5],x[6],x[7],x[10],x[11],x[22],x[23],x[24],x[25],x[26],x[27]]
  x_g = [x[2],x[4],x[5],x[6],x[7],x[12],x[13],x[28],x[29],x[30],x[31],x[32],x[33]]
  x_b = [x[3],x[4],x[5],x[6],x[7],x[14],x[15],x[34],x[35],x[36],x[37],x[38],x[39]]

  model_k, z = transit_model(x_k,variable,variable)
  model_r, z = transit_model(x_r,variable,variable)
  model_g, z = transit_model(x_g,variable,variable)
  model_b, z = transit_model(x_b,variable,variable)

  dif_k = (data_k - model_k)/error_k
  dif_r = (data_r - model_r)/error_r
  dif_g = (data_g - model_g)/error_g
  dif_b = (data_b - model_b)/error_b

  chi2 = sum(dif_k**2.0 + dif_r**2.0 + dif_g**2.0 +dif_b**2.0)

  if plotting == True:
    print x
    plot(variable,model_k,'b-')
    plot(variable,data_k,'ro')
    show()
    plot(variable,dif_k,'ro')
    show()

  print chi2
  return chi2
#  return [chi2]*50
 
def transit_model(x,hjd,airmass,fix=[]):

  p0 = abs(x[0])
  gamma = x[5:7]

  if len(fix) > 0:
    x = fix

  tt = x[1]
  per = x[2]
  inc = x[3]
  ars = x[4]
  airmass_2 = x[7]
  norm_1 = x[8]
  norm_2 = x[9]
  t_10 = x[10] 
  t_11 = x[11]
  t_20 = x[12]
  t_21 = x[13]
  poly_p = x[14:]

  poly_p_1 = x[14:14+len(poly_p)/2]
  poly_p_2 = x[14+len(poly_p)/2:]

  #poly_p_2[1] = poly_p_1[1]
  #poly_p_2[2] = poly_p_1[2]

  per = 1.7497798
  inc = 85.35
  ars = 7.38
#  gamma = [0,0]
  tt = -1.25135820e-04

  tt += 2456892.54327

  gamma[0] = abs(gamma[0])
  gamma[1] = abs(gamma[1])
  
  if gamma[1] > gamma[0]:
    gamma[1] = gamma[0]

  #print airmass_2
  #print norm_1, norm_2
  #print t_1, t_2

  z = transit.t2z(tt, per, inc, hjd, ars)

  divider = 2456895

  poly_m_1 = poly_model(poly_p_1,hjd[hjd < divider] - min(hjd[hjd < divider]))
  poly_m_2 = poly_model(poly_p_2,hjd[hjd > divider] - min(hjd[hjd > divider]))

  poly_t_1 = poly_model([t_10,t_11],hjd[hjd < divider] - min(hjd[hjd < divider]))
  poly_t_2 = poly_model([t_20,t_21],hjd[hjd > divider] - min(hjd[hjd > divider]))

  t_dep = array(list(poly_t_1) + list(poly_t_2))

  extinction = array(list(poly_m_1) + list(poly_m_2))

  poly_m = 1.0 + 10.0**(0.4*(airmass + airmass_2*airmass**2.0)*extinction) 
  
  print x
  model = transit.occultquad(z, p0, gamma)

  detrended = t_dep*model / poly_m

  for i in range(0,len(hjd)):
    if hjd[i] < divider:
      detrended[i] = detrended[i]*norm_1
    else:
      detrended[i] = detrended[i]*norm_2

#  detrended = detrended*sum(model)/sum(detrended)

  return detrended, z


def detrend_poly(full_data,full_variable,full_errors,order,transit_window=False,time=False,repeats=1,clip=3.0,istime=False):

  prior = [median(full_data)] + [1.0]*(order)

  if istime != False:
    not_transit = logical_or((time<transit_window[0]),(time>transit_window[1]))
    data = full_data[not_transit].copy()
    variable = full_variable[not_transit].copy()
    error = full_errors[not_transit].copy()
  else:
    data = full_data.copy()
    variable = full_variable.copy()
    error = full_errors.copy()

  final, jacob = leastsq(polyfunc, prior, args=(data,variable))

  for i in range(0,repeats):
    final, jacob = leastsq(polyfunc, prior, args=(data,variable))
    final_model = poly_model(final,variable)
    dif = abs((final_model - data)/(error))
    data = data[dif < clip].copy()
    variable = variable[dif < clip].copy()
    error = error[dif < clip].copy()
    prior = final.copy()

  final_model = poly_model(final,full_variable)
  detrended = full_data / final_model

#  detrended = detrended*sum(full_data)/sum(detrended)


  print final

  plot(full_variable,full_data,'ro')
  plot(full_variable,final_model,'b-')
  show()

  return detrended

def fit_line(full_data,full_errors,full_wvl,order,repeats=1,clip=5.0):

  line_strength = -1.0
  line_center = 6684
  clipped = [(full_wvl > (100 + min(full_wvl))) & (full_wvl < (max(full_wvl - 100)))]
  line_center = full_wvl[clipped][argmin(full_data[clipped])]
  line_width = 10.0

  prior = [line_strength] + [line_center]  + [line_width] + [median(full_data)] + [1.0]*(order)

  data = full_data.copy()
  wvl = full_wvl.copy()
  error = full_errors.copy()

  for i in range(0,repeats):
    final, jacob = leastsq(linefunc, prior, args=(data,wvl))
    final_model = linefunc(final,data,wvl,return_model=True)
    dif = abs((final_model - data)/(error*final_model))
    data = data[dif < clip].copy()
    wvl = wvl[dif < clip].copy()
    error = error[dif < clip].copy()
    prior = final.copy()


  final_model = linefunc(final,full_data,full_wvl,return_model=True)

  #plot(full_wvl,full_data,'ro')
  #plot(full_wvl,final_model,'b-')
  #show()

  fit_cen = final[1]

  return fit_cen

def linefunc(x,data,variable,return_model=False):
  continuum = poly_model(x[3:],variable)
  absorbtion = gauss_model([(abs(x[0])*-1),x[1],x[2],1.0],variable)
  model = continuum*absorbtion
  dif = data - model

  if return_model == True:
    return model

  return dif

def polyfunc(x,data,variable):
  model = poly_model(x,variable)
  dif = data - model
  return dif

def poly_model(x,variable):

  model = array([0.0]*len(variable))
  for i in range(0,len(x)):
    model += x[i]*variable**i
  return model

def gauss_model(x,variable):
  a = x[0]
  b = x[1]
  c = x[2]
  d = x[3]  
  model = a*exp(-((variable-b)**2.0)/(2.0*(c**2))) + d
  return model

def align_lightcurves(flux_table,flux_frac_errors,wvl_table,reg_min,reg_max):

  lightcurve = []
  error = []

  centers = []
  order = 2
  for i in range(0,len(flux_table)):
    reg = [(wvl_table[i] > reg_min) & (wvl_table[i] < reg_max)]
    center = fit_line(flux_table[i][reg],flux_frac_errors[i][reg],wvl_table[i][reg],order,repeats=3.0)
    centers += [center]

  center_shifts = centers[0] - array(centers)

  for i in range(0,len(wvl_table)):
    wvl_table[i] += center_shifts[i]

  return wvl_table, centers

def create_lightcurve(flux_table,flux_frac_errors,wvl_table,reg_min,reg_max,hjd):


  divider = 2456895
  lightcurve = []
  error = []

  first_reg = [(wvl_table[0] > reg_min) & (wvl_table[0] < reg_max)]

  clean_flux_first, clean_errors_first = clean_light_curve(flux_table[0][first_reg],flux_frac_errors[0][first_reg],wvl_table[0][first_reg])

  spec_1_1 = [wvl_table[hjd < divider][0],flux_table[hjd < divider][0]]

  spec_1_2 = [wvl_table[hjd > divider][0],flux_table[hjd > divider][0]]

  wvl_list = []

  curve_list = []

  ends = []

  for i in range(0,len(flux_table)):
    spec_2 = [wvl_table[i],flux_table[i]]
    if hjd[i] < divider:
      dif = cross_correlate(spec_1_1,spec_2,low_lim=reg_min-200,high_lim=reg_max+200)
    else:
      dif = cross_correlate(spec_1_2,spec_2,low_lim=reg_min-200,high_lim=reg_max+200)

    wvls = wvl_table[i]

    reg = [(wvls > reg_min) & (wvls < reg_max)]

    clean_flux, clean_errors = clean_light_curve(flux_table[i][reg],flux_frac_errors[i][reg],wvls[reg])

    indexes = arange(0,len(wvls))
    high_step = abs(wvls[max(indexes[reg])] - wvls[max(indexes[reg]) + 1])
    high_end = abs(max(wvls[reg]) - reg_max)
    high_flux = (high_end/high_step)*flux_table[i][max(indexes[reg]) + 1]

    low_step = abs(wvls[min(indexes[reg])] - wvls[min(indexes[reg]) - 1])
    low_end = abs(min(wvls[reg]) - reg_min)
    low_flux = (low_end/low_step)*flux_table[i][min(indexes[reg]) - 1]

    #plot(wvls[reg],clean_flux)

    wvl_list += [wvls[reg]]
    curve_list += [(clean_flux.copy())]

    new_error = []
    index = arange(0,len(clean_flux))
    for x in index:
      new_error += [(clean_errors[x]*clean_flux[x])**2.0]
    new_error = array(new_error)

    ends += [high_flux + low_flux]

    lightcurve += [sum(clean_flux) + high_flux + low_flux]

    error += [sqrt(sum(new_error))]
  
  ends = array(ends)

  wvl_list = array(wvl_list)
  curve_list = array(curve_list)

  first_night = curve_list[hjd<divider].copy()
  second_night = curve_list[hjd>divider].copy()

  ends_first = ends[hjd<divider].copy()
  ends_second = ends[hjd>divider].copy()

  wvls_first = wvl_list[hjd<divider].copy()
  wvls_second = wvl_list[hjd>divider].copy()

  first_night = remove_outliers(first_night)

  for i in range(0,len(first_night)):
    plot(wvls_first[i],first_night[i]/median(first_night[i]))
  savefig((str(int((reg_min+reg_max)/2.0))+'night1.pdf'))
  close()

  second_night = remove_outliers(second_night)

  for i in range(0,len(second_night)):
    plot(wvls_second[i],second_night[i]/median(second_night[i]))
  savefig((str(int((reg_min+reg_max)/2.0))+'night2.pdf'))
  close()

  lightcurve = []
  for i in range(0,len(first_night)):
    lightcurve += [sum(first_night[i]) + ends_first[i]]
  for i in range(0,len(second_night)):
    lightcurve += [sum(second_night[i])+ ends_second[i]]

  #show()
  #plot(lightcurve)
  #show()

  return array(lightcurve), array(error)

def remove_outliers(raw_flux_array):

  avg_spec = []
  norm_spec = []

  array_lens = []

  for i in range(0,len(raw_flux_array)):
    array_lens += [len(raw_flux_array[i])]

  min_len = min(array_lens)

  flux_array = raw_flux_array.copy()

  for i in range(0,len(flux_array)):
    flux_array[i] = raw_flux_array[i][:min_len]

  for i in range(0,len(flux_array)):
    avg_spec += [median(flux_array[i])]
    norm_spec += [flux_array[i]/avg_spec[i]]
  norm_spec = array(norm_spec)
  std_spec = []

  med_spec = median(norm_spec,axis=0)
  std_spec = std(norm_spec,axis=0)

  for i in range(0,len(flux_array)):  
    deviation = abs((norm_spec[i] - med_spec)/std_spec)
    indexes = arange(0,len(flux_array[:][i]))
    bad_indexes = indexes[deviation>3]
    should_be = med_spec*avg_spec[i]
    #plot(should_be)
    #plot(flux_array[:][i])
    for x in bad_indexes:
      raw_flux_array[:][i][x] = should_be[x]
    #plot(flux_array[:][i])
    #show()

  return raw_flux_array

def clean_light_curve(flux,errors,wvl,boxlen=20,clip=3):

  clean_flux = flux.copy()
  clean_errors = errors.copy()

  return clean_flux, clean_errors

  boxr = boxlen/2

  for i in range(boxr,len(flux)-boxr):
    seg = array(list(flux[i-boxr:i]) + list(flux[i+1:i+boxr+1]))

    x=[1.0,wvl[i],14,0.0]
    g_m = gauss_model(x,wvl[i-boxr:i+boxr+1])
    g_m = g_m/(sum(g_m[0:boxr]) + sum(g_m[boxr+1:]))
    g_m = array(list(g_m[0:boxr]) + list(g_m[boxr+1:]))
    g_m = g_m*seg
    seg_mean = sum(g_m)

    sig = abs(flux[i] -seg_mean)/(flux[i]*errors[i])
    if sig > clip:
      clean_flux[i] = seg_mean
#      clean_flux[i] = 0
      clean_errors[i] = (std(seg)/sqrt(len(seg)))/seg_mean
  
  return clean_flux, clean_errors
      

def load_data(filelist,cache=False):

  nfilelist = glob('*')

  if cache:
    if cache in nfilelist:
      cached_data = fits.open(cache)
      app1_table = cached_data[0].data
      app2_table = cached_data[1].data
      airmass = cached_data[2].data
      hjd = cached_data[3].data
      return hjd, airmass, app1_table, app2_table

  flux_table = []
  wvl_table = []
  hjd = []
  airmass = []
  app1_table = []
  app2_table = []

  i = 0

  gain = 0.92
  #because we're in slow imaging mode, see http://www.ing.iac.es/Engineering/detectors/g3_ultra_auxcam.html

  ra = 348.494833333
  dec = 8.76127777778

#  for file in filelist[:10]:
  for file in filelist:
    if i > 0:
      raw = file.lstrip('t_app1') + '.fit'
      print raw
      jd = fits.open(raw)[0].header['jd']
      hjd_line = helio_jd(jd-2400000,ra,dec) + 2400000
      hjd += [hjd_line]
      airmass += [fits.open(raw)[0].header['AIRMASS']]
      app1_name = 't_app1_' + file.lstrip('t_app1')
      app2_name = 't_app2_' + file.lstrip('t_app1')
      app1 = loadtxt(app1_name)
      app2 = loadtxt(app2_name)

      err1_name = 't_err1_' + file.lstrip('t_app1')
      err2_name = 't_err2_' + file.lstrip('t_app1')
      err1 = loadtxt(err1_name)
      err2 = loadtxt(err2_name)

      app1_electrons = app1[:,1]*gain
      app2_electrons = app2[:,1]*gain

      app1_error_electrons = err1[:,1]*gain
      app2_error_electrons = err2[:,1]*gain

      app1_table += [[app1[:,0],app1[:,1]*gain,err1[:,1]/app1[:,1],sqrt((1.0/app1_electrons))]]
      app2_table += [[app2[:,0],app2[:,1]*gain,err2[:,1]/app2[:,1],sqrt((1.0/app2_electrons))]]
#      plot(app1[:,0],app1_electrons)
#      plot(app2[:,0],app2_electrons)

    i += 1

  hjd = array(hjd)
  airmass = array(airmass)
  app1_table = array(app1_table)
  app2_table = array(app2_table)      

  if cache:
    hdu1 = fits.PrimaryHDU(app1_table)
    hdu2 = fits.ImageHDU(app2_table)
    hdu3 = fits.ImageHDU(airmass)
    hdu4 = fits.ImageHDU(hjd)
    hdulist = fits.HDUList([hdu1,hdu2,hdu3,hdu4])
    hdulist.writeto(cache)

    #  show()
  return hjd, airmass, app1_table, app2_table


def cross_correlate(spec_1,spec_2,low_lim=4000,high_lim=8400):
#  to first order, the deltas should be the same... find an offset.

  # first, let's make sure the spectra are normalised with respect to each other

  flux_1 = spec_1[1][(spec_1[0] > low_lim) & (spec_1[0] < high_lim)]/median(spec_1[1][(spec_1[0] > low_lim) & (spec_1[0] < high_lim)])
  flux_2 = spec_2[1][(spec_2[0] > low_lim) & (spec_2[0] < high_lim)]/median(spec_2[1][(spec_2[0] > low_lim) & (spec_2[0] < high_lim)])

  wvl_1 = spec_1[0][(spec_1[0] > low_lim) & (spec_1[0] < high_lim)]
  wvl_2 = spec_2[0][(spec_2[0] > low_lim) & (spec_2[0] < high_lim)]

  dif_a = []
  shift_a = []

  clip = 20

  for x in range(-(clip-1),clip-1):
    dif = 0.0

    #plot(flux_1[clip:-clip])
    #plot(flux_2[clip+x:(-clip)+x])
    #show()
    dif += sum(abs(flux_1[clip:-clip]/median(flux_1[clip:-clip]) - flux_2[clip+x:(-clip)+x]/median(flux_2[clip+x:(-clip)+x])))
    dif_a += [dif]
    shift_a += [x]

  #plot(shift_a,dif_a)
  #show()

  best_offset = shift_a[argmin(dif_a)]

  scale = median(diff(wvl_1))

  final_dif = best_offset*scale

  #plot(wvl_1,flux_1)
  #plot(wvl_2-final_dif,flux_2)
  #show()
  return final_dif

def plot_ratio(wvl_table,flux_table):
  mean_flux = median(flux_table, axis=0)

  for i in range(0,len(flux_table)):
    plot(wvl_table[i],(flux_table[i]-mean_flux)/mean_flux)

  show()

def binning(in_array,bins):
  
  result = []
  for i in range(1,floor(len(in_array)/bins)):
    sub = in_array[(i-1)*bins:i*bins]
    sub = sub[isreal(sub)]
    result += [median(sub)]
  return array(result)

main()
