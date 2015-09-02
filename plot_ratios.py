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
from fprism import fprism

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


  width = 200
  zero = 4000

  for i in range(0,45/2):
#  for i in range(1,3):
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

  #regs_cen = [4100.0, 4300.0, 4500.0, 4700.0, 4900.0, 5100.0, 5300.0, 5500.0, 5700.0, 5900.0, 6100.0, 6300.0, 6500.0, 6700.0, 6900.0, 7100.0, 7300.0, 7500.0, 7700.0, 7900.0, 8100.0, 8300.0]
  #depth = [0.16198369194043519, 0.16458595810730817, 0.16929305832975097, 0.17642224828554715, 0.17075893639442172, 0.167184126635751, 0.16846019261522283, 0.16727633655857316, 0.16641163122933766, 0.16638606107343043, 0.16447163801392337, 0.1659595191674049, 0.16569349545243778, 0.16508181566668687, 0.1618419776306903, 0.16800972006214931, 0.16499007006574393, 0.16390080131374507, 0.16378675590487962, 0.16373546785608806, 0.16426892324320638, 0.1659484303341503]
  #depth_error = [0.0023644369502292739, 0.0012981277894602787, 0.00051952022649596314, 0.00065108044491394418, 0.00058937211517202329, 0.00056269171653920702, 0.00093709800469547516, 0.00066771796566806208, 0.00038964069641791744, 0.00039738705896492655, 0.00029115619246652681, 0.00068737957314824133, 0.00074604218959434439, 0.00064428493025280129, 0.0012670859402076544, 0.00031447980363070549, 0.00035619921095150715, 0.00040754518645238537, 0.00097456179847171549, 0.0021260522983182024, 0.00093681996848496679, 1.1151778311738591e-05]

    errorbar(regs_cen,depth,yerr=depth_error,xerr=100,fmt='ko')
    ylabel('Transit depth')
    xlabel('Wavelength (A)')

    ylim([0.15,0.18])

    #show()
    savefig('spectrum.pdf')

    results_array = array(results)

    print 'results time'

    name_list = ['depth','airmass_2','norm_1','norm_2','t_11','t21','splotlang','spotlat','spotsize','spotflux','gamma_0','gamma_1','p10','p11','p12','p20','p21','p22']

    close()

    for i in range(0,len(results_array[0])):
      plot(regs_cen,results_array[:,i],'o')
      savefig(str(name_list[i]))
      close()

  print regs_cen
  print depth
  print depth_error



def plot_lightcurve(lightcurve,frac_error,x,airmass,hjd,name=False):
  # first remove the polynomial from both transit and lightcurve.

  first_min = min(hjd)
  second_min = min(hjd[abs(hjd - 2456892.54327) > 3.0])

  t = hjd.copy()
  for i in range(0,len(t)):
    if (abs(t[i] - 2456892.54327) > 3.0):
      t[i] = ((t[i] - second_min)*(24*60)) + 500
    else:
      t[i] = (t[i] - first_min)*(24*60)

  start_spot = 0
  end_spot = 30
  not_spot = logical_or((t<start_spot),(t>end_spot))
  lightcurve = lightcurve[not_spot]
  frac_error = frac_error[not_spot]
  airmass = airmass[not_spot]
  hjd = hjd[not_spot]
  t = t[not_spot]

  start_spot = 450
  end_spot = 530
  not_spot = logical_or((t<start_spot),(t>end_spot))
  lightcurve = lightcurve[not_spot]
  frac_error = frac_error[not_spot]
  airmass = airmass[not_spot]
  hjd = hjd[not_spot]
  t = t[not_spot]

  #start_spot = 146
  #end_spot = 187
  #not_spot = logical_or((t<start_spot),(t>end_spot))
  #lightcurve = lightcurve[not_spot]
  #frac_error = frac_error[not_spot]
  #airmass = airmass[not_spot]
  #hjd = hjd[not_spot]
  #t = t[not_spot]

  #start_spot = 50
  #end_spot = 500
  #not_spot = logical_or((t<start_spot),(t>end_spot))
  #lightcurve = lightcurve[not_spot]
  #frac_error = frac_error[not_spot]
  #airmass = airmass[not_spot]
  #hjd = hjd[not_spot]
  #t = t[not_spot]

  #start_spot = 676
  #end_spot = 710
  #not_spot = logical_or((t<start_spot),(t>end_spot))
  #lightcurve = lightcurve[not_spot]
  #frac_error = frac_error[not_spot]
  #airmass = airmass[not_spot]
  #hjd = hjd[not_spot]
  #t = t[not_spot]

  start_spot = 40
  end_spot = 500
  not_spot = logical_or((t<start_spot),(t>end_spot))
  lightcurve = lightcurve[not_spot]
  frac_error = frac_error[not_spot]
  airmass = airmass[not_spot]
  hjd = hjd[not_spot]
  t = t[not_spot]

  poly_coefs = x.copy()
  poly_coefs[0] = 0

  poly_m, z = transit_model(poly_coefs,hjd,airmass)

  transit_m, z = transit_model(x,hjd,airmass,SPOT_DATA=False)

  norm_m = transit_m/poly_m

  norm_lightcurve = lightcurve / poly_m

  f, axarr = plt.subplots(2, sharex=True)

  residuals = norm_lightcurve - norm_m

  zero = norm_m - norm_m

  dummy = arange(0,len(norm_lightcurve))

  dummy = t

  outfile = open('wasp52-spot.dat','w')

  for i in range(0,len(hjd)):
    outline = str(hjd[i]-2400000.5)+' '+str(norm_lightcurve[i])+' '+str(frac_error[i]*norm_lightcurve[i])+'\n'
    outfile.write(outline)

  outfile.close()

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

  #start_spot = 146
  #end_spot = 187
  #not_spot = logical_or((t<start_spot),(t>end_spot))
  #b_lightcurve = b_lightcurve[not_spot]
  #b_frac_errors = b_frac_errors[not_spot]
  #airmass = airmass[not_spot]
  #hjd = hjd[not_spot]
  #t = t[not_spot]

  #start_spot = 50
  #end_spot = 500
  #not_spot = logical_or((t<start_spot),(t>end_spot))
  #b_lightcurve = b_lightcurve[not_spot]
  #b_frac_errors = b_frac_errors[not_spot]
  #airmass = airmass[not_spot]
  #hjd = hjd[not_spot]
  #t = t[not_spot]

  #start_spot = 676
  #end_spot = 710
  #not_spot = logical_or((t<start_spot),(t>end_spot))
  #b_lightcurve = b_lightcurve[not_spot]
  #b_frac_errors = b_frac_errors[not_spot]
  #airmass = airmass[not_spot]
  #hjd = hjd[not_spot]
  #t = t[not_spot]

  start_spot = 40
  end_spot = 500
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

  start_spot = 450
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
  gamma_2 = 0.5
  gamma_3 = 0.5

  x = [p0,tt,per,inc,ars,gamma_0,gamma_1,gamma_2,gamma_3,-0.3,1.00,1.00,1.0,0.002,1.0,0.002]

  poly_prior = [-10.0,10.0]

  x += poly_prior*2

  #x = [1.62858863e-01,-1.25135820e-04,1.74977980e+00,8.53500000e+01,7.38000000e+00,0.34,0.66,-0.0625,-0.0907,-1.97555202e-01,9.72396652e-01,9.60859487e-01,1.03377789e+00,-8.28656790e-03,1.04746041e+00,-1.96077433e-02,-5.61442377e+00,-5.27598514e+00,-6.11726671e+00,-4.77924359e+00]
  #b_transit, b_final, b_final_errors = fit_transit(hjd,airmass,b_lightcurve,b_frac_errors,x,log=True)
  #b_transit, b_final, b_final_errors = fit_transit(hjd,airmass,b_lightcurve,b_frac_errors,b_final,log=False,fix=b_final)

  spotlong = 0
  spotlat = 50
  spotsize = 10
  spotflux = 0.9

  x = [1.62858863e-01,-1.97555202e-01,9.72396652e-01,9.60859487e-01,-8.28656790e-03,-1.96077433e-02,spotlong,spotlat,spotsize,spotflux,0.34,0.66,-5.27598514,-6.1172667,-4.77924359,-5.27598514e+00,-6.11726671e+00,-4.77924359e+00]

#  x = [1.60260218e-01,-2.29551791e-01,1.00531878e+00,1.00775710e+00,-5.65343590e-03,6.97034020e-02,-1.00000000e+00,5.00000000e+01,0.90000000e+01,9.00000000e-01,-2.08496158e-01,8.73535223e-01,-6.13166846e+00,-2.60665710e+00,-3.68217786e+00,-5.27299970e+00,-6.11855418e+00,-4.77921417e+00]

#  x = [1.76952892e-01,-3.00075664e-01,1.00474131e+00,9.99598130e-01,-8.28656790e-03,-1.96077433e-02,1.58317370e+00,-1.03777274e+00,-5.63915533e+00,-1.34156135e+01,-5.84934418e+00,-3.66706891e+01,-3.14032082e+00,-3.74213491e+00]
  b_transit, b_final, b_final_errors = fit_transit(hjd,airmass,b_lightcurve,b_frac_errors,x,log=True)
#  b_transit, b_final, b_final_errors = fit_transit(hjd,airmass,b_lightcurve,b_frac_errors,b_final,log=False,fix=b_final)

  print 'finale!'

  print b_final

  return b_transit, b_lightcurve_raw, b_frac_errors_raw, b_final, b_final_errors

def fit_transit(hjd,airmass,flux,error,prior,log=False,fix=[]):

  first_min = min(hjd)
  second_min = min(hjd[abs(hjd - 2456892.54327) > 3.0])

  t = hjd.copy()
  for i in range(0,len(t)):
    if (abs(t[i] - 2456892.54327) > 3.0):
      t[i] = ((t[i] - second_min)*(24*60)) + 500
    else:
      t[i] = (t[i] - first_min)*(24*60)

  flux_all = flux.copy()
  error_all = error.copy()
  airmass_all = airmass.copy()
  hjd_all = hjd.copy()
  t_all = t.copy()

  start_spot = 116
  end_spot = 157
  not_spot = logical_or((t<start_spot),(t>end_spot))
  flux = flux[not_spot]
  error = error[not_spot]
  airmass = airmass[not_spot]
  hjd = hjd[not_spot]
  t = t[not_spot]

  clip = 50.0
  repeats = 2.0
  for i in range(0,repeats):
    limits = ((0.05,0.25),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(0,1),(0,1),(None,None),(None,None),(None,None),(None,None))
    f = minimize(transit_func, prior, method='L-BFGS-B',bounds=limits,args=(flux,error,hjd,airmass,log,fix,False,False))
    x = f['x']
    prior = x
    model, z = transit_model(x,hjd,airmass)

    x2 = x.copy()
    x2[0] = 0.0
    model2, z = transit_model(x2,hjd,airmass)

    #plot(t,model)
    #plot(t,model2)
    #errorbar(t,flux,yerr=error,fmt='r.')
    #show()

    divider = 2456895

    #plot(airmass[t > divider],model[t > divider])
    #errorbar(airmass[t > divider],flux[t > divider],yerr=error[t > divider],fmt='r.')
    #show()

    #plot(airmass[t < divider],model[t < divider])
    #errorbar(airmass[t < divider],flux[t < divider],yerr=error[t < divider],fmt='r.')
    #show()

    sig = abs((model - flux)/error)
    print sig[sig > clip]
    index = [sig < clip]
    flux = flux[index]
    error = error[index]
    hjd = hjd[index]
    airmass = airmass[index]

#  return model, x, model

  #flux_all = flux.copy()
  #error_all = error.copy()
  #airmass_all = airmass.copy()
  #hjd_all = hjd.copy()
  #t_all = t.copy()

  log = False
  fix = x
  prior = x
  #f = minimize(transit_func, prior, method='BFGS',args=(flux_all,error_all,hjd_all,airmass_all,log,fix,True,True))
  f = minimize(transit_func, prior, method='BFGS',options={'eps':0.000001},args=(flux_all,error_all,hjd_all,airmass_all,log,fix,False,False))
  x = f['x']
  jacob = f['jac']
  hess = f['hess_inv']

  errors = []

  for i in range(0,len(hess)):
    errors += [sqrt(hess[i][i])]

  errors = array(errors)

  print errors

  #f, jacob  = leastsq(transit_func, prior, args=(flux,error,t,airmass))
  #x = f

  model, z = transit_model(x,hjd,airmass)

  print 'white depth =',x[0]

  return model, x, errors

def transit_func(x,data,error,variable,airmass,log,fix,plotting=False,spot_data=False):

  model, z = transit_model(x,variable,airmass,fix,spot_data)

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

def transit_model(x,hjd,airmass,fix=[],SPOT_DATA=False):

  spotlong = x[6]
  spotlat = x[7]
  spotsize = abs(x[8])
  spotflux = abs(x[9])

  if spotsize > 20:
    spotsize = 20

  if spotflux > 1.0:
    spotflux = 1.0

  print SPOT_DATA

  print x

  p0 = x[0]

  #p0 = abs(x[0])
  #if p0 > 0.25:
    #p0 = 0.25

  gamma = x[10:12]

  if len(fix) > 0:
    x = fix

  #if abs(gamma[0]) > 1.0:
    #gamma[0] = 2.0*(gamma[0]/abs(gamma[0]))

  #if abs(gamma[0] + gamma[1]) >= 1.0:
    #gamma[1] = 1.0 - gamma[0]


  airmass_2 = x[1]
  norm_1 = x[2]
  norm_2 = x[3]
  t_11 = x[4]
  t_21 = x[5]

  poly_p = x[12:]

  poly_p_1 = x[12:12+len(poly_p)/2]
  poly_p_2 = x[12+len(poly_p)/2:]

  #poly_p_2[1] = poly_p_1[1]
  #poly_p_2[2] = poly_p_1[2]

  per = 1.7497798
  inc = 85.35
  ars = 7.38
#  gamma = [0,0]
  tt = -1.25135820e-04

  tt += 2456892.54327

  t_10 = 1.0
  t_20 = 1.0


  #print airmass_2
  #print norm_1, norm_2
  #print t_1, t_2

  print gamma

  divider = 2456895

  poly_m_1 = poly_model(poly_p_1,hjd[hjd < divider] - min(hjd[hjd < divider]))
  poly_m_2 = poly_model(poly_p_2,hjd[hjd > divider] - min(hjd[hjd > divider]))

  poly_t_1 = poly_model([t_10,t_11],hjd[hjd < divider] - min(hjd[hjd < divider]))
  poly_t_2 = poly_model([t_20,t_21],hjd[hjd > divider] - min(hjd[hjd > divider]))

  t_dep = array(list(poly_t_1) + list(poly_t_2))

  extinction = array(list(poly_m_1) + list(poly_m_2))

  poly_m = 1.0 + 10.0**(0.4*(airmass + airmass_2*airmass**2.0)*extinction) 
  
  z = transit.t2z(tt, per, inc, hjd, ars)
  #model = transit.occultquad(z, p0, gamma)
#  model = transit.occultnonlin(z, p0, gamma)

  radiustotal = (1.0 + p0)/ars

  phase_offset = 0

  INPUT = [p0,radiustotal,gamma[0],gamma[1],inc,phase_offset,spotlong,spotlat,spotsize,spotflux]

  if SPOT_DATA == True:
    plotting = False
    model = fprism((hjd-tt)/(per),INPUT,SPOT_DATA,plotting=plotting)
    if plotting == True:
      quit()
  else:
    model = transit.occultquad(z, p0, gamma)

  #plot(hjd,model2,'b-')
  #plot(hjd,model,'ro')
  #show()

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

  raw_path = '/storage/astro2/phrmat/WHT_08_2014_data/combined_large_fixed'

#  for file in filelist[:10]:
  for file in filelist:
    if i > 0:
      raw = raw_path +'/'+file.lstrip('t_app1') + '.fit'
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
