# -*- coding: utf-8 -*-
from fprism import fprism
from numpy import *

# night 1
p0 = 0.16598202
radiustotal = 0.15896937
gamma_0 = 0.65794940
gamma_1 = -0.31492361
inc = 85.119340
phase = -0.028844446
phase = 0
spotlong = -1.1158468
spotlat = 44.913135
spotsize = 10.939604
spotflux = 0.79823888
fname = 'wasp52-spot_n1.dat'

# night 2
p0 = 0.16306021
radiustotal = 0.15634212
gamma_0 = 0.94426940
gamma_1 = -0.41811691
inc = 85.604770
phase = -0.014684627
phase = 0
spotlong = 59.382099
spotlat = 48.924085
spotsize = 22.062314
spotflux = 0.94007068
fname = 'wasp52-spot_n2.dat'


period = 1.7497798

input = [p0,radiustotal,gamma_0,gamma_1,inc,phase,spotlong,spotlat,spotsize,spotflux]

tt = -1.25135820e-04

tt += 56892.04327

data = loadtxt(fname)
data_x = (data[:,0]-tt)/period

fprism(data_x, input, spot_data=True, type_data='FLUX',plotting=True)