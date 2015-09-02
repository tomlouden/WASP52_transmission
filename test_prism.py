# -*- coding: utf-8 -*-
from fprism import fprism
from pylab import *
from numpy import *

radiusratio = 0.164

ars = 7.38

radiustotal = (1.0 + 0.164)/ars

u1 = 0.1

u2 = 0.1

inclination = 85.35

midpoint = 2e7

spotlong = 30

spotlat = 50

spotsize = 10

spotflux = 0.8

INPUT = [radiusratio,radiustotal,u1,u2,inclination,midpoint,spotlong,spotlat,spotsize,spotflux]

SPOT_DATA = True

data_x = linspace(2e7-0.05,2e7+0.05,1000)

plotting = True

flux_out = fprism(data_x,INPUT,SPOT_DATA,plotting=True)