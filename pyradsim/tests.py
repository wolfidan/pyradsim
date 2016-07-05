# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:02:31 2016

@author: wolfensb
"""
from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, GammaPSD
from pytmatrix import orientation, radar, tmatrix_aux, refractive

from tictoc import *

def drop_ar(D_eq):
    if D_eq < 0.7:
        return 1.0;
    elif D_eq < 1.5:
        return 1.173 - 0.5165*D_eq + 0.4698*D_eq**2 - 0.1317*D_eq**3 - \
            8.5e-3*D_eq**4
    else:
        return 1.065 - 6.25e-2*D_eq - 3.99e-3*D_eq**2 + 7.66e-4*D_eq**3 - \
            4.095e-5*D_eq**4 

scatterer = Scatterer(wavelength=tmatrix_aux.wl_C, m=refractive.m_w_10C[tmatrix_aux.wl_C])
scatterer.psd_integrator = PSDIntegrator()
scatterer.psd_integrator.axis_ratio_func = lambda D: 1.0/drop_ar(D)
scatterer.psd_integrator.D_max = 10.0
scatterer.psd_integrator.geometries = (tmatrix_aux.geom_horiz_back, tmatrix_aux.geom_horiz_forw)
scatterer.or_pdf = orientation.gaussian_pdf(20.0)
scatterer.orient = orientation.orient_averaged_fixed
tic()
scatterer.psd_integrator.init_scatter_table(scatterer,verbose=True)
toc()
#scatterer.psd = GammaPSD(D0=2.0, Nw=1e3, mu=4)
#radar.refl(scatterer)

#tic()
#scatterer = Scatterer(radius=2.0, wavelength=6.5, m=complex(1.5,0.5), axis_ratio=1.0/0.6)
#scatterer.or_pdf = orientation.gaussian_pdf(20.0)
#scatterer.orient = orientation.orient_averaged_fixed
#for i in range(1000):
#    scatterer.radius = 2.0+i/1000.
#    scatterer.m = complex(1.5,0.5)+i/1000.
#    scatterer.axis_ratio=1.0/(0.6+i/10000.)
#    scatterer.get_SZ()
#toc()