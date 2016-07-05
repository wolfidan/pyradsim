# -*- coding: utf-8 -*-
"""
Created on Thu May 26 15:38:06 2016

@author: wolfensb
"""
import numpy as np
from pyradsim.utilities import Lambda
AR_MODELS = {}

def thurai_2007(D):
    ar = np.zeros((len(D),))
    
    ar[D<0.7] = 1.0
    mid_diam = np.logical_and(D<1.5,D>=0.7)
    ar[mid_diam] = 1.173 - 0.5165*D[mid_diam] + 0.4698*D[mid_diam]**2 - 0.1317*D[mid_diam]**3 - \
            8.5e-3*D[mid_diam]**4
    ar[D>=1.5] = 1.065 - 6.25e-2*D[D>=1.5] - 3.99e-3*D[D>=1.5]**2 + 7.66e-4*D[D>=1.5]**3 - \
            4.095e-5*D[D>=1.5]**4 
        
    return 1./ar

def brandes_2002(D):
    ar = 0.9951 + 0.02510*D - 0.03644*D**2 +  0.005030*D**3 \
         - 0.0002492*D**4
    return 1./ar
    
#def beard_1990(D):
#    D = D*0.1
#    ar = 1.01668 -0.98055*D - 2.52686*D**2 +3.75061*D**3+1.68692*D**4
#    return 1./ar

def andsager_1999(D):
    D = D*0.1
    
    # Ratios using empirical data, valid from 1.1 mm to 4.4 mm.
    ar = 1.012 - 0.144*D - 1.03*D**2
  
    # Ratios using theoretical equilibrium axis ratios, valid outside 
    # the 1 mm to 4.4mm range.
    ar_eq = 1.0048 + 0.0057*D - 2.628*D**2 + 3.682*D**3 - 1.677*D**4
  
    # Use theory outside empirical range.
    ar[D < 0.11] = ar_eq[D < 0.11]
    ar[D > 0.44] = ar_eq[D  > 0.44]
    return 1./ar

            
           
AR_MODELS['Thurai_2007'] = Lambda(lambda D: thurai_2007(D),'Aspect-ratio model for rain (Thurai, 2007)')
AR_MODELS['Andsager_1999'] = Lambda(lambda D: andsager_1999(D),'Aspect-ratio model for rain (Andsager, 1999)')
AR_MODELS['Brandes_2002'] = Lambda(lambda D: brandes_2002(D),'Aspect-ratio model for rain (Branded, 2002)')
