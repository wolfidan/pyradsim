# -*- coding: utf-8 -*-
"""
Created on Mon May 23 17:09:27 2016

@author: wolfensb
"""

import json
import numpy as np
import warnings
import copy
import matplotlib.pyplot as plt
from collections import OrderedDict

from scipy.special import gamma
from scipy.interpolate import PchipInterpolator, interp1d
from  pyradsim.utilities import Lambda
from  pyradsim.constants_cosmo import *

def create_PSD(name,*args):
    if name == 'NormalizedGammaPSD':
        return NormalizedGammaPSD(*args)
    elif name == 'UnnormalizedGammaPSD':
        return UnnormalizedGammaPSD(*args)
    elif name == 'ExponentialPSD':
        return ExponentialPSD(*args)
    elif name == 'BinnedPSD':
        return BinnedPSD(*args)
    elif name == 'COSMO_1mom_graupel':
        return COSMO_1mom_graupel(*args)
    elif name == 'COSMO_1mom_snow':
        return COSMO_1mom_snow(*args)
    elif name == 'COSMO_1mom_rain':
        return COSMO_1mom_rain(*args)
    else:
        raise NameError('Invalid psd model name!')
        
class PSD(object):
    def __init__(self,psd_type, psd_func,dmin = 0,dmax = 10):
        self.func = psd_func
        self.type = psd_type
        self.dmin = dmin
        self.dmax = dmax
    def __call__(self,*inputs):
        return self.func(*inputs)
        
    def moment(self,order,nbins_diam = 1024):
        dbins = np.linspace(self.dmin,self.dmax,nbins_diam)
        dD = dbins[1] - dbins[0]
        return np.trapz(dbins**order * self.psd_func(dbins), dx = dD)

    def get_dic(self):
        if self.type == 'BinnedPSD':
            dic = OrderedDict([('PSD type',self.type),('Bin edges',str(self.bin_edges)),
                   ('Bin values',str(self.bin_values))])
        else:
            dic = OrderedDict([('PSD type',self.type),('Function',str(self.func.expression)),
                   ('Dmin',str(self.dmin)),('Dmax',str(self.dmax))])
        return dic
    
    def __mul__(self,scalar):
        obj = copy.deepcopy(self)
        obj.func = self.func*scalar
        return obj
    def __rdiv__(self,scalar):
        obj = copy.deepcopy(self)
        obj.func = self.func/scalar
        return obj
    def __add__(self,scalar):
        obj = copy.deepcopy(self)
        obj.func = self.func+scalar
        return obj
    def __sub__(self,scalar):
        obj = copy.deepcopy(self)
        obj.func = self.func-scalar
        return obj
        
    def __str__(self):
        msg = ''

        dic = self.get_dic()

        msg = json.dumps(dic,indent=3) +'\n'
                 
#            d = np.linspace(self.dmin,self.dmax,100)
#            plt.plot(d,self.psd_func(d))
#            plt.xlabel('diameter')
#            plt.ylabel('number')
        return msg
           
            
            
'''
Classical PSD forms
'''

class NormalizedGammaPSD(PSD):
    def __init__(self,Nw,D0,mu,dmin = 0,dmax = 10):
        self.type = 'NormalizedGammaPSD'
        f =  6/3.67**4 * (3.67 + mu)**(mu+4)/(gamma(mu+4))
        Lamb = (3.67+mu)/D0
        expression = str(Nw) +' * '+str(f)+' * (D/'+str(D0)+')**'+str(mu)+' * np.exp(-'+str(Lamb)+' * D)'
        psd_func = Lambda(lambda D, T: Nw * f * (D/D0)**mu * np.exp(-Lamb * D),expression)
        super(NormalizedGammaPSD,self).__init__('NormalizedGammaPSD',psd_func,dmin,dmax)

class UnnormalizedGammaPSD(PSD):
    def __init__(self,N0,Lamb,mu,dmin,dmax):
        self.type = 'UnnormalizedGammaPSD'
        expression = str(N0) +' * D**'+str(mu)+' * np.exp(-'+str(Lamb)+' * D)'
        psd_func = Lambda(lambda D, T: N0 * D**mu * np.exp(-Lamb * D),expression)
        super(UnnormalizedGammaPSD,self).__init__('UnnormalizedGammaPSD',psd_func,dmin,dmax)

class ExponentialPSD(PSD):
    def __init__(self,N0,Lamb,dmin = 0,dmax = 10):
        self.type = 'ExponentialPSD'

        expression = str(N0)+' * np.exp(-'+str(Lamb)+'*D)'
        psd_func = Lambda(lambda D, T: N0 * np.exp(-Lamb*D),expression)
        
        super(ExponentialPSD,self).__init__('ExponentialPSD',psd_func,dmin,dmax)    

class BinnedPSD(PSD):
    def __init__(self,bin_edges,bin_values,interp_method = 'linear'):
        """ interp_method can be either : 
                  - rect
                  - linear 
                  - pchip
        """
        self.type = 'BinnedPSD'
        self.bin_edges = bin_edges
        self.bin_values = bin_values
        
        if interp_method not in ['rect','linear','pchip']:
            msg = "Invalid interp_method argument, please provide" \
                          " a valid method: 'rect', 'linear' or 'pchip' "
            warnings.warn(msg)
            warnings.warn("Using default = 'linear'")
            interp_method = 'linear'
            
        expression = 'Binned psd'
        if interp_method == 'rect':
            psd_func = Lambda(lambda D, T: interp1d(bin_edges, bin_values, kind = 'zero')(D),expression)
        elif interp_method ==  'linear':
            psd_func = Lambda(lambda D, T: interp1d(bin_edges, bin_values, kind = 'linear')(D),expression)
        elif interp_method ==  'pchip':
            psd_func = Lambda(lambda D, T: PchipInterpolator(bin_edges, bin_values)(D),expression)
            
        super(BinnedPSD,self).__init__('BinnedPSD',psd_func,np.min(bin_edges),np.max(bin_edges))    
        
    
'''
COSMO 1-moment PSD
'''

class COSMO_1mom_graupel(PSD):
    def __init__(self,Q,dmin = 0,dmax = 10):
        Lamb = (LAMBDA_FACTOR_G/Q)**(1./(4.+MU_G))
        expression = str(N0_G) +' * D**'+str(MU_G)+' * np.exp(-'+str(Lamb)+' * D)'
        psd_func = Lambda(lambda D,T: N0_G * D**MU_G * np.exp(-Lamb * D),expression)
        super(COSMO_1mom_graupel,self).__init__('COSMO_1mom_graupel',psd_func,dmin,dmax)

class COSMO_1mom_rain(PSD):
    def __init__(self,Q,dmin = 0,dmax = 10):
        Lamb = (LAMBDA_FACTOR_R/Q)**(1./(4.+MU_R))
        expression = str(N0_R) +' * D**'+str(MU_R)+' * np.exp(-'+str(Lamb)+' * D)'
        psd_func = Lambda(lambda D,T: N0_R * D**MU_R * np.exp(-Lamb * D),expression)
        super(COSMO_1mom_rain,self).__init__('COSMO_1mom_rain',psd_func,dmin,dmax)

class COSMO_1mom_snow(PSD):
    def __init__(self,Q,T,dmin = 0,dmax = 10):
        N0_S = lambda T: 13.5*(5.65*10**5*np.exp(-0.107*(T-273.15)))/1000 # mm^-1 m^-3
        Lamb = (LAMBDA_FACTOR_R/Q)**(1./(4.+MU_R))
        expression = str(N0_S(T)) +' * D**'+str(MU_S)+' * np.exp(-'+str(Lamb)+' * D)'
        psd_func = Lambda(lambda D,T: N0_S(T) * D**MU_S * np.exp(-Lamb * D),expression)
        super(COSMO_1mom_snow,self).__init__('COSMO_1mom_rain',psd_func,dmin,dmax)
        
if __name__ == '__main__':
    a = create_PSD('ExponentialPSD',1,1)
    print(a)
    d = np.linspace(1,10,100)
    N=a(d,200)
#    plt.figure()
#    plt.plot(d,a(d))
#    plt.hold(True)
#    a = a-1000
#    plt.plot(d,a(d))
    
