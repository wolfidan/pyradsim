# -*- coding: utf-8 -*-
"""
Created on Mon May 23 17:54:25 2016

@author: wolfensb
"""

import numpy as np
import multiprocessing as mp
from collections import OrderedDict
import json

from  pyradsim.tmatrix import  flatten_matrices, create_scatterer, compute_pol_var
import  pyradsim.constants as constants

MATH_FUNCTIONS = ['cos','sin','exp','log','log10','tan']

global GEOM
global SCATT

from tictoc import tic,toc

class Hydrometeor(object):
    def __init__(self,name,hydro_specs, frequency, temperature, elevation_angle, nbins_d):
        self.name = name
        
        self.psd = hydro_specs['psd']
        self.aspect_ratio = hydro_specs['aspect_ratio']
        self.canting_angle_std = hydro_specs['canting_angle_std']
        self.permittivity = hydro_specs['permittivity']
        
        self.theta = elevation_angle
        self.frequency = frequency
        self.temperature = temperature
        self.nbins_d = nbins_d
        
        self._scatter_signature = None
#        self._integ_signature = None
        
        self._S = None
        self._Z = None
        self._integ_S = None
        self._integ_Z = None
    
    def update(self,hydro_specs, frequency, temperature, elevation_angle, nbins_d):
        self.psd = hydro_specs['psd']
        self.aspect_ratio = hydro_specs['aspect_ratio']
        self.canting_angle_std = hydro_specs['canting_angle_std']
        self.permittivity = hydro_specs['permittivity']
        
        self.theta = elevation_angle
        self.frequency = frequency
        self.temperature = temperature
        self.nbins_d = nbins_d
        
    def __str__(self):
        dic = OrderedDict([('Hydrometeor name',self.name),('PSD',self.psd.get_dic()),
                           ('Aspect ratio',self.aspect_ratio.expression),
                            ('Canting angle std',self.canting_angle_std),
                            ('Permittitivity',self.permittivity.expression)])
        msg = json.dumps(dic,indent=3) +'\n'
        return msg
        
    def _set_scatter_signature(self):
        self._scatter_signature = (self.psd.dmin,self.psd.dmax,self.aspect_ratio.expression,
                                   self.canting_angle_std,self.permittivity.expression,
                                   self.theta,self.temperature,self.frequency,self.nbins_d)
        
#    def _set_integ_signature(self):
#        self._integ_signature = (self.psd.dmin,self.psd.dmax,self.aspect_ratio.expression,
#                                   self.canting_angle_std,self.permittivity.expression,
#                                   self.theta,self.temperature,self.frequency,self.nbins_d)
    def get_pol_vars(self):
        # Scattering
        S, Z = self.get_SZ()
            
        # Integration on DSD
        integ_S, integ_Z = self.get_SZ_integrated()
     
        pol = compute_pol_var(integ_S,integ_Z,self.frequency)  
        
        return pol
  
    def get_SZ_integrated(self):
        
        # Get amplitude and phase matrices
        S, Z = self.get_SZ()
        
#        outdated_psd = False
#        if self._integ_signature is None:
#            outdated_psd = True
#        elif self.integ_signature != (self.psd.dmin,self.psd.dmax,self.psd.func.expression,self.nbins_d):
#            outdated_psd = True
#        
#        if outdated_psd:
        list_D = np.linspace(self.psd.dmin,self.psd.dmax,self.nbins_d)
        
        N = self.psd(list_D,self.temperature)

        integ_S=np.trapz(S*N[:,None],dx=list_D[1]-list_D[0],axis=0)# multiply by integration step
        integ_Z=np.trapz(Z*N[:,None],dx=list_D[1]-list_D[0],axis=0) # multiply by integration step
        
        self._integ_S = integ_S
        self._integ_Z = integ_Z
        

        return self._integ_S, self._integ_Z

    def get_SZ(self):
        
        # Check if outdated
        outdated_scatter = False
        if self._scatter_signature is None:
            outdated_scatter = True
        elif self._scatter_signature != (self.psd.dmin,self.psd.dmax,self.aspect_ratio.expression,
                                   self.canting_angle_std,self.permittivity.expression,
                                   self.theta,self.temperature,self.frequency,self.nbins_d):
            outdated_scatter = True
        if outdated_scatter:
            wavelength=constants.C/(self.frequency*1E09)*1000 # in mm

            orientation_std = self.canting_angle_std
            
            list_D = np.linspace(self.psd.dmin,self.psd.dmax,self.nbins_d)
            
            # Aspect ratio and dielectric constant
            all_ar = self.aspect_ratio(list_D)
            all_m = self.permittivity(self.temperature,self.frequency,list_D)

            if not isinstance(all_m,np.ndarray):
                all_m = np.ones(list_D.shape) * all_m
                
            global SCATT
            # Initialize SCATTerer
            SCATT = create_scatterer(wavelength,orientation_std)
            
            global GEOM
            
            # Define GEOMetries
            GEOM = {} 
            GEOM['back'] = (90-self.theta, 180-(90-self.theta), 0., 180, 0.0,0.0) # Backward
            GEOM['forw'] =(90-self.theta, 90-self.theta, 0., 0.0, 0.0,0.0) # Forward
            

            # Initialize computing pool
            pool = mp.Pool(processes = mp.cpu_count())
            
            SZ=list(pool.map(get_SZ_by_D,zip(list_D,all_m,all_ar)))

            SZ = flatten_matrices(SZ)
    
            self._S = SZ[0]
            self._Z = SZ[1]
            
            self._set_scatter_signature()
            
        return self._S, self._Z
    

def get_SZ_by_D(params):
    SCATT.m =  params[1]

    SCATT.axis_ratio = params[2]
    SCATT.radius = params[0]/2.0
    
    SCATT.set_geometry(GEOM['back'])
    (S_back, Z_back) = SCATT.get_SZ_orient()
    
    SCATT.set_geometry(GEOM['forw'])
    (S_forw, Z_forw) = SCATT.get_SZ_orient()

    return (S_forw,Z_back)
    
    
if __name__ == '__main__':
    from parse import *

    config = {}
    config['permittivity'] = parse_permittivity([[0.1,0.9],[1+1j,2+2j]])

    config['psd'] = parse_psd(['ExponentialPSD',0.1,10])
#    
    config['aspect_ratio'] = parse_aspect_ratio(0.8)
    config['canting_angle_std'] = 16
    rain = Hydrometeor('rain',config,5.6,298,10,1024)
    print(rain)
#    
    tic()
    [S,Z] = rain.get_SZ()
    toc()

    pol = rain.get_pol_vars()
