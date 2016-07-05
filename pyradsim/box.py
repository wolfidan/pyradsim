# -*- coding: utf-8 -*-
"""
Created on Mon May 30 13:44:43 2016

@author: wolfensb
"""

import json
import numpy as np

from pyradsim.hydrometeor import Hydrometeor
from pyradsim.tmatrix import compute_pol_var

class Box(object):
    def __init__(self,name,box,config):
        self.name = name
        self.radar = box['radar']
        self.atmosphere = box['atmosphere']
        self.geometry = box['geometry']
        self.weight = box['weight']
        self.hydrometeors = []
        self.config = config
        
        for k in box['hydrometeors']:
            self.hydrometeors.append(Hydrometeor(k,box['hydrometeors'][k],
                 self.radar['frequency'],self.atmosphere['T'],
                 self.geometry['elevation_angle'],self.config['nbins_d']))   
    
    def update(self,box,config):
        self.radar = box['radar']
        self.atmosphere = box['atmosphere']
        self.geometry = box['geometry']
        self.weight = box['weight']
        self.config = config
        for h in self.hydrometeors:
            h.update(box['hydrometeors'][h.name],
                 self.radar['frequency'],self.atmosphere['T'],
                 self.geometry['elevation_angle'],self.config['nbins_d'])   
  
    def get_ensemble_SZ(self):
        ensemble_S = np.zeros((4,),dtype=complex) # Amplitutde matrix
        ensemble_Z = np.zeros((16,)) # Phase matrix
        
        for h in self.hydrometeors:      
            print('Simulating scattering of hydrometeor '+h.name)
            # Compute scattering matrices
            integ_S, integ_Z = h.get_SZ_integrated()

            ensemble_S += integ_S
            ensemble_Z += integ_Z
        return ensemble_S, ensemble_Z
        
    def get_pol_vars(self):
        ensemble_S, ensemble_Z = self.get_ensemble_SZ()
        pol = compute_pol_var(ensemble_S,ensemble_Z,self.radar['frequency'])       
        return pol
            
    def __str__(self):
        msg = '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
        msg += 'BOX: '+self.name+'\n'
        msg += '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n'
        msg += 'Hydrometeors\n'
        msg += '-------------------\n'
        for h in self.hydrometeors:
            msg += h.__str__()
        msg += 'Radar\n'
        msg += '-------------------\n'
        msg += json.dumps(self.radar,indent=3) +'\n\n'
        msg += 'Geometry\n'
        msg += '-------------------\n'
        msg += json.dumps(self.geometry,indent=3) +'\n\n'
        msg += 'Atmosphere\n'
        msg += '-------------------\n'
        msg += json.dumps(self.atmosphere,indent=3) +'\n\n'
        msg += 'Weight\n'
        msg += '-------------------\n'
        msg += json.dumps(self.weight,indent=3) +'\n'
        
        return msg
        
        
