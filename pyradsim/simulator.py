# -*- coding: utf-8 -*-
"""
Created on Mon May 30 13:48:49 2016

@author: wolfensb
"""

import numpy as np
import copy
import itertools
import numbers
# These tags will be read separately, their arguments needs to be specified 
# in brackets after the yaml input for example

# azimuth: 150 %sens(100,1,2)


VALID_TAGS = ['%sens']

from pyradsim.parse import parse
from pyradsim.box import Box
from pyradsim.tmatrix import compute_pol_var
from pyradsim.utilities import set_from_dict, get_from_dict,InfoArray

class Simulator(object):    
    def __new__(self,box_file, config_file):
        config = parse('config_file',config_file)
        boxes, tags_params = parse('box_file',box_file,VALID_TAGS)
        if tags_params['%sens'] and config['sens_analysis'] != None:
            return SensSimulator(boxes,config,tags_params)
        else:
            return SingleSimulator(boxes,config)
            
class SingleSimulator(object):    
    def __init__(self,boxes, config):
        
        self.config = config
        self.boxes = []
        
        self._boxes_dic = boxes
        # Create list of all boxes
        for k in boxes.keys():
            self.boxes.append(Box(k,boxes[k],self.config))
            
    # Getter and setters for the dictionary that contains all params from the YAML
    @property
    def boxes_dic(self):
        return copy.deepcopy(self._boxes_dic)
    @boxes_dic.setter
    def boxes_dic(self,dic):
        self._boxes_dic = dic
        self.update()
        
    def update(self):
        # TODO proper check
        for k,box in enumerate(self._boxes_dic):
            self.boxes[k].update(self._boxes_dic[box],self.config)
            
    def get_pol_vars(self):
        boxes_pol_vars = {}
        for box in self.boxes:
            print('Simulating scattering of box: '+box.name)
            # Compute scattering matrices
            box_pol_vars = box.get_pol_vars()
            boxes_pol_vars[box.name] = box_pol_vars
        
        return boxes_pol_vars
        
            
    def get_integrated_pol_vars(self):
        ensemble_S = np.zeros((4,),dtype=complex) # Amplitutde matrix
        ensemble_Z = np.zeros((16,)) # Phase matrix
        
        # Check if radar frequency is the same in every box
        all_freq = [b.radar['frequency'] for b in self.boxes]  
        
        if len(np.unique(all_freq)) != 1:
            print('Integration of polarimetric variables over multiple boxes',
                  'is only possible when the radar frequency is the same in every box')
            print('Would you like to use the frequency of the first box ('+str(self.boxes[0].radar['frequency'])+' GHz) for all boxes?')
            userInput = input("[Y/N]")
            if userInput == 'Y':
                freq = self.boxes[0].radar['frequency']
            elif userInput == 'N':
                while True:
                  try:
                     userInput = float(input('Please specify the frequency to use for all boxes (in GHz):'))       
                  except ValueError:
                     print("Not a number!")
                     continue
                  else:
                     break 

                freq = float(userInput)
        else:
            freq = all_freq[0]
            
        for box in self.boxes:
            print('Simulating scattering of box: '+box.name)
            # Compute scattering matrices
            box_S, box_Z = box.get_ensemble_SZ()

            ensemble_S += box_S
            ensemble_Z += box_Z
            
        return compute_pol_var(ensemble_S,ensemble_Z,freq)
        
    def __str__(self):
        msg = ''
        for box in self.boxes:
            msg += box.__str__()
            
        return msg
        
class SensSimulator(SingleSimulator): # For sensitivity analysis
    def __init__(self,boxes, config, tags):
        self.sens_params = tags['%sens']
        super(SensSimulator,self).__init__(boxes,config)
        
    def get_pol_vars(self):
        print('Estimation of polarimetric variables with sensitivty analysis')
        print('-------------------------------------------------------------')
        type_sens_analysis = self.config['sens_analysis']
    
        initial_box_config = self.boxes_dic
        
        if type_sens_analysis == 'serial':
            dic_sens = {}
            for key,param in zip(self.sens_params[0],self.sens_params[1]):
                
                key_str = '/'.join(key)
                dic_sens[key_str] = {}
                print('Sensitivity analysis on param:'+str(key))
                working_copy = copy.deepcopy(initial_box_config)
                
                n_pts = int(param[2])
                range_var = np.linspace(param[0],param[1],n_pts)/100
                current_val = get_from_dict(initial_box_config,key)
                
                sens_values = []
                sens_param_name = key_str.split('/')[-1]
                if not isinstance(current_val, numbers.Number):
                    sens_param_name += ' factor'
                    
                for i,r in enumerate(range_var):
                    print('Step '+str(i+1)+'/'+str(len(range_var)))
                    if isinstance(current_val, numbers.Number):
                        sens_values.append(r * current_val)
                    else:
                        sens_values.append(r)
                        
                    # Modify working copy of boxes dict and assign it to current
                    # simulator
                    set_from_dict(working_copy ,key,current_val*r)
                    self.boxes_dic = working_copy
                                        
                    # Simulate polarimetric variables
                    pol_vars = super(SensSimulator,self).get_pol_vars()
                    
                    for k in pol_vars.keys():
                        if k not in dic_sens[key_str].keys():
                            dic_sens[key_str][k] = {}
                            for k2 in pol_vars[k].keys():
                                dic_sens[key_str][k][k2] = []

                    for k in pol_vars.keys():
                        for k2 in pol_vars[k].keys():
                            dic_sens[key_str][k][k2].append(pol_vars[k][k2])
                            
                # Convert lists to InfoArray
                for k in dic_sens[key_str].keys():
                    for k2 in dic_sens[key_str][k].keys():
                        dic_sens[key_str][k][k2] = InfoArray(dic_sens[key_str][k][k2],
                                        [sens_param_name],[sens_values])
                                           
        return dic_sens
        
    def get_integrated_pol_vars(self):
        print('Estimation of polarimetric variables with sensitivty analysis')
        print('-------------------------------------------------------------')
        type_sens_analysis = self.config['sens_analysis']
        print(type_sens_analysis)
        initial_box_config = self.boxes_dic
        
        if type_sens_analysis == 'serial':
            # Initialize output
            dic_sens = {}
            for key,param in zip(self.sens_params[0],self.sens_params[1]):
                
                key_str = '/'.join(key)
                dic_sens[key_str] = {}
                print('Sensitivity analysis on param:'+str(key))
                working_copy = copy.deepcopy(initial_box_config)
                
                n_pts = int(param[2])
                range_var = np.linspace(param[0],param[1],n_pts)/100
                current_val = get_from_dict(initial_box_config,key)
                
                sens_values = []
                sens_param_name = key_str.split('/')[-1]
                if not isinstance(current_val, numbers.Number):
                    sens_param_name += ' factor'
                    
                for i,r in enumerate(range_var):
                    if isinstance(current_val, numbers.Number):
                        sens_values.append(r * current_val)
                    else:
                        sens_values.append(r)
                        
                    print('Step '+str(i+1)+'/'+str(len(range_var)))
                    # Modify working copy of boxes dict and assign it to current
                    # simulator
                    set_from_dict(working_copy ,key,current_val*r)
                    self.boxes_dic = working_copy
                                        
                    # Simulate polarimetric variables
                    pol_vars = super(SensSimulator,self).get_integrated_pol_vars()
                    
                    for k in pol_vars.keys():
                        if k not in dic_sens[key_str].keys():
                            dic_sens[key_str][k] = []

                    for k in pol_vars.keys():
                            dic_sens[key_str][k].append(pol_vars[k])
                            
                # Convert lists to InfoArray
                for k in dic_sens[key_str].keys():
                    dic_sens[key_str][k] = InfoArray(dic_sens[key_str][k],[sens_param_name],
                                           [sens_values])
                                           
        elif type_sens_analysis == 'parallel': # ==> This will be very time-consuming
            dic_sens = {}
            # Generate all sequences of values
            dimensions = [p[2] for p in self.sens_params[1]]
            sens_values = [np.linspace(p[0],p[1],p[2]) for p in self.sens_params[1]]           
            combinations = list(itertools.product(*sens_values))
            
            for c in combinations:
                print(c)
#                for key,param in zip(self.sens_params[0],self.sens_params[1]):
#                    
#                    key_str = '/'.join(key)
#                    dic_sens[key_str] = {}
#                    print('Sensitivity analysis on param:'+str(key))
#                    working_copy = copy.deepcopy(initial_box_config)
#                    
#                    n_pts = int(param[2])
#                    range_var = np.linspace(param[0],param[1],n_pts)/100
#                    current_val = get_from_dict(initial_box_config,key)
#                    
#                    sens_values = []
                
#                    sens_param_name = key_str.split('/')[-1]
#                    if not isinstance(current_val, numbers.Number):
#                        sens_param_name += ' factor'
#                        
#                    for i,r in enumerate(range_var):
#                        if isinstance(current_val, numbers.Number):
#                            sens_values.append(r * current_val)
#                        else:
#                            sens_values.append(r)
#                            
#                        print('Step '+str(i+1)+'/'+str(len(range_var)))
#                        # Modify working copy of boxes dict and assign it to current
#                        # simulator
#                        set_from_dict(working_copy ,key,current_val*r)
#                        self.boxes_dic = working_copy
#                                            
#                        # Simulate polarimetric variables
#                        pol_vars = super(SensSimulator,self).get_integrated_pol_vars()
#                        
#                        for k in pol_vars.keys():
#                            if k not in dic_sens[key_str].keys():
#                                dic_sens[key_str][k] = []
#    
#                        for k in pol_vars.keys():
#                                dic_sens[key_str][k].append(pol_vars[k])
#                                
#                    # Convert lists to InfoArray
#                    for k in dic_sens[key_str].keys():
#                        dic_sens[key_str][k] = InfoArray(dic_sens[key_str][k],[sens_param_name],
#                                               [sens_values])
                
#        return dic_sens

if __name__ == '__main__':
    from pympler import muppy
    s = Simulator('box_file_example.yml','configuration_file_example.yml')
    out = s.get_integrated_pol_vars()
    all_objects = muppy.get_objects()
    from pympler import summary
    sum1 = summary.summarize(all_objects)
    summary.print_(sum1)
    
    