# -*- coding: utf-8 -*-
"""
Created on Thu May 26 16:21:20 2016

@author: wolfensb
"""

import re
import functools
import collections
import numpy as np

class Float(float):
    def __new__(cls, arg = 0.0):
        return float.__new__(cls, arg)
        
class Lambda(object):
    def __init__(self,lambda_f, expression):
        self.expression = expression
        self.lambda_f = lambda_f
    def __call__(self,*inputs):
        return self.lambda_f(*inputs)
    def __str__(self):
        return self.expression
    def __mul__(self,scalar):
        return Lambda(lambda *args: self.__call__(*args)*scalar,str(scalar)+'*'+self.expression)
    def __rdiv__(self,scalar):
        return Lambda(lambda *args: self.__call__(*args)/scalar,str(scalar)+'/'+self.expression)        
    def __add__(self,scalar):
        return Lambda(lambda *args: self.__call__(*args)+scalar,str(scalar)+'+'+self.expression) 
    def __sub__(self,scalar):
        return Lambda(lambda *args: self.__call__(*args)-scalar,str(scalar)+'-'+self.expression)         
    
class InfoArray(np.ndarray):
    def __new__(cls, input_array,attr_names,attr_data):
        self = np.asarray(input_array).view(cls)
        if not isinstance(attr_names,list):
            attr_names = [attr_names]
        if not isinstance(attr_data,list):
            attr_data = [attr_data]
            
        if len(attr_names) != len(attr_data):
            print('Number of attributes names does not correspond with number of',
                  ' attributes vectors, ignoring all attributes...')
            return self
        
        dim_array = self.shape
        dim_attr = np.array([len(a) for a in attr_data])
   
        if not (dim_attr==dim_array).all():
            print('Length of attributes vectors does not correspond with dim of ',
                  ' data array, ignoring all attributes...')
            return self            
        
        self.attributes = {}
        for i,attr in enumerate(attr_names):
            self.attributes[attr] = np.array(attr_data[i])
        
        return self
        
    def __str__(self):
        strr = 'Data : \n' + super().__str__() +'\n'
        strr +='Attributes: \n'
        for attr in self.attributes.keys():
            strr += attr + ' : ' + str(self.attributes[attr]) +'\n'
        return strr
        
def create_defaults():
    content = '''
hydrometeors:
    permittivity: 8.5871375786139676+1.6977965395728176j
    canting_angle_std: 10
    psd: [COSMO_1mom_rain,0.0001] # Marshall-Palmer with R = 2 mm/hr
    psd_range: [0.1,20]
    aspect_ratio: 0.9
geometry:
    elevation_angle: 0
    azimuth: 0
    position: [1,1,1]
    size: [1,1,1]
atmosphere:
   T: 283
   P: 1015 # Not used currently...
radar:
   frequency: 5.6 # In GHZ
weight: 1
'''
    text_file = open("./.defaults.yml", "w")
    text_file.write(content)
    text_file.close()
    return
    
def get_from_dict(dataDict, mapList):
    return functools.reduce(lambda d, k: d[k], mapList, dataDict)

def set_from_dict(dic, keys, value):
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value
    
def remove_tags(lines,tags):
    new_lines = lines
    for i,l in enumerate(lines):
        for tag in tags:
            new_lines[i] =  re.sub(tag+'(.*)','',l)
    return new_lines

def remove_keys_with_none_values(data):
    new_data = {}
    for k, v in data.items():
        if isinstance(v, dict):
            v = remove_keys_with_none_values(v)
        if not v in (u'', None, {}):
            new_data[k] = v
    return new_data

def flatten_dic(d,parent_key='', sep='/'):
    flat = flatten(d,parent_key,sep)
    keys = [k.split(sep) for k in flat.keys()]
    values = [v for v in flat.values()]
    return keys, values
    
def flatten(d, parent_key='', sep='/'):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


if __name__ == '__main__':
    l = Lambda(lambda D,T: D**2*T**2, 'lambda D: D**2')
    v = l*2
    print(v(3,2),l(3,2))
