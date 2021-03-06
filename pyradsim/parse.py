
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 15:56:53 2016

@author: wolfensb
"""

import collections
import numbers
import functools
import numpy as np
import re
import yaml
import os
import warnings

from pyradsim.utilities import remove_tags, flatten_dic,remove_keys_with_none_values
from pyradsim.permittivity_models import PERMITTITIVITY_MODELS
from pyradsim.aspect_ratio_models import AR_MODELS
from pyradsim.psd import PSD, BinnedPSD, create_PSD
from pyradsim.utilities import Lambda, create_defaults

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))

yaml.add_representer(collections.OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)


MATH_FUNCTIONS = ['cos','sin','exp','log','log10','tan','abs']
        
def parse(input_type,*inputs):
    if input_type == 'psd':
        return parse_psd(*inputs)
    elif input_type == 'aspect_ratio':
        return parse_aspect_ratio(*inputs)
    elif input_type == 'vol_fraction':
        return parse_vol_fraction(*inputs)
    elif input_type == 'permittivity':
        return parse_permittivity(*inputs)
    elif input_type == 'config_file':
        return parse_config_file(*inputs)
    elif input_type == 'box_file':
        return parse_box_file(*inputs)      
  
  
def parse_box_file(filename, tags = ''):
    print('Parsing box file...\n')
    
    ###########################################################################
    # Check if tags are present in the file
    tags_params = None
    if tags != '':
        tags_params = {}
        # Start by getting the tags from the file
        for tag in tags:
            tags_params[tag] = parse_tag(filename, tag)
        
        # Now remove all tags from the file and write a copy to /tmp/
        lines = open(filename).readlines()
        new_lines = remove_tags(lines,tags)

        tmp_filename = '/tmp/box_file.yml'
        with open(tmp_filename,'w') as tmp_file:
            tmp_file.writelines(new_lines)
        filename = tmp_filename
    ###########################################################################
        
    # Read default file
    defaults_file = './.defaults.yml'

    try:
        with open(defaults_file, 'r') as defaults_file:
                defaults = yaml.safe_load(defaults_file)
    except:
        warnings.warn('Could not find or read default file, creating a new one...')
        create_defaults()
        with open(defaults_file, 'r') as defaults_file:
            defaults = yaml.safe_load(defaults_file)
    
    # Read specified box file
    with open(filename, 'r') as filename:
        boxes = yaml.load(filename)


    try:       
        # Parse the content of the file
        parsed_boxes = boxes
        
        for box_name in boxes.keys(): # Loop on boxes  
            print("Box: '"+box_name+"'")
            print('===\n')
            for cat in defaults.keys(): # Loop on category
                if cat == 'weight':
                    if 'weight' not in boxes[box_name].keys():
                        print("Could not read weight in box '",box_name,"'")
                        print("Assigning default value: ",defaults['weight'])
                        print('---\n')
                        parsed_boxes[box_name]['weight'] = defaults['weight']
                # Check for hydrometeor parameters
                elif cat == 'hydrometeors':
                    for hydro_name in boxes[box_name][cat].keys(): # Loop on hydrometeors
                        for hydro_param in defaults['hydrometeors'].keys(): # Loop on hydrometeors parameters
                            if hydro_param in ['psd','permittivity','aspect_ratio']:
                                try:
                                    parsed_boxes[box_name][cat][hydro_name][hydro_param] = parse(hydro_param,boxes[box_name][cat][hydro_name][hydro_param])
                                except:
                                    print("Could not read '"+hydro_param+"' for hydrometeor ",hydro_name," in box '",box_name,"'")
                                    print("Assigning default value: ",defaults[cat][hydro_param])
                                    print('---\n')
                                    parsed_boxes[box_name][cat][hydro_name][hydro_param] = parse(hydro_param,defaults[cat][hydro_param]) 
                            else:
                                # The other parameters don't require specific parsing
                                if hydro_param not in parsed_boxes[box_name][cat][hydro_name].keys():
                                    print("Could not read '"+hydro_param+"' for hydrometeor ",hydro_name," in box '",box_name,"'")
                                    print("Assigning default value: ",defaults[cat][hydro_param])
                                    print('---\n')
                                    parsed_boxes[box_name][cat][hydro_name][hydro_param] = defaults[cat][hydro_param]

                        # Finally assign diameter range to psd
                        parsed_boxes[box_name][cat][hydro_name]['psd'].dmin = parsed_boxes[box_name][cat][hydro_name]['psd_range'][0]
                        parsed_boxes[box_name][cat][hydro_name]['psd'].dmax = parsed_boxes[box_name][cat][hydro_name]['psd_range'][1]

                else:                   
                    # Check for all other parameters
                    for param in defaults[cat].keys():
                        if param not in parsed_boxes[box_name][cat].keys():
                            print("Could not read '"+cat+'/'+param+"' in box '",box_name,"'")
                            print("Assigning default value: ",defaults[cat][param])
                            print('---\n')
                            parsed_boxes[box_name][cat][param] = defaults[cat][param]    
    except:
        raise
        raise RuntimeError('Could not parse content of box file, please make sure you respect the required synthax, check the wiki and the given example') 
    print('Box file was succesfully read!\n')
    print('====================')

    return parsed_boxes, tags_params


def parse_config_file(filename):
    print('Parsing configuration file...\n')
    if not os.path.exists(filename):
        raise RuntimeError('File could not be found!')
    try:
        with open(filename, 'r') as filename:
            config = yaml.load(filename)
    except:
        raise RuntimeError('Could not parse content of config file, please make sure you respect the required synthax, check the wiki and the given example') 
    print('Configuration file was succesfully read!\n')
    print('====================')
    return config

        
def parse_math_expression(string):
    for func in MATH_FUNCTIONS:
        if func in string:
            string = string.replace(func,'np.'+func)
    return string
    
def parse_psd(psd):
    flag = 0
    try:
        if isinstance(psd,str):
            if 'D' in psd or 'T' in psd:
                flag = 1
                psd_corr = parse_math_expression(psd)
                out = PSD('custom',Lambda(eval('lambda D, T:'+psd_corr),psd_corr))
                # I guess we should try to run the function here...   
        elif isinstance(psd,list):
            if isinstance(psd[0],list) and isinstance(psd[1],list):
                # Binned psd
                flag = 1
                try:
                    out = BinnedPSD(psd[0],psd[1],psd[2])
                except:
                    out = BinnedPSD(psd[0],psd[1])
            elif isinstance(psd[0],str):
                try:
                    flag = 1
                    out = create_PSD(psd[0],*psd[1:])
                    print(out)
                except:
                    flag == 0
                    raise
    except:
        flag = 0

    # Final output        
    if flag == 1:
        return out
    else:
        print('Specified PSD model is invalid, enter either a valid PSD model or', 
              ' a functional form that depends on D and/or T: e.g. 0.8*D**1.08')
        return
        
def parse_aspect_ratio(aspect_ratio):
    flag = 1
    try:
        if isinstance(aspect_ratio,str):
            if aspect_ratio in AR_MODELS.keys():
                out = AR_MODELS[aspect_ratio]
            elif 'D' in aspect_ratio:
                aspect_ratio_corr = parse_math_expression(aspect_ratio)
                expression = 'lambda D:'+ aspect_ratio_corr
                out = Lambda(eval(expression),expression)
                # I guess we should try to run the function here...
            else:
                flag = 0            
        elif isinstance(aspect_ratio,numbers.Number):
            expression = 'lambda D:(D>0)*'+str(aspect_ratio)
            out = Lambda(eval(expression),expression)
        else:
            flag = 0
    except:
        flag = 0

    # Final output        
    if flag == 1:
        return out
    else:
        print('Specified aspect ratio is invalid, enter either a constant number', 
              ' or a functional form that depends on D: e.g. 0.8*D**1.08')
        return

def parse_vol_fraction(vol_fraction):
    flag = 1
    try:
        if isinstance(vol_fraction,str):
            if 'D' in vol_fraction or 'T' in vol_fraction:
                vol_fraction_corr = parse_math_expression(vol_fraction)
                expression = 'lambda D, T:'+vol_fraction_corr
                out = Lambda(eval(expression),expression)
                # I guess we should try to run the function here...
            else:
                flag = 0            
        elif isinstance(vol_fraction,numbers.Number):
            expression = 'lambda D, T:(D>0)*'+str(vol_fraction)
            out = Lambda(eval(expression),expression)
        else:
            flag = 0
    except:
        flag = 0

    # Final output        
    if flag == 1:
        return out
    else:
        print('Specified volumetric fraction is invalid, enter either a constant number', 
              ' or a functional form that depends on D and/or T:',
              ' e.g. 0.8*(D**1.08)*(T**0.002)')
        return
        
def parse_permittivity(diel_cst):
    try: # Try converting the input to a complex number
        diel_cst = complex(diel_cst)
    except:
        pass
    
    flag = 1
    try:
        # Case 1 : input is string
        ############################
        if isinstance(diel_cst,str):

            if diel_cst == 'water':
                expression = "lambda T, F, D: PERMITTITIVITY_MODELS['water'](T,F)"
                out = Lambda(eval(expression),"Dielectric model for liq. water taken from 'A model for the complex permittivity of"+
                "water at frequencies below 1 THz' (Liebe,1991)")
            elif diel_cst == 'ice':
                expression = "lambda T, F, D: PERMITTITIVITY_MODELS['ice'](T,F)"
                out = Lambda(eval(expression),"Dielectric model for liq. water taken from 'A model for the complex permittivity of"+
                "ice at frequencies below 1 THz' (Hufford,1991)")
            elif 'D' in diel_cst or 'T' in diel_cst or 'F' in diel_cst:
                diel_cst_checked = parse_math_expression(diel_cst)
                expression = 'lambda T, F, D:'+str(diel_cst_checked)
                out = Lambda(eval(expression),expression)
            else:
                flag = 0                
             # I guess we should try to run the function here...
                
        # Case 2 : input is number
        ############################        
        elif isinstance(diel_cst,numbers.Number):
            expression = 'lambda T, F, D:(D>0)*'+str(diel_cst)
            out = Lambda(eval(expression),expression)
   
       # Case 3 : input is list
       ############################
        elif isinstance(diel_cst,list):        
            # When input is a list, it is assumed to be a list of two lists.
            # The first list contains the dielectric constants (any valid format)
            # and the second the vol. fraction relations (can be constants or 
            # functions of D and T)
            
            if len(diel_cst) == 2 and len(diel_cst[0]) == len(diel_cst[1]):
                # First we read all the volumetric fractions
                vol_frac_func = []
                diel_func = []
                for component in zip(diel_cst[0],diel_cst[1]):
                    
                    vol_frac_func.append(parse_vol_fraction(component[0]))
                    diel_func.append(parse_permittivity(component[1]))

                def func(vol_frac_func,T, F,D):
                    frac = [v(D,T) for v in vol_frac_func]
                    m = [m(T,F,D) for m in diel_func]
                    return PERMITTITIVITY_MODELS['mixture'](frac, m)
                expression = "Mixture model taken from 'Radar Backscattering by Inhomogeneous Precipitation Particles'"
                expression+=" (Bohren and Battan, 1980), with "
                expression += 'vol fractions = '+' / '.join([v.expression for v in vol_frac_func])
                expression += ' and dielectric constants = '+' / '.join([d.expression for d in diel_func])
                out = Lambda(functools.partial(func, vol_frac_func),expression)      
                
            else:
                flag = 0
        else:
            flag = 0
    except:
        raise
        flag = 0
    
    # Final output
    if flag == 1:
        return out
    else:
        print('Specified volumetric fraction is invalid, enter either a constant number', 
              ' or a functional form that depends on D and/or T and/or F:',
              ' e.g. 0.8*(D**1.08)*(T**0.002)*(F**0.002)')
        return


def parse_tag(filename,tag):
    f = open(filename)
    lines = f.readlines()
    lines = [x for x in lines if not re.match(r'\n', x.strip(' '))]
    
    levels = [len(l)-len(l.lstrip()) for l in lines]
    ranks = sorted(list(set(levels)))
    levels = [ranks.index(lev) for lev in levels]
    
    levels.insert(0,-1)
    lines.insert(0,':')

    out = recursive_parse(lines,levels,tag)
    
    # Remove empty values from dic
    out = remove_keys_with_none_values(out)
    if out == {}:
        out = None
    else:
        return flatten_dic(out)
        
def recursive_parse(lines,levels,tag):
    current_level = levels[0]
    
    if len(lines) == 1:
        data = lines[0][lines[0].index(':')+1:].strip()
        if tag in data:
            try:
                data_after_tag = data[data.index(tag):] 
                idx_startag = data_after_tag.index('(')
                idx_endtag = data_after_tag.index(')')
                return parse_string(data_after_tag[idx_startag+1:idx_endtag])
            except:
                print('Invalid tag input at line :'+lines[0])
                print('Tag must be followed by left bracket (',
                    'and closed with right bracket')

    else:
        odict = collections.OrderedDict()
        start = -1
        end = -1
        for i,l in enumerate(lines):
            if levels[i] == current_level + 1 and start == -1:
                start = i
                key = lines[i][0:lines[i].index(':')].strip()
            elif (levels[i] == current_level + 1 and start != -1):
                end = i
                odict[key] = recursive_parse(lines[start:end],levels[start:end],tag)
                start = end
                key = lines[i][0:lines[i].index(':')].strip()
            if i == len(lines) - 1:
                end = i + 1
                odict[key] = recursive_parse(lines[start:end],levels[start:end],tag)
        
        return odict

def parse_string(string):
    out = string
    try:
        out = string.strip('()[]')
        out = out.split(',')
        out = [float(s) for s in out]
    except:
        print('Could not extract numbers from string')
    return out
    
if __name__ == '__main__':

    filee = 'boxes.yml'
#    a = yaml_tag_parser(filee,'%sens')
    box = parse_box_file(filee,['%sens'])
#    filename = '/tmp/box_file.yml'
#    ooo = open(filename)
#    boxes = yaml.load(ooo)
#    config = parse('config_file','configuration.yml')
#    box = parse('box_file','boxes_test.yml')
