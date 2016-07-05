# -*- coding: utf-8 -*-
"""
Created on Mon May 30 13:36:03 2016

@author: wolfensb
"""

from pytmatrix import orientation
from pytmatrix.tmatrix import Scatterer
import numpy as np

import pyradsim.constants as constants
 

from pyradsim.utilities import Float


def create_scatterer(wavelength,orientation_std):
    scatt = Scatterer(radius = 1.0, wavelength = wavelength)
    scatt.or_pdf = orientation.gaussian_pdf(std=orientation_std)
    scatt.orient = orientation.orient_averaged_fixed
    return scatt
    
def flatten_matrices(list_matrices):
    arr_S=np.zeros((len(list_matrices),4),dtype='c8')
    arr_Z=np.zeros((len(list_matrices),16),dtype='float32')
    
    for i,mat_D in enumerate(list_matrices):
            # Backward
            # Z
            arr_Z[i,:]=mat_D[1].ravel()
            # Forward
            # S
            arr_S[i,:]=mat_D[0].ravel()
            
    return arr_S,arr_Z
       
       
def compute_pol_var(S,Z,F):
    # Transform SZ to two S matrices and one Z matrices
    wavelength=constants.C/(F*1E09)*1000 # in mm
    siz = S.shape
    
    if len(siz) == 1:
        Z = Z.reshape(4,4)
        S = S.reshape(2,2)
    elif len(siz) == 2:
        Z = Z.reshape(siz[0],4,4)
        Z = np.swapaxes(Z,0,2)
        S = S.reshape(siz[0],2,2)
        S = np.swapaxes(S,0,2)
    elif len(siz) == 3:
        Z = Z.reshape(siz[0],siz[1],4,4)
        Z = np.swapaxes(Z,0,2)
        Z = np.swapaxes(Z,1,3)
        S = S.reshape(siz[0],siz[1],2,2)       
        S = np.swapaxes(S,0,2)
        S = np.swapaxes(S,1,3)
        
    # Horizontal reflectivity
    radar_xsect_h=2*np.pi*(Z[0,0] - Z[0,1] - Z[1,0] + Z[1,1])
    zh=wavelength**4/(np.pi**5*constants.KW)*radar_xsect_h
    
    # Vertical reflectivity
    radar_xsect_v=2*np.pi*(Z[0,0] + Z[0,1] + Z[1,0] + Z[1,1])
    zv=wavelength**4/(np.pi**5*constants.KW)*radar_xsect_v 
    
    # Differential reflectivity
    zdr=radar_xsect_h/radar_xsect_v
    
    # Differential phase shift
    kdp=1e-3 * (180.0/np.pi) * wavelength * (S[1,1]-S[0,0]).real
    
    # Backscattering differential phase
    delta_hv = np.arctan2(Z[2,3] - Z[3,2], -Z[2,2] - Z[3,3])
    
    # Attenuation
    ext_xsect_h = 2 * wavelength * S[1,1].imag
    ext_xsect_v = 2 * wavelength * S[0,0].imag
    ah= 4.343e-3 * ext_xsect_h
    av= 4.343e-3 * ext_xsect_v
    
    # Copolar correlation coeff.
    a = (Z[2,2] + Z[3,3])**2 + (Z[3,2] - Z[2,3])**2
    b = (Z[0,0] - Z[0,1] - Z[1,0] + Z[1,1])
    c = (Z[0,0] + Z[0,1] + Z[1,0] + Z[1,1])
    rhohv = np.sqrt(a / (b*c))
    
    
    # Create output dictionary
    
    out = {}    
    
    zh = Float(zh)
    zh.units = 'mm^6*m^-3'
    out['Zh'] = zh
    ##
    zv = Float(zv)
    zv.units = 'mm^6*m^-3'
    out['Zv'] = zv  
    ##
    zdr = Float(zdr)
    zdr.units = '-'
    zdr.name = 'Specific attenuation at horizontal pol'
    out['Zdr'] = zdr      
    ##
    kdp = Float(kdp)
    kdp.units = 'deg*km^-1'
    kdp.name = 'Specific differential phase shift on propagation'
    out['Kdp'] = kdp
    ##
    delta_hv = Float(delta_hv)
    delta_hv.units = 'deg*km^-1'
    delta_hv.name = 'Phase shift on backscattering'
    out['Delta_hv'] = kdp
    ##
    ah = Float(ah)
    ah.units = 'dB/km'
    ah.name = 'Specific attenuation at horizontal polarization'
    out['Ah'] = ah
    ##
    av = Float(av)
    av.units = 'dB/km'
    av.name = 'Specific attenuation at vertical polarization'
    out['Av'] = av    
    ##
    rhohv = Float(rhohv)
    rhohv.units = '-'
    rhohv.name = 'Copolar correlation coefficient'
    out['Rhohv'] = rhohv       
    return out
    
