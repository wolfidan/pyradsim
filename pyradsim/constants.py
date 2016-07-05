# -*- coding: utf-8 -*-
"""
Created on Wed May 25 12:51:32 2016

@author: wolfensb
"""

import numpy as np
##############################################################################
# Parameters
##############################################################################

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
# Numerical parameters
EPS = np.finfo(np.float32).eps

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
# Physical parameters
C=299792458 # m/s
RHO_W=1000./(1000**3) # kg/mm3 at 10°C
RHO_I=916./(1000**3)  # kg/mm3 at 0°C
M_AIR = 1 # dielectric constant of air
KW = 0.93 # Dielectric factor for water at weather radar frequencies