See the wiki for info about installation and use


Structure of the code
-------


METHODS
-------
aspect_ratio_models.py : contains definition of aspect ratio models taken from 
                         litterature
                         
permittivity_models.py : contains definition of permittivity models taken from 
                         litterature

utilities.py           : a set of useful methods used throughout the code

parse.py               : a set of parsing methods to read the user specified box
                         and configuration files
                         
tmatrix.py             : a set of methods to compute scattering and phase matrix
                         of individual particles
                         
CONSTANTS
---------
constants.py           : definition of physical and numerical constants

constants_cosmo.py     : definition of COSMO specific constants

CLASSES
-------
hydrometeor.py         : definition of the hydrometeor class which defines a 
                         hydrometeor which all its parameters (PSD,perm.,etc)
            
box.py                 : definition of the box class which defines a 
                         box which all its parameters (hydrometeors, temp, etc)
                         
psd.py                 : definition of the psd class

simulator.py           : the MAIN class, defines the simulator class, which requires
                         to give a configuration and boxes file, does all the job
                         by recursively calling the specific methods of the box and
                         hydrometeor classes