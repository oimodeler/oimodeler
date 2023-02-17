# -*- coding: utf-8 -*-
"""
Custom model components added by the community
"""



"""
If you add python files with custom components in the oimCustomComponents 
directory add the import here so that all component will be available through
 the customComponent module
"""

#TODO: find a way to automatically add all class that derives from oimComponent 
#from all filesin the oimCustomComponents directory

#TODO: import crashes randomly with ModuleNotFoundError: No module named 'oimFastRotator' 

from .oimBox import oimBox
from .oimFastRotator import oimFastRotator
from .oimSpiral import oimSpiral
from .oimRadialPowerLaw import oimRadialPowerLaw, oimAsymmetricRadialPowerLaw
