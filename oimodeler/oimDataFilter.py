# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 13:21:23 2022

@author: Ame
"""

import numpy as np
from astropy.io import fits
import os
import oimodeler as oim
from oimodeler import oimParam,_standardParameters
#import numbers


###############################################################################
 
class oimDataFilter(object):    
    """
    Base class for data filter
    """
    
    name = "Generic Filter"
    shortname = "Gen filter"
    description = "This is the class from which all filters derived"
    
    def __init__(self,data=None,**kwargs): 
        
       self.params={}
       self.params["dest"]=oimParam(name="dest",value=None,description="Destination of the filter")
       
       self._eval(**kwargs)
      
    def _eval(self,**kwargs):
        for key, value in kwargs.items():
            if key in self.params.keys(): 
                   self.params[key].value=value
        
    def _filteringFunction(self):
        pass
    
    
###############################################################################


class oimDataFilteringStack(object):
    """
    class for data filter stack
    """
    def __init__(self,data=None,**kwargs): 
        pass
    
    def applyFilters(self,):
        pass