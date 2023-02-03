# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:27:15 2022

@author: Ame
"""
import numpy as np
from astropy import units as units

from ..oimComponent import oimComponentImage
from ..oimParam import oimParam, _standardParameters


class oimSpiral(oimComponentImage):
    
    name="Spiral component"
    shorname="Sp"
    
    #Set elliptic to True to use  elong keyword for a change of variable
    #to "flatten" objects
    elliptic=True
    
    def __init__(self,**kwargs):
        super(). __init__(**kwargs)
        
        #Component parameters. Note that as it inherits from the oimComponentImage class it already has
        #x,y,f and dim as parameters
        self.params["fwhm"]=oimParam(**_standardParameters["fwhm"])
        self.params["P"]=oimParam(name="P",
               value=1,description="Period in mas",unit=units.mas)
        self.params["width"]=oimParam(name="width",
               value=0.01,description="Width as filling factor",unit=units.one)
       
        self._pixSize=0.05*units.mas.to(units.rad)
        
        self._t = np.array([0]) # constant value <=> static model
        self._wl = np.array([0])  # constant value <=> achromatic model
        
        #Finally evalutating paramters as for all other components
        self._eval(**kwargs)
    
    def _imageFunction(self,xx,yy,wl,t):
        
        # As xx and yy are transformed coordinates, r and phi takes into account 
        #the ellipticity and orientation using the pa and elong keywords
        r=np.sqrt(xx**2+yy**2)  
        phi=np.arctan2(yy,xx)
        
        p=self.params["P"](wl,t)
        sig=self.params["fwhm"](wl,t)/2.35
        w=self.params["width"](wl,t)
        
        im=1 + np.cos(-phi-2*np.pi*np.log(r/p+1))
        im=(im<2*w)*np.exp(-r**2/(2*sig**2))
        return im
    


