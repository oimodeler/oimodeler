# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:39:36 2021

@author: Ama
"""

import numpy as np
from oimodeler import oimParam, oimComponentRadialProfile,_standardParameters
from astropy import units as units
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import astropy.units as u
import matplotlib.cm as cm
from scipy.interpolate import interp1d
import os

mas2rad=u.mas.to(u.rad)
dim = 256

class oimRadialRing(oimComponentRadialProfile):
    name="A ring defined by a radial intensity profile in r^p"
    shortname = "RadialRing"
    
    elliptic=True
    
    def __init__(self,**kwargs): 
         # First call the __init__ function from the parent class
         super().__init__(**kwargs)
         # Then define the component parameters  
         self.params["din"]=oimParam(_standardParameters["din"]) 
         self.params["dout"]=oimParam(_standardParameters["dout"])
         self.params["p"]=oimParam(name="p", value=0,description="Power-law exponent")
         
         self._dim=dim
         
         self._t = np.array([0,1e+99]) # constant value <=> static model
         self._wl = None#np.array([0.5,1])*1e-6
         
         # Finally call the _eval function that allow the parameters to be processed
         self._eval(**kwargs)
         
    
    def _radialProfileFunction(self,r,wl,t):
        #r=np.arange(0,dim)*pixSize*mas2rad
        din=self.params["din"](wl,t)
        dout=self.params["dout"](wl,t)
        p=self.params["p"](wl,t)
        #r=np.linspace(din/2.,dout/2.,dim)
        rin=din/2.
        rout=dout/2.
        I_in=1.
        I=np.nan_to_num(np.logical_and(r>rin,r<rout).astype(int)*I_in*np.divide(r,rin)**p,nan=0)
        #I=np.nan_to_num((r>rin).astype(float)*np.exp(-0.692*np.divide(r-rin,rout)),nan=0)
    #I=np.nan_to_num((r>r0).astype(float)*np.exp(-0.692*np.divide(r-r0,fwhm)),nan=0)
        #I=I_in*(r/rin)**p
        #I=(r<=self.params["d"].value*self.params["d"].unit.to(units.rad)/2).astype(float)
        return I
    
    @property
    def _r(self):
        if False:
            dout=self.params["dout"](self._wl,self._t)
        else:
            dout=self.params["dout"](1e99)
        rmax=4*dout/2.
        r=np.linspace(0,1, self._dim)*rmax
        return r
        #din=self.params["din"](self._wl,self._t)* mas2rad
        #dout=self.params["dout"](self._wl,self._t)* mas2rad
        #r=np.linspace(din/2.,dout/2.,self._dim)
        #return r

    @_r.setter
    def _r(self, r):
        pass
    
    # def _imageFunction(self,xx,yy,wl,t):
    #     r2=(xx**2+yy**2)
    #     I_in=1.
    #     Ir=I_in*(np.sqrt(r2)/self.params["din"](wl,t)/2)**self.params["p"](wl,t)
    #     return  ((r2<=(self.params["dout"](wl,t)/2)**2) & 
    #              (r2>=(self.params["din"](wl,t)/2)**2)).astype(float)*Ir


############################################################################### 
