# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:27:12 2026

@author: ame
"""

import numpy as np
import astropy.units as u
from ..oimComponent import oimComponentImage
from ..oimParam import oimParam,_standardParameters

###############################################################################
class oimInnerRim(oimComponentImage):
    name = "Inner rim"
    shortname = "InRim"

    elliptic = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        # essential to be able to have chromatic flux on this component <=> Temperature of inner rim
        self.interpFillValue = None # extrapolation of image in x,y,t and wl

        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["h"] = oimParam(
            name="h", value=1, description="height", unit=u.mas)
        self.params["incl"] =  oimParam(
            name="incl", value=0, description="inclination angle", unit=u.deg)

        self._t = np.array([0])
        self._wl = np.array([0])

        self._eval(**kwargs)

    def _internalImage(self):
        
        dim  = self.params["dim"]()
        incl = self.params["incl"]()
        h0   = self.params["h"]()
        d   = self.params["d"]()
        h    = h0*np.tan(np.deg2rad(incl))/d
        
        xy = np.linspace(-1,1,dim)
        xx,yy = np.meshgrid(xy,xy)        
        yy2 = yy/np.cos(np.deg2rad(incl))


        r1 = xx**2+(yy2+h)**2
        r2 = xx**2+(yy2)**2
        im = (r2<0.25) & (1-(r1<0.25)) # 0.7071 => sqrt(0.5)
        
        im = im[np.newaxis, np.newaxis, :, :]
        im = im/np.sum(im)
        self.getPixelSize()
        
        return im

    def getPixelSize(self,mas=False):
        dim  = self.params["dim"]()
        d    = self.params["d"]()*self.params["d"].unit.to(u.rad)
        self._pixSize = 2*d/dim
        fact=u.rad.to(u.mas)*float(mas)+float(not(mas))

        return  self._pixSize*fact
    
#%%