# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 16:09:56 2025

@author: ame
"""

import oimodeler as oim
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from scipy.special import j0,j1
from oimodeler.oimComponent import oimComponentFourier
from oimodeler.oimParam import oimParam, _standardParameters
import time
from pathlib import Path
from pprint import pprint
import matplotlib.colors as colors

class oimCustomSkewedRing(oimComponentFourier):
    """Infinitesimal Ring component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    """
    name = " Ring with Custom intensity profile and skwedeness"
    shortname = "CSRing"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["n"] = oimParam(name="n",value=10,mini=1,maxi=100,free=False,description="number of rings")
        self.params["fwhm"]=oimParam(**_standardParameters["fwhm"])
        self.params["din"] = oimParam(**_standardParameters["din"])
        self.params["din2"] = oimParam(**_standardParameters["din"])
        self.params["din2"].name="din2"
        self.params["dout"] = oimParam(**_standardParameters["dout"])
        self.params["dout"].free=False
        self.params["skwIn"] = oimParam(**_standardParameters["skw"])
        self.params["skwIn"].name="skwIn"
        self.params["skwIn2"] = oimParam(**_standardParameters["skw"])
        self.params["skwIn2"].name="skwIn2"
        self.params["skwOut"] = oimParam(**_standardParameters["skw"])
        self.params["skwOut"].name="skwOut"           
        self.params["skwPa"] = oimParam(**_standardParameters["skwPa"])       
        self._eval(**kwargs)

    

        
    def _visFunction(self, xp, yp, rho, wl, t):
        
        phi = (self.params["skwPa"](wl, t)-self.params["pa"](wl, t)) * \
            self.params["skwPa"].unit.to(u.rad) + np.arctan2(yp, xp)
            
            
        n = self.params["n"].value
        fwhm = self.params["fwhm"].value#(wl,t)
        din = self.params["din"].value#(wl,t)
        din2 = self.params["din2"].value#(wl,t)
        dout = self.params["dout"].value#(wl,t)
       
        
        rs=np.linspace(din/2,dout/2,n)
        Ir=np.exp(-(rs-din2/2)/fwhm)
        Ir=np.clip(Ir,0,1)
    
        
        ftot=0
        for i in range(n):
            ftot+=Ir[i]*rs[i]
        
        res = 1j*xp*0
        
        skwIn = self.params["skwIn"].value#(wl,t)
        skwIn2 = self.params["skwIn2"].value#(wl,t)
        skwOut = self.params["skwOut"].value#(wl,t)
        
        skw=np.interp(rs,[din/2,din2/2,dout/2],[skwIn,skwIn2,skwOut])


        for i in range(n):
            
            xx = np.pi*rs[i]*2*self.params["din"].unit.to(u.rad)*rho
            #xx = np.pi*self.params["d"](wl, t)*self.params["d"].unit.to(u.rad)*rho
            res+= np.nan_to_num(j0(xx)-1j*np.sin(phi)*j1(xx)*skw[i])*Ir[i]*rs[i]/ftot

            #res+=j0(xx)*Ir[i]*rs[i]/ftot
        return res
#%%
ring = oimCustomSkewedRing(n=100,din=10,din2=25,dout=50,fwhm=10,elong=1.5,pa=90,skwIn=0,skwIn2=0.5,skwOut=0)
m= oim.oimModel(ring)
dim = 256
fov=60
pix=fov/dim
m.showModel(dim,pix,fromFT=True,normPow=1)