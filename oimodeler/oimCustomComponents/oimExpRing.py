# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:27:15 2022

@author: Ame
"""

import numpy as np

from ..oimComponent import oimComponentRadialProfile,oimComponentFourier
from ..oimParam import oimParam
import astropy.units as u
from scipy.special import  j1


class oimExpRing(oimComponentFourier):
    """Uniform Disk component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the disk (in mas). The default is 0.
    """

    name = "Exponential Ring"
    shortname = "ExpRing"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"]     =  oimParam(base="d")
        self.params["fwhm"]  = oimParam(base="fwhm")
        self.params["dim"]   = oimParam(base="dim")   
        self.params["nfwhm"] = oimParam(name="nfwhm", value=8, mini=1,
            maxi=np.inf, free=False,description="Extension in number of fwhm")       
               
        self._eval(**kwargs)
        
    def _visFunction(self, ucoord, vcoord, rho, wl, t):
        
        
        dim = self.params["dim"].value
        d0 = self.params["d"](wl, t)* self.params["d"].unit.to(u.rad)
        fwhm = self.params["fwhm"](wl, t) * self.params["fwhm"].unit.to(u.rad)
        nfwhm = self.params["nfwhm"].value

        ftot=0
        res=rho*0
        for i in range(dim):
            
            xi    = d0/2 + fwhm*nfwhm/dim*i
            xmidi = xi + fwhm*nfwhm/dim/2
            yi    = np.exp(-0.692*np.divide(xi-d0/2,fwhm))
            yip   = np.exp(-0.692*np.divide(xi+fwhm*nfwhm/dim-d0/2,fwhm))

            xx = (np.pi*  2*xmidi * rho)

            if i!=dim-1:
                fi = np.pi*(yi-yip)*xmidi**2
            else:
                fi =  np.pi*(yi)*xmidi**2
            ftot+=fi
            #print(f"{i} \t {xmidi:.3e} {yi:.3e} {yip:.3e} {fi:.3e}")

           
            res+=np.nan_to_num(np.divide(2 * j1(xx), xx), nan=1)*fi
        #print(ftot)
        xx = (np.pi*  d0 * rho)
        f0 = np.pi*(d0/2)**2
        res-=(np.nan_to_num(np.divide(2 * j1(xx), xx), nan=1)*f0)
        ftot-=f0
        #print(f0)
        res/=ftot
        return res
