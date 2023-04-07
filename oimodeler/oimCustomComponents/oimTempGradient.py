# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:39:36 2021

@author: Ama
"""

import numpy as np
from oimodeler import oimParam, oimComponentRadialProfile,_standardParameters
from astropy import units as units
from astropy.modeling import models
from scipy.constants import c,h,k
from astorpy.constants import M_sun
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import astropy.units as u
import matplotlib.cm as cm
from scipy.interpolate import interp1d
import os

mas2rad=u.mas.to(u.rad)
dim=256

class oimTempGradient(oimComponentRadialProfile):
    name="A ring defined by a radial temperature profile in r^q and a radial dust surface density profile in r^p."
    shortname = "TempGradient"
    
    elliptic=True
    
    def __init__(self,**kwargs): 
         # First call the __init__ function from the parent class
         super().__init__(**kwargs)
         # Then define the component parameters  
         self.params["rin"]=oimParam(name="rin", value=0,description="Inner radius of the disk [au]")
         self.params["rout"]=oimParam(name="rout", value=0,description="Outer radius of the disk [au]")
         self.params["Tin"]=oimParam(name="Tin", value=0,description="Inner radius temperature [K]")
         self.params["Mdust"]=oimParam(name="Mdust", value=0,description="Mass of the dusty disk [M_sun]")
         self.params["q"]=oimParam(name="q", value=0,description="Power-law exponent for the temperature profile")
         self.params["p"]=oimParam(name="p", value=0,description="Power-law exponent for the dust surface density profile")
         self.params["kappa_abs"]=oimParam(name="kappa_abs", value=0,description="Dust mass absorption coefficient [cm2.g-1]")
         self.params["dist"]=oimParam(name="dist", value=0,description="Distance of the star [pc]")
         self._dim=dim
         
         self._t = np.array([0,1e+99]) # constant value <=> static model
         self._wl = None#np.array([0.5,1])*1e-6
         
         # Finally call the _eval function that allow the parameters to be processed
         self._eval(**kwargs)
         
    
    def _radialProfileFunction(self,r,wl,t):
        #r=np.arange(0,dim)*pixSize*mas2rad
        dist=self.params["dist"](wl,t)
        rin=self.params["rin"](wl,t)
        rout=self.params["rout"](wl,t)
        q=self.params["q"](wl,t)
        p=self.params["p"](wl,t)
        Tin=self.params["Tin"](wl,t)
        Mdust=self.params["Mdust"](wl,t)*M_sun*1e+3
        kappa_abs=self.params["kappa_abs"](wl,t)
        
        rin_mas=1e+3*rin/dist
        rout_mas=1e+3*rout/dist
        #r=np.linspace(din/2.,dout/2.,dim)
        
        #Temperature radial profile
        T_profile=Tin*np.divide(r,rin_mas)**q
        
        #Surface density radial profile
        rin_cm=rin*u.au.to(units.cm)
        rout_cm=rout*u.au.to(units.cm)
        if p==-2:
            sigma_in=Mdust/(2.*np.pi*np.log(rout_cm/rin_cm)*rin_cm**2)
        else:
            f=((rout_cm/rin_cm)**(p+2)-1)/(p+2)
            sigma_in=Mdust/(2.*np.pi*f*rin_cm**2)
        sigma_profile=sigma_in*np.divide(r,rin_mas)**p
        
        # BB emission radialprofile
        BB_profile=(2*h*c/wl**5)*np.divide(1.,np.exp((h*c)/(wl*k*T_profile))-1)
        
        # Emissivity factor
        emiss_factor=1-np.exp(-sigma_profile*kappa_abs)
        
        # Final intensity profile
        I=np.nan_to_num(np.logical_and(r>rin_mas,r<rout_mas).astype(int)*BB_profile*emiss_factor,nan=0)
        #I=np.nan_to_num((r>rin).astype(float)*np.exp(-0.692*np.divide(r-rin,rout)),nan=0)
    #I=np.nan_to_num((r>r0).astype(float)*np.exp(-0.692*np.divide(r-r0,fwhm)),nan=0)
        #I=I_in*(r/rin)**p
        #I=(r<=self.params["d"].value*self.params["d"].unit.to(units.rad)/2).astype(float)
        return I
    
    @property
    def _r(self):
        if False:
            rout=self.params["rout"](self._wl,self._t)
        else:
            rout=self.params["rout"](1e99)
        rmax=4*rout
        r=np.linspace(0,1, self._dim)*1e+3*rmax/self.params["dist"](self._wl,self._t)
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
