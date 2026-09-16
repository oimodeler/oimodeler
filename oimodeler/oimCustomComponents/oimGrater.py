# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 10:04:23 2026

@author: ame
"""



import numpy as np
import astropy.units as u

from pyGrater.stargrains import Grain, Star
from pyGrater.SED import SED
from pyGrater.image import Image
from pyGrater.density import two_power_law
from pyGrater.size_distributions import power_law_distribution
from pyGrater.phase_functions import isotropic
import pyGrater


from oimodeler.oimComponent import oimComponentImage
from oimodeler.oimParam import oimParam
from oimodeler.oimUtils import oimAckWarning


num_size_bins = 400         # Number of size bins for integrat


class oimGrater (oimComponentImage):
    name = "Grater disk model"
    shortname = "grater"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        text="The Grater model was created by J.C. Augereau & P. Priollet"
        oimAckWarning(oimGrater,text)
        
        self.grain = None
        self.star = None
        self.normalizeImage=False
        
        self.params["r0"]=oimParam(name="r0",value=10,
                                description="Disk reference radius",
                                unit=u.au,free=True,mini=0,maxi=np.inf)
        self.params["h0"]=oimParam(name="h0",value=0.1,
                                description="scale height at the reference radius",
                                unit=u.au,free=True,mini=0,maxi=np.inf)
        self.params["alphaIn"]=oimParam(name="alphaIn",value=10,
                                description="Inner power-law index",
                                unit=u.one,free=True,mini=-np.inf,maxi=np.inf)
        self.params["alphaOut"]=oimParam(name="alphaOut",value=-4,
                                description="Outer power-law index",
                                unit=u.one,free=True,mini=-np.inf,maxi=np.inf)
        self.params["gamma"]=oimParam(name="gamma",value=2,
                                description="Vertical profile exponent",
                                unit=u.one,free=True,mini=-np.inf,maxi=np.inf)
        self.params["beta"]=oimParam(name="beta",value=2,
                                description="Scale height flaring exponent",
                                unit=u.one,free=True,mini=-np.inf,maxi=np.inf)
        self.params["incl"]=oimParam(name="incl",value=45,
                                description="Inclination Angle",
                                unit=u.deg,free=True,mini=0,maxi=90)
        self.params["omega"]=oimParam(name="omega",value=60,
                                description="Longitude of ascending node",
                                unit=u.deg,free=True,mini=0,maxi=90)
        self.params["fov"]=oimParam(name="fov",value=40,
                                description="field of view",
                                unit=u.au,free=True,mini=0,maxi=np.inf)
        self.params["Mtot"]=oimParam(name="Mtot",value= 1e-8,
                                description="disk total mass",
                                unit=u.Msun,free=True,mini=0,maxi=100)
        self.params["dist"]=oimParam(name="dist",value=108,
                                description="Distance",
                                unit=u.pc,free=False,mini=0,maxi=np.inf)
        self.params["grainSizeMin"]=oimParam(name="grainSizeMin",value=0.01,
                                description="Minimum grain size",
                                unit=u.um,free=False,mini=0,maxi=np.inf)        
        self.params["grainSizeMax"]=oimParam(name="grainSizeMax",value=3000,
                                description="Maximum grain size",
                                unit=u.um,free=False,mini=0,maxi=np.inf)        
        self.params["grainPow"]=oimParam(name="grainPow",value=3.5,
                                description="Power-law index for size distribution",
                                unit=u.one,free=False,mini=0,maxi=np.inf)           

        self._wl = np.array([0.55, 0.9, 1.65, 2.15, 3.5, 8, 10, 12])*1e-6
        self._t= np.array([0])
        self._eval(**kwargs)
    
    def setDataPath(self,path):
        pyGrater.set_data_path(path)
    
    def setGrain(self,**kwargs):
        self.grain = Grain(**kwargs)
        self._createImgObj()
        
        
    def setStar(self,**kwargs):
        self.star = Star(**kwargs)
        self._createImgObj()
        
        
    def _createImgObj(self):
        if (self.grain!=None) & (self.star!=None):
            self.imgObj = Image(self.grain, self.star, two_power_law, 
                            power_law_distribution,isotropic,self._wl*1e6)
            
            
    def getSED(self,wl):
        sedObj = SED(self.grain, self.star, two_power_law, 
                     power_law_distribution, wl*1e6, N_distances=400)
        
        image_params = {
            'r0'       : self.params["r0"].value,
            'alphain'  : self.params["alphaIn"].value, 
            'alphaout' : self.params["alphaOut"].value,
            'h0'       : self.params["h0"].value,
            'beta'     : self.params["beta"].value,
            'gamma'    : self.params["gamma"].value,
            'itilt'    : self.params["incl"].value,
            'PA'       : self.params["pa"].value, 
            'omega'    : self.params["omega"].value, 
            'nx'       : self.params["dim"].value,
            'ny'       : self.params["dim"].value,
            'FOV_AU'   : self.params["fov"].value,
            'M_tot'    : self.params["Mtot"].value,
            'a_min'    : self.params["grainSizeMin"].value*self.params["grainSizeMin"].unit.to(u.m),
            'a_max'    : self.params["grainSizeMax"].value*self.params["grainSizeMax"].unit.to(u.m),
            'kappa'    : self.params["grainPow"].value,
            'N_sizes_integral': 400
        }
        
        sed_therm, sed_sca = sedObj.get_SED(keep_separate_fluxes=True, **image_params)
        return sed_therm + sed_sca
    
    
    def setInternalWavelengths(self,wl):
        self._wl = wl
        self._createImgObj()
        
    def _internalImage(self):

        image_params = {
            'r0'       : self.params["r0"].value,
            'alphain'  : self.params["alphaIn"].value, 
            'alphaout' : self.params["alphaOut"].value,
            'h0'       : self.params["h0"].value,
            'beta'     : self.params["beta"].value,
            'gamma'    : self.params["gamma"].value,
            'itilt'    : self.params["incl"].value,
            'PA'       : self.params["pa"].value, 
            'omega'    : self.params["omega"].value, 
            'nx'       : self.params["dim"].value,
            'ny'       : self.params["dim"].value,
            'FOV_AU'   : self.params["fov"].value,
            'M_tot'    : self.params["Mtot"].value,
            'a_min'    : self.params["grainSizeMin"].value*self.params["grainSizeMin"].unit.to(u.m),
            'a_max'    : self.params["grainSizeMax"].value*self.params["grainSizeMax"].unit.to(u.m),
            'kappa'    : self.params["grainPow"].value,
            'N_sizes_integral': 400
        }
                
        im0 = self.imgObj.get_image(**image_params)
        im = im0[np.newaxis, :, :, :]
        self.getPixelSize()
        
        return im

    def getPixelSize(self,mas=False):
         dist = self.params["dist"].value*self.params["dist"].unit.to(u.m)
         dim   = self.params["dim"].value
         fov  = self.params["fov"].value*self.params["fov"].unit.to(u.m)
         
         fovrad  = fov/dist
         
         self._pixSize = fovrad/dim
         fact=u.rad.to(u.mas)*float(mas)+float(not(mas))
         return self._pixSize*fact