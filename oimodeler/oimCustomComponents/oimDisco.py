# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:27:12 2026

@author: ame
"""

import numpy as np
import astropy.units as u
from .disco import model as discomodel
from ..oimComponent import oimComponentImage
from ..oimParam import oimParam
from ..oimUtils import oimAckWarning


###############################################################################

class oimDisco (oimComponentImage):
    name = "DISCO gaseous disk model"
    shortname = "Disco"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        text="The DISCO model was created by A. Domiciano de Souza "\
                    "based on the work of Rodrigo G. Vieira" 
        
        oimAckWarning(oimDisco,text)

        
        self.params["Rstar"]=oimParam(name="Rstar",value=10,
                                description="Stellar Radius",
                                unit=u.Rsun,free=True,mini=0,maxi=100)
        self.params["Mstar"]=oimParam(name="Mstar",value=5,
                                description="Stellar Mass",
                                unit=u.Msun,free=True,mini=0,maxi=100)        
        self.params["Teff"]=oimParam(name="Teff",value=20000,
                                description="Stellar Effective Temperature",
                                unit=u.K,free=True,mini=2000,maxi=100000)
        self.params["incl"]=oimParam(name="incl",value=45,
                                description="Inclination Angle",
                                unit=u.deg,free=True,mini=0,maxi=90)
        self.params["Rd"]=oimParam(name="Rd",value=50,
                                description="Disk outer radius",
                                unit=u.Rsun,free=True,mini=0,maxi=1000)
        self.params["Td0"]=oimParam(name="Td0",value=12000.,
                                description="Temperature at the disk basis",
                                unit=u.K,free=True,mini=0,maxi=100000)
        self.params["powTd"]=oimParam(name="powTd",value=0.75,
                                description="Power law coefficient for disk temperature",
                                unit=u.one,free=True,mini=0,maxi=1)
        self.params["rho0"]=oimParam(name="rho0",value=1.0e-8,
                                description="Density at disk basis",
                                unit=u.kg/u.m**3,free=True,mini=0,maxi=1)
        self.params["powRho"]=oimParam(name="powRho",value=-3.5,
                                description="Power law coefficient for disk density",
                                unit=u.one,free=True,mini=0,maxi=1)
        self.params["powHd"]=oimParam(name="powHd",value=1.5,
                                description="Power-law coefficient of disk flaring",
                                unit=u.one,free=False,mini=0,maxi=1)
        self.params["ionfrac"]=oimParam(name="ionfrac",value=1,
                                description="Ionization fraction",
                                unit=u.one,free=False,mini=0,maxi=1)
        self.params["dist"]=oimParam(name="dist",value=100,
                                description="Distance",
                                unit=u.pc,free=True,mini=0,maxi=1)

        self._t = np.array([0])
        #self._wl = np.array([0.1,0.4,0.7,1.0,1.5,2.0,2.5,3.,5.,8.,12.])*1e-6
        self._wl = np.logspace(-7,-4.8,num=30)
        self._eval(**kwargs)
        

        self._eval(**kwargs)
        
        self.outunit=u.Jy

    def _internalImage(self):
        star={"model":"BBsphere",
            "Req_sun"  : self.params['Rstar'].value*self.params['Rstar'].unit.to(u.Rsun),
            "Teq"      : self.params['Teff'].value, #K
            "Mstar_sun": self.params['Mstar'].value*self.params['Mstar'].unit.to(u.Msun),
            }
        
        maps={}
        im, xx, yy=discomodel(
           input_wlen= self._wl*u.m.to(u.micron), 
           star=star, 
           Rd=self.params["Rd"].value,
           Td0=self.params["Td0"].value, 
           powTd=self.params["powTd"].value, 
           rho0=self.params["rho0"].value,
           powrho=self.params['powRho'].value, 
           powHd=self.params['powHd'].value, 
           ionfrac=self.params['ionfrac'].value,
           incl=self.params["incl"].value,
           N=self.params["dim"].value,    
           maps=maps,)

        dist = self.params["dist"].value*self.params["dist"].unit.to(u.m)
        Rd   = self.params["Rd"].value*self.params["Rd"].unit.to(u.m)
        fovrad  = 2*Rd/dist
        # make a nt,nwl,dim,dim hcube (even if t and/or wl are not relevent)
        im = im[np.newaxis, :, :, :]
        
        # computing the pixelSize based on the internal image size and the polar diameter
        self._pixSize = fovrad/self.params["dim"].value
        
        # integrating flux on pixel <=> sum(im)==total flux
        A = self._pixSize * self._pixSize 
        im=im*A*(u.W/u.m**3).to(u.Unit(self.outunit),equivalencies=u.spectral_density(
                            self._wl[np.newaxis,:,np.newaxis, np.newaxis]*u.m))
        return im
    
    def getPixelSize(self,mas=False):
        dist = self.params["dist"].value*self.params["dist"].unit.to(u.m)
        Rd   = self.params["Rd"].value*self.params["Rd"].unit.to(u.m)
        fovrad  = 2*Rd/dist
        self._pixSize = fovrad/self.params["dim"].value
        fact=u.rad.to(u.mas)*float(mas)+float(not(mas))
        return self._pixSize*fact
        