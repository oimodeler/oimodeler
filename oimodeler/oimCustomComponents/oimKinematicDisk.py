# -*- coding: utf-8 -*-
"""
<<<<<<< Updated upstream
Created on Fri Oct 21 12:27:15 2022
=======
Created on Thu Apr 13 14:04:37 2023
>>>>>>> Stashed changes

@author: Ame
"""
import numpy as np
from astropy import units as units
from astropy import constants as cst
from ..oimComponent import oimComponentImage
from ..oimParam import oimParam
from ..oimOptions import standard_parameters

angToSize=lambda x: (units.AU,units.mas,lambda y:1000*y/x.to(units.parsec).value,lambda y:y/1000*x.to(units.parsec).value)


class oimKinematicDisk(oimComponentImage):
    name = "kinematic disk component"
    shorname = "kinDisk"
    elliptic = False

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)
        # The many parameters of the rotating disk model
        self.dim = oimParam(**standard_parameters["dim"])
        self.fov = oimParam(name="fov",value=30,description="Field of view in stellar diameters",unit=units.one,free=False)
        self.wl0 = oimParam(name="wl0",value=2.1656e-6,description="central wavelength of the line",unit=units.m,free=False)
        self.dwl = oimParam(name="dwl",value=0.9e-10,description="pixel size in wavelength ",unit=units.m,free=False)
        self.res = oimParam(name="R",value=1.8e-10,description="spectral resolution",unit=units.m,free=False)        
        self.nwl = oimParam(name="nwl",value=51,description="number of wavelengths",unit=units.one,free=False)        
        self.Rstar = oimParam(name="Rstar",value=5,description="Stellar Radius",unit=units.R_sun)
        self.dist = oimParam(name="dist",value=100,description="inclination angle",unit=units.pc)
        self.fwhmCont = oimParam(name="fwhmCont",value=5,description="FWHM of the disk in the continuum in DStar",unit=units.one)
        self.fwhmLine = oimParam(name="fwhmLine",value=2,description="FWHM of the disk in the line in DStar",unit=units.one)
        self.fluxDiskCont = oimParam(name="fluxDiskCont",value=0.5,description="Flux of the disk in the continuum",unit=units.one)
        self.incl = oimParam(name="incl",value=45,description="inclination angle",unit=units.deg)
        self.EW = oimParam(name="EW",value=40,description="Equivalent width of the line",unit=units.AA)
        self.vrot = oimParam(name="vrot",value=400,description="rotational velocity at the photosphere",unit=units.km/units.s)
        self.beta = oimParam(name="beta",value=-0.5,description="exponent of the rotational law",unit=units.one)
        self.v0 = oimParam(name="v0",value=0,description="expension velocity at the photosphere",unit=units.km/units.s)
        self.vinf = oimParam(name="vinf",value=0,description="expension velocity at the infinity ",unit=units.km/units.s)
        self.gamma = oimParam(name="vinf",value=0.86,description="exponent of the expansion velocity law",unit=units.one)

        self._t = np.array([0]) # constant value <=> static model

        # will be set in the _internalImage function
        # self._wl = None
        self._eval(**kwargs)

    def _internalImage(self):
        dim=self.dim.value
        fov=self.fov.value
        rstar=self.Rstar.value*self.Rstar.unit
        dist=self.dist.value*self.dist.unit
        wl0=self.wl0.value
        dwl=self.dwl.value        
        nwl=self.nwl.value        
        res=self.res.value
        incl = self.incl.value*self.incl.unit.to(units.rad)
        fwhmLine=self.fwhmLine.value
        fwhmCont=self.fwhmCont.value
        Fcont=self.fluxDiskCont.value
        EW=self.EW.value*self.EW.unit.to(units.m)
        beta = self.beta.value
        vrot = self.vrot.value
        vinf = self.vinf.value
        v0 = self.v0.value
        gamma = self.gamma.value

        Rstar2mas = rstar.to(units.mas,equivalencies=[angToSize(dist)]).value

        # NOTE: We define the pixelSize in rad from the fov and dim
        self._pixSize=fov*2*Rstar2mas/dim*units.mas.to(units.rad)

        # NOTE: intrisinct wl table is computed from the parameters wl0, nwl and dwl
        self._wl = np.linspace(wl0-dwl*(nwl//2),wl0+dwl*(nwl//2),num=nwl)

        #internal gird is used return 1D vectors for x, y  (=xy) 
        #Then 2D  xx and yy are computed
        _,_,x,y=self._getInternalGrid()
        xx,yy=np.meshgrid(x,y)
        xx=xx/Rstar2mas
        yy=yy/Rstar2mas      
        yp=yy/np.cos(incl)
        #r_incl is the radius projected with the inclination angle used of the disk
        r_incl=np.sqrt(xx*xx+yp*yp)
        # r is used for the central spherical star
        r=np.sqrt(xx*xx+yy*yy)

        rin = 1
        mask=r>rin
        yp=yy/np.cos(incl)

        #Velocity map computation
        phi=(np.arctan2(xx,yp))
        vphi=vrot*r_incl**beta
        vr=np.abs(v0+(vinf-v0)*(1.-1/r_incl))**gamma
        vmap = ((vphi*np.sin(phi)-vr*np.cos(phi))*np.sin(incl))*mask

        # 3 Intensity maps : star, disk in continuum and disk in the line
        sigmaLine=2*fwhmLine/2.3548
        sigmaCont=2*fwhmCont/2.3548
        mapEnvL=np.exp((-r_incl**2)/(2*sigmaLine**2.))*(r_incl>1)
        mapEnvL=mapEnvL/np.sum(mapEnvL)
        mapEnvC = np.exp((-r_incl**2)/(2*sigmaCont**2.))
        mapEnvC=mapEnvC/np.sum(mapEnvC)
        mapStar = (r <= 1)
        mapStar = mapStar/np.sum(mapStar)

        #Continuum map with Star + disk
        mapC= (1-Fcont)*mapStar+Fcont*mapEnvC

        #conversion wl to velocity
        c=cst.c.to(units.Unit("km/s")).value
        v=(self._wl-wl0)/wl0*c
        resv=res/wl0*c

        #Some normalization cst
        C0=2*(resv/2.3548)**2.
        C1=mapEnvL*EW
        C2=np.sqrt(2*np.pi)/2.3548*res
        C=C1/C2

        mapL=np.ndarray([1,nwl,dim,dim])

        # NOTE: The loop is surprisingly faster than the matrix
        # It computes narrow band images through the emission line
        for iwl in range(nwl):
            mapL[0, iwl, :, :] = mapC + np.exp(-(vmap-v[iwl])**2./C0)*C
        return mapL
