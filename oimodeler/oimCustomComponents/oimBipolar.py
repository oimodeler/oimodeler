# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 11:52:59 2026

@author: ame
"""


import numpy as np
from oimodeler.oimComponent import oimComponentImage
from oimodeler.oimParam import oimParam, _standardParameters
from astropy import units as units


#%%
mas=units.rad.to(units.mas)
c=3e5

def bipolar(wl0,res,dim,incl,vpole,veq,alpha,beta,F0,time,dr,dist):

    vmax=np.max([vpole,veq])
    t_seconds = time*units.day.to(units.s)
    size_max_km = vmax * t_seconds *1.2
    pix_km = size_max_km / dim *2
    d_km=dist*units.pc.to(units.km)
    
    pix_rad = pix_km/d_km
    

    incl=incl*units.deg.to(units.rad)
    
    vrangle=np.arange(181)*units.deg.to(units.rad)
    vr=vpole+(veq-vpole)*np.sin(vrangle)**alpha
    
    xyz = (np.arange(dim)-dim//2)*pix_km
    ones=np.ones(dim)
    xm=xyz[:,np.newaxis,np.newaxis]*ones[np.newaxis,:,np.newaxis]*ones[np.newaxis,np.newaxis,:]
    ym=np.swapaxes(xm,0,1)
    zm=np.swapaxes(xm,0,2)

    r=np.sqrt(xm*xm+ym*ym+zm*zm)
    rxym=np.sqrt(xm*xm+ym*ym)
    
    thetam=np.arctan2(rxym,zm)
    #phim=np.arctan2(xm,ym)
    
    xp=xm
    yp=ym*np.cos(incl)-zm*np.sin(incl)
    zp=ym*np.sin(incl)+zm*np.cos(incl)
    
      
    r=np.sqrt(xp*xp+yp*yp+zp*zp)
    rxy=np.sqrt(xp*xp+yp*yp)
      
    theta=np.arctan2(rxy,zp)
    #phi=np.arctan2(xp,yp)
      
      
    vri=np.interp(theta,vrangle,vr)
       
    mapp=(np.abs(vri*t_seconds-r) < dr*pix_km)*1./r**beta
    mapp=np.nan_to_num(mapp,0)
      
    vproj=vri*np.cos(thetam)
    
    
    mapt=np.sum(mapp,axis=2)
    
    dimv=int(1.2*res*vmax/3e5)*2+1
    dimv = np.max([dimv,3])
    
    dv=1.5*2.*vpole/(dimv-1)
    v=-1.5*vpole+np.arange(dimv)*dv
    
    
    wl = wl0 * (1 + v/c)
      
    resv=c/res
    
    C2=np.sqrt(2*np.pi)/2.3548*res
    
    kinmap0=np.zeros((dim,dim,dimv))
    
    for iv in range(dimv):
        #mapi= (np.abs(vproj-v[iv]) <= dv/2.)*mapp
        mapi=np.exp(-(vproj-v[iv])**2./(2*(resv/2.3548)**2.))*mapp/C2
        #Vmapi=mapi*(zp>0)
        kinmap0[:,:,iv]=np.sum(mapi,axis=2)
    
    kinmap0=kinmap0/np.sum(kinmap0)
    
    kinmap = kinmap0*F0+mapt[:,:,np.newaxis]/np.sum(mapt)
    
    return kinmap,wl,pix_rad
#%%

class oimBipolar(oimComponentImage):
    name = "Bipolar nebula "
    shorname = "bipolar"
    elliptic = False
    
    def __init__(self, **kwargs):
        super(). __init__(**kwargs)
       
        self.params["dim"]=oimParam(**_standardParameters["dim"])
        self.params["wl0"]=oimParam(name="wl0",value=2.1656e-6,description="central wavelength of the line",unit=units.m,free=False)
        self.params["incl"]=oimParam(name="incl",value=45,description="inclination angle",unit=units.deg)
        self.params["res"]=oimParam(name="R",value=1.8e-10,description="spectral resolution",unit=units.m,free=False) 
        self.params["veq"]=oimParam(name="veq",value=500,description="Expansion velocity at the equator",unit=units.km/units.s)
        self.params["vpole"]=oimParam(name="vpole",value=2000,description="Expansion velocity at the vpole",unit=units.km/units.s)
        self.params["alpha"]=oimParam(name="alpha",value=1,description="exponent of the expansion velocity ",unit=units.one)
        self.params["beta"]=oimParam(name="beta",value=-3,description="exponent of radial intensity profile ",unit=units.one)
       
        self.params["time"]=oimParam(name="time",value=1,description="Time (in days) aftert the explosion",unit=units.day)
       
        self.params["dr"]=oimParam(name="dr",value=5,description="width of the nebula (in pixel)",unit=units.one,free=False)        
        self.params["dist"]=oimParam(name="dist",value=100,description="distance in parsec",unit=units.pc,free=False)
        self.params["EW"]=oimParam(name="EW",value=40,description="Equivalent width of the line",unit=units.AA)

    
        self._t = np.array([0]) # constant value <=> static model

        # will be set in the _internalImage function
        # self._wl = None
        self._eval(**kwargs)
       

    def _internalImage(self):
        
        wl0   = self.params["wl0"].value
        dim   = self.params["dim"].value
        res   = self.params["res"].value
        incl  = self.params["incl"].value
        veq   = self.params["veq"].value
        vpole = self.params["vpole"].value
        alpha = self.params["alpha"].value
        beta = self.params["beta"].value
        dist  = self.params["dist"].value
        EW    = self.params["EW"].value
        time     = self.params["time"].value
        dr    = self.params["dr"].value
     
        im,wl,pix_rad=bipolar(wl0,res,dim,incl,vpole,veq,alpha,beta,EW,time,dr,dist)
        im,wl,pix_rad=bipolar(wl0,res,dim,incl,vpole,veq,alpha,beta,EW,time,dr,dist)
        

        im = np.swapaxes(im,2,0)[np.newaxis,:,:,:]
        
        self._pixSize=pix_rad
        self._wl = wl

        return im
