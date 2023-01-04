# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:30:21 2022

@author: Ame
"""
import numpy as np
from astropy import units as units
from oimodeler import oimParam, oimComponentImage
from fastRotator import fastRotator

class oimFastRotator(oimComponentImage):
    name="Fast Rotator"
    shortname="FRot"
    def __init__(self,**kwargs):
        super(). __init__(**kwargs)
        
        #Component parameters. Note that as it inherits from the oimComponentImage class it already has
        #x,y,f and dim as parameters
        self.params["incl"]=oimParam(name="incl",value=0,description="Inclination angle",unit=units.deg)
        self.params["rot"]=oimParam(name="rot",value=0,description="Rotation Rate",unit=units.one)     
        self.params["Tpole"]=oimParam(name="Tpole",value=20000,description="Polar Temperature",unit=units.K)
        self.params["dpole"]=oimParam(name="dplot",value=1,description="Polar diameter",unit=units.mas)
        self.params["beta"]=oimParam(name="beta",value=0.25,description="Gravity Darkening Exponent",unit=units.one)
       
        # constant value <=> static model
        self._t = np.array([0]) 
        
        #The component is chromatic. Here we set a fixed array of reference wavelengths. This can be
        #modified later as, in our case the model is recomputed at each call to the fastRotator function
        self._wl = np.linspace(0.5e-6,15e-6,num=10)  
        
        #Finally evalutating paramters as for all other components
        self._eval(**kwargs)
    
    def _internalImage(self):
        dim=self.params["dim"].value
        incl=self.params["incl"].value        
        rot=self.params["rot"].value
        Tpole=self.params["Tpole"].value
        dpole=self.params["dpole"].value
        beta=self.params["beta"].value  
       
        im=fastRotator(dim,1.5,incl,rot,Tpole,self._wl,beta=beta)
        
        #make a nt,nwl,dim,dim hcube (even if t and/or wl are not relevent)
        im=np.tile(np.moveaxis(im,-1,0)[None,:,:,:],(1,1,1,1))
       
        # computing the pixelSize based on the internal image size and the polar diameter
        self._pixSize=1.5*dpole/dim*units.mas.to(units.rad)
        
        return im
