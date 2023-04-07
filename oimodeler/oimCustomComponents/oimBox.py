import astropy.units as u
import numpy as np

from ..oimComponent import oimComponentFourier
from ..oimParam import oimParam


class oimBox(oimComponentFourier):
    name="2D Box"
    shortname = "BOX"
    def __init__(self,**kwargs): 
         super().__init__(**kwargs)
         self.params["dx"]=oimParam(name="dx", value=1,description="Size in x",unit=u.mas)
         self.params["dy"]=oimParam(name="dy", value=1,description="Size in y",unit=u.mas)    
         self._eval(**kwargs)
    
  
    def _visFunction(self,ucoord,vcoord,rho,wl,t):
        x=self.params["dx"](wl,t)*self.params["dx"].unit.to(u.rad)*ucoord
        y=self.params["dy"](wl,t)*self.params["dy"].unit.to(u.rad)*vcoord
        return np.sinc(x)*np.sinc(y) 
  
    def _imageFunction(self,xx,yy,wl,t):
        return ((np.abs(xx)<=self.params["dx"](wl,t)/2) &
                (np.abs(yy)<=self.params["dy"](wl,t)/2)).astype(float)
