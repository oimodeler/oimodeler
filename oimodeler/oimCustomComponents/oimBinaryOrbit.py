# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 14:38:05 2026

@author: ame
"""

import numpy as np

from scipy import special
from ..oimComponent import oimComponentFourier
from ..oimBasicFourierComponents import oimPt
from ..oimParam import oimParam,_standardParameters
import astropy.units as u

def Rotate_2D(Coord,theta):
    rotMat = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    Coord1 = np.matmul(rotMat,Coord);
    return Coord1

def RotateOrbit(XY_Coord, i, omega, OMEGA):

     Coord3 = XY_Coord
     Coord1 = XY_Coord
     Coord1 = Rotate_2D(XY_Coord,np.pi/2-omega)
     Coord2 = Coord1
     Coord2[0] = Coord1[0] * np.cos(i);
     Coord3 = Rotate_2D(Coord2,-OMEGA);
     return Coord3 

def keplerCoordinates(phase, eccentricity, majorAxis):
    rho = majorAxis * (1 - eccentricity * np.cos(phase));
    
    #rho = majorAxis *(1-eccentricity**2) / (1 + eccentricity * np.cos(phase));
    
    theta = 2*np.arctan(np.tan(phase/2.0) * np.sqrt((1+eccentricity)/(1-eccentricity)));
    X = rho*np.cos(theta);
    Y = rho*np.sin(theta);   
    XY = np.array([X,Y])  
    return XY

def binaryCoordinates(time, eccentricity, majorAxis, period):
    omega = 2*np.pi/period
    order = 20
    phi = omega * time  
    for k in range(1,order):
        phi = phi + 2 * special.jn(float(k),k * eccentricity)/float(k) * np.sin(float(k) * omega * time)
    #phi = (phi-np.pi) % (2*np.pi) + np.pi
    #phi = phi % (2*np.pi) 
    XY = keplerCoordinates(phi,eccentricity,majorAxis)    
    return XY

def binaryProjectedCoordinates(time, t0, eccentricity, majorAxis, period, i, omega, OMEGA):  
    CoordOrbit = binaryCoordinates(time-t0,eccentricity,majorAxis,period)
    CoordP2 = RotateOrbit(CoordOrbit,i,omega,OMEGA)
    return CoordP2


def orbit2Mass(period,majorAxis):
    G=6.673e-11
    omega = 2.0*np.pi/period;
    totalMass = omega**2 * majorAxis**3 / G;
    return totalMass

def radialVelocity(time,T,t0,e,omega,Ka,V0):
    ninterpol=10000
    EEi=np.linspace(0,2*np.pi,ninterpol)
    ti=T/(2*np.pi)*(EEi-e*np.sin(EEi))
    tr=(t0-time) % T
    EE=np.interp(tr,ti,EEi)
    nu=2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(EE/2.))
    return V0-Ka*(np.cos(omega+nu)+e*np.cos(omega))



class oimBinaryOrbit(oimComponentFourier):
    
    _params = {}
    @property
    def params(self):
        
        res=self._params.copy()
        
        for key in self.primary.params:
            if not(key in ["x","y"]):
                res[f"primary_{key}"]=self.primary.params[key]
        for key in self.secondary.params:
            if not(key in ["x","y"]):
                res[f"secondary_{key}"]=self.secondary.params[key]                
        
        return res
    
    def __init__(self, **kwargs):
        #super().__init__(**kwargs)
        
        self._params["x"] = oimParam(**_standardParameters["x"])
        self._params["y"] = oimParam(**_standardParameters["y"])
        self._params["f"] = oimParam(**_standardParameters["f"])
        self._params["e"]=oimParam(name="e",value=0,description="Eccentricity",unit=u.one,free=True,mini=0,maxi=1)
        self._params["a"]=oimParam(name="a",value=0,description="Semi-major axis",unit=u.mas,free=True,mini=0)        
        self._params["T"]=oimParam(name="T",value=30,description="Period",unit=u.day,free=True)        
        self._params["T0"]=oimParam(name="T0",value=30,description="Epoch of periastron passage (MJD)",unit=u.day,free=True,mini=0)
        self._params["i"]=oimParam(name="i",value=30,description="Inclination angle",unit=u.deg,free=True,mini=-90,maxi=90)
        self._params["O"]=oimParam(name="O",value=30,description="Longitude of the ascending node ",unit=u.deg,free=True,mini=-180,maxi=180)
        self._params["o"]=oimParam(name="o",value=30,description="Argument of periapsis",unit=u.deg,free=True,mini=-180,maxi=180)     
        
        self.primary = oimPt()
        self.secondary = oimPt()
        self._eval(**kwargs)
        
    def _getSeparation(self,t,mas=False):
        e = self.params["e"].value
        a = self.params["a"].value*self.params["a"].unit.to(u.rad)
        T = self.params["T"].value*self.params["T"].unit.to(u.day)
        T0 = self.params["T0"].value
        i = self.params["i"].value*self.params["i"].unit.to(u.rad)
        O = self.params["O"].value*self.params["O"].unit.to(u.rad)
        o = self.params["o"].value*self.params["o"].unit.to(u.rad)
        x,y = binaryProjectedCoordinates(t,T0,e,a,T,i,o,O)
        if mas:
            x=x*u.rad.to(u.mas)
            y=y*u.rad.to(u.mas)            
        return x,y

    def _visFunction(self, ucoord, vcoord, rho, wl, t):
        x,y = self._getSeparation(t) #in rad
        
        transFact = np.exp(-2 * 1j * np.pi * (ucoord * x + vcoord * y))
 
        p = self.primary.params["f"](wl,t)*\
            self.primary._visFunction(ucoord, vcoord, rho, wl, t) 
        s = self.secondary.params["f"](wl,t)*\
            self.secondary._visFunction(ucoord, vcoord, rho, wl, t)*transFact

        return p+s
#