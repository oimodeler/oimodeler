# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:39:36 2021

@author: Ame
"""

from scipy import integrate

import numpy as np
import oimodeler as oim
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import astropy.units as u
import time

mas2rad=u.mas.to(u.rad)
dim = 64

pixSize=0.1#mas
D=10 #mas


class myDisk1D(oim.oimComponentRadialProfile):
    name="A disk define by a radial profile"
    shortname = "myDisk"
    
    def __init__(self,**kwargs): 
         # First call the __init__ function from the parent class
         super().__init__(**kwargs)
         # Then define the component parameters  
         self.params["d"]=oim.oimParam(**(oim._standardParameters["d"]))
         # Finally call the _eval function that allow the parameters to be processed
         self._eval(**kwargs)
         
    
    def _radialProfileFunction(self):
        r=np.arange(0,dim)*pixSize*mas2rad
        I=(r<=self.params["d"].value*self.params["d"].unit.to(u.rad)/2).astype(float)
        return r,I




wl=0.65e-6
B=np.linspace(0.01,100,num=2000)
spf=B/wl


n=100

spf0=1.22*wl/D/mas2rad/wl


# Creating UD model for radial profile
c=myDisk1D(d=D)
modelHT=oim.oimModel([c])


# Creating UD model directly in Fourier
ud=oim.oimUD(d=D)
modelFT=oim.oimModel([ud])


start_time = time.time()
for i in range(n):
    vH=np.abs(modelHT.getComplexCoherentFlux(spf,spf*0))
print("Hankel {:.1f}ms".format((time.time() - start_time)/n*1000))    


start_time = time.time()
for i in range(n):
    vFT=np.abs(modelFT.getComplexCoherentFlux(spf,spf*0))
    vFT=vFT/vFT[0]
print("Analytical FT {:.1f}ms".format((time.time() - start_time)/n*1000))



fig,ax=plt.subplots(1,1,figsize=(10,7))

ax.plot(B/wl*mas2rad,vH,label="Radial Profile")
ax.plot(B/wl*mas2rad,vFT,label="FT from Formula")

ax.set_xlabel("B/$\\lambda$ (cycles/mas)")
ax.set_ylabel("Visbility")
ax.plot([1.22*wl/D/mas2rad/wl*mas2rad]*2,[0,1],color="k",linestyle="--")
plt.legend()

 

