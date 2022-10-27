# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:27:15 2022

@author: Ame
"""

import numpy as np
import matplotlib.pyplot as plt
import oimodeler as oim
from astropy import units as units
from datetime import datetime

class oimSpiral(oim.oimComponentImage):
    
    name="Spiral component"
    shorname="Sp"
    
    #Set elliptic to True to use  elong keyword for a change of variable
    #to "flatten" objects
    elliptic=True
    
    def __init__(self,**kwargs):
        super(). __init__(**kwargs)
        
        #Component parameters. Note that as it inherits from the oimComponentImage class it already has
        #x,y,f and dim as parameters
        self.params["fwhm"]=oim.oimParam(**oim._standardParameters["fwhm"])
        self.params["P"]=oim.oimParam(name="P",
               value=1,description="Period in mas",unit=units.mas)
        self.params["width"]=oim.oimParam(name="width",
               value=0.01,description="Width as filling factor",unit=units.one)
       
        self._pixSize=0.05*units.mas.to(units.rad)
        
        self._t = np.array([0]) # constant value <=> static model
        self._wl = np.array([0])  # constant value <=> achromatic model
        
        #Finally evalutating paramters as for all other components
        self._eval(**kwargs)
    
    def _imageFunction(self,xx,yy,wl,t):
        
        # As xx and yy are transformed coordinates, r and phi takes into account 
        #the ellipticity and orientation using the pa and elong keywords
        r=np.sqrt(xx**2+yy**2)  
        phi=np.arctan2(yy,xx)
        
        p=self.params["P"](wl,t)
        sig=self.params["fwhm"](wl,t)/2.35
        w=self.params["width"](wl,t)
        
        im=1.+np.cos(-phi-2*np.pi*np.log(r/p+1))
        im=(im<2*w)*np.exp(-r**2/(2*sig**2))
        return im
    


#%%
dim=256

#Uncomment this to use FFTW backend instead of numpy backend for FFT computation
oim.oimOptions['FTBackend']=oim.FFTWBackend

ud=oim.oimUD(d=2)#,f=oim.oimInterpWl(wl=[0.5e-6,2e-6],value=[1,0]))
c=oimSpiral(dim=dim,fwhm=5,P=0.1,width=0.2,pa=30,elong=2,x=10)
m=oim.oimModel(c,ud)


#%%
#Computing and plotting visibilities for various baselines and walvelengths

nB=5000
nwl=5
wl=np.linspace(0.5e-6,2e-6,num=nwl)

B=np.linspace(0,100,num=nB//2)
Bx=np.append(B,B*0)
By=np.append(B*0,B)

Bx_arr=np.tile(Bx[None,:], (nwl, 1)).flatten()
By_arr=np.tile(By[None,:], (nwl,  1)).flatten()
wl_arr=np.tile(wl[:,None], (1, nB)).flatten()

spfx_arr=Bx_arr/wl_arr
spfy_arr=By_arr/wl_arr

n0=1
t0 = datetime.now()
for i in range(n0):
    vc=m.getComplexCoherentFlux(spfx_arr,spfy_arr,wl_arr)
dt=(datetime.now() - t0).total_seconds()*1000/n0 
print("Vcompl dim={} =>{:.1f}ms".format(dim,dt))

v=np.abs(vc.reshape(nwl,nB))
fig,ax=plt.subplots(1,1)

label=["East-West Baselines","North-South Baselines"]
for iwl in range(nwl):
    cwl=iwl/(nwl-1)
    ax.plot(B/wl[iwl],v[iwl,:nB//2]/v[iwl,0],label=label[0],color=plt.cm.viridis(cwl))
    ax.plot(B/wl[iwl],v[iwl,nB//2:]/v[iwl,0],alpha=0.5,label=label[1],color=plt.cm.viridis(cwl))  
    label=[None,None]
ax.legend()

#%%
fig,ax=plt.subplots(1,2)
m.showModel(256,0.1,swapAxes=True,fromFT=False,normPow=1,axe=ax[0],colorbar=False)
m.showModel(256,0.1,swapAxes=True,fromFT=True,normPow=1,axe=ax[1],colorbar=False)
ax[1].get_yaxis().set_visible(False)
ax[0].set_title("Direct Image")
ax[1].set_title("From FFT")

