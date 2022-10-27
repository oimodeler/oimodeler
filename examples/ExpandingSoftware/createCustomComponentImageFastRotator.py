# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:30:21 2022

@author: Ame
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from astropy import units as units
import oimodeler as oim
from fastRotator import fastRotator

import os

path = os.path.dirname(oim.__file__)


class oimFastRotator(oim.oimComponentImage):
    name="Fast Rotator"
    shortname="FRot"
    def __init__(self,**kwargs):
        super(). __init__(**kwargs)
        
        #Component parameters. Note that as it inherits from the oimComponentImage class it already has
        #x,y,f and dim as parameters
        self.params["incl"]=oim.oimParam(name="incl",value=0,description="Inclination angle",unit=units.deg)
        self.params["rot"]=oim.oimParam(name="rot",value=0,description="Rotation Rate",unit=units.one)     
        self.params["Tpole"]=oim.oimParam(name="Tpole",value=20000,description="Polar Temperature",unit=units.K)
        self.params["dpole"]=oim.oimParam(name="dplot",value=1,description="Polar diameter",unit=units.mas)
        self.params["beta"]=oim.oimParam(name="beta",value=0.25,description="Gravity Darkening Exponent",unit=units.one)
       
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

#%% Creating a model

c=oimFastRotator(dpole=5,dim=128,incl=-70,rot=0.99,Tpole=20000,beta=0.25)
m=oim.oimModel(c)

#%% Plotting the model image

m.showModel(512,0.025,wl=[1e-6,10e-6 ],legend=True,normalize=True,
        savefig=os.path.join(path,os.pardir,"images","customCompImageFastRotator.png")))

#%% Computing and plotting visibilities for various baselines and walvelengths

nB=1000
nwl=20
wl=np.linspace(1e-6,2e-6,num=nwl)

B=np.linspace(0,100,num=nB//2)
Bx=np.append(B,B*0) # 1st half of B array are baseline in the East-West orientation
By=np.append(B*0,B) # 2nd half are baseline in the North-South orientation

Bx_arr=np.tile(Bx[None,:], (nwl, 1)).flatten()
By_arr=np.tile(By[None,:], (nwl,  1)).flatten()
wl_arr=np.tile(wl[:,None], (1, nB)).flatten()

spfx_arr=Bx_arr/wl_arr
spfy_arr=By_arr/wl_arr

vc=m.getComplexCoherentFlux(spfx_arr,spfy_arr,wl_arr)
v=np.abs(vc.reshape(nwl,nB))

fig,ax=plt.subplots(1,2,figsize=(15,5))
titles=["East-West Baselines","North-South Baselines"]
for iwl in range(nwl):
    cwl=iwl/(nwl-1)
    ax[0].plot(B/wl[iwl]/units.rad.to(units.mas),v[iwl,:nB//2],
            color=plt.cm.plasma(cwl))
    ax[1].plot(B/wl[iwl]/units.rad.to(units.mas),v[iwl,nB//2:],
           color=plt.cm.plasma(cwl))  

for i in range(2):
    ax[i].set_title(titles[i])
    ax[i].set_xlabel("B/$\lambda$ (cycles/rad)")
ax[0].set_ylabel("Visibility")    
ax[1].get_yaxis().set_visible(False)   

norm = colors.Normalize(vmin=np.min(wl)*1e6,vmax=np.max(wl)*1e6)
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax,label="$\\lambda$ ($\\mu$m)")
fig.savefig(os.path.join(path,os.pardir,"images","customCompImageFastRotatorVis.png"))


#%% Modifying the fast rotator parameters and  creating a new model with UD component
c.params['f'].value=0.9
c.params['pa'].value=45
ud=oim.oimUD(d=1,f=0.1,y=10)
m2=oim.oimModel(c,ud)

#%% Plotting the new model image

m2.showModel(512,0.06,wl=[1e-6,10e-6],legend=True, normalize=True,normPow=0.5,
   savefig=os.path.join(path,os.pardir,"images","customCompImageFastRotator2.png"))

#%% Computing and plotting visibilities for various baselines and walvelengths

vc=m2.getComplexCoherentFlux(spfx_arr,spfy_arr,wl_arr)
v=np.abs(vc.reshape(nwl,nB))

fig,ax=plt.subplots(1,2,figsize=(15,5))
titles=["East-West Baselines","North-South Baselines"]
for iwl in range(nwl):
    cwl=iwl/(nwl-1)
    ax[0].plot(B/wl[iwl]/units.rad.to(units.mas),v[iwl,:nB//2],
            color=plt.cm.plasma(cwl))
    ax[1].plot(B/wl[iwl]/units.rad.to(units.mas),v[iwl,nB//2:],
           color=plt.cm.plasma(cwl))  

for i in range(2):
    ax[i].set_title(titles[i])
    ax[i].set_xlabel("B/$\lambda$ (cycles/rad)")
ax[0].set_ylabel("Visibility")    
ax[1].get_yaxis().set_visible(False)   

norm = colors.Normalize(vmin=np.min(wl)*1e6,vmax=np.max(wl)*1e6)
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax,label="$\\lambda$ ($\\mu$m)")
fig.savefig(os.path.join(path,os.pardir,"images","customCompImageFastRotatorVis.png"))

