# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:30:21 2022

@author: Ame
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import oimodeler as oim
from fastRotator import fastRotator
from astropy import units as units
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


#%%
dim=128

#Uncomment this to use FFTW backend instead of numpy backend for FFT computation
oim.oimOptions['FTBackend']=oim.FFTWBackend

#Padding factor, by default 8x de size of non-zero element of the image can be 
#reduced to redcue FFT computation time
oim.oimOptions['FTpaddingFactor']=8

c=oimFastRotator(dpole=5,dim=dim,incl=70,rot=0.99,Tpole=20000,x=0,y=0,pa=180,f=1,beta=0.25)
m=oim.oimModel(c)

#%%
#Computing and plotting visibilities for various baselines and walvelengths

nB=5000
nwl=20
wl=np.linspace(0.5e-6,2e-6,num=nwl)

B=np.linspace(1,100,num=nB//2)
Bx=np.append(B,B*0) # 1st half of B array are baseline in the East-West orientation
By=np.append(B*0,B) # 2nd half are baseline in the North-South orientation

Bx_arr=np.tile(Bx[None,:], (nwl, 1)).flatten()
By_arr=np.tile(By[None,:], (nwl,  1)).flatten()
wl_arr=np.tile(wl[:,None], (1, nB)).flatten()

spfx_arr=Bx_arr/wl_arr
spfy_arr=By_arr/wl_arr



vc=m.getComplexCoherentFlux(spfx_arr,spfy_arr,wl_arr)
v=np.abs(vc.reshape(nwl,nB))
fig,ax=plt.subplots(1,1)
label=["East-West Baselines","North-South Baselines"]
for iwl in range(nwl):
    cwl=iwl/(nwl-1)
    ax.plot(B/wl[iwl]/units.rad.to(units.mas),v[iwl,:nB//2],
            label=label[0],color=plt.cm.plasma(cwl))
    ax.plot(B/wl[iwl]/units.rad.to(units.mas),v[iwl,nB//2:],
            alpha=0.1,label=label[1],color=plt.cm.plasma(cwl))  
    label=[None,None]
ax.legend()
ax.set_xlabel("B/$\lambda$ (cycles/mas)" )
ax.set_ylabel("Visibility" )

norm = colors.Normalize(vmin=np.min(wl)*1e6,vmax=np.max(wl)*1e6)
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax,label="$\\lambda$ ($\\mu$m)")
fig.savefig(os.path.join(path,os.pardir,"images","customCompImageFastRotatorVis.png"))


#%%,

fig,ax,im=m.showModel(512,0.025,wl=[wl[0],wl[-1]],legend=True, 
   savefig=os.path.join(path,os.pardir,"images","customCompImageFastRotator.png"))


#%%
c.params['f'].value=0.9
c.params['pa'].value=30
ud=oim.oimUD(d=1,f=0.1,y=10)
m2=oim.oimModel(c,ud)

vc=m2.getComplexCoherentFlux(spfx_arr,spfy_arr,wl_arr)
v=np.abs(vc.reshape(nwl,nB))

fig,ax=plt.subplots(1,1)
label=["East-West Baselines","North-South Baselines"]
for iwl in range(nwl):
    cwl=iwl/(nwl-1)
    ax.plot(B/wl[iwl]/units.rad.to(units.mas),v[iwl,:nB//2],
            label=label[0],color=plt.cm.plasma(cwl))
    ax.plot(B/wl[iwl]/units.rad.to(units.mas),v[iwl,nB//2:],
            alpha=0.1,label=label[1],color=plt.cm.plasma(cwl))  
    label=[None,None]
ax.legend()
ax.set_xlabel("B/$\lambda$ (cycles/mas)" )
ax.set_ylabel("Visibility" )

norm = colors.Normalize(vmin=np.min(wl)*1e6,vmax=np.max(wl)*1e6)
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax,label="$\\lambda$ ($\\mu$m)")
fig.savefig(os.path.join(path,os.pardir,"images","customCompImageFastRotatorVis2.png"))

#%%,

fig,ax,im=m2.showModel(512,0.1,wl=[wl[0],wl[-1]],legend=True, normalize=True,normPow=0.5,
   savefig=os.path.join(path,os.pardir,"images","customCompImageFastRotator2.png"))



