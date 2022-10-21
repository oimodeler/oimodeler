# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:30:21 2022

@author: Ame
"""
import numpy as np
import matplotlib.pyplot as plt
import oimodeler as oim
from fastRotator import fastRotator
from astropy import units as units
from datetime import datetime

class oimFastRotator(oim.oimComponentImage):
    def __init__(self,**kwargs):
        super(). __init__(**kwargs)
        
        #Component parameters. Note that as it inherits from the oimComponentImage class it already has
        #x,y,f and dim as parameters
        self.params["incl"]=oim.oimParam(name="incl",value=0,description="Inclination angle",unit=units.deg)
        self.params["rot"]=oim.oimParam(name="rot",value=0,description="Rotation Rate",unit=units.one)     
        self.params["Tpole"]=oim.oimParam(name="Tpole",value=20000,description="Polar Temperature",unit=units.K)
        self.params["dpole"]=oim.oimParam(name="dplot",value=1,description="Polar diameter",unit=units.mas)
        self.params["beta"]=oim.oimParam(name="beta",value=0.25,description="Gravity Darkening Exponent",unit=units.one)
        self._t = np.array([0]) # constant value <=> static model
        
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
dim=32

#Uncomment this to use FFTW backend instead of numpy backend for FFT computation
#oim.oimOptions['FTBackend']=oim.FFTWBackend

c=oimFastRotator(dpole=5,dim=dim,incl=70,rot=0.99,Tpole=20000,x=0,y=0,pa=0,f=1,beta=0.25)
m=oim.oimModel(c)

#%%
#Computing and plotting visibilities for various baselines and walvelengths

nB=5000
nwl=5
wl=np.linspace(0.5e-6,2e-6,num=nwl)

B=np.linspace(1,100,num=nB//2)
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
    ax.plot(B/wl[iwl],v[iwl,:nB//2],label=label[0],color=plt.cm.viridis(cwl))
    ax.plot(B/wl[iwl],v[iwl,nB//2:],alpha=0.5,label=label[1],color=plt.cm.viridis(cwl))  
    label=[None,None]
ax.legend()

#%%,
fig,axe,im=m.showModel(256,0.1,wl=wl,swapAxes=True,legend=True)


