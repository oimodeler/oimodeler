# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 12:03:58 2023

@author: Ame
"""

import oimodeler as oim
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import os
from astropy.io import fits
from datetime import datetime


oim.oimOptions['FTpaddingFactor']=2
oim.oimOptions['FTBackend']=oim.FFTWBackend

path = os.path.dirname(oim.__file__)
pathData=os.path.join(path,os.pardir,"examples","BasicExamples")
filename=os.path.join(pathData,"KinematicsBeDiskModel.fits")

im=fits.open(filename)

#%% creating the model

c=oim.oimComponentFitsImage(im) # load image from an opened astropy.io.primaryHDU
#c=oim.oimComponentFitsImage(filename) # load image from a valid filename of fits file

m=oim.oimModel(c,c2)

#%% Plotting the model image
wl0=2.1661e-6
Dwl=100e-10
nwl=7
wl=np.linspace(wl0-Dwl/2,wl0+Dwl/2,num=nwl)
m.showModel(256,0.07,wl=wl,legend=True,normalize=True,fromFT=False,normPow=0.5)
#%% Computing and plotting visibilities for various baselines and walvelengths

nB=1000
nwl=51
wl=np.linspace(wl0-Dwl/2,wl0+Dwl/2,num=nwl)


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
v=v/np.tile(v[:,0][:,None],(1,nB))

fig,ax=plt.subplots(1,2,figsize=(15,5))
titles=["East-West Baselines","North-South Baselines"]


for iB in range(nB):
    cB=(iB % (nB//2))/(nB//2-1)
    ax[2*iB//nB].plot(wl*1e9,v[:,iB],
            color=plt.cm.plasma(cB))
    

for i in range(2):
    ax[i].set_title(titles[i])
    ax[i].set_xlabel("$\lambda$ (nm)")
ax[0].set_ylabel("Visibility")    
ax[1].get_yaxis().set_visible(False)   

        
norm = colors.Normalize(vmin=np.min(B),vmax=np.max(B))
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax,label="B (m)")


#fig.savefig(os.path.join(path,os.pardir,"images","customCompImageFastRotatorVis.png"))

