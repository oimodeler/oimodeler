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

dim = 1024
pixSize=0.5#mas


D=10 #mas


class myDisk2D(oim.oimComponentImage):
    name="A disk define by a 2D image"
    shortname = "myDisk"
    
    def __init__(self,**kwargs): 
         # First call the __init__ function from the parent class
         super().__init__(**kwargs)
         # Then define the component parameters  
         self.params["d"]=oim.oimParam(**(oim._standardParameters["d"]))
        
         
         # Finally call the _eval function that allow the parameters to be processed
         self._eval(**kwargs)
         
    def _imageFunction(self,xx,yy,wl,t):
        
        return  ((xx**2+yy**2)<=(self.params["d"](wl,t)/2)**2).astype(float)




wl=0.65e-6
B=np.linspace(0.01,100,num=2000)
spf=B/wl


n=10 

spf0=1.22*wl/D/mas2rad/wl


# Creating UD model for radial profile
c=myDisk2D(d=D,dim=dim,pixSize=pixSize)
modelFFT=oim.oimModel([c])


# Creating UD model directly in Fourier
ud=oim.oimUD(d=D)
modelFT=oim.oimModel([ud])


start_time = time.time()
for i in range(n):
    vFFT=np.abs(modelFFT.getComplexCoherentFlux(spf,spf*0))
    vFFT=vFFT/vFFT[0]
print("FFT {:.1f}ms".format((time.time() - start_time)/n*1000))    


start_time = time.time()
for i in range(n):
    vFT=np.abs(modelFT.getComplexCoherentFlux(spf,spf*0))
    vFT=vFT/vFT[0]
print("Analytical FT {:.1f}ms".format((time.time() - start_time)/n*1000))



fig,ax=plt.subplots(1,1,figsize=(10,7))

ax.plot(B/wl*mas2rad,vFFT,label="Direct image + FFT")
ax.plot(B/wl*mas2rad,vFT,label="FT from Formula")

ax.set_xlabel("B/$\\lambda$ (cycles/mas)")
ax.set_ylabel("Visbility")
ax.plot([1.22*wl/D/mas2rad/wl*mas2rad]*2,[0,1],color="k",linestyle="--")
plt.legend()




#%%
"""
Creating images of the models from their formula in direct space or using the 
FT formula
"""

normPow=0.5  # exponent for normalization of the image 

#dim=128   #Number of pixel for the image
pix=pixSize  #size of the pixel in the image in mas

fig,ax=plt.subplots(1,3,figsize=(10,5))


im0=modelFT.getImage(dim,pix)
im0=im0/np.sum(im0)
ax[0].imshow(im0,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
               norm=colors.PowerNorm(gamma=normPow))
ax[0].text(0,0.9*dim/2*pix,"FFt based model",va="top",ha="center",
           color="w",fontsize=14)
    
#get the image for the Fourier space formula
im1=c.getInternalImage(dim,pix)
im1=im1/np.sum(im1)
ax[1].imshow(im1,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
               norm=colors.PowerNorm(gamma=normPow))
ax[1].text(0,0.9*dim/2*pix,"Image based model :\n internal image",va="top",ha="center",
           color="w",fontsize=14)

# Plot the difference between the two to check
im2=modelFFT.getImage(dim,pix)
im2=im2/np.sum(im2)
ax[2].imshow(im2,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
               norm=colors.PowerNorm(gamma=normPow))
ax[2].text(0,0.9*dim/2*pix,"Image based model :\n Rescaled image",va="top",ha="center",
           color="w",fontsize=14)
    
ax[0].set_ylabel("$\\delta$(mas)")
for i in range(3):
    ax[i].set_xlabel("$\\alpha$(mas)")
    
 

