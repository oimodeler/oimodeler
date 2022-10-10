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

dim = 128
pixSize=0.2#mas
fact=8

D=12.8 #mas


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

interp='linear'
# Creating UD model for radial profile
c1=myDisk2D(d=D,dim=dim*fact,pixSize=pixSize)
c2=myDisk2D(d=D,dim=dim,pixSize=pixSize,interpMethod='linear')
c3=myDisk2D(d=D,dim=dim,pixSize=pixSize,interpMethod=interp)
modelFFT1=oim.oimModel([c1])
modelFFT2=oim.oimModel([c2])
modelFFT3=oim.oimModel([c3])
# Creating UD model directly in Fourier
ud=oim.oimUD(d=D)
modelFT=oim.oimModel([ud])


start_time = time.time()
for i in range(n):
    vFFT1=np.abs(modelFFT1.getComplexCoherentFlux(spf,spf*0))
    vFFT1=vFFT1/vFFT1[0]    
dt1=(time.time() - start_time)/n*1000  
print("FFT linear {:.1f}ms".format(dt1))    


start_time = time.time()
for i in range(n):
    vFFT2=np.abs(modelFFT2.getComplexCoherentFlux(spf,spf*0))
    vFFT2=vFFT2/vFFT2[0]   
dt2=(time.time() - start_time)/n*1000    
print("FFT linear {:.1f}ms".format(dt2))    

start_time = time.time()
for i in range(n):
    vFFT3=np.abs(modelFFT3.getComplexCoherentFlux(spf,spf*0))
    vFFT3=vFFT3/vFFT3[0]   
dt3=(time.time() - start_time)/n*1000    
print("FFT {} {:.1f}ms".format(interp,dt3)) 

start_time = time.time()
for i in range(n):
    vFT=np.abs(modelFT.getComplexCoherentFlux(spf,spf*0))
    vFT=vFT/vFT[0] 
dt0=(time.time() - start_time)/n*1000
print("Analytical FT {:.1f}ms".format(dt0))


#%%
fig,ax=plt.subplots(2,1,figsize=(10,7),sharex=True)


ax[0].plot(B/wl*mas2rad,vFT,label="FT from Formula  $\Rightarrow${:.1f}ms".format(dt0),color="k",lw=4)
ax[0].plot(B/wl*mas2rad,vFFT1,label="Image+FFT (linear interp. + zeropad) x{} $\Rightarrow${:.1f}ms".format(fact,dt1),color="tab:blue")
ax[0].plot(B/wl*mas2rad,vFFT2,label="Image+FFT (linear interp.)  $\Rightarrow${:.1f}ms".format(dt2),color="tab:green")
#ax[0].plot(B/wl*mas2rad,vFFT3,label="Image+FFT ({} interp.)  $\Rightarrow${:.1f}ms".format(interp,dt3),color="tab:red")

ax[1].plot(B/wl*mas2rad,(vFFT1-vFT)/vFT,color="tab:blue")
ax[1].plot(B/wl*mas2rad,(vFFT2-vFT)/vFT,color="tab:green")
#ax[1].plot(B/wl*mas2rad,(vFFT3-vFT)/vFT,color="tab:red")

#ax[0].set_yscale('log')
ax[0].set_yscale('log')

ax[1].set_xlabel("B/$\\lambda$ (cycles/mas)")
ax[0].set_ylabel("Visbility")
ax[0].plot([1.22*wl/D/mas2rad/wl*mas2rad]*2,[0,1],color="k",linestyle="--")
ax[0].legend()
ax[1].set_ylabel("Residuals")
ax[1].set_ylim(-1,1)


#%%
#Creating images of the models from their formula in direct space or using the 
#FT formula


normPow=0.5  # exponent for normalization of the image 

#dim=128   #Number of pixel for the image
pix=pixSize  #size of the pixel in the image in mas

fig,ax=plt.subplots(1,3,figsize=(10,5))


modelFT.showModel(dim,pix*fact,axe=ax[0],colorbar=False,fromFT=True)
modelFFT1.showModel(dim,pix*fact,axe=ax[1],colorbar=False,fromFT=True)
modelFFT2.showModel(dim,pix*fact,axe=ax[2],colorbar=False,fromFT=True)

#%%
fig,ax=plt.subplots(1,3,figsize=(10,5))

im0=modelFT.getImage(dim,pix)
im0=im0/np.sum(im0)
ax[0].imshow(im0,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
               norm=colors.PowerNorm(gamma=normPow))
ax[0].text(0,0.9*dim/2*pix,"FFt based model",va="top",ha="center",
           color="w",fontsize=14)
    
#get the image for the Fourier space formula
im1=c1.getInternalImage(dim,pix)
im1=im1/np.sum(im1)
ax[1].imshow(im1,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
               norm=colors.PowerNorm(gamma=normPow))
ax[1].text(0,0.9*dim/2*pix,"Image based model :\n internal image",va="top",ha="center",
           color="w",fontsize=14)

# Plot the difference between the two to check
im2=modelFFT1.getImage(dim,pix)
im2=im2/np.sum(im2)
ax[2].imshow(im2,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
               norm=colors.PowerNorm(gamma=normPow))
ax[2].text(0,0.9*dim/2*pix,"Image based model :\n Rescaled image",va="top",ha="center",
           color="w",fontsize=14)
    
ax[0].set_ylabel("$\\delta$(mas)")
for i in range(3):
    ax[i].set_xlabel("$\\alpha$(mas)")

