import oimodeler as oim
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import astropy.units as u

mas2rad=u.mas.to(u.rad)

# Create some spatial frequencies (Baselines from 0 to 400m at 2.1 microns)
lam=2.1e-6
B=np.linspace(0.01,400,num=1000)
spf=B/lam
spf0=spf*0


# Some components with parameter chromaticity
ud=oim.oimUD(d=1,f=oim.oimInterpWl([2e-6,2.5e-6],[1,0.5]))
eg=oim.oimEGauss(fwhm=4,pa=50,elong=2,f=0.4)
skwer=oim.ESKRing(d=15,elong=2,skw=1,f=0.5)

# linking parameters together
# 1st solution : replacing one param by another
skwer.params["elong"]=eg.params["elong"]
# 2nd solution using the oimParamLinker class
skwer.params["pa"]=oim.oimParamLinker(eg.params["pa"])
# This allows to add an offset or a factor to the linked param
skwer.params["skwPa"]=oim.oimParamLinker(eg.params["pa"],"add",180)


# Building the model from the components 
model=oim.oimModel([ud, eg, skwer])

# Getting all model parameters:
# Duplicated or linked parameters are shown only once
params=model.getParameters()

# Getting only the free parameters:
# Note that by default x & y are not free
fparams=model.getFreeParameters()




fig,ax=plt.subplots(2,1,figsize=(10,7))

#Compute visibilities for the models and spatial frequencies
v1=np.abs(model.getComplexCoherentFlux(spf,spf0)) # East-West baselines
v2=np.abs(model.getComplexCoherentFlux(spf0,spf)) # North-South baselines
Ftot=v1[0]
plt.plot(B/lam*mas2rad,v1/Ftot)
plt.plot(B/lam*mas2rad,v2/Ftot,ls=":")
plt.xlabel("B/$\\lambda$ (cycle/mas)")
plt.ylabel("Visibility")  

plt.margins(0,0)




#%%
# Create images of the model from their formula in direct space or using 
# the FT formula

normPow=0.1
dim=1024   #Number of pixel for the image
pix=0.1  #size of the pixel in the image in mas
fig,ax=plt.subplots(1,1,figsize=(9,9))

#get the image for the direct space formula
im=model.getImage(dim,pix,fromFT=True)
im=im/np.sum(im)
ax.imshow(im,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
               norm=colors.PowerNorm(gamma=normPow))
ax.set_xlabel("$\\alpha$(mas)")
ax.set_ylabel("$\\delta$(mas)")