# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 12:03:58 2023

@author: Ame
"""
import os
from pprint import pprint as print

import numpy as np
import oimodeler as oim
from astropy.io import fits
from matplotlib import pyplot as plt

#You can change FFT option, for instance reduce the standard zero-padding factor 
#from 8 to 2 or use FFTW backend instead of the standard numpy FFT module
#oim.oimOptions['FTpaddingFactor']=2
#oim.oimOptions['FTBackend']=oim.FFTWBackend


path = os.path.dirname(oim.__file__)
pathData=os.path.join(path,os.pardir,"examples","BasicExamples")
filename=os.path.join(pathData,"BeDISCO.fits")

im=fits.open(filename)

#%% creating the model

im=fits.open(filename)
c=oim.oimComponentFitsImage(im) # load image from an opened astropy.io.primaryHDU
#c=oim.oimComponentFitsImage(filename) # load image from a valid filename of fits file

m=oim.oimModel(c)

#%% Plotting the model image
m.showModel(512,0.05,legend=True,normalize=True,normPow=1,cmap="hot",figsize=(7,5.5),
     savefig=os.path.join(path,os.pardir,"images","FitsImage_Disco_image.png"))

#%%
#Create some spatial frequencies (Baselines from 0 to 120m at 1.5 microns)
wl=1.5e-6
nB=1000
B=np.linspace(0,120,num=nB)

spfx=np.append(B,B*0)/wl # 1st half of B array are baseline in the East-West orientation
spfy=np.append(B*0,B)/wl # 2nd half are baseline in the North-South orientation



ccf = m.getComplexCoherentFlux(spfx,spfy)
v = np.abs(ccf)
v=v/v.max()


plt.figure()
plt.plot(B , v[0:nB],label="East-West")
plt.plot(B , v[nB:],label="North-South")
plt.xlabel("B (m)")
plt.ylabel("Visbility")
plt.legend()
plt.margins(0)

plt.savefig(os.path.join(path,os.pardir,"images",
                         "FitsImage_Disco_visibility.png"))


#%%
print(m.getParameters())

#%% scaling and rotating

c.params['pa'].value=40
c.params['scale'].value=0.8

m.showModel(512,0.05,legend=True,normalize=True,normPow=1,cmap="hot",figsize=(7,5.5),
     savefig=os.path.join(path,os.pardir,"images","FitsImage_Disco_image2.png"))

#%%Adding a companion

c2=oim.oimUD(x=20,d=1,f=0.03)
m2=oim.oimModel(c,c2)

m2.showModel(512,0.1,legend=True,normalize=True,fromFT=True,normPow=1,
            cmap="hot",savefig=os.path.join(path,os.pardir,"images",
                 "FitsImage_Disco_image3.png"))


#%%

ccf = m2.getComplexCoherentFlux(spfx,spfy)
v = np.abs(ccf)
v=v/v.max()

plt.figure()
plt.plot(B , v[0:nB],label="East-West")
plt.plot(B , v[nB:],label="North-South")
plt.xlabel("B (m)")
plt.ylabel("Visbility")
plt.legend()
plt.margins(0)

plt.savefig(os.path.join(path,os.pardir,"images",
                         "FitsImage_Disco_visibility2.png"))
