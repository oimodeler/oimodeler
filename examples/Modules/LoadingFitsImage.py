# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 13:52:46 2025

@author: ame
"""

from pathlib import Path
from pprint import pprint
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim
from astropy.io import fits

path = Path(__file__).parent.parent.parent

# NOTE: Change these path if you want to save the products at another location
save_dir = path / "images"
product_dir = path / "data"
if not save_dir.exists():
    save_dir.mkdir(parents=True)
    
#%% The oimFitImageComponent

file_name = path / "data" / "BeDISCO.fits"
im = fits.open(file_name)

# load image from an opened astropy.io.primaryHDU
cdisco = oim.oimComponentFitsImage(im)

# load image from a valid filename of fits file
# cdisco=oim.oimComponentFitsImage(filename)
mdisco = oim.oimModel(cdisco)


# %% Plotting the model image
mdisco.showModel(512, 0.1, legend=True, normalize=True, normPow=1, cmap="hot", 
            figsize=(7, 5.5),savefig=save_dir / "FitsImage_Disco_image.png")

#%%
im_disco = cdisco._internalImage()
print(im_disco.shape)

#%%

pixSize = cdisco._pixSize*u.rad.to(u.mas)
dim = im_disco.shape[-1]
mdisco.showModel(dim, pixSize, legend=True, normalize=True, normPow=1, cmap="hot", 
            figsize=(7, 5.5),savefig=save_dir / "FitsImage_Disco_internal_image.png")

#%% Create some spatial frequencies (Baselines from 0 to 120m at 1.5 microns)
wl, nB = 1.5e-6, 1000
B = np.linspace(0, 120, num=nB)

# 1st half of B array are baseline in the East-West orientation
spfx = np.append(B, B*0)/wl
# 2nd half are baseline in the North-South orientation
spfy = np.append(B*0, B)/wl


ccf = mdisco.getComplexCoherentFlux(spfx, spfy)
v = np.abs(ccf)
v = v/v.max()

plt.figure()
plt.plot(B, v[0:nB], label="East-West")
plt.plot(B, v[nB:], label="North-South")
plt.xlabel("B (m)")
plt.ylabel("Visbility")
plt.legend()
plt.margins(0)

plt.savefig(save_dir / "FitsImage_Disco_visibility.png")
#plt.close()

# %%
pprint(mdisco.getParameters())

# %% scaling and rotating
cdisco.params['pa'].value = 40
cdisco.params['scale'].value = 0.8

mdisco.showModel(512, 0.05, legend=True, normalize=True, normPow=1, cmap="hot", 
            figsize=(7, 5.5),savefig=save_dir / "FitsImage_Disco_image2.png")

# %%Adding a companion

cud = oim.oimUD(x=20, d=1, f=0.03)
mdisco_ud = oim.oimModel(cdisco, cud)

mdisco_ud.showModel(512, 0.1, legend=True, normalize=True, fromFT=True, normPow=1,
             cmap="hot", savefig=save_dir / "FitsImage_Disco_image3.png")

# %%
ccf = mdisco_ud.getComplexCoherentFlux(spfx, spfy)
v = np.abs(ccf)
v = v/v.max()

plt.figure()
plt.plot(B, v[0:nB], label="East-West")
plt.plot(B, v[nB:], label="North-South")
plt.xlabel("B (m)")
plt.ylabel("Visbility")
plt.legend()
plt.margins(0)
plt.savefig(save_dir / "FitsImage_Disco_visibility2.png")

#%%
