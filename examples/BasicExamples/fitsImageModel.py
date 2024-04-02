# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 12:03:58 2023

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import numpy as np
import oimodeler as oim
from astropy.io import fits
from matplotlib import pyplot as plt


# NOTE: You can change FFT option, for instance reduce the standard
# zero-padding factor from 8 to 2 or use FFTW backend instead of the
# standard numpy FFT module 
oim.oimOptions.ft.padding = 8
#oim.oimOptions.ft.backend.active = oim.FFTWBackend

path = Path(__file__).parent.parent.parent
file_name = path / "examples" / "BasicExamples" / "BeDISCO.fits"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# %% creating the model
im = fits.open(file_name)

# load image from an opened astropy.io.primaryHDU
c = oim.oimComponentFitsImage(im)

# c=oim.oimComponentFitsImage(filename) # load image from a valid filename of fits file
m = oim.oimModel(c)

# %% Plotting the model image
m.showModel(512, 0.05, legend=True, normalize=True, normPow=1, cmap="hot", 
            figsize=(7, 5.5),savefig=save_dir / "FitsImage_Disco_image.png")

# %%
# Create some spatial frequencies (Baselines from 0 to 120m at 1.5 microns)
wl, nB = 1.5e-6, 1000
B = np.linspace(0, 120, num=nB)

# 1st half of B array are baseline in the East-West orientation
spfx = np.append(B, B*0)/wl
# 2nd half are baseline in the North-South orientation
spfy = np.append(B*0, B)/wl


ccf = m.getComplexCoherentFlux(spfx, spfy)
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
plt.close()

# %%
pprint(m.getParameters())

# %% scaling and rotating
c.params['pa'].value = 40
c.params['scale'].value = 0.8

m.showModel(512, 0.05, legend=True, normalize=True, normPow=1, cmap="hot", 
            figsize=(7, 5.5),savefig=save_dir / "FitsImage_Disco_image2.png")

# %%Adding a companion

c2 = oim.oimUD(x=20, d=1, f=0.03)
m2 = oim.oimModel(c, c2)

m2.showModel(512, 0.1, legend=True, normalize=True, fromFT=True, normPow=1,
             cmap="hot", savefig=save_dir / "FitsImage_Disco_image3.png")

# %%
ccf = m2.getComplexCoherentFlux(spfx, spfy)
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
