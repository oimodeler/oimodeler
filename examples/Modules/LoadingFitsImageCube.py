# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 12:03:58 2023

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import oimodeler as oim
from matplotlib import pyplot as plt

# NOTE: You can change FFT option, for instance reduce the standard
# zero-padding factor from 8 to 2 or use FFTW backend instead of the
# standard numpy FFT module oim.oimOptions.ft.padding = 2
# oim.oimOptions.ft.backend.active = oim.FFTWBackend

path = Path(__file__).parent.parent.parent

# NOTE: Change these path if you want to save the products at another location
save_dir = path / "images"
product_dir = path / "data"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

file_name = product_dir / "IMAGES" / "KinematicsBeDiskModel.fits"

# NOTE: Change this path if you want to save the products at a certain location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# %% creating the model
c = oim.oimComponentFitsImage(file_name)
m = oim.oimModel(c)

#%% accessing the internal image and  wavelength table
print(c._image.shape)
print(c._wl)

# %% Plotting the model image
wl0, dwl, nwl = 2.1661e-6, 60e-10, 5
wl = np.linspace(wl0-dwl/2, wl0+dwl/2, num=nwl)
fig, ax, im  = m.showModel(256, 0.04, wl=wl, legend=True, normPow=0.4, colorbar=False,
            figsize=(2, 2.5), savefig=save_dir / "FitsImageCube_BeDiskKinematicsModel_images.png")


# %% Computing and plotting visibilities for various baselines and walvelengths
nB = 1000
nwl = 51

wl = np.linspace(wl0-dwl/2, wl0+dwl/2, num=nwl)
B = np.linspace(0, 100, num=nB//2)

# The 1st half of B array are baseline in the East-West orientation
# and the 2nd half are baseline in the North-South orientation
Bx = np.append(B, B*0) 
By = np.append(B*0, B)  

#creating the spatial frequencies and wls arrays nB x nwl by matrix multiplication
spfx = (Bx[np.newaxis,:]/wl[:,np.newaxis]).flatten()
spfy = (By[np.newaxis,:]/wl[:,np.newaxis]).flatten()
wls  = (wl[:,np.newaxis]*np.ones([1,nB])).flatten()

#computing the complex coherent flux and visbility
vc = m.getComplexCoherentFlux(spfx, spfy, wls)
v = np.abs(vc.reshape(nwl, nB))
v = v/np.tile(v[:, 0][:, None], (1, nB))

#plotting the results
fig, ax = plt.subplots(1, 2, figsize=(8, 4))
titles = ["East-West Baselines", "North-South Baselines"]

for iB in range(nB):
    cB = (iB % (nB//2))/(nB//2-1)
    ax[2*iB//nB].plot(wl*1e9, v[:, iB],
                      color=plt.cm.plasma(cB))

for i in range(2):
    ax[i].set_title(titles[i])
    ax[i].set_xlabel(r"$\lambda$ (nm)")
ax[0].set_ylabel("Visibility")
ax[1].get_yaxis().set_visible(False)

norm = colors.Normalize(vmin=np.min(B), vmax=np.max(B))
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax, label="B (m)")

fig.savefig(save_dir / "FitsImageCube_BeDiskKinematicsModel_visibility.png")
