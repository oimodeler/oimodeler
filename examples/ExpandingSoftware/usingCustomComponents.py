# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:30:21 2022

@author: Ame
"""
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim
from astropy import units as units
from oimodeler.oimCustomComponents import oimFastRotator, oimSpiral, oimBox


path = Path().resolve().parent.parent

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# %% Creating a model
fastRot = oimFastRotator(dpole=5, dim=128, incl=-70, rot=0.99,
                         Tpole=20000, beta=0.25, f=0.9, pa=30, x=10, y=5)
spiral = oimSpiral(dim=128, fwhm=2.5, P=0.1, width=0.1,
                   pa=30, elong=2, x=10, y=-5, f=0.2)
gauss = oim.oimGauss(y=0, x=-5, fwhm=oim.oimInterp("wl", wl=[1e-6, 2e-6], values=[0.5, 4]),
                     f=oim.oimInterp("wl", wl=[1e-6, 2e-6], values=[0, 0.1]))
box = oimBox(dx=2, dy=oim.oimInterp(
    "wl", wl=[1e-6, 2e-6], values=[0.5, 2]), f=0.1, pa=20)
m = oim.oimModel(fastRot, spiral, gauss, box)

# %% Computing and  visibilities for various baselines and walvelengths
nB = 1000
nwl = 20
wl = np.linspace(1e-6, 2e-6, num=nwl)

B = np.linspace(0, 100, num=nB//2)

# 1st half of B array are baseline in the East-West orientation
Bx = np.append(B, B*0)
By = np.append(B*0, B)  # 2nd half are baseline in the North-South orientation

Bx_arr = np.tile(Bx[None, :], (nwl, 1)).flatten()
By_arr = np.tile(By[None, :], (nwl,  1)).flatten()
wl_arr = np.tile(wl[:, None], (1, nB)).flatten()

spfx_arr = Bx_arr/wl_arr
spfy_arr = By_arr/wl_arr

vc = m.getComplexCoherentFlux(spfx_arr, spfy_arr, wl_arr)
v = np.abs(vc.reshape(nwl, nB))

# %%  Plot an image and the visibility side by side ( for main page of documentation)
fig, ax = plt.subplot_mosaic([['upper left', 'upper middle', 'upper right'],
                              ['lower', 'lower', 'lower']], figsize=(8, 6))
ax = np.array(list(ax.values()))

m.showModel(512, 0.06, wl=[1e-6, 1.5e-6, 2e-6], legend=True, normalize=True,
            normPow=0.2, axe=ax[:3], colorbar=False, cmap=plt.cm.plasma)

labels = ["East-West Baselines", "North-South Baselines"]
for iwl in range(nwl):
    cwl = iwl/(nwl-1)
    ax[3].plot(B/wl[iwl]/units.rad.to(units.mas), v[iwl, :nB//2],
               color=plt.cm.plasma(cwl), label=labels[0])
    ax[3].plot(B/wl[iwl]/units.rad.to(units.mas), v[iwl, nB//2:],
               color=plt.cm.plasma(cwl), alpha=0.1, label=labels[1])
    labels = [None, None]

ax[3].set_xlabel("B/$\lambda$ (cycles/rad)")
ax[3].set_ylabel("Visibility")
ax[3].legend()

plt.subplots_adjust(left=0.10, bottom=0.15, right=0.95,
                    top=0.95, wspace=0.25, hspace=0.25)

norm = colors.Normalize(vmin=np.min(wl)*1e6, vmax=np.max(wl)*1e6)
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax[3], label="$\\lambda$ ($\\mu$m)")
fig.savefig(save_dir / "usingCustomComponents.png")
