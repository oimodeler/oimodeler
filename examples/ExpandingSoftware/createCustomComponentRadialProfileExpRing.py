# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:39:36 2021

@author: Ame
"""
import astropy.units as u
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


mas2rad = u.mas.to(u.rad)
dim = 256

###############################################################################


class oimExpRing(oim.oimComponentRadialProfile):
    name = "A ring with a descreasing exponential profil"
    shortname = "ExpR"

    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oim.oimParam(**(oim._standardParameters["d"]))
        self.params["fwhm"] = oim.oimParam(**(oim._standardParameters["fwhm"]))

        self._dim = dim

        self._t = np.array([0])  # constant value <=> static model
        self._wl = None  # np.array([0.5,1])*1e-6
        # self._r = np.arange(0, self._dim)*pixSize

        # Finally call the _eval function that allow the parameters to be processed
        self._eval(**kwargs)

    def _radialProfileFunction(self, r, wl, t):

        r0 = self.params["d"](wl, t)/2
        fwhm = self.params["fwhm"](wl, t)
        I = np.nan_to_num((r > r0).astype(float) * np.exp(-0.692*np.divide(r-r0, fwhm)), nan=0)

        return I

    @property
    def _r(self):
        if False:
            fwhm_max = np.max(self.params["fwhm"](self._wl, self._t))
            r0_max = np.max(self.params["d"](self._wl, self._t))/2
        else:
            fwhm_max = self.params["fwhm"](1e99)
            r0_max = self.params["d"](1e99)
        rmax = r0_max+8*fwhm_max
        return np.linspace(0, 1, self._dim)*rmax

    @_r.setter
    def _r(self, r):
        pass


# %% Creating UD model for radial profile
ring = oimExpRing(d=5, fwhm=oim.oimInterp(
    "wl", wl=[0.5e-6, 1e-6], values=[0.2, 4]), elong=2, pa=30)
m = oim.oimModel(ring)

# %% Computing and plotting visibilities for various baselines and walvelengths

nB = 1000
nwl = 100
wl = np.linspace(0.5e-6, 1e-6, num=nwl)

B = np.linspace(0, 100, num=nB//2)
# 1st half of B array are baseline in the East-West orientation
Bx = np.append(B, B*0)
By = np.append(B*0, B)  # 2nd half are baseline in the North-South orientation

Bx_arr = np.tile(Bx[None, :], (nwl, 1)).flatten()
By_arr = np.tile(By[None, :], (nwl,  1)).flatten()
wl_arr = np.tile(wl[:, None], (1, nB)).flatten()

spfx_arr = Bx_arr/wl_arr
spfy_arr = By_arr/wl_arr

v = np.abs(m.getComplexCoherentFlux(
    spfx_arr, spfy_arr, wl_arr)).reshape(nwl, nB)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
title = ["East-West Baselines\n", "North-South Baselines\n"]

for iwl in range(nwl):
    cwl = iwl/(nwl-1)
    ax[0].plot(B/wl[iwl]/u.rad.to(u.mas),
               v[iwl, :nB//2], color=plt.cm.plasma(cwl))
    ax[1].plot(B/wl[iwl]/u.rad.to(u.mas),
               v[iwl, nB//2:], color=plt.cm.plasma(cwl))

ax[0].set_title("East-West Baselines")
ax[1].set_title("North-South Baselines")

norm = colors.Normalize(vmin=np.min(wl)*1e6, vmax=np.max(wl)*1e6)
sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax, label="$\\lambda$ ($\\mu$m)")

# %% plttnig images of the models from direct formula and from FT
pix = 0.15
dim = 512
fig, ax, im = m.showModel(dim, pix, wl=[wl[0], wl[nwl//2-1], wl[-1]], legend=True, normalize=True)
fig.suptitle("Direct images")
figFT, axFT, imFT = m.showModel(dim, pix, wl=[wl[0], wl[nwl//2-1], wl[-1]], legend=True, normalize=True, fromFT=True)
figFT.suptitle("Images computed from inverse FT")
