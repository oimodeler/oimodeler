# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 12:21:51 2022

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import oimodeler as oim


path = Path(__file__).parent.parent.parent

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


# %% Function to simplify plotting of parameter an visibility
def plotParamAndVis(B, wl, t, model, param, ax=None, colorbar=True):
    nB = B.size

    if t is None:
        n = wl.size
        x = wl*1e6
        y = param(wl, 0)
        xlabel = r"$\lambda$ ($\mu$m)"
    else:
        n = t.size
        x = t-60000
        y = param(0, t)
        xlabel = "MJD - 60000 (days)"

    Bx_arr = np.tile(B[None, :], (n, 1)).flatten()
    By_arr = Bx_arr*0

    if t is None:
        t_arr = None
        wl_arr = np.tile(wl[:, None], (1, nB)).flatten()
        spfx_arr = Bx_arr/wl_arr
        spfy_arr = By_arr/wl_arr
    else:
        t_arr = np.tile(t[:, None], (1, nB)).flatten()
        spfx_arr = Bx_arr/wl
        spfy_arr = By_arr/wl
        wl_arr = None

    v = np.abs(model.getComplexCoherentFlux(
        spfx_arr, spfy_arr, wl=wl_arr, t=t_arr).reshape(n, nB))

    if ax is None:
        fig, ax = plt.subplots(2, 1)
    else:
        fig = ax.flatten()[0].get_figure()

    ax[0].plot(x, y, color="r")

    ax[0].set_ylabel("{} (mas)".format(param.name))
    ax[0].get_xaxis().set_visible(False)

    for iB in range(1, nB):
        ax[1].plot(x, v[:, iB]/v[:, 0], color=plt.cm.plasma(iB/(nB-1)))

    ax[1].set_xlabel(xlabel)
    ax[1].set_ylabel("Visibility")

    if colorbar == True:
        norm = colors.Normalize(vmin=np.min(B[1:]), vmax=np.max(B))
        sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
        fig.colorbar(sm, ax=ax, label="Baseline Length (m)")

    return fig, ax, v


# %%
nB = 200
B = np.linspace(0, 60, num=nB)

nwl = 1000
nt = 1000

# %%
c1 = oim.oimUD(d=oim.oimInterp('GaussWl', val0=2,
               value=4, x0=2.1656e-6, fwhm=1e-8))
m1 = oim.oimModel(c1)

wl = np.linspace(2.1e-6, 2.3e-6, num=nwl)

fig, ax, im = plotParamAndVis(B, wl, None, m1, c1.params['d'])
fig.suptitle(
    "Gaussian interpolator in $\lambda$ on a uniform disk diameter", fontsize=10)
plt.savefig(save_dir / "interp1.png")

pprint(c1.params['d'].params)
pprint(c1.params['d'].x0)
pprint(m1.getParameters())

# %%
c2 = oim.oimUD(f=0.5, d=oim.oimInterp("mGaussWl", val0=2, values=[4, 0, 0],
                                      x0=[2.05e-6, 2.1656e-6, 2.3e-6],
                                      fwhm=[2e-8, 2e-8, 1e-7]))
pt = oim.oimPt(f=0.5)
m2 = oim.oimModel(c2, pt)

c2.params['d'].values[1] = oim.oimParamLinker(
    c2.params['d'].values[0], "*", 3)
c2.params['d'].values[2] = oim.oimParamLinker(
    c2.params['d'].values[0], "+", -1)

wl = np.linspace(1.9e-6, 2.4e-6, num=nwl)

fig, ax, im = plotParamAndVis(B, wl, None, m2, c2.params['d'])
fig.suptitle(
    "Multiple Gaussian interpolator in $\lambda$ on a uniform disk diameter", fontsize=10)
plt.savefig(save_dir / "interp2.png")


# %%
c3 = oim.oimGauss(fwhm=oim.oimInterp(
    "cosTime", T0=60000, P=1, values=[1, 3], x0=0.8))
m3 = oim.oimModel(c3)

t = np.linspace(60000, 60006, num=nt)
wl = 2.2e-6

fig, ax, im = plotParamAndVis(B, wl, t, m3, c3.params['fwhm'])
fig.suptitle(
    "Assym. Cosine interpolator in Time on a Gaussian fwhm", fontsize=10)
plt.savefig(save_dir / "interp3.png")

# %%
c4 = oim.oimIRing(d=oim.oimInterp("wl", wl=[2e-6, 2.4e-6, 2.7e-6, 3e-6], values=[2, 6, 5, 6],
                                  kind="linear", extrapolate=True))
m4 = oim.oimModel(c4)
wl = np.linspace(1.8e-6, 3.2e-6, num=nwl)
fig, ax = plt.subplots(2, 6, figsize=(18, 6), sharex=True, sharey="row")

plotParamAndVis(B, wl, None, m4, c4.params['d'], ax=ax[:, 0], colorbar=False)
c4.params['d'].extrapolate = False
plotParamAndVis(B, wl, None, m4, c4.params['d'], ax=ax[:, 1], colorbar=False)

c4.params['d'].extrapolate = True
c4.params['d'].kind = "quadratic"
plotParamAndVis(B, wl, None, m4, c4.params['d'], ax=ax[:, 2], colorbar=False)
c4.params['d'].extrapolate = False
plotParamAndVis(B, wl, None, m4, c4.params['d'], ax=ax[:, 3], colorbar=False)

c4.params['d'].extrapolate = True
c4.params['d'].kind = "cubic"
plotParamAndVis(B, wl, None, m4, c4.params['d'], ax=ax[:, 4], colorbar=False)
c4.params['d'].extrapolate = False
plotParamAndVis(B, wl, None, m4, c4.params['d'], ax=ax[:, 5], colorbar=False)

plt.subplots_adjust(left=0.05, bottom=0.1, right=0.99, top=0.9,
                    wspace=0.05, hspace=0.05)
for i in range(1, 6):
    ax[0, i].get_yaxis().set_visible(False)
    ax[1, i].get_yaxis().set_visible(False)

fig.suptitle("Linear, Quadratic and Cubic interpolators (with extrapolation"
             r" or fixed values outside the range) in $\lambda$ on a uniform"
             " disk diameter", fontsize=18)
plt.savefig(save_dir / "interp4.png")
# %%
c5 = oim.oimUD(d=oim.oimInterp('polyTime', coeffs=[1, 3.5, -0.5], x0=60000))
m5 = oim.oimModel(c5)

wl = 2.2e-6
t = np.linspace(60000, 60006, num=nt)

fig, ax, im = plotParamAndVis(B, wl, t, m5, c5.params['d'])
fig.suptitle(
    "Polynomial interpolator in Time on a uniform disk diameter", fontsize=10)
plt.savefig(save_dir / "interp5.png")


#%%
star1 = oim.oimUD(f=oim.oimInterp('polyTime', coeffs=[1, 3.5, -0.5], x0=60000))
