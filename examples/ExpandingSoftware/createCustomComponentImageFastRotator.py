# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:30:21 2022

@author: Ame
"""
import os

import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim
from astropy import units as units


path = os.path.dirname(oim.__file__)


# %%
def fastRotator(dim0, size, incl, rot, Tpole, lam, beta=0.25):
    """"""
    h = 6.63e-34
    c = 3e8
    kb = 1.38e-23

    a = 2./3*(rot)**0.4+1e-9
    K = np.sin(1./3.)*np.pi

    K1 = h*c/kb
    nlam = np.size(lam)
    incl = np.deg2rad(incl)

    x0 = np.linspace(-size, size, num=dim0)
    idx = np.where(np.abs(x0) <= 1.5)
    x = np.take(x0, idx)
    dim = np.size(x)
    unit = np.ones(dim)
    x = np.outer(x, unit)
    x = np.einsum('ij, k->ijk', x, unit)

    y = np.swapaxes(x, 0, 1)
    z = np.swapaxes(x, 0, 2)

    yp = y*np.cos(incl)+z*np.sin(incl)
    zp = y*np.sin(incl)-z*np.cos(incl)

    r = np.sqrt(x**2+yp**2+zp**2)

    theta = np.arccos(zp/r)

    x0 = (1.5*a)**1.5*np.sin(1e-99)
    r0 = a*np.sin(1/3.)*np.arcsin(x0)/(1.0/3.*x0)

    x2 = (1.5*a)**1.5*np.sin(theta)
    rin = a*np.sin(1/3.)*np.arcsin(x2)/(1.0/3.*x2)

    rhoin = rin*np.sin(theta)/a/K

    dr = (rin/r0-r) >= 0
    # dr=(rin/(r0*1.5)*1.5-r)>=0

    Teff = Tpole*(np.abs(1-rhoin*a)**beta)
    # TODO : implement a correct limb-darkening law
    # limb=np.abs((np.cos(np.arctan2(np.sqrt(x**2+y**2),-np.abs(z)))))*0+1

    if nlam == 1:
        flx = 1./(np.exp(K1/(lam*Teff))-1)

        im = np.zeros([dim, dim])

        for iz in range(dim):
            im = im*(im != 0)+(im == 0) * \
                dr[:, :, iz]*flx[:, :, iz]  # *limb[:,:,iz]

        im = np.rot90(im)

        tot = np.sum(im)
        im = im/tot
        im0 = np.zeros([dim0, dim0])

        im0[dim0//2-dim//2:dim0//2+dim//2, dim0//2-dim//2:dim0//2+dim//2] = im

        return im0

    else:
        unit = np.zeros(nlam)+1
        dr = np.einsum('ijk, l->ijkl', dr, unit)
        flx = 1./(np.exp(K1/np.einsum('ijk, l->ijkl', Teff, lam))-1)

        im = np.zeros([dim, dim, nlam])

        for iz in range(dim):
            im = im*(im != 0)+dr[:, :, iz, :]*flx[:, :, iz, :]*(im == 0)

        im = np.rot90(im)

        tot = np.sum(im, axis=(0, 1))

        for ilam in range(nlam):
            im[:, :, ilam] = im[:, :, ilam]/tot[ilam]

        im0 = np.zeros([dim0, dim0, nlam])
        im0[dim0//2-dim//2:dim0//2+dim//2, dim0//2-dim//2:dim0//2+dim//2, :] = im
        return im0


# %%
class oimFastRotator(oim.oimComponentImage):
    name = "Fast Rotator"
    shortname = "FRot"

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)

        # Component parameters. Note that as it inherits from the oimComponentImage class it already has
        # x,y,f and dim as parameters
        self.params["incl"] = oim.oimParam(
            name="incl", value=0, description="Inclination angle", unit=units.deg)
        self.params["rot"] = oim.oimParam(
            name="rot", value=0, description="Rotation Rate", unit=units.one)
        self.params["Tpole"] = oim.oimParam(
            name="Tpole", value=20000, description="Polar Temperature", unit=units.K)
        self.params["dpole"] = oim.oimParam(
            name="dplot", value=1, description="Polar diameter", unit=units.mas)
        self.params["beta"] = oim.oimParam(
            name="beta", value=0.25, description="Gravity Darkening Exponent", unit=units.one)

        # constant value <=> static model
        self._t = np.array([0])

        # The component is chromatic. Here we set a fixed array of reference wavelengths. This can be
        # modified later as, in our case the model is recomputed at each call to the fastRotator function
        self._wl = np.linspace(0.5e-6, 15e-6, num=10)

        # Finally evalutating paramters as for all other components
        self._eval(**kwargs)

    def _internalImage(self):
        dim = self.params["dim"].value
        incl = self.params["incl"].value
        rot = self.params["rot"].value
        Tpole = self.params["Tpole"].value
        dpole = self.params["dpole"].value
        beta = self.params["beta"].value

        im = fastRotator(dim, 1.5, incl, rot, Tpole, self._wl, beta=beta)

        # make a nt,nwl,dim,dim hcube (even if t and/or wl are not relevent)
        im = np.tile(np.moveaxis(im, -1, 0)[None, :, :, :], (1, 1, 1, 1))

        # computing the pixelSize based on the internal image size and the polar diameter
        self._pixSize = 1.5*dpole/dim*units.mas.to(units.rad)

        return im

# %% Creating a model
c = oimFastRotator(dpole=5, dim=128, incl=-70,
                   rot=0.99, Tpole=20000, beta=0.25)
m = oim.oimModel(c)

# %% Plotting the model image
m.showModel(512, 0.025, wl=[1e-6, 10e-6], legend=True, normalize=True,
            savefig=os.path.join(path, os.pardir, "images", "customCompImageFastRotator.png"))

# %% Computing and plotting visibilities for various baselines and walvelengths
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

fig, ax = plt.subplots(1, 2, figsize=(15, 5))
titles = ["East-West Baselines", "North-South Baselines"]
for iwl in range(nwl):
    cwl = iwl/(nwl-1)
    ax[0].plot(B/wl[iwl]/units.rad.to(units.mas), v[iwl, :nB//2],
               color=plt.cm.plasma(cwl))
    ax[1].plot(B/wl[iwl]/units.rad.to(units.mas), v[iwl, nB//2:],
               color=plt.cm.plasma(cwl))

for i in range(2):
    ax[i].set_title(titles[i])
    ax[i].set_xlabel("B/$\lambda$ (cycles/rad)")
ax[0].set_ylabel("Visibility")
ax[1].get_yaxis().set_visible(False)

norm = colors.Normalize(vmin=np.min(wl)*1e6, vmax=np.max(wl)*1e6)
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax, label="$\\lambda$ ($\\mu$m)")
fig.savefig(os.path.join(path, os.pardir, "images",
            "customCompImageFastRotatorVis.png"))


# %% Modifying the fast rotator parameters and  creating a new model with UD component
c.params['f'].value = 0.9
c.params['pa'].value = 30
ud = oim.oimUD(d=1, f=0.1, y=10)
m2 = oim.oimModel(c, ud)

# %% Plotting the new model image
m2.showModel(512, 0.06, wl=[1e-6, 10e-6], legend=True, normalize=True, normPow=0.5,
             savefig=os.path.join(path, os.pardir, "images", "customCompImageFastRotator2.png"))

# %% Computing and plotting visibilities for various baselines and walvelengths
vc = m2.getComplexCoherentFlux(spfx_arr, spfy_arr, wl_arr)
v = np.abs(vc.reshape(nwl, nB))

fig, ax = plt.subplots(1, 2, figsize=(15, 5))
titles = ["East-West Baselines", "North-South Baselines"]
for iwl in range(nwl):
    cwl = iwl/(nwl-1)
    ax[0].plot(B/wl[iwl]/units.rad.to(units.mas), v[iwl, :nB//2],
               color=plt.cm.plasma(cwl))
    ax[1].plot(B/wl[iwl]/units.rad.to(units.mas), v[iwl, nB//2:],
               color=plt.cm.plasma(cwl))

for i in range(2):
    ax[i].set_title(titles[i])
    ax[i].set_xlabel("B/$\lambda$ (cycles/rad)")
ax[0].set_ylabel("Visibility")
ax[1].get_yaxis().set_visible(False)

norm = colors.Normalize(vmin=np.min(wl)*1e6, vmax=np.max(wl)*1e6)
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax, label="$\\lambda$ ($\\mu$m)")
fig.savefig(os.path.join(path, os.pardir, "images",
            "customCompImageFastRotatorVis2.png"))


# %%  Plot an image and the visibility side by side ( for main page of documentation)
fig, ax = plt.subplot_mosaic([['upper left', 'upper right'],
                              ['lower', 'lower']], figsize=(8, 9))
ax = np.array(list(ax.values()))

m2.showModel(512, 0.06, wl=[1e-6, 2e-6], legend=True, normalize=True,
             normPow=1, axe=ax[:2], colorbar=False, cmap=plt.cm.plasma)

labels = ["East-West Baselines", "North-South Baselines"]
for iwl in range(nwl):
    cwl = iwl/(nwl-1)
    ax[2].plot(B/wl[iwl]/units.rad.to(units.mas), v[iwl, :nB//2],
               color=plt.cm.plasma(cwl), label=labels[0])
    ax[2].plot(B/wl[iwl]/units.rad.to(units.mas), v[iwl, nB//2:],
               color=plt.cm.plasma(cwl), alpha=0.1, label=labels[1])
    labels = [None, None]

ax[2].set_xlabel("B/$\lambda$ (cycles/rad)")
ax[2].set_ylabel("Visibility")
ax[2].legend()

plt.subplots_adjust(left=0.10, bottom=0.15, right=0.95,
                    top=0.95, wspace=0.25, hspace=0.25)

norm = colors.Normalize(vmin=np.min(wl)*1e6, vmax=np.max(wl)*1e6)
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sm, ax=ax[2], label="$\\lambda$ ($\\mu$m)")
fig.savefig(os.path.join(path, os.pardir, "images",
            "customCompImageFastRotatorImageAndVis.png"))
