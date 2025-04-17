# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 10:22:47 2025

@author: ame
"""

from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim
from astropy.io import fits
import time

path = Path(__file__).parent.parent.parent

# NOTE: Change these path if you want to save the products at another location
save_dir = path / "images"
product_dir = path / "data"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


#%% getting the list of all radial-profile-based components currently available.
print(oim.listComponents(componentType="radial"))


#%% Creating a flattened exponential ring component

c = oim.oimExpRing(d=10,fwhm=1,elong=1.5,pa=90)
m = oim.oimModel(c)

#%% plotting it's image
fig, ax, im = m.showModel(256,0.5,figsize=(5,4))

fig.savefig(save_dir / "radialProfile_image_exp.png")

#%%

c1 = oim.oimIRing(d=10,elong=1.5,pa=90)
m1 = oim.oimModel(c1)

c2 = oim.oimRing(din=10,dout=13,elong=1.5,pa=90)
m2 = oim.oimModel(c2)

fig, ax = plt.subplots(1,3,figsize=(15,5))

m.showModel(256,0.15,figsize=(5,4),axe=ax[0],colorbar=False)
m1.showModel(256,0.15,figsize=(5,4),axe=ax[1],fromFT=True,colorbar=False)
m2.showModel(256,0.15,figsize=(5,4),axe=ax[2],fromFT=True,colorbar=False)

fig.savefig(save_dir / "radialProfile_image_comp.png")


#%%
wl = 2.1e-6
B = np.linspace(0, 100, num=10000)
spf = B/wl

start = time.time()
ccf = m.getComplexCoherentFlux(spf, spf*0)
v = np.abs(ccf/ccf[0])
dt = (time.time() - start)*1000


start = time.time()
ccf1 = m1.getComplexCoherentFlux(spf, spf*0)
v1 = np.abs(ccf1/ccf1[0])
dt1 = (time.time() - start)*1000

start = time.time()
ccf2 = m2.getComplexCoherentFlux(spf, spf*0)
v2 = np.abs(ccf2/ccf2[0])
dt2 = (time.time() - start)*1000

plt.figure()
plt.plot(B, v, label=f"Exponential Ring ({dt:.1f}ms)")
plt.plot(B, v1, label=f"Infinitesimal Ring ({dt1:.1f}ms)")
plt.plot(B, v2, label=f"Uniform Ring ({dt2:.1f}ms)")
plt.xlabel("B (m)")
plt.ylabel("Visbility")
plt.legend()
plt.margins(0)

plt.savefig(save_dir / "radialProfile_visi_comp.png")

#%%