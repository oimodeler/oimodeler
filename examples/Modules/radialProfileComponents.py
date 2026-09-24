# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 10:22:47 2025

@author: ame
"""

import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import oimodeler as oim

path = Path(__file__).parent.parent.parent

# NOTE: Change these path if you want to save the products at another location
save_dir = path / "images"
product_dir = path / "data"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


# NOTE: Get the list of all radial-profile-based components currently available.
print(oim.listComponents(componentType="radial"))

# NOTE: Create a flattened exponential ring component
c = oim.oimRadialExpRing(d=10, fwhm=1, elong=1.5, pa=90)
m = oim.oimModel(c)

# NOTE: Plot it's image
fig, ax, im = m.showModel(256, 0.5, figsize=(5, 4))
fig.savefig(save_dir / "radialProfile_image_exp.png")


# %%

c1 = oim.oimRadialPowRing(din=10, dout=100,p=-3, elong=1.5, pa=90,dim=64)
m1 = oim.oimModel(c1)

c2 = oim.oimIRing(d=10, elong=1.5, pa=90)
m2 = oim.oimModel(c2)

c3 = oim.oimRing(din=10, dout=25, elong=1.5, pa=90)
m3 = oim.oimModel(c3)

fig, ax = plt.subplots(1, 4, figsize=(13, 4))

dim=256
pix=0.3
m.showModel(256, 0.3, figsize=(5, 4), axe=ax[0], colorbar=False,normPow=1)
m1.showModel(256, 0.3, figsize=(5, 4), axe=ax[1], colorbar=False,normPow=1)
m2.showModel(256, 0.3, figsize=(5, 4), axe=ax[2], fromFT=True, colorbar=False,normPow=1)
m3.showModel(256, 0.3, figsize=(5, 4), axe=ax[3], fromFT=True, colorbar=False,normPow=1)

cs=[c,c1,c2,c3]
for i in range(4):
    if i!=0:
        ax[i].get_yaxis().set_visible(False)
    ax[i].text(0,dim*pix/2.5,cs[i].name,color="w",ha="center",fontsize=15)
fig.tight_layout()

fig.savefig(save_dir / "radialProfile_image_comp.png")

# %%
wl = 2.1e-6
B = np.linspace(0, 100, num=10000)

fig, ax = plt.subplots()

ms = [m,m1,m2,m3]
for i in range(4):
    ms[i].plotVis(B,wl,axe=ax,label=cs[i].name,addTimeToLabel=True)
plt.legend()
plt.savefig(save_dir / "radialProfile_visi_comp.png")
#%%
spf= B/wl

dts=[]
for i in range(4):
    start = time.time()
    ccf = ms[i].getComplexCoherentFlux(spf, spf * 0)
    dts.append((time.time() - start) * 1000)
print(dts)

