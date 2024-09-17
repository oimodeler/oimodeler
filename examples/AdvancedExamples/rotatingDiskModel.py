# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:00:30 2023

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


# as the model as no sharp outer edge, no zero-padding is needed
oim.oimOptions.ft.padding = 1
#oim.oimOptions.ft.backend.active = oim.FFTWBackend

# Path to the AMBER oifits files: the classical Be star Alpha Col
path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "RealData" / "AMBER" / "AlphaCol"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

files = list(data_dir.glob("*.fits"))

# %%
dim, nwl, wl0 = 256, 51, 2.1656e-6
c = oim.oimKinematicDisk(dim=dim, 
                         fov=20, 
                         incl=45, 
                         Rstar=5.8, 
                         dist=80, 
                         fwhmCont=2.2,
                         fluxDiskCont=0.25, 
                         EW=10.4, 
                         fwhmLine=5.7, 
                         nwl=nwl,
                         vrot=360, 
                         beta=-0.5, 
                         pa=-83.2, 
                         wl0=2.1656e-6, 
                         dwl=0.9e-10,
                         res=1.8e-10)
m = oim.oimModel(c)

# %% Showing 5 images in narrow bands through the disk
dwl, nwl = 32e-10, 5
wl = np.linspace(wl0-dwl/2, wl0+dwl/2, num=nwl)
fig0, _, _ = m.showModel(dim, 0.05, wl=wl, legend=True, normalize=True,
            fromFT=False, normPow=0.2, colorbar=False,
            savefig=save_dir / "ExampleRotatingDiskModel_Model_images.png")

# %% Simulating data and comparing data and model of the disk around the
# classical Be star Alpha Col
sim = oim.oimSimulator(files, m)
pprint(f"chi2: {sim.chi2r}")

#%%
fig=plt.figure(FigureClass=oim.oimWlTemplatePlots,figsize=(12, 7))
fig.autoShape(sim.data.data,shape=[["VIS2DATA",None],["VISPHI","T3PHI"]])
fig.set_xunit("micron")
fig.plot(sim.data.data,plotFunction=plt.Axes.errorbar,
         plotFunctionkwarg=dict(color="tab:red",alpha=0.5))
fig.plot(sim.simulatedData.data,plotFunctionkwarg=dict(color="tab:blue"))
fig.set_ylim(["VISPHI","T3PHI"],-25,25)
fig.set_ylim(["VIS2DATA"],0,1.2)
fig.set_xlim(2.1624,2.1688)
fig.set_legends(0.5,0.1,"$BASELINE$ $LENGTH$m $PA$$^o$",["VIS2DATA","VISPHI"],
                fontsize=12,ha="center")
fig.set_legends(0.5,0.1,"$BASELINE$",["T3PHI"],fontsize=12,ha="center")
title = f"$\\alpha$ Col AMBER data + Kinematic disk model: chi$^2_r$= {sim.chi2r:.1f}"
fig.suptitle(title)

fig.tight_layout()

plt.savefig(save_dir / "ExampleRotatingDiskModel_data_model.png")

