# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 11:53:04 2022

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import numpy as np
import oimodeler as oim


path = Path(__file__).parent.parent.parent
data_dir = path / "examples" / "testData" / "ASPRO_CHROMATIC_SKWDISK"

# NOTE: Change these paths if you want to save the products at another location
save_dir = path / "images"
product_dir = path / "examples" / "testData" / "IMAGES"
if not save_dir.exists():
    save_dir.mkdir(parents=True)
if not product_dir.exists():
    product_dir.mkdir(parents=True)

# %% Model creation
star = oim.oimUD(d=1, f=oim.oimInterp(
    "wl", wl=[3e-6, 4e-6], values=[0.5, 0.1]))
ring = oim.oimESKRing(din=8, dout=oim.oimInterp(
    "wl", wl=[3e-6, 4e-6], values=[9, 14]), elong=1.5, skw=0.8, pa=50)
ring.params['f'] = oim.oimParamNorm(star.params['f'])
ring.params["skwPa"] = oim.oimParamLinker(ring.params["pa"], "add", 90)
model = oim.oimModel(star, ring)

# %% Showing images at the reference walvnegths (3 and 4 microns) and also at the  3.5microns
model.showModel(256, 0.06, wl=np.linspace(3., 4, num=3)*1e-6, legend=True, normalize=True, normPow=1, fromFT=True,
                savefig=save_dir / "chromaticModelFitImageInit.png")

# %% Exporting a 50 wl image cube for ASPRO
img = model.getImage(256, 0.06, toFits=True, fromFT=True,
                     wl=np.linspace(3, 4, num=50)*1e-6)
img.writeto(product_dir / "skwDisk.fits", overwrite=True)

# %% Load the simulated data from ASPRO and apply some filter to keep only
# VIS2DATA and T3PHI for model fitting
files = list(data_dir.glob("*.fits"))
data = oim.oimData(files)
f1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX"])
f2 = oim.oimDataTypeFilter(targets="all", dataType=["T3AMP"])
data.setFilter(oim.oimDataFilter([f1, f2]))

# %%
sim = oim.oimSimulator(data, model)
pprint(sim.chi2r)
sim.plot(["VIS2DATA", "T3PHI"])

# %% Specifying the parameter space
params = model.getFreeParameters()
pprint(params)

params['c1_UD_f_interp1'].set(min=0.0, max=1)
params['c1_UD_f_interp2'].set(min=-0.0, max=1)
params['c1_UD_d'].set(min=0, max=5, free=True)
params['c2_SKER_pa'].set(min=0., max=180)
params['c2_SKER_elong'].set(min=1, max=3)
params['c2_SKER_din'].set(min=5, max=20.)
params['c2_SKER_skw'].set(min=0, max=1.5)
params['c2_SKER_dout_interp1'].set(min=5., max=30.)
params['c2_SKER_dout_interp2'].set(min=5., max=30.)

# %% perfoming the model-fitting
fit = oim.oimFitterEmcee(data, model, nwalkers=50)
fit.prepare(init="random")
fit.run(nsteps=3000, progress=True)

# %%
figWalkers, axeWalkers = fit.walkersPlot(
    savefig=save_dir / "chromaticModelFitWalkers.png")

# %%
figCorner, axeCorner = fit.cornerPlot(discard=2500, chi2limfact=3,
                                      savefig=save_dir / "chromaticModelFitCorner.png")

# %%
best, err_l, err_u, err = fit.getResults(mode='best', discard=1500)

figSim, axSim = fit.simulator.plot(["VIS2DATA", "T3PHI"],
                                   savefig=save_dir / "chromaticModelFitVisCP.png")
pprint(f"Chi2r = {fit.simulator.chi2r}")

# %%
figImg, axImg, im = model.showModel(256, 0.06, wl=[3e-6, 4e-6], normPow=1, legend=True, normalize=True, fromFT=True,
                                    savefig=save_dir / "chromaticModelFitImageFinal.png")
