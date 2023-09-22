# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:59:15 2022

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import oimodeler as oim


# Path to a fake MATISSE-L-band binary observation (3 oifits) created with ASPRO
path = Path(__file__).parent.parent.parent
data_dir = path / "examples" / "testData" / "ASPRO_MATISSE2"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# TODO: After pathlib change of all `oimodeler` modules, remove str casting.
files = list(map(str, data_dir.glob("*.fits")))

# %%
ud = oim.oimUD(d=3, f=0.5, x=5, y=-5)
pt = oim.oimPt(f=1)

# %%
pprint(ud)
pprint(ud.params['d'])
pprint(ud.params['x'])

# %%
ud.params['d'].set(min=0.01, max=20)
ud.params['x'].set(min=-50, max=50, free=True)
ud.params['y'].set(min=-50, max=50, free=True)
ud.params['f'].set(min=0., max=10.)
pt.params['f'].free = False

# %%
model = oim.oimModel(ud, pt)

# %%
pprint(model.getParameters())
pprint(model.getFreeParameters())

# %%
sim = oim.oimSimulator(data=files, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)

# %%
pprint(f"Chi2r = {sim.chi2r}")


fig0, ax0 = sim.plot(["VIS2DATA", "T3PHI"],
                     savefig=save_dir / "gettingStarted_model0.png")

# %%

fit = oim.oimFitterEmcee(files, model, nwalkers=20)
fit.prepare(init="random")
pprint(fit.initialParams)
fit.run(nsteps=2000, progress=True)

# %%
figWalkers, axeWalkers = fit.walkersPlot(
    savefig=save_dir / "gettingStarted_Walkers.png")

figCorner, axeCorner = fit.cornerPlot(discard=1000,
                                      savefig=save_dir / "gettingStarted_corner.png")

# %%
median, err_l, err_u, err = fit.getResults(mode='median', discard=1000)

figSim, axSim = fit.simulator.plot(["VIS2DATA", "T3PHI"],
                                   savefig=save_dir / "gettingStarted_modelFinal.png")
pprint(f"Chi2r = {fit.simulator.chi2r}")

figImg, axImg, im = model.showModel(512, 0.1, normPow=0.1, figsize=(6, 4.8),
                                    savefig=save_dir / "gettingStarted_modelImage.png")
