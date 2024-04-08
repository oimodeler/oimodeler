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
data_dir = path / "data" / "ASPRO_MATISSE2"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

files = list(data_dir.glob("*.fits"))

# Building a oimodeler model with the same parameters
ud = oim.oimUD(d=3, f=0.5)
pt = oim.oimPt(f=1)
model = oim.oimModel([ud, pt])

# Setting limits of the parameter space but also setting x,y of the UD as a
# free parameters and f of the pt as fixed (to 1)
ud.params['d'].set(min=0.01, max=20)
ud.params['x'].set(min=-50, max=50, free=True)
ud.params['y'].set(min=-50, max=50, free=True)
ud.params['f'].set(min=0., max=10.)
pt.params['f'].free = False
pprint(model.getFreeParameters())

# Create a new fitter with 32 walkers and the list of oifits files and the model
fit = oim.oimFitterEmcee(files, model, nwalkers=32, ncores=6)
# pprint(fit._logProbability([0,10,1,5]))


# Prepare the fitter.Here we ste the intial positions of all walkers to
# the current parameters values of our model.
fit.prepare(init="random")

# pprinting the initial values of the walkers
pprint("Initial values of the free parameters for the {} walkers".format(
    fit.params["nwalkers"].value))
pprint(fit.initialParams)

# run a 1000 steps fit with fixed starting inital and 1000 steps
fit.run(nsteps=2000, progress=True)

# %%
sampler = fit.sampler
chain = fit.sampler.chain
lnprob = fit.sampler.lnprobability

# %%
figWalkers, axeWalkers = fit.walkersPlot(cmap="plasma_r",
                                         savefig=save_dir / "exampleOimFitterEmceeWalkers.png")
figCorner, axeCorner = fit.cornerPlot(discard=1000,
                                      savefig=save_dir / "exampleOimFitterEmceeCorner.png")

# %%
median, err_l, err_u, err = fit.getResults(mode='median', discard=1000, chi2limfact=20)

# %%
fig0, ax0 = fit.simulator.plot(["VIS2DATA", "VISAMP", "VISPHI", "T3AMP", "T3PHI"],
                               savefig=save_dir / "ExampleOimFitterEmcee_fittedData.png")
