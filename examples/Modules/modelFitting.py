# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:59:15 2022

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import oimodeler as oim

#%%
# Path tothe binary star  Beta Ari observed with MIRCX (1 files)
path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "RealData" "/MIRCX" / "Beta Ari"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

file = data_dir / "MIRCX_L2.2023Oct14._bet_Ari.MIRCX_IDL.nn.AVG10m.fits"

# Building a oimodeler model with the same parameters
ud1 = oim.oimUD()
ud2 = oim.oimUD()
model = oim.oimModel([ud1, ud2])
#model.normalizeFlux()
# Setting limits of the parameter space but also setting x,y of the UD as a
# free parameters and f of the pt as fixed (to 1)
ud1.params["d"].set(min=0, max=2)
ud1.params["f"].set(min=0.8, max=1)
ud2.params["d"].set(min=0, max=2)
ud2.params["x"].set(min=-100, max=100, free=True)
ud2.params["y"].set(min=-100, max=100, free=True)
model.normalizeFlux()
pprint(model.getFreeParameters())


data=oim.oimData(file)
filt=oim.oimRemoveArrayFilter(arr="OI_VIS")
data.setFilter(filt)
#%%
# Create a new fitter with 32 walkers and the list of oifits files and the model
fit = oim.oimFitterEmcee(data, model, nwalkers=20,dataTypes=["VIS2DATA","T3PHI"])
# pprint(fit._logProbability([0,10,1,5]))

# Prepare the fitter. Here we set the intial positions of all walkers to
# the current parameters values of our model.
fit.prepare(init="random")

# pprinting the initial values of the walkers
pprint(fit.initialParams)

# run a 1000 steps fit with fixed starting inital and 1000 steps
fit.run(nsteps=40000, progress=True)
#%%
# Get results from the fit (updates the class internal logic)
best, err_l, err_u, err = fit.getResults(mode="best", discard=30000, chi2limfact=3)

# %%
sampler = fit.sampler
chain = fit.sampler.chain
lnprob = fit.sampler.lnprobability

# %%
figWalkers, axeWalkers = fit.walkersPlot(chi2limfact=20)
figCorner, axeCorner = fit.cornerPlot(discard=30000, chi2limfact=3)
#%%
fit2 = oim.oimFitterEmcee(data, model, nwalkers=20,dataTypes=["VIS2DATA","T3PHI"])
# pprint(fit._logProbability([0,10,1,5]))

# Prepare the fitter. Here we set the intial positions of all walkers to
# the current parameters values of our model.
fit2.prepare(init="gaussian")
fit2.run(nsteps=5000, progress=True)
# %%
figWalkers, axeWalkers = fit2.walkersPlot(chi2limfact=5)
figCorner, axeCorner = fit2.cornerPlot(discard=2000)

# %%
fig0, ax0 = fit.simulator.plot(["VIS2DATA", "T3PHI"])#kwargsData=dict(marker=".",ls=""),
                              # savefig=save_dir / f"Example{class_name}_fittedData.png")
#%%
#%% residual plot
fig2, ax2 = fit.simulator.plot_residuals(["VIS2DATA", "T3PHI"],levels=[1,2])
                     #savefig=save_dir / "ExampleOimSimulator_residuals_plot.png")
                     
#%%
fig2, ax2 = fit.simulator.plotWithResiduals(["VIS2DATA", "T3PHI"],levels=[1,2,3],
                                            xunit="cycle/mas",
                                            kwargsData=dict(color="byBaseline"
                                                            ,marker="."))
                     #savefig=save_dir / "ExampleOimSimulator_residuals_plot.png")
ax2[0].legend(fontsize=6)    
#%%

lmfit = oim.oimFitterMinimize(data, model,dataTypes=["VIS2DATA", "T3PHI"])
lmfit.prepare()
lmfit.run()

