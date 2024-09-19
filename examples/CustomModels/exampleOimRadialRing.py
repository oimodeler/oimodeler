from pathlib import Path

import matplotlib.pyplot as plt
import oimodeler as oim


# NOTE: Load the simulated data from ASPRO and apply several filters to
# keep only VIS2DATA and T3PHI, and a narrower wavelength range, for model fitting
path = Path(__file__).parent.parent.parent
files = list((path / "data" / "ASPRO_MATISSE").glob("*.fits"))
data = oim.oimData(files)
f1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX"])
f2 = oim.oimDataTypeFilter(targets="all", dataType=["T3AMP"])
f3 = oim.oimWavelengthRangeFilter(targets="all", wlRange=[3.2e-6, 3.8e-6])
data.setFilter(oim.oimDataFilter([f1, f2, f3]))

# NOTE: Model creation
star=oim.oimPt(f=0.6)
rr = oim.oimRadialRing(dim=128, din=2, dout=5, p=0.5, f=0.8)
model = oim.oimModel([star, rr])

# NOTE: Simulate and plot the initial model observables and compute the associated reduced Chi2
sim = oim.oimSimulator(data=data, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)
fig0, ax0 = sim.plot(["VIS2DATA", "T3PHI"])
print(f"Chi2r = {sim.chi2r}")

# NOTE: Specifying the parameter space
star.params["f"].set(min=0,max=1)
rr.params["din"].set(min=0, max=5)
rr.params["dout"].set(min=5.01, max=25)
rr.params["p"].set(min=-1, max=1)
rr.params["f"] = oim.oimParamNorm(star.params["f"])

# NOTE: Perfoming the model-fitting
fit = oim.oimFitterEmcee(data, model, nwalkers=25)
fit.prepare(init="random")
fit.run(nsteps=1000, progress=True)

# NOTE: Plot the walkers path and make the corner plot
figWalkers, axeWalkers = fit.walkersPlot()
figCorner, axeCorner=fit.cornerPlot(discard=200)
figSim, axSim=fit.simulator.plot(["VIS2DATA", "T3PHI"])

# NOTE: Get the best-fit reduced chi2 and best-fit values of the free parameters (+ their errors)
best, err_l, err_u, err = fit.getResults(mode="best", discard=10)
print(f"Chi2r = {fit.simulator.chi2r}")

# NOTE: Plotting images of the model
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
model.showModel(128, 0.15, swapAxes=True, fromFT=False,
                normPow=0.2, axe=ax[0], colorbar=False)
model.showModel(128, 0.15, swapAxes=True, fromFT=True,
                normPow=0.2, axe=ax[1], colorbar=False)
ax[1].get_yaxis().set_visible(False)
ax[0].set_title("Direct Image")
ax[1].set_title("From FFT")
plt.show()
