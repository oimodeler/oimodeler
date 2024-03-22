from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import oimodeler as oim


# NOTE: Load the simulated data from ASPRO and apply some filter to
# keep only VIS2DATA and T3PHI for model fitting
oimodeler_dir = Path(oim.__file__).parent
files = list((oimodeler_dir / ".." / "examples" /
         "testData" / "PIONIER" / "nChannels3").glob("*.fits"))
data = oim.oimData(files)
f1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX", "OI_T3"])
data.setFilter(f1)

# # NOTE: This trend may need different computation for different wavelength grids.
wl, ks = oim.compute_photometric_slope(data, 7500)
ks = oim.oimInterp('multiParam', values=ks, wl=wl)

# NOTE: Specifying the parameter space
shgl = oim.oimStarHaloGaussLorentz(dim=128, fwhm=1, fs=0.4, fc=0.58,
                                   flor=0.4, ks=ks, kc=-3, wl0=1.68e-6)

shgl.params["kc"].set(min=-10, max=10)
shgl.params["fwhm"].set(min=0, max=32)
shgl.params["elong"].set(min=1, max=10)
shgl.params["pa"].set(min=0, max=180)
shgl.params["f"].free = False


# NOTE: Model creation
model = oim.oimModel([shgl])
params = model.getParameters()
free_params = model.getParameters(free=True)

print("Free Parameters: ")
pprint(free_params)

sim = oim.oimSimulator(data=data, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)

# NOTE: Perfoming the model-fitting
fit = oim.oimFitterEmcee(data, model, nwalkers=35)
fit.prepare(init="random")
fit.run(nsteps=int(2.5e4), progress=True)

discard = 1000
figWalkers, axeWalkers = fit.walkersPlot()
fig0, ax0 = sim.plot(["VIS2DATA"])
figCorner, axeCorner = fit.cornerPlot(discard=discard)

best, err_l, err_u, err = fit.getResults(mode="median", discard=discard)

# NOTE: Plotting images of the model
# fig, ax = plt.subplots(1, 2, figsize=(10, 5))
# model.showModel(128, 0.15, swapAxes=True, fromFT=False,
#                 normPow=1, axe=ax[0], colorbar=False)
# model.showModel(128, 0.15, swapAxes=True, fromFT=True,
#                 normPow=1, axe=ax[1], colorbar=False)
# ax[1].get_yaxis().set_visible(False)
# ax[0].set_title("Direct Image")
# ax[1].set_title("From FFT")
plt.show()
