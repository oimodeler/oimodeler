from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


# NOTE: Change this path if you want to save the products at another location
path = Path(__file__).parent.parent.parent
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


# def _logProbability(self, theta: np.ndarray) -> float:
#     """The log probability for dynesty so it keeps (fs + fc) <= 1"""
#     theta[-2] *= 1-theta[-3]
#     for iparam, parami in enumerate(self.freeParams.values()):
#         parami.value = theta[iparam]

#     self.simulator.compute(computeChi2=True, dataTypes=self.dataTypes)
#     return -0.5 * self.simulator.chi2r


def _logProbability(self, theta: np.ndarray) -> float:
    """The log probability for emcee so it keeps (fs + fc) <= 1"""
    keys_free_params = list(self.freeParams.keys())
    indices = tuple([index for index, element in enumerate(keys_free_params)
                     if "fs" in element or "fc" in element or "fh" in element])

    for iparam, parami in enumerate(self.freeParams.values()):
        parami.value = theta[iparam]

    if sum(theta[np.array(indices)]) > 1:
        return -np.inf

    for val, key in zip(theta, self.freeParams):
        low, up = self.limits[key]
        if not low < val < up:
            return -np.inf

    self.simulator.compute(computeChi2=True, dataTypes=self.dataTypes)
    return -0.5 * self.simulator.chi2r


# NOTE: Load the simulated data from ASPRO and apply some filter to
# keep only VIS2DATA for model fitting
path = Path(oim.__file__).parent.parent
files = list((path / "data" / "RealData" / "PIONIER" / "nChannels3").glob("*.fits"))
data = oim.oimData(files)
f1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX", "OI_T3"])
data.setFilter(f1)

# # NOTE: The calculation of the photometric slope from the star's effective temperature
wl, ks = oim.compute_photometric_slope(data, 7500)
ks = oim.oimInterp("wl", values=ks, wl=wl, kind="linear", extrapolate=False)

# NOTE: Specifying the parameter space
shgl = oim.oimStarHaloGaussLorentz(fwhm=2.143038610475213,
                                   fs=0.44, fc=0.54, fh=0.02,
                                   flor=0.41, kc=-3.64, ks=ks,
                                   pa=162.72001381715378,
                                   elong=2, wl0=1.68e-6)


shgl.params["kc"].set(min=-10, max=10)
shgl.params["fwhm"].set(min=0, max=32)
shgl.params["pa"].set(min=0, max=180)
shgl.params["elong"].set(min=1, max=10)
shgl.params["f"].free = False


# NOTE: Model creation
model = oim.oimModel([shgl])
params = model.getParameters()

print("Free Parameters:")
pprint(model.getParameters(free=True))

sim = oim.oimSimulator(data=data, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)
fig0, ax0 = sim.plot(["VIS2DATA"], savefig=save_dir / f"ExampleOim{shgl.shortname}_prefit.png")

# NOTE: Perfect parameter chi_sqr is a bit higher than Lazareff+2017 (chi_sqr = 1.53)
print("Pre-fit Chi2r (with Lazareff+2017 best-fit params):", sim.chi2r)

# NOTE: Perfoming the model-fitting
fit = oim.oimFitterEmcee(data, model, nwalkers=50)

# NOTE: Overwrite the existing _logProbability method of the class "fit"
fit._logProbability = lambda theta: _logProbability(fit, theta)

fit.prepare(init="random")
fit.run(nsteps=20000, progress=True)

best, err_l, err_u, err = fit.getResults(mode="median")
figWalkers, axeWalkers = fit.walkersPlot(savefig=save_dir / f"ExampleOim{shgl.shortname}_walkers.png")
figCorner, axeCorner = fit.cornerPlot(savefig=save_dir / f"ExampleOim{shgl.shortname}_corner.png")
fig0, ax0 = fit.simulator.plot(["VIS2DATA"], savefig=save_dir / f"ExampleOim{shgl.shortname}_fit.png")
print("Chi2r:", fit.simulator.chi2r)
plt.close()

print("Free Parameters:")
pprint(fit.simulator.model.getParameters(free=True))

# NOTE: Plotting images of the model
model.showModel(512, 0.02, swapAxes=True, fromFT=False,
                wl=1.68e-6, normPow=0.5, colorbar=False)
plt.savefig(save_dir / f"ExampleOim{shgl.shortname}_images.png")

model.showModel(512, 0.02, swapAxes=True, fromFT=True,
                wl=1.68e-6, normPow=0.5, colorbar=False)
plt.savefig(save_dir / f"ExampleOim{shgl.shortname}_images_from_ft.png")
