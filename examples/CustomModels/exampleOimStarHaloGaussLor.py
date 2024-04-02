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


def _logProbability(self, theta: np.ndarray) -> float:
    """The log probability so it keeps fs+fc<=1"""
    theta[-2] *= 1-theta[-3]
    for iparam, parami in enumerate(self.freeParams.values()):
        parami.value = theta[iparam]

    self.simulator.compute(computeChi2=True, dataTypes=self.dataTypes)
    return -0.5 * self.simulator.chi2r


# NOTE: Load the simulated data from ASPRO and apply some filter to
# keep only VIS2DATA for model fitting
oimodeler_dir = Path(oim.__file__).parent
files = list((oimodeler_dir / ".." / "examples" /
         "data" / "RealData" / "PIONIER" / "nChannels3").glob("*.fits"))
data = oim.oimData(files)
f1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX", "OI_T3"])
data.setFilter(f1)

# # NOTE: The calculation of the photometric slope from the star's effective temperature
wl, ks = oim.compute_photometric_slope(data, 7500)
ks = oim.oimInterp('multiParam', values=ks, wl=wl)

# NOTE: Specifying the parameter space
shgl = oim.oimStarHaloGaussLorentz(dim=128, fwhm=2*1.0715, fs=0.44, fc=0.54,
                                   pa=162, elong=2, flor=0.41, ks=ks,
                                   kc=-3.64, wl0=1.68e-6)

shgl.params["kc"].set(min=-10, max=10)
shgl.params["fwhm"].set(min=0, max=32)
shgl.params["elong"].set(min=1, max=10)
shgl.params["f"].free = False


# NOTE: Model creation
model = oim.oimModel([shgl])
params = model.getParameters()
free_params = model.getParameters(free=True)

print("Free Parameters: ")
pprint(free_params)

sim = oim.oimSimulator(data=data, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)
fig0, ax0 = sim.plot(["VIS2DATA"], savefig=save_dir / f"ExampleOim{shgl.shortname}_prefit.png")

# NOTE: Perfect parameter chi_sqr is a bit higher than Lazareff+2017 (chi_sqr = 1.53)
print("Pre-fit Chi2r (with Lazareff+2017 best-fit params): ", sim.chi2r)

# NOTE: Perfoming the model-fitting
fit = oim.oimFitterDynesty(data, model, nlive=1500)

# NOTE: Overwrite the existing _logProbability method of the class "fit"
fit._logProbability = lambda theta: _logProbability(fit, theta)

fit.prepare(init="random")
fit.run(dlogz=0.010, progress=True)

best, err_l, err_u, err = fit.getResults(mode="median")
figWalkers, axeWalkers = fit.walkersPlot(savefig=save_dir / f"ExampleOim{shgl.shortname}_walkers.png")
figCorner, axeCorner = fit.cornerPlot(savefig=save_dir / f"ExampleOim{shgl.shortname}_corner.png")
fig0, ax0 = fit.simulator.plot(["VIS2DATA"], savefig=save_dir / f"ExampleOim{shgl.shortname}_fit.png")
print("Chi2r: ", fit.simulator.chi2r)
plt.close()

print("Free Parameters: ")
pprint(fit.simulator.model.getParameters(free=True))

# NOTE: Plotting images of the model
model.showModel(128, 0.15, swapAxes=True, fromFT=False,
                normPow=1, colorbar=False)
plt.savefig(save_dir / f"ExampleOim{shgl.shortname}_images.png")
