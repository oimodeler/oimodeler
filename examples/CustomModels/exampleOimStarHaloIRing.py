from pathlib import Path
from pprint import pprint

import astropy.units as u
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
# keep only VIS2DATA and T3PHI for model fitting
oimodeler_dir = Path(oim.__file__).parent
files = list((oimodeler_dir / ".." / "examples" /
         "testData" / "PIONIER" / "nChannels3").glob("*.fits"))
data = oim.oimData(files)
f1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX"])
f2 = oim.oimDataTypeFilter(targets="all", dataType=["T3AMP"])
data.setFilter(oim.oimDataFilter([f1, f2]))

# # NOTE: The calculation of the photometric slope from the star's effective temperature
wl, ks = oim.compute_photometric_slope(data, 7500)
ks = oim.oimInterp('multiParam', values=ks, wl=wl)

# NOTE: Specifying the parameter space (best-fit parameters from Lazareff+2017)
shglr = oim.oimStarHaloIRing(fs=0.45, fc=0.53,
                             pa=2.84*u.rad.to(u.deg), elong=1/0.5,
                             rin=0.6381143066747668, fwhm=2*0.8607925145549679,
                             flor=0.58, ks=ks, a=0.04, phi=180,
                             kc=-4.18, wl0=1.68*u.um.to(u.m))

shglr.params["kc"].set(min=-10, max=10)
shglr.params["rin"].set(min=0, max=32)
shglr.params["fwhm"].set(min=0, max=32)
shglr.params["elong"].set(min=1, max=10)
shglr.params["pa"].set(min=10, max=180)
shglr.params["phi"].set(min=10, max=180)
shglr.params["f"].free = False


# NOTE: Model creation
model = oim.oimModel([shglr])
params = model.getParameters()
free_params = model.getParameters(free=True)

print("Free Parameters: ")
pprint(free_params)

sim = oim.oimSimulator(data=data, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)
fig0, ax0 = sim.plot(["VIS2DATA", "T3PHI"], savefig=save_dir / f"ExampleOim{shglr.shortname}_prefit.png")

# NOTE: Perfect parameter chi_sqr is a bit higher than Lazareff+2017 (chi_sqr = 1.53)
print("Pre-fit Chi2r (with Lazareff+2017 best-fit params): ", sim.chi2r)

# NOTE: Perfoming the model-fitting
fit = oim.oimFitterDynesty(data, model, method="dynamic")

# NOTE: Overwrite the existing _logProbability method of the class "fit"
fit._logProbability = lambda theta: _logProbability(fit, theta)

# fit.prepare(nlive=2000)
# fit.run(progress=True)

# best, err_l, err_u, err = fit.getResults(mode="median")
# figWalkers, axeWalkers = fit.walkersPlot(savefig=save_dir / f"ExampleOim{shglr.shortname}_walkers.png")
# figCorner, axeCorner = fit.cornerPlot(savefig=save_dir / f"ExampleOim{shglr.shortname}_corner.png")
# fig0, ax0 = fit.simulator.plot(["VIS2DATA", "T3PHI"], savefig=save_dir / f"ExampleOim{shglr.shortname}_fit.png")
# print("Chi2r: ", fit.simulator.chi2r)

# print("Free Parameters: ")
# pprint(fit.simulator.model.getParameters(free=True))

# NOTE: Plotting images of the model
model.showModel(128, 0.15, swapAxes=True, fromFT=False,
                normPow=1, colorbar=False)
plt.savefig(save_dir / f"ExampleOim{shglr.shortname}_images.png")
