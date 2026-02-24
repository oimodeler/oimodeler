from pathlib import Path
from pprint import pprint

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

import oimodeler as oim


def _logProbability(self, theta: np.ndarray) -> float:
    """The log probability for emcee so it keeps (fs + fc) <= 1"""
    keys_free_params = list(self.freeParams.keys())
    indices = tuple(
        [
            index
            for index, element in enumerate(keys_free_params)
            if "fs" in element or "fc" in element or "fh" in element
        ]
    )

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


# NOTE: Change this path if you want to save the products at another location
path = Path(__file__).parent.parent.parent
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# NOTE: Load some PIONIER archive data and apply some filters to it
# to keep only VIS2DATA and T3PHI
# TODO: Examples files don't fit very well to this model
path = Path(oim.__file__).parent.parent
files = list(
    (path / "data" / "RealData" / "PIONIER" / "HD142527").glob("*.fits")
)
data = oim.oimData(files)
f1 = oim.oimWavelengthRangeFilter(targets="all", wlRange=[1.58e-6, 1.78e-6])
f2 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX"])
f3 = oim.oimDataTypeFilter(targets="all", dataType=["T3AMP"])
data.setFilter(oim.oimDataFilter([f1, f2, f3]))

# NOTE: The calculation of the photometric slope from the star's effective temperature
wl, ks = oim.compute_photometric_slope(data, 6500)
ks = oim.oimInterp("wl", values=ks, wl=wl, kind="linear", extrapolate=False)

# NOTE: Specifying the parameter space (best-fit parameters from Lazareff+2017)
shglr = oim.oimStarHaloIRing(
    fs=0.42,
    fc=0.55,
    fh=0.03,
    flor=0.42,
    pa=0.09 * u.rad.to(u.deg),
    elong=1 / 0.84,
    la=0.08,
    lkr=0.54,
    ks=ks,
    kc=-3.67,
    skw=0.8275868534480233,
    skwPa=46.46880071438582,
    wl0=1.68e-06,
)

shglr.params["elong"].set(min=1, max=3)
shglr.params["pa"].set(min=-180, max=180)
shglr.params["fs"].set(min=0.3, max=0.5)
shglr.params["fc"].set(min=0.4, max=0.6)
shglr.params["fh"] = oim.oimParamNorm([shglr.params["fs"], shglr.params["fc"]])
shglr.params["f"].free = False

# NOTE: Model creation
model = oim.oimModel([shglr])
params = model.getParameters()
free_params = model.getParameters(free=True)

print("Free Parameters: ")
pprint(free_params)

sim = oim.oimSimulator(data=data, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)
fig0, ax0 = sim.plot(
    ["VIS2DATA", "T3PHI"],
    savefig=save_dir / f"ExampleOim{shglr.shortname}_prefit.png",
)

# NOTE: Perfect parameter chi_sqr is a bit higher than Lazareff+2017 (chi_sqr = 1.53)
print("Pre-fit Chi2r (with Lazareff+2017 best-fit params): ", sim.chi2r)

# NOTE: Perfoming the model-fitting
fit = oim.oimFitterEmcee(data, model, nwalkers=50)

# NOTE: Overwrite the existing _logProbability method of the class "fit"
fit._logProbability = lambda theta: _logProbability(fit, theta)

nsteps = 10000
discard = int(nsteps * 0.1)
fit.prepare(init="random")
fit.run(nsteps=nsteps, progress=True)

best, err_l, err_u, err = fit.getResults(mode="median")
figWalkers, axeWalkers = fit.walkersPlot(
    savefig=save_dir / f"ExampleOim{shglr.shortname}_walkers.png"
)
figCorner, axeCorner = fit.cornerPlot(
    savefig=save_dir / f"ExampleOim{shglr.shortname}_corner.png"
)
fig0, ax0 = fit.simulator.plot(
    ["VIS2DATA", "T3PHI"],
    savefig=save_dir / f"ExampleOim{shglr.shortname}_fit.png",
)
print("Chi2r: ", fit.simulator.chi2r)

print("Free Parameters: ")
pprint(fit.simulator.model.getParameters(free=True))

# NOTE: Plotting images of the model
model.showModel(
    512,
    0.02,
    swapAxes=True,
    fromFT=False,
    wl=1.68e-6,
    normPow=0.5,
    colorbar=False,
)
plt.savefig(save_dir / f"ExampleOim{shglr.shortname}_images.png")

model.showModel(
    512,
    0.02,
    swapAxes=True,
    fromFT=True,
    wl=1.68e-6,
    normPow=0.5,
    colorbar=False,
)
plt.savefig(save_dir / f"ExampleOim{shglr.shortname}_images_from_ft.png")
