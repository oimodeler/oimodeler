from pathlib import Path
from pprint import pprint

import astropy.units as u
import matplotlib.pyplot as plt
import oimodeler as oim


# NOTE: Change this path if you want to save the products at another location
path = Path(__file__).parent.parent.parent
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


# NOTE: Load the simulated data from ASPRO and apply some filter to
# keep only VIS2DATA for model fitting
# TODO: Examples files don't fit very well to this model
path = Path(oim.__file__).parent.parent
files = list((path / "data" / "RealData" / "PIONIER" / "HD142527").glob("*.fits"))
data = oim.oimData(files)
f1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX", "OI_T3"])
data.setFilter(f1)

# NOTE: Specifying the parameter space
gl = oim.oimGaussLorentz(hlr=10**0.06, flor=0.43,
                         pa=0.12*u.rad.to(u.deg), elong=1/0.83)

gl.params["hlr"].set(min=0, max=32)
gl.params["elong"].set(min=1, max=4)
gl.params["f"].free = False

# NOTE: Model creation
model = oim.oimModel([gl])
params = model.getParameters()
free_params = model.getParameters(free=True)

print("Free Parameters: ")
pprint(free_params)

sim = oim.oimSimulator(data=data, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)
fig0, ax0 = sim.plot(["VIS2DATA"], savefig=save_dir / f"ExampleOim{gl.shortname}_prefit.png")

print("Pre-fit Chi2r: ", sim.chi2r)

# NOTE: Perfoming the model-fitting
fit = oim.oimFitterEmcee(data, model, nwalkers=32)

fit.prepare(init="random")
fit.run(nsteps=2000, progress=True)

best, err_l, err_u, err = fit.getResults(mode="median")
figWalkers, axeWalkers = fit.walkersPlot(savefig=save_dir / f"ExampleOim{gl.shortname}_walkers.png")
figCorner, axeCorner = fit.cornerPlot(savefig=save_dir / f"ExampleOim{gl.shortname}_corner.png")
fig0, ax0 = fit.simulator.plot(["VIS2DATA"], savefig=save_dir / f"ExampleOim{gl.shortname}_fit.png")
print("Chi2r: ", fit.simulator.chi2r)
plt.close()

print("Free Parameters: ")
pprint(fit.simulator.model.getParameters(free=True))

# NOTE: Plotting images of the model
model.showModel(512, 0.02, swapAxes=True, fromFT=False, normPow=0.5, colorbar=False)
plt.savefig(save_dir / f"ExampleOim{gl.shortname}_images.png")

model.showModel(512, 0.02, swapAxes=True, fromFT=True, normPow=0.5, colorbar=False)
plt.savefig(save_dir / f"ExampleOim{gl.shortname}_images_from_ft.png")
