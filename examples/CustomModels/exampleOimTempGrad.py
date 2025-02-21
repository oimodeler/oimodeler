from pathlib import Path

import matplotlib.pyplot as plt
import oimodeler as oim

# NOTE: Load the simulated data from ASPRO and apply some filter to
# keep only VIS2DATA and T3PHI for model fitting
path = Path(__file__).parent.parent.parent
files = list((path / "data" / "ASPRO_MATISSE2").glob("*.fits"))

wavelengths = [8e-6, 10e-6]

data = oim.oimData(files)
f1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX"])
f2 = oim.oimDataTypeFilter(targets="all", dataType=["T3AMP"])
data.setFilter(oim.oimDataFilter([f1, f2]))

# NOTE: Grid can be changed from 'linear' to 'logarithmic'
oim.oimOptions.model.grid.type = "logarithmic"

# NOTE: The padding of the 1D-grid can be set (Multiplies to outer radius at
# grid creation). Default is 'None/1' and doesn't pad
oim.oimOptions.ft.padding = 4

# NOTE: An interpolator that contains different values for different wavelengths
kappa_abs = oim.oimInterp("wl", wl=wavelengths, values=[400e1, 300e-5])

# NOTE: Define a point source with the flux from the star (blackbody)
# interpolator
s = oim.oimPt(
    f=oim.oimInterp("starWl", temp=7500, dist=146, lum=9, wl=wavelengths)
)

# NOTE: Define the temperature gradient
tg = oim.oimTempGrad(
    dim=128,
    dist=146,
    kappa_abs=kappa_abs,
    temp0=1500,
    rin=0.5,
    rout=5,
    p=0.5,
    q=0.5,
    dust_mass=1e-3,
    pa=30,
    elong=1,
)

# NOTE: Specifying the parameter space.
tg.params["rin"].set(min=0, max=20)
tg.params["rout"].set(min=3, max=20)
tg.params["q"].set(min=0, max=1)
tg.params["p"].set(min=0, max=1)
tg.params["dust_mass"].set(min=0, max=3)
tg.params["elong"].set(min=1, max=50)
tg.params["pa"].set(min=0, max=360)
tg.params["f"].free = False

# NOTE: Set the wavelengths to be fitted
# WARN: Needs to be the same for all components.
s._wl = tg._wl = wavelengths

# NOTE: Model creation.
model = oim.oimModel([s, tg])

sim = oim.oimSimulator(data=data, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)
print("Chi2r = {}".format(sim.chi2r))

# # NOTE: Perfoming the model-fitting
fit = oim.oimFitterEmcee(data, model, nwalkers=25)
fit.prepare(init="random")
fit.run(nsteps=100, progress=True)

figWalkers, axeWalkers = fit.walkersPlot()
fig0, ax0 = sim.plot(["VIS2DATA", "T3PHI"])

best, err_l, err_u, err = fit.getResults(mode="best", discard=10)

# NOTE: Plotting images of the model
tg._wl = None
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
model.showModel(
    128,
    0.35,
    wl=8e-6,
    swapAxes=True,
    fromFT=False,
    normPow=0.3,
    axe=ax[0],
    colorbar=False,
)
model.showModel(
    128,
    0.35,
    wl=8e-6,
    swapAxes=True,
    fromFT=True,
    normPow=0.3,
    axe=ax[1],
    colorbar=False,
)
ax[1].get_yaxis().set_visible(False)
ax[0].set_title("Direct Image")
ax[1].set_title("From FFT")
plt.show()
