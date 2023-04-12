from pathlib import Path

import matplotlib.pyplot as plt
import oimodeler as oim


# NOTE: Load the simulated data from ASPRO and apply some filter to
# keep only VIS2DATA and T3PHI for model fitting
oimodeler_dir = Path(oim.__file__).parent
files = (oimodeler_dir / ".." / "examples" /
         "testData" / "ASPRO_MATISSE2").glob("*.fits")
files = list(map(str, files))
data = oim.oimData(files)
f1 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX"])
f2 = oim.oimDataTypeFilter(targets="all", dataType=["T3AMP"])
data.setFilter(oim.oimDataFilter([f1, f2]))

# NOTE: Specifying the parameter space
atg = oim.oimAsymTempGradient(dim=128, dist=140, kappa_abs=276,
                              Tin=1500, rin=0.5, rout=3, p=0.5, q=0.5,
                              a=0.5, phi=45, Mdust=0.11, pa=30, elong=2, f=0.8)

atg.params["rin"].set(min=0, max=20)
atg.params["rout"].set(min=10, max=20)
atg.params["q"].set(min=0, max=1)
atg.params["p"].set(min=0, max=1)
atg.params["Mdust"].set(min=0, max=3)
atg.params["a"].set(min=0, max=1)
atg.params["phi"].set(min=0, max=360)
atg.params["elong"].set(min=1, max=50)
atg.params["pa"].set(min=0, max=360)
atg.params["f"].free = False
atg.params["pixSize"].free = False

# NOTE: Model creation
model = oim.oimModel([atg])

sim = oim.oimSimulator(data=data, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)

# NOTE: Perfoming the model-fitting
fit = oim.oimFitterEmcee(data, model, nwalkers=25)
fit.prepare(init="random")
fit.run(nsteps=100, progress=True)

figWalkers, axeWalkers = fit.walkersPlot()
fig0, ax0 = sim.plot(["VIS2DATA", "T3PHI"])

best, err_l, err_u, err = fit.getResults(mode='best', discard=10)

# NOTE: Plotting images of the model
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
model.showModel(128, 0.15, swapAxes=True, fromFT=False,
                normPow=1, axe=ax[0], colorbar=False)
model.showModel(128, 0.15, swapAxes=True, fromFT=True,
                normPow=1, axe=ax[1], colorbar=False)
ax[1].get_yaxis().set_visible(False)
ax[0].set_title("Direct Image")
ax[1].set_title("From FFT")
plt.show()
