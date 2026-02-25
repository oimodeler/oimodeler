from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants as const

import oimodeler as oim

# NOTE: Load MATISSE data of AS209 (YSO) into an oimData object
path = Path(__file__).parent.parent.parent
files = sorted((path / "data" / "AS209_MATISSE").glob("*.fits"))
data = oim.oimData(files)

# NOTE: Apply filters
f1 = oim.oimWavelengthRangeFilter(targets=[0, 2], wlRange=[3.2e-6, 3.8e-6])
filt_bin_L = oim.oimWavelengthBinningFilter(
    targets=[0, 2], bin=5, normalizeError=False
)
filt_bin_N = oim.oimWavelengthBinningFilter(
    targets=1, bin=7, normalizeError=False
)
data.setFilter(oim.oimDataFilter([f1,filt_bin_L,filt_bin_N]))
wave_data = np.unique(data.vect_wl)

# NOTE: plot the unfiltered and filtered data (VISAMP)
fig = plt.figure()
ax = plt.subplot(projection="oimAxes")
data.useFilter = False
ax.oiplot(data, "SPAFREQ", "VISAMP", lw=3, alpha=0.2, label="unfiltered")
data.useFilter = True
ax.oiplot(data, "SPAFREQ", "VISAMP", label="filtered")
ax.set_ylabel("Visibility", fontsize=22)
ax.set_xlabel("Spatial frequency (rad$^{-1}$)", fontsize=22)
ax.tick_params(axis="both", labelsize="20")
ax.legend(fontsize=16)
ax.autolim()

fig = plt.figure()
ax = plt.subplot(projection="oimAxes")
data.useFilter = False
ax.oiplot(
    data,
    "EFF_WAVE",
    "FLUXDATA",
    errorbar=False,
    xunit="micron",
    lw=3,
    alpha=0.1,
    label="unfiltered",
)
data.useFilter = True
ax.oiplot(
    data,
    "EFF_WAVE",
    "FLUXDATA",
    errorbar=False,
    xunit="micron",
    label="filtered",
)
ax.set_ylabel("Flux density (Jy)", fontsize=22)
ax.set_xlabel(r"Wavelength ($\mu$m)", fontsize=22)
ax.tick_params(axis="both", labelsize="20")
ax.legend(fontsize=16)
ax.autolim()


# NOTE: Plot the unfiltered and filtered data (FLUXDATA)
# figFlux,axFlux =plt.subplots(figsize=(15,5),subplot_kw=dict(projection="oimAxes"))
# data.useFilter = False
# data.plot("EFF_WAVE","FLUXDATA",color="byFile",xunit="micron",errorbar=False,axe=axFlux, lw=3, alpha=0.1,label="unfiltered")
# data.useFilter = True
# data.plot("EFF_WAVE","FLUXDATA",color="byFile",xunit="micron",errorbar=False,axe=axFlux,label="filtered")
# axFlux.legend()
# axFlux.set_xlim(3,13)
# axFlux.set_ylabel('Flux density (Jy)')

# %%
# NOTE: Load a SED from a text file
path = Path(__file__).parent.parent.parent
sed_file = list((path / "data" / "AS209_MATISSE").glob("*.dat"))
sed_data = np.loadtxt(sed_file[0])
wl = sed_data[:, 0] * 1e-6  # in m
dwl = 1e-9  # dummy value since dwl is currently not used in oimodeler
flx = sed_data[:, 1]  # in W.m-2
err_flx = sed_data[:, 2]  # in W.m-2
flx_dens_Jy = (flx / wl) * wl**2 / const.c * 1e26  # conversion in Jy
err_flx_dens_Jy = (err_flx / wl) * wl**2 / const.c * 1e26  # conversion in Jy

# NOTE: create an oimFluxData object containing only an OI_FLUX table
oitarget = data.data[0]["OI_TARGET"].copy()
SEDFluxData = oim.oimFluxData(oitarget, wl, dwl, flx_dens_Jy, err_flx_dens_Jy)

# NOTE: (OPTIONAL) Uncomment the following lines if you want to include the SED data in your "data" object and plot it
# data.addData(SEDFluxData)
# fig = plt.figure()
# ax = plt.subplot(projection='oimAxes')
# ax.oiplot(data.data[0:3], "EFF_WAVE", "FLUXDATA", color="tab:blue", errorbar=True, xunit="micron", label="MATISSE spectra")
# ax.oiplot(data.data[3], "EFF_WAVE", "FLUXDATA", color="tab:red", errorbar=False, xunit="micron", label="SED",marker='o',linestyle='none')
# ax.set_ylabel('Flux density (Jy)',fontsize=25)
# ax.set_xlabel('Wavelength ($\mu$m)',fontsize=25)
# ax.set_xscale("log")
# ax.tick_params(axis='both',labelsize='25')
# ax.legend(fontsize=18)
# ax.autolim()

# NOTE: Display info about the data
data.info()

# %%

# NOTE: Radial grid of the model can be set to 'linear' or 'logarithmic'
oim.oimOptions.model.grid.type = "logarithmic"

# NOTE: The padding of the 1D-grid can be set (Multiplies to outer radius at
# grid creation). Default is 'None/1' and doesn't pad
oim.oimOptions.ft.padding = 1

# NOTE: Define a point source for the central start with a flux defined with a blackbody interpolator
s = oim.oimPt(
    f=oim.oimInterp("starWl", temp=4600, dist=121, lum=2.23, wl=None)
)

# NOTE: Sets the dust opacities
opac_file = (
    path / "data" / "FSCMa_MATISSE" / "dustkappa_olivine_graphite_1_20.inp"
)
op_wl, op = np.loadtxt(opac_file, usecols=[0, 1], unpack=True)
plt.figure()
plt.plot(op_wl,op)
plt.xlim(3,13)
plt.show()

# NOTE: Define a first instance of the temperature gradient model for the circumstellar emission
tg = oim.oimTempGrad(
    dim=128,
    dist=121.0,
    temp0=350,
    rin=0.1,
    rout=30,
    p=-1,
    q=-0.5,
    dust_mass=7e-9,
    pa=76,
    elong=1.5,
    kappa_abs=oim.oimInterp("wl", wl=op_wl * 1e-6, values=op),
)

# NOTE: Set the number of spatial frequency elements at which the Hankel transform will be computed.
# If not specified, the Hankel Transform will be computed for all the spatial frequencies of the data (usually slower).
tg.precision = 128

# NOTE: (OPTIONAL) Uncomment the following two lines if you want to set manually the model wavelengths to be fitted (to speed up fitting process).
# tg._wl = np.array([3.0e-6,3.5e-6,4.0e-6,8.0e-6,10.e-6,13e-6])
# s._wl = tg._wl

# NOTE: Model creation
model = oim.oimModel([s, tg])

# NOTE: Simulate the initial model observables and compute the associated reduced Chi2
sim = oim.oimSimulator(data=data, model=model)
sim.compute(computeChi2=True, computeSimulatedData=True)
print("Chi2r = {}".format(sim.chi2r))

# NOTE: plot the model observables and data (without SED)
fig0, ax0 = sim.plot(
    [
        "VISAMP",
    ]
)
fig3 = plt.figure()
ax3 = plt.subplot(projection="oimAxes")
ax3.oiplot(
    data.data,
    "EFF_WAVE",
    "FLUXDATA",
    errorbar=True,
    xunit="micron",
    kwargs_error={"alpha": 0.3},
    color="grey",
    label="MATISSE spectra",
)
ax3.oiplot(
    sim.simulatedData.data,
    "EFF_WAVE",
    "FLUXDATA",
    xunit="micron",
    color="red",
    label="Model spectra",
)
ax3.set_ylabel("Flux density (Jy)", fontsize=25)
ax3.set_xlabel(r"Wavelength ($\mu$m)", fontsize=25)
ax3.tick_params(axis="both", labelsize="25")
ax3.legend(fontsize=18)

# NOTE: (OPTIONAL) Uncomment the following lines if you want to plot the model observables and data (with SED)
# fig0, ax0 = sim.plot(["VISAMP",])
# fig3= plt.figure()
# ax3 = plt.subplot(projection='oimAxes')
# ax3.oiplot(data.data[0:3],"EFF_WAVE","FLUXDATA",
#             errorbar=True,xunit='micron',kwargs_error={"alpha":0.3},color='grey',label='MATISSE spectra')
# ax3.oiplot(sim.simulatedData.data[0:3],"EFF_WAVE","FLUXDATA",xunit='micron',color='red',label='Model spectra')
# ax3.oiplot(data.data[3],"EFF_WAVE","FLUXDATA",
#             errorbar=False,xunit='micron',kwargs_error={"alpha":0.3},color='grey',linestyle='none',marker='o',label='Measured SED')
# ax3.oiplot(sim.simulatedData.data[3],"EFF_WAVE","FLUXDATA",xunit='micron',color='red',linestyle='none',marker='o',label='Model SED')
# ax3.set_xscale('log')
# ax3.set_xlim(0.3,20)
# ax3.set_ylabel('Flux density (Jy)',fontsize=25)
# ax3.set_xlabel('Wavelength ($\mu$m)',fontsize=25)
# ax3.tick_params(axis='both',labelsize='25')
# ax3.legend(fontsize=18)

# NOTE: Plot of the model images at two wavelengths
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
model.showModel(
    512,
    0.15,
    wl=[3.5e-6],
    fromFT=False,
    normPow=0.3,
    axe=ax[0],
    colorbar=False,
)
model.showModel(
    512,
    0.15,
    wl=[10.5e-6],
    fromFT=False,
    normPow=0.3,
    axe=ax[1],
    colorbar=False,
)
ax[1].get_yaxis().set_visible(False)
ax[0].set_title(r"$\lambda = 3.5~\mu$m")
ax[1].set_title(r"$\lambda = 10.5~\mu$m")
plt.show()

#%%
# NOTE: Specifying the parameter space for the fit
tg.params["rin"].free = False
tg.params["rout"].free = False
tg.params["temp0"].free = False
tg.params["q"].set(min=-1, max=0)
tg.params["p"].set(min=-1.5, max=-0.5)
tg.params["dust_mass"].set(min=6e-9, max=6e-8)
tg.params["elong"].set(min=1, max=5)
tg.params["pa"].set(min=0, max=180)

# NOTE: Perfoming the model-fitting on VISAMP and fluxes
fit = oim.oimFitterEmcee(
    data, model, nwalkers=32, dataTypes=["VISAMP", "FLUXDATA"]
)
fit.prepare(init="random")
print(fit.initialParams)
fit.run(nsteps=5000, progress=True)

# NOTE: Plot the walkers path and make the corner plot
figWalkers, axeWalkers = fit.walkersPlot()
figCorner, axeCorner = fit.cornerPlot(discard=2000, chi2limfact=3)

# NOTE: Get the best-fit values (here the mode of the parameters posterior distribution) and print them
best, err_l, err_u, err = fit.getResults(
    mode="best", discard=2000, chi2limfact=3
)
fit.printResults(mode="best", format=".2e", discard=2000)

# NOTE: Plot the data and best-fit model (without SED)
fig0, ax0 = fit.simulator.plot(
    ["VISAMP"], kwargsSimulatedData=dict(ls="none", marker=".", lw=3)
)
fig1, ax1 = fit.simulator.plot(
    ["FLUXDATA"],
    xaxis="EFF_WAVE",
    xunit="micron",
    kwargsData=dict(color="blue"),
    kwargsSimulatedData=dict(ls="none", marker=".", lw=3, color="red"),
)


#%%
# NOTE: (OPTIONAL) Uncomment the following lines if you want to plot the data and best-fit model (with SED)
# fig3= plt.figure()
# ax3 = plt.subplot(projection='oimAxes')
# ax3.oiplot(data.data[0:3],"EFF_WAVE","FLUXDATA",
#             errorbar=True,xunit='micron',kwargs_error={"alpha":0.3},color='grey',label='MATISSE spectra')
# ax3.oiplot(fit.simulator.simulatedData.data[0:3],"EFF_WAVE","FLUXDATA",xunit='micron',color='red',label='Model spectra')
# ax3.oiplot(data.data[3],"EFF_WAVE","FLUXDATA",
#             errorbar=False,xunit='micron',kwargs_error={"alpha":0.3},color='grey',linestyle='none',marker='o',label='Measured SED')
# ax3.oiplot(fit.simulator.simulatedData.data[3],"EFF_WAVE","FLUXDATA",xunit='micron',color='red',linestyle='none',marker='o',label='Model SED')
# ax3.set_xscale('log')
# ax3.set_xlim(0.3,20)
# ax3.set_ylabel('Flux density (Jy)',fontsize=25)
# ax3.set_xlabel('Wavelength ($\mu$m)',fontsize=25)
# ax3.tick_params(axis='both',labelsize='25')
# ax3.legend(fontsize=18)

# NOTE: Setup model
tg_bestfit = oim.oimTempGrad(
    dim=128,
    dist=121.0,
    temp0=350,
    rin=0.1,
    rout=30,
    p=best[4],
    q=best[3],
    dust_mass=best[2],
    pa=best[1],
    elong=best[0],
    kappa_abs=oim.oimInterp("wl", wl=op_wl * 1e-6, values=op),
)

model_bestfit = oim.oimModel([s, tg_bestfit])

# NOTE: Best-fit model images at two wavelengths
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
model_bestfit.showModel(
    512,
    0.15,
    wl=[3.5e-6],
    fromFT=False,
    normPow=0.3,
    axe=ax[0],
    colorbar=False,
)
model_bestfit.showModel(
    512,
    0.15,
    wl=[10.5e-6],
    fromFT=False,
    normPow=0.3,
    axe=ax[1],
    colorbar=False,
)
ax[1].get_yaxis().set_visible(False)
ax[0].set_title(r"$\lambda = 3.5~\mu$m")
ax[1].set_title(r"$\lambda = 10.5~\mu$m")
plt.show()
