# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 14:07:54 2023

@author: Ame
"""

from datetime import datetime
from pathlib import Path
from pprint import pprint

import numpy as np

import oimodeler as oim

# NOTE: Set the FFT parameters and backend. As the kinematicsDisk model as no
# outer sharp edge we decide not to zero-pad the images during the FFT process
# no oder to save time and without noticeable effect on the simulated visibilities
oim.oimOptions.ft.padding = 1

# NOTE: In this example we will use an old VLTI/AMBER observation of a classical
# Be star Alpha Col published in Cochetti et al. 2019
path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "AMBER_AlphaCol"
files = list(data_dir.glob("*.fits"))

# NOTE: Create a keplerian rotating disk model
nwl = 81
wl0 = 2.1656e-6
c = oim.oimKinematicDisk(
    dim=64,
    fov=20,
    incl=45,
    Rstar=5.8,
    dist=80,
    fwhmCont=2.2,
    fluxDiskCont=0.25,
    EW=10.4,
    fwhmLine=5.7,
    nwl=nwl,
    vrot=360,
    beta=-0.5,
    pa=-83.2,
    wl0=2.1656e-6,
    dwl=0.9e-10,
    res=1.8e-10,
)
m = oim.oimModel(c)

# NOTE: Define the parameters space and free/fixed parameters
c.params["f"].free = False
c.params["Rstar"].free = False
c.params["dist"].free = False
c.params["v0"].free = False
c.params["vinf"].free = False
c.params["gamma"].free = False
c.params["fwhmCont"].free = False
c.params["beta"].free = False

c.params["pa"].set(min=-180, max=0, error=10)
c.params["fwhmLine"].set(min=2, max=10, error=2)
c.params["fluxDiskCont"].set(min=0, max=0.6, error=0.1)
c.params["incl"].set(min=20, max=70, error=10)
c.params["EW"].set(min=2, max=12, error=1)
c.params["vrot"].set(min=100, max=600, error=20)

# NOTE: Filter and modify the data before fitting.
data = oim.oimData(files)
for datai in data.data:
    datai["OI_VIS"].data["VISPHIERR"] *= 0.5
    datai["OI_VIS2"].data["VIS2ERR"] *= 0.5
    datai["OI_T3"].data["T3PHIERR"] *= 0.5

# NOTE: Only fit the Br Gamma line within a 20 angstrom band
dlam = 2e-9
f1 = oim.oimWavelengthRangeFilter(
    targets="all", wlRange=([(wl0 - dlam), (wl0 + dlam)])
)
filters = oim.oimDataFilter([f1])
data.setFilter(filters)

# NOTE: Create a simple simulator to compare the data vs. model
sim = oim.oimSimulator(data, m)
t0 = datetime.now()
sim.compute(computeChi2=True)
dt = (datetime.now() - t0).total_seconds() * 1000
pprint("compute chi2: {:.1f}ms/model".format(dt))
pprint("chi2:{} ".format(sim.chi2r))

# NOTE: Create the fitter and prepare the parameter space
fit = oim.oimFitterEmcee(data, m, nwalkers=12)

c.params["pa"].set(min=-180, max=0, error=10)
c.params["fwhmLine"].set(min=2, max=10, error=2)
c.params["fluxDiskCont"].set(min=0, max=0.6, error=0.1)
c.params["incl"].set(min=20, max=70, error=10)
c.params["EW"].set(min=2, max=12, error=1)
c.params["vrot"].set(min=100, max=600, error=20)

pprint(m.getFreeParameters())

fit = oim.oimFitterEmcee(data, m, nwalkers=12)
fit.prepare(init="random")

# NOTE: Run the mcmc fit
fit.run(nsteps=2000, progress=True)

# NOTE: Plot the walkers and parameter probability distributions
fit.walkersPlot(chi2limfact=3)
fit.cornerPlot(chi2limfact=5, discard=1000)

# NOTE: Generate plots of the data and model per baseline
sim.data.setFilter()
fig = sim.plotWlTemplate(
    [["VIS2DATA"], ["VISPHI", "T3PHI"]], xunit="micron", figsize=(12, 7)
)
fig.set_xlim((wl0 - 4e-9) * 1e6, (wl0 + 4e-9) * 1e6)
fig.set_legends("$LENGTH$m $PA$$^o$", "VIS2DATA", fontsize=8)
fig.set_legends(
    0.02, 0.9, "$BASELINE$", ["VIS2DATA", "VISPHI", "T3PHI"], fontweight=200
)
title = f"$\\alpha$ Col AMBER data + Kinematics disk model chi$^2_r$={sim.chi2r:.1f}"
fig.suptitle(title)
fig.tight_layout()


# NOTE: Plot images of the best model
wl0 = 2.1656e-6
Dwl = 32e-10
nwl = 5
wl = np.linspace(wl0 - Dwl / 2, wl0 + Dwl / 2, num=nwl)
mydim = 256
t0 = datetime.now()
dim = c.params["dim"].value
c.params["dim"].value = mydim
fig, ax, im = m.showModel(
    mydim,
    0.025 / mydim * 512,
    wl=wl,
    legend=True,
    normalize=True,
    fromFT=False,
    normPow=0.2,
    colorbar=False,
    cmap="plasma",
)
c.params["dim"].value = dim
