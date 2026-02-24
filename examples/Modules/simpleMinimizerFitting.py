# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 14:10:18 2025

@author: ame
"""

from pathlib import Path
from pprint import pprint

import oimodeler as oim

# NOTE: Path tothe binary star  Beta Ari observed with MIRCX (1 files)
path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "RealData" "/PIONIER" / "canopus"

files = list(data_dir.glob("PION*.fits"))
data = oim.oimData(files)


# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# NOTE: Build an oimodeler model with the same parameters
pldd = oim.oimPowerLawLDD(d=8, a=0)
model = oim.oimModel(pldd)
model.normalizeFlux()

pprint(model.getFreeParameters())

# %%
lmfit = oim.oimFitterMinimize(data, model, dataTypes=["VIS2DATA", "T3PHI"])
lmfit.prepare()
lmfit.run()
lmfit.printResults()

# %%
fig, ax = lmfit.simulator.plotWithResiduals(
    ["VIS2DATA", "T3PHI"],
    xunit="cycle/mas",
    kwargsData=dict(color="byConfiguration"),
)

ax[0].set_yscale("log")
ax[0].set_ylim(1e-4, 1)
fig.savefig(save_dir / "simpleMinimizerFitting_plotwithresiduals.png")

# %%
pldd.params["d"].value = 15

lmfit.prepare()
lmfit.run()
lmfit.printResults()

# %%
fig2, ax2 = lmfit.simulator.plotWithResiduals(
    ["VIS2DATA", "T3PHI"],
    xunit="cycle/mas",
    kwargsData=dict(color="byConfiguration"),
)

ax2[0].set_yscale("log")
ax2[0].set_ylim(1e-4, 1)
fig2.savefig(save_dir / "simpleMinimizerFitting_plotwithresiduals_bad.png")


# %%
lmfit = oim.oimFitterMinimize(
    data, model, dataTypes=["VIS2DATA", "T3PHI"], mezthod="BFGS"
)

