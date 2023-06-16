# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""
from datetime import datetime
from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim
from tqdm import tqdm


# Path to a fake MATISSE-L-band binary observation (3 oifits) created with ASPRO
path = Path().resolve().parent.parent
data_dir = path / "examples" / "testData" / "ASPRO_MATISSE"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# TODO: After pathlib change of all `oimodeler` modules, remove str casting.
files0 = list(map(str, data_dir.glob("*.fits")))

text = ["Complex Corr Flux only", "Complex Corr Flux + Chi2",
        "Complex Corr Flux + Sim. Data", "Full Computation "]
computeChi2 = [False, True, False, True]
computeSimulatedData = [False, False, True, True]

ntest, ndata = len(text), 30
dt = np.ndarray([ntest, ndata])

start_time0 = datetime.now()
ud = oim.oimUD(d=20, f=4)
pt = oim.oimPt(f=6)
model = oim.oimModel([ud, pt])

for idata in tqdm(range(ndata)):
    files = files0*(idata+1)
    sim = oim.oimSimulator(data=files, model=model)

    if idata == 0:
        x0 = np.size(sim.data.vect_u)

    navg = 10
    for itype in range(ntest):
        start_time = datetime.now()
        for i in range(navg):
            sim.compute(computeChi2=computeChi2[itype],
                        computeSimulatedData=computeSimulatedData[itype],
                        checkSimulatedData=False)
        end_time = datetime.now()
        dt[itype, idata] = (end_time - start_time).total_seconds() * 1000/navg

end_time0 = datetime.now()
pprint(f"Full Computation time {(end_time0 - start_time0).total_seconds():.3f}s")

# %%
x = np.linspace(1, ndata, ndata)*x0

col = plt.rcParams['axes.prop_cycle'].by_key()['color']
for itype in range(ntest):
    y = np.poly1d(np.polyfit(x, dt[itype, :], 1))(x)
    plt.plot(x, dt[itype, :], label=text[itype],
             marker="o", ls="", color=col[itype])
    plt.plot(x, y, color=col[itype])
plt.legend()
plt.xlabel("Number of data points")
plt.ylabel("Computation time (ms)")
plt.xlim(0, np.max(x))
plt.ylim(0, np.max(dt))
txt = ""

for c in model.components:
    txt += c.__str__().split("x")[0]
    txt += "+ "
txt = txt[:-3]

plt.title(f"Computation time for a {txt} model")
plt.savefig(save_dir / "oimodel_test_simulator_speed.png")
