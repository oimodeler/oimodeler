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
path = Path(__file__).parent.parent.parent
data_dir = path / "examples" / "testData" / "ASPRO_MATISSE"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

files0 = list(data_dir.glob("*.fits"))

text = ["Complex Corr Flux only", "Complex Corr Flux + Chi2"]
computeChi2 = [False, True]
computeSimulatedData = [False, False]

ndata = 50
dt = np.ndarray([2, ndata])
start_time0 = datetime.now()

ud = oim.oimUD(d=20, f=4)
pt = oim.oimPt(f=6)
model = oim.oimModel([ud, pt])


for idata in tqdm(range(ndata)):
    files = files0*(idata+1)
    sim = oim.oimSimulator(data=files, model=model)

    if idata == 0:
        x0 = np.size(sim.data.vect_u)

    navg = 100
    for itype in range(2):
        start_time = datetime.now()
        for i in range(navg):
            sim.compute(computeChi2=computeChi2[itype],
                        computeSimulatedData=computeSimulatedData[itype])
        end_time = datetime.now()
        dt[itype, idata] = (end_time - start_time).total_seconds() * 1000/navg

end_time0 = datetime.now()
pprint(f"Full Computation time {(end_time0 - start_time0).total_seconds():.3f}s")

# %% Plot Time Ratio between CorrFlux and Chi2
plt.figure()
x = np.linspace(1, ndata, ndata)*len(files0)

r = (dt[1, :]-dt[0, :])/dt[0, :]
y = np.poly1d(np.polyfit(x, r, 1))(x)
plt.plot(x, r, marker="o", ls="")
plt.plot(x, y)
plt.xlim(0, np.max(x))
plt.ylim(0, 3)
plt.xlabel("Number of OIFits files")
plt.ylabel("Computation time ratio ")
txt = ""


for c in model.components:
    txt += c.__str__().split("x")[0]
    txt += "+ "
txt = txt[:-3]
plt.title("Computation time Ratio \n dt($\chi^2$)/dt(F$_{corr}$)")
plt.savefig(save_dir / "oimodel_test_simulator_speed_ratio.png")
