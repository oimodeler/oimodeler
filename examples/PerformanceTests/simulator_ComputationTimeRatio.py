# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim
from datetime import datetime
from tqdm import tqdm


path = os.path.dirname(oim.__file__)
pathData = os.path.join(path, os.pardir, "examples",
                        "testData", "ASPRO_MATISSE")
files0 = [os.path.abspath(os.path.join(pathData, fi))
          for fi in os.listdir(pathData) if ".fits" in fi]


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
print('Full Computation time {:.3f}s'.format(
    (end_time0 - start_time0).total_seconds()))

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
filename = os.path.join(path, os.pardir, "images",
                        "oimodel_test_simulator_speed_ratio.png")
plt.savefig(filename)
