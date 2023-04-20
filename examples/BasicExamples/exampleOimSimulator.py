# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import oimodeler as oim


# Path to a fake MATISSE-L-band binary observation (3 oifits) created with ASPRO
path = Path(oim.__file__).parent.parent
pathData = path / Path() / "examples" / "testData" / "ASPRO_MATISSE2"

# TODO: After pathlib change of all `oimodeler` modules, remove str here
files = list(map(str, pathData.glob("*.fits")))

# Building a oimodeler model with the same parameters
ud = oim.oimUD(d=3, f=1, x=10, y=20)
pt = oim.oimPt(f=0.5)
model = oim.oimModel([ud, pt])

# Creating the simulator with the filename list and the model
sim = oim.oimSimulator(data=files, model=model)

# Preparing data (building vectors of coordinates and structure of data types)
# This is automatically called when creating the simulator.
# sim.data.prepareData()

# Computing the complex corr flux from the model at the data spatial freq
# with option to compute chi2 and final simulated data in oifits format
sim.compute(computeChi2=True, computeSimulatedData=True)
# Printing data model chi2r
pprint(f"Chi2r = {sim.chi2r}")


# %%
fig0, ax0 = sim.plot(["VIS2DATA", "VISAMP", "VISPHI", "T3AMP", "T3PHI"],
                     savefig=path / Path().parent / "images" / "ExampleOimSimulator_model0.png")
