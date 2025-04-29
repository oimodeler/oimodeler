# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 11:02:12 2025

@author: ame
"""

from pathlib import Path
from pprint import pprint

import oimodeler as oim

# Path tothe binary star  Beta Ari observed with MIRCX (1 files)
path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "RealData" "/MIRCX" / "Beta Ari"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


file = data_dir / "MIRCX_L2.2023Oct14._bet_Ari.MIRCX_IDL.nn.AVG10m.fits"

# Building a oimodeler model with the same parameters
ud1 = oim.oimUD(d=1, f=0.8)
ud2 = oim.oimUD(d=0.8, f=0.2, x=5, y=15)
model = oim.oimModel([ud1, ud2])

# Creating the simulator with the filename list and the model
sim = oim.oimSimulator(data=file, model=model)


#accessing the data inside our simulator
sim.data.info()


# Computing the complex corr flux from the model at the data spatial freq
# with option to compute chi2 and final simulated data in oifits format
sim.compute(computeChi2=True, computeSimulatedData=True)

# Accessing the simulated data
sim.simulatedData.info()

#plotting the data model comparison
fig0, ax0 = sim.plot(["VIS2DATA", "T3PHI"],
                     savefig=save_dir / "ExampleOimSimulator_plot.png")


# Printing data model chi2r
pprint(f"Chi2r = {sim.chi2r}")

#Using the dataTypes option of the oimSimulator compute method allows to 
#compute chi2 only on certain data types such as square-visibility and closure phase
sim.compute(computeChi2=True, dataTypes=["VIS2DATA","T3PHI"])
pprint(f"Chi2r = {sim.chi2r}")

#%% Wavelength template plots (per baseline plot as funciton of wl)

fig1 = sim.plotWlTemplate([["VIS2DATA"],["T3PHI"]],xunit="micron",figsize=(22,3))
#fig.set_legends("$LENGTH$m $PA$$^o$","VIS2DATA")
fig1.set_legends(0.5,0.8,"$BASELINE$",["VIS2DATA","T3PHI"],fontsize=10,ha="center")
fig1.tight_layout()
fig1.savefig(save_dir / "ExampleOimSimulator_WlTemplatePlot.png")

#%% residual plot

fig2, ax2 = sim.plot_residuals(["VIS2DATA", "T3PHI"],
                     savefig=save_dir / "ExampleOimSimulator_residuals_plot.png")

