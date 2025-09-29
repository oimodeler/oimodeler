# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 14:02:23 2025

@author: ame
"""

import oimodeler as oim
from pathlib import Path
import matplotlib.pyplot as plt


dir0=Path(__file__).parents[2]
dirdata= dir0 / "data" / "realData"


dirplot = dir0 / "images"

fdata = list((dirdata / "PIONIER"/ "canopus" ).glob("*.fits"))


data = oim.oimData(fdata)

fig,ax = plt.subplots(1,2,subplot_kw=dict(projection='oimAxes'),figsize=(15,6))
data.uvplot(color="byBaseline",axe=ax[0])
data.plot("SPAFREQ","VIS2DATA",color="byBaseline",axe=ax[1],xunit="cycle/mas",errorbar=True)
ax[1].set_yscale("log")
fig.savefig(dirdata / "gridFitting_uv_and_v2_plot.png")

#%%
ud = oim.oimUD(d=45)
mud = oim.oimModel(ud)
mud.normalizeFlux()



grid = oim.oimFitterRegularGrid(data,mud,dataTypes=["VIS2DATA","T3PHI"])
grid.prepare(params=[ud.params["d"]],min=[6],max=[9],steps=[0.05])
grid.run()
fig_grid, ax_grid = grid.plotMap()
ax_grid.set_yscale("log")
fig_grid.savefig(dirplot / "gridFitting_grid1D.png")


grid.printResults()
#%%

miniz = oim.oimFitterMinimize(data,mud,dataTypes=["VIS2DATA","T3PHI"])
miniz.prepare()
miniz.run()
miniz.printResults()

#%%
fig_sim, ax_sim = grid.simulator.plot(["VIS2DATA","T3PHI"])
ax_sim[0].set_yscale("log")
ax_sim[0].set_ylim(1e-4,1)
fig_sim.savefig(dirplot / "gridFitting_simulator.png")



#%%

pld = oim.oimPowerLawLDD()
mpldd=oim.oimModel(pld)
mpldd.normalizeFlux()

#%%

grid2 = oim.oimFitterRegularGrid(data,mpldd,dataTypes=["VIS2DATA","T3PHI"])
grid2.prepare(params=[pld.params["d"],pld.params["a"]],min=[6,0],max=[9,2],steps=[0.05,0.05])
grid2.run()
#%%
import matplotlib.colors as colors

fig2_grid2, ax_grid2 = grid2.plotMap(plotContour=True,
                                     contour_kwargs=dict(levels=[2,4]),
                                     norm=colors.LogNorm(),cmap="plasma")
fig2_grid2.savefig(dirplot / "gridFitting_grid2D.png")

#%%
grid2.printResults()

#%%
miniz2 = oim.oimFitterMinimize(data,mpldd,dataTypes=["VIS2DATA","T3PHI"])
miniz2.prepare()
miniz2.run()
miniz2.printResults()

#%%

fig_sim2, ax_sim2 = miniz2.simulator.plot(["VIS2DATA","T3PHI"])
ax_sim2[0].set_yscale("log")
ax_sim2[0].set_ylim(1e-4,1)
fig_sim2.savefig(dirplot / "gridFitting_simulator2.png")
#%%
