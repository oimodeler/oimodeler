# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:00:30 2023

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


# as the model as no sharp outer edge, no zero-padding is needed
oim.oimOptions['FTpaddingFactor'] = 1
oim.oimOptions['FTBackend'] = oim.FFTWBackend

# Path to the AMBER oifits files: the classical Be star Alpha Col
path = Path(oim.__file__).parent.parent
pathData = path / Path().parent / "examples" / "testData" / "RealData" / "AMBER" / "AlphaCol"

# TODO: After pathlib change of all `oimodeler` modules, remove str here
files = list(map(str, pathData.glob("*.fits")))

# %%
dim, nwl, wl0 = 256, 51, 2.1656e-6
c = oim.oimKinematicDisk(dim=dim, fov=20, incl=45, Rstar=5.8, dist=80, fwhmCont=2.2,
                         fluxDiskCont=0.25, EW=10.4, fwhmLine=5.7, nwl=nwl,
                         vrot=360, beta=-0.5, pa=-83.2, wl0=2.1656e-6, dwl=0.9e-10,
                         res=1.8e-10)
m = oim.oimModel(c)

# %% Showing 5 images in narrow bands through the disk
dwl, nwl = 32e-10, 5
wl = np.linspace(wl0-dwl/2, wl0+dwl/2, num=nwl)
m.showModel(dim, 0.05, wl=wl, legend=True, normalize=True,
            fromFT=False, normPow=0.5, colorbar=False, cmap="hot")

# %% Simulating data and comparing data and model of the disk around the
# classical Be star Alpha Col
sim = oim.oimSimulator(files, m)
pprint(f"chi2: {sim.chi2r}")

# %% Plotting the data/model comparison (V² phi diff. and CP) for the 6 baselines
fig, ax = plt.subplots(nrows=4, ncols=4, sharex=True, figsize=(12, 7))
for k in range(2):
    lam = sim.data.data[k]['OI_WAVELENGTH'].data['EFF_WAVE']
    B, PA = oim.getBaselineLengthAndPA(sim.data.data[k])
    for i in range(3):
        ax[0+2*k][i].errorbar(lam*1e6, sim.data.data[k]['OI_VIS2'].data['VIS2DATA'][i, :],
                              sim.data.data[k]['OI_VIS2'].data['VIS2ERR'][i, :], color="tab:red", alpha=0.5)
        ax[0+2*k][i].plot(lam*1e6, sim.simulatedData.data[k]['OI_VIS2'].data['VIS2DATA'][i, :],
                          color="tab:blue")
        ax[0+2*k][i].set_ylim(0.3, 1.1)
        ax[0+2*k][i].text(wl0*1e6, 0.05, "B={0:.0f}m PA={1:.0f}$^o$".format(
            B[i], PA[i]), horizontalalignment="center")
        ax[0+2*k][i].set_ylim(0, 1.2)

        ax[1+2*k][i].errorbar(lam*1e6, sim.data.data[k]['OI_VIS'].data['VISPHI'][i, :],
                              sim.data.data[k]['OI_VIS'].data['VISPHIERR'][i, :], color="tab:red", alpha=0.5)
        ax[1+2*k][i].plot(lam*1e6, sim.simulatedData.data[k]
                          ['OI_VIS'].data['VISPHI'][i, :], color="tab:blue")
        ax[1+2*k][i].set_ylim(-30, 30)

        if k == 1:
            ax[1+2*k][i].set_xlabel("$\\lambda$ ($\mu$m)")

    ax[0+2*k][0].set_ylabel("V$^2$")
    ax[1+2*k][0].set_ylabel("$\\varphi$ ($^o$)")

    ax[0+2*k][3].axis('off')

    ax[1+2*k][3].errorbar(lam*1e6, sim.data.data[k]['OI_T3'].data['T3PHI'][0, :],
                          sim.data.data[k]['OI_T3'].data['T3PHIERR'][0, :], color="tab:red", alpha=0.5)
    ax[1+2*k][3].plot(lam*1e6, sim.simulatedData.data[k]
                      ['OI_T3'].data['T3PHI'][0, :], color="tab:blue")
    ax[1+2*k][3].set_ylim(-30, 30)
    ax[1+2*k][3].yaxis.set_label_position("right")
    ax[1+2*k][3].set_ylabel("CP  ($^o$)")

title = f"$\\alpha$ Col AMBER data + Kinematic disk model: chi$^2_r$= {sim.chi2r:.1f}"
fig.suptitle(title)
ax[0, 0].set_xlim(1e6*(wl0-dwl), 1e6*(wl0+dwl))
fig.tight_layout()
plt.show()
