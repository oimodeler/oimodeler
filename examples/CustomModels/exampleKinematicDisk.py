# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 14:07:54 2023

@author: Ame
"""
from datetime import datetime
from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim

# %%
oim.oimOptions['FTpaddingFactor'] = 1
oim.oimOptions['FTBackend'] = oim.FFTWBackend

path = Path(oim.__file__).parent.parent
data_path = path / "tempTests"
oifits = [data_path / "ALPHACOL_2010-01-09T00_58.fits",
          data_path / "ALPHACOL_2010-01-20T10_36.fits"]


def errorFill(axe, X, Y, dY, color="k",
              alpha=1, zorder=0, line=True, smooth=3):
    Ys = Y
    if smooth != 1:
        ker = np.ones(smooth)/smooth
        Ys = np.convolve(Y, ker, mode="same")  # [smooth//2:-smooth//2-1]
    XX = np.concatenate([X, np.flip(X)])
    YY = np.concatenate([Ys-dY, np.flip(Ys+dY)])
    axe.fill(XX, YY, zorder=zorder, alpha=alpha, color=color)
    if line:
        axe.plot(X, Y, zorder=zorder+1, color=color)


# %%
nwl, wl0 = 81, 2.1656e-6
c = oim.oimKinematicDisk(dim=64, fov=20, incl=45, Rstar=5.8, dist=80, fwhmCont=2.2,
                         fluxDiskCont=0.25, EW=10.4, fwhmLine=5.7, nwl=nwl,
                         vrot=360, beta=-0.5, pa=-83.2, wl0=2.1656e-6, dwl=0.9e-10,
                         res=1.8e-10)
m = oim.oimModel(c)

# %%
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
# c.params["beta"].set(mini=-1, maxi=0,error=0.1)
# pprint(m.getFreeParameters())

# %%
dlam = 2e-9
f1 = oim.oimWavelengthRangeFilter(
    targets="all", wlRange=([(wl0-dlam), (wl0+dlam)]))
filters = oim.oimDataFilter([f1])

# TODO: After pathlib change of all `oimodeler` modules, remove str here
data = oim.oimData(str(oifits))

for datai in data.data:
    datai["OI_VIS"].data["VISPHIERR"] *= 0.5
    datai["OI_VIS2"].data["VIS2ERR"] *= 0.5
    datai["OI_T3"].data["T3PHIERR"] *= 0.5

data.setFilter(filters)

# %%
sim = oim.oimSimulator(data, m)
t0 = datetime.now()
n = 10
for i in range(n):
    sim.compute(computeChi2=True)
dt = (datetime.now() - t0).total_seconds()*1000/n
pprint(f"compute chi2: {dt:.1f}ms/model")
pprint(f"chi2: {sim.chi2r}")

# %%
fit = oim.oimFitterEmcee(data, m, nwalkers=12)

fit.prepare(init="gaussian")
fit.run(nsteps=2000, progress=True)

# %%
fit.walkersPlot(chi2limfact=3)
fit.cornerPlot(chi2limfact=3)

# %%
data.setFilter(filters)
sim.compute(computeChi2=True, computeSimulatedData=True)
pprint(f"chi2:{sim.chi2r} ")

data.setFilter()
# sim.prepareData()
plotDlam = True
dlam, alpha, alpha2, lw = 8e-9, 0.3, 0.7, 2
wlrange = [(wl0-dlam)*1e6, (wl0-dlam)*1e6]

fig, ax = plt.subplots(nrows=4, ncols=4, sharex=True, figsize=(12, 7))
for k in range(2):
    lam = sim.data.data[k]['OI_WAVELENGTH'].data['EFF_WAVE']
    lamsim = sim.simulatedData.data[k]['OI_WAVELENGTH'].data['EFF_WAVE']
    B, PA = oim.getBaselineLengthAndPA(sim.data.data[k])
    for i in range(3):
        errorFill(ax[0+2*k][i], lam*1e6, sim.data.data[k]['OI_VIS2'].data['VIS2DATA'][i, :],
                  sim.data.data[k]['OI_VIS2'].data['VIS2ERR'][i, :], color="tab:red", alpha=alpha)
        ax[0+2*k][i].plot(lamsim*1e6, sim.simulatedData.data[k]['OI_VIS2'].data['VIS2DATA'][i, :],
                          color="tab:blue", linewidth=lw, alpha=alpha2)
        ax[0+2*k][i].set_ylim(0.3, 1.1)
        ax[0+2*k][i].text(wl0*1e6, 0.05, "B={0:.0f}m PA={1:.0f}$^o$".format(
            B[i], PA[i]), horizontalalignment="center")
        if plotDlam:
            ax[0+2*k][i].plot([(wl0-dlam)*1e6, (wl0-dlam)*1e6],
                              [0, 1.2], ls=":", color="k")
            ax[0+2*k][i].plot([(wl0+dlam)*1e6, (wl0+dlam)*1e6],
                              [0, 1.2], ls=":", color="k")

        ax[0+2*k][i].set_ylim(0, 1.2)

        errorFill(ax[1+2*k][i], lam*1e6, sim.data.data[k]['OI_VIS'].data['VISPHI'][i, :],
                  sim.data.data[k]['OI_VIS'].data['VISPHIERR'][i, :], color="tab:red", alpha=alpha)
        ax[1+2*k][i].plot(lamsim*1e6, sim.simulatedData.data[k]['OI_VIS'].data['VISPHI']
                          [i, :], color="tab:blue", linewidth=lw, alpha=alpha2)
        ax[1+2*k][i].set_ylim(-30, 30)
        if plotDlam:
            ax[1+2*k][i].plot([(wl0-dlam)*1e6, (wl0-dlam)*1e6],
                              [-100, 100], ls=":", color="k")
            ax[1+2*k][i].plot([(wl0+dlam)*1e6, (wl0+dlam)*1e6],
                              [-100, 100], ls=":", color="k")


    # ax[0+2*k][0].set_xlim([(wl0-Dlam)*1e6,(lam0+Dlam)*1e6])
    # sp=np.mean(sim.data.data[k]["AMBER_SPECTRUM"].data["SPECTRUM"],axis=1)
    # sperr=np.mean(sim.data.data[k]["AMBER_SPECTRUM"].data["SPECTRUM_ERROR"],axis=1)
    # sperr/=np.mean(sp)/10
    # sp/=np.mean(sp)

    # ax[0+2*k][3].plot(lam*1e6,sp,color="tab:red",linewidth= 3,alpha=alpha)
#    errorFill(ax[0+2*k][3],lam*1e6,sp,sperr,color="tab:red",alpha=alpha)
#    ax[0+2*k][3].plot(lamodel*1e6,spmodel-0.07,color="tab:blue",linewidth=2)
    ax[0+2*k][3].axis('off')
    #######################################################################

    errorFill(ax[1+2*k][3], lam*1e6, sim.data.data[k]['OI_T3'].data['T3PHI'][0, :],
              sim.data.data[k]['OI_T3'].data['T3PHIERR'][0, :], color="tab:red", alpha=alpha)
    ax[1+2*k][3].plot(lamsim*1e6, sim.simulatedData.data[k]['OI_T3'].data['T3PHI']
                      [0, :], color="tab:blue", linewidth=lw, alpha=alpha2)
    ax[1+2*k][3].set_ylim(-30, 30)

    if plotDlam:
        ax[1+2*k][3].plot([(wl0-dlam)*1e6, (wl0-dlam)*1e6],
                          [-100, 100], ls=":", color="k")
        ax[1+2*k][3].plot([(wl0+dlam)*1e6, (wl0+dlam)*1e6],
                          [-100, 100], ls=":", color="k")

title = f"$\\alpha$ Col AMBER data + Kinematics disk model chi$^2_r$={sim.chi2r:.1f}"
fig.suptitle(title)
ax[0][0].set_xlim((wl0-dlam)*1e6, (wl0+dlam)*1e6)
# %%
wl0, dwl, nwl = 2.1656e-6, 32e-10, 5
wl = np.linspace(wl0-dwl/2, wl0+dwl/2, num=nwl)
mydim = 256
t0 = datetime.now()
dim = c.params["dim"].value
c.params["dim"].value = mydim
m.showModel(mydim, 0.025/mydim*512, wl=wl, legend=True, normalize=True,
            fromFT=False, normPow=0.2, colorbar=False, cmap="hot")
c.params["dim"].value = dim
