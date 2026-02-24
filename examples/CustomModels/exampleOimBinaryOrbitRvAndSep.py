# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 09:06:35 2026

@author: ame
"""

from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.time import Time

import oimodeler as oim

plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True

dir0 = Path(__file__).parent.parent.parent / "data" / "NonInterferometricData"

fname_pos = dir0 / "delsco_position.dat"
data = ascii.read(fname_pos)
mjd = np.array(data["mjd"])
x = np.array(data["x"])
y = np.array(data["y"])
err_x = np.array(data["err_x"])
err_y = np.array(data["err_y"])
source = np.array(data["source"])

fname_rv = dir0 / "delsco_rv.dat"
data_rv = ascii.read(fname_rv)
mjdRv = Time(
    np.array((data_rv["col1"] - 544.5) / 365.0 + 2000.0), format="byear"
).mjd
rv = np.array(data_rv["col2"])

# %%
t0 = Time(2000.6941, format="byear").mjd
e = 0.94
a = 98.7403
T = 10.811 * u.yr.to(u.day)
i = 30.2
o = 0.7
O = 174.0
v0 = -6.7
Ka = 23.9

orb = oim.oimBinaryOrbit(e=e, a=a, T=T, T0=t0, i=i, o=o, O=O, V0=v0, Ka=Ka)
morb = oim.oimModel(orb)

# NOTE: Define the parameter space
orb.params["a"].set(min=0, max=200)
orb.params["T"].set(min=0, max=365 * 20)
orb.params["T0"].set(min=51000, max=52000)
orb.params["primary_f"].free = False
orb.params["secondary_f"].free = False
orb.params["O"].set(min=0, max=180)
orb.params["o"].set(min=-180, max=180)
orb.params["i"].set(min=0, max=90)
orb.params["Ka"].set(min=0, max=100, free=True)
orb.params["V0"].set(min=-50, max=50, free=True)


# NOTE: Define an external constraint of the prior
def sepRvPrior(whatever):
    x_mod, y_mod = orb.getSeparation(mjd, mas=True)
    dist2 = ((x_mod - x) / err_x) ** 2 + ((y_mod - y) / err_y) ** 2
    chi2r = np.sum(dist2 / len(mjd))
    chi2r_sep = np.sum(dist2 / len(mjd))
    rv_model = orb.getPrimaryRadialVelocity(mjdRv)
    chi2r_rv = np.sum((rv - rv_model) ** 2) / len(mjdRv)
    chi2r = chi2r_sep + 10 * chi2r_rv
    return chi2r


# %%
sim = oim.oimSimulator(oim.oimData(), morb, cprior=sepRvPrior)
fit = oim.oimFitterEmcee(sim, nwalkers=50)
fit.prepare()

# %%
fit.run(nsteps=40000, progress=True)

# %%
fit.walkersPlot(chi2limfact=4, ncolors=8)

# %%
fit.cornerPlot(discard=35000, chi2limfact=2, fontsize=6)

# %%
fit.printResults(discard=15000, chi2limfact=2)

nt = 10000
tmod = np.linspace(t0 - T / 2, t0 + T / 2, nt)
x_mod0, y_mod0 = orb.getSeparation(tmod, mas=True)

x_mod, y_mod = orb.getSeparation(mjd, mas=True)
fig, ax = plt.subplots()
ax.grid(which="major", lw=1, alpha=0.2)
ax.grid(which="minor", lw=0.5, alpha=0.2)
ax.plot(x_mod0, y_mod0, c="b", label="Model Orbit")

ax.axis("equal")

ax.scatter(x_mod, y_mod, color="grey", label="model", marker=".")
for i in range(len(x)):
    ax.plot([x[i], x_mod[i]], [y[i], y_mod[i]], color="grey", lw=0.5)
ax.scatter(x, y, color="r", label="data", marker=".", zorder=10)
ax.scatter(0, 0, marker="+", color="k", s=50, lw=1, zorder=10)
ax.legend()

ax.set_xlabel("$\\alpha$ (mas)")
ax.set_ylabel("$\\delta$ (mas)")
# %%

fi2, ax2 = plt.subplots()
ax2.grid(which="major", lw=1, alpha=0.2)
ax2.grid(which="minor", lw=0.5, alpha=0.2)

nt = 10000
tmod = np.linspace(t0 - T / 8, t0 + T / 8, nt)
rv_mod = orb.getPrimaryRadialVelocity(tmod)
ax2.plot(Time(tmod, format="mjd").byear, rv_mod, color="b", label="model Rv")
ax2.scatter(
    Time(mjdRv, format="mjd").byear, rv, color="r", label="data", marker="."
)
ax2.legend()
ax2.set_xlim(2000, 2001.4)
ax2.set_xlabel("Time (yr)")
ax2.set_ylabel("Radial Velocity (km/s)")
