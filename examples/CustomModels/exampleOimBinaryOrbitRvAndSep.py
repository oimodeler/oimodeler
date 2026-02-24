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


path = Path(__file__).parent.parent.parent

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)
    
#%% Loading the separation and radial velocity data from Meilland et al. 2011

fname_pos= path / "data" / "NonInterferometricData" /  "delsco_position.dat"

data = ascii.read(fname_pos)

mjd = np.array(data["mjd"])
x = np.array(data["x"])
y = np.array(data["y"])
err_x = np.array(data["err_x"])
err_y = np.array(data["err_y"])
source = np.array(data["source"])


fname_rv= path / "data" / "NonInterferometricData" / "delsco_rv.dat"
data_rv=ascii.read(fname_rv)
mjdRv = Time(data_rv["yr"],format="byear").mjd
rv  =   np.array(data_rv["rv"])
#%% 
T0=Time(2000.6941,format="byear").mjd
T = 10.811*u.yr.to(u.day)

orb = oim.oimBinaryOrbit(e  = 0.94,   # Eccentricity
                         a  = 98.74,  # semi-major axis (mas)
                         T0 = T0,     # Time Periastron passage (MJD by default or decimal year)
                         T  = T,      # Period (in days by default or any compatible astropy unit if specified)
                         i = 30.2,    # inclination angle (deg)
                         O = 174,     # Longitude of ascending node (deg)
                         o =0.7,      # Argument of periastron
                         Ka = 23.9,   # Radial Velocity semi-amplitude of the first component (km/s by default)
                         V0 = -6.7    # Systemic velocity  (km/s by default)
                         )

morb = oim.oimModel(orb)

#%% Define the parameter space

orb.params["a"].set(min=80,max=100)
orb.params["T"].set(min=8*365.25,max=12*365.25)
orb.params["T0"].set(min=51000,max=52000)
orb.params["primary_f"].free=False
orb.params["secondary_f"].free=False
orb.params["O"].set(min=150,max=180)
orb.params["o"].set(min=-20,max=20)
orb.params["i"].set(min=0,max=90)    
orb.params["Ka"].set(min=20,max=25,free=True)
orb.params["V0"].set(min=-10,max=0,free=True)

#%% defining the external constraint prior 

def sepRvPrior(whatever):
    x_mod, y_mod = orb.getSeparation(mjd, mas=True)
    dist2 = ((x_mod - x) / err_x) ** 2 + ((y_mod - y) / err_y) ** 2
    chi2r = np.sum(dist2 / len(mjd))
    chi2r_sep = np.sum(dist2 / len(mjd))
    rv_model = orb.getPrimaryRadialVelocity(mjdRv)
    chi2r_rv = np.sum((rv - rv_model) ** 2) / len(mjdRv)
    chi2r = chi2r_sep + 10 * chi2r_rv
    return chi2r

#%% Setting a fitter without interferometric data but with the prior function
fit=oim.oimFitterEmcee(oim.oimData(),morb,nwalkers=20)
fit.simulator.cprior = sepRvPrior
fit.prepare()
#%% Fitting Rv and separation using the simulator prior
fit.run(nsteps=20000,progress=True)
#%% Plotting results of the fit
fit.walkersPlot(chi2limfact=2,ncolors=8)
fit.cornerPlot(discard=15000,chi2limfact=2,fontsize=6)
fit.printResults(discard=15000,chi2limfact=2)

#%% Plotting the fitted orbit

nt = 10000
tmod = np.linspace(T0 - T/2, T0 + T/2, nt)
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
#%% Plotting the fitted Radial Velocity


fi2, ax2 = plt.subplots()
ax2.grid(which="major", lw=1, alpha=0.2)
ax2.grid(which="minor", lw=0.5, alpha=0.2)

nt = 10000
tmod = np.linspace(T0 - T/8, T0 + T/8, nt)
rv_mod = orb.getPrimaryRadialVelocity(tmod)
ax2.plot(Time(tmod, format="mjd").byear, rv_mod, color="b", label="model Rv")
ax2.scatter(Time(mjdRv, format="mjd").byear, 
            rv, color="r", label="data", marker=".")

ax2.legend()
ax2.set_xlim(2000, 2001.4)
ax2.set_xlabel("Time (yr)")
ax2.set_ylabel("Radial Velocity (km/s)")
