# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 09:06:35 2026

@author: ame

Here we present advanced features of the new oimBinaryOrbit component.The more
basic features are detailed in the oimBinaryOrbit.py example in the same folder.

For this example, we are using some data from Meilland et al. (2011) paper:
    
    The binary Be star δ Scorpii at high spectral and spatial resolution.
       I. Disk geometry and kinematics before the 2011 periastron 
       
Our dataset consist in :
    
- delsco_rv.dat : RV around the 2001 periastron from Miroshnichenko (2001)

- delsco_position.dat : separation measurtement (in mas) from litterature
      source => 0 : 4th Catalog of Interferometric Measurements of Binary Star (Hartkopf 2001)
             => 1 : from Tycner et al. (2011) derived  from NPOI data 
             => 2 : Meilland et al. (2011) derived from AMBER data
- 2 AMBER oifits files from Meilland et al. 2011 
    - DelSco_K_2010-04-15T07h12_MR.fits
    - DelSco_K_2010-04-15T07h55_MR.fits
    
In this example we show how to setup priors on the binary separation and radial
velocity and use them to perform fit with and without interferometric data.

"""

from pathlib import Path
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

#%% Setting the model. We use Meilland et al. 2011 orbit

T0=Time(2000.6941,format="byear").mjd
T = 10.811*u.yr.to(u.day)

orb = oim.oimBinaryOrbit(e  = 0.94,            # Eccentricity
                         a  = 98.74,           # semi-major axis (mas)
                         T0 = 2000.6941*u.yr,  # Time Periastron passage (MJD by default or decimal year)
                         T  = 10.811*u.yr,     # Period (days by default or any compatible astropy unit if specified)
                         i  = 30.2,            # inclination angle (deg)
                         O  = 174,             # Longitude of ascending node (deg)
                         o  = 0.7,             # Argument of periastron
                         Ka = 23.9,            # Radial Velocity semi-amplitude of the first component (km/s by default)
                         V0 = -6.7             # Systemic velocity  (km/s by default)
                         )

morb = oim.oimModel(orb)

#%% We define the parameter space around value found by Meilland et al. (2011)
# note that you can search in a larger parameter space but convergence will be slower

orb.params["a"].set(min=80,max=100)
orb.params["T"].set(min=8,max=12)
orb.params["T0"].set(min=2000,max=2001)
orb.params["primary_f"].free=False
orb.params["secondary_f"].free=False
orb.params["O"].set(min=150,max=180)
orb.params["o"].set(min=-20,max=20)
orb.params["i"].set(min=0,max=90)    
orb.params["Ka"].set(min=20,max=25,free=True)
orb.params["V0"].set(min=-10,max=0,free=True)

#%% defining the external constraint prior which consist of the weighted sum
# of the chi2r on the RV and the separation

def sepRvPrior(whatever):
    x_mod, y_mod = orb.getSeparation(mjd, mas=True)
    dist2 = ((x_mod - x) / err_x) ** 2 + ((y_mod - y) / err_y) ** 2
    chi2r_sep = np.sum(dist2 / len(mjd))
    rv_model = orb.getPrimaryRadialVelocity(mjdRv)
    chi2r_rv = np.sum((rv - rv_model) ** 2) / len(mjdRv)
    chi2r = chi2r_sep + 10 * chi2r_rv # here we put a factor 10 on the RV to put more weight on it 
    return chi2r

#%% Setting a mcmc fitter without interferometric data but with the prior function

data_empty = oim.oimData() # No interferometric data

fit=oim.oimFitterEmcee(data_empty, morb, nwalkers=20)
fit.simulator.cprior = sepRvPrior # the prior is set on the simulator
fit.prepare()

#%% Fitting Rv and separation using the simulator prior
fit.run(nsteps=20000,progress=True)

#%% Plotting results of the fit
fit.walkersPlot(chi2limfact=2,ncolors=8,savefig=save_dir / "ExampleBinary_delsco_walker.png")
fit.cornerPlot(discard=15000,chi2limfact=2,savefig=save_dir / "ExampleBinary_delsco_corner.png")
fit.printResults(discard=15000,chi2limfact=2)

#%% Plotting the best-fit orbit with the separation measurements

nt = 10000
tmod = np.linspace(T0 - T/2, T0 + T/2, nt)
x_mod0, y_mod0 = orb.getSeparation(tmod, mas=True)

x_mod, y_mod = orb.getSeparation(mjd, mas=True)
fig, ax = plt.subplots()
ax.grid(which="major", lw=1, alpha=0.2)
ax.grid(which="minor", lw=0.5, alpha=0.2)
ax.plot(x_mod0, y_mod0, c="pink", label="Best fit-Orbit")

ax.axis("equal")

ax.scatter(x_mod, y_mod, color="grey", label="at epochs of Obs", marker=".",s=20)
for i in range(len(x)):
    ax.plot([x[i], x_mod[i]], [y[i], y_mod[i]], color="grey", lw=0.5)
    

data_text = ["Hartkopf+ (2001)","Tycner+ (2011)","Meilland+ (2011)"]
col=["r","g","b"]
for i in range(3):
    idx = np.where(source == i)
    ax.scatter(x[idx], y[idx], color=col[i], label=data_text[i], marker=".", zorder=10)

ax.scatter(0, 0, marker="+", color="k", s=50, lw=1, zorder=10)

p,h = ax.get_legend_handles_labels()

legend = ax.legend(p[0:2],h[0:2],title="Model", title_fontsize=13,loc=1)
ax.legend(p[2:],h[2:],title="Data", title_fontsize=13,loc=4)
ax.add_artist(legend)

ax.set_xlabel("$\\alpha$ (mas)")
ax.set_ylabel("$\\delta$ (mas)")
plt.savefig(save_dir / "ExampleBinary_delsco_fitted_orbit.png")
#%% Plotting the best-fit Radial Velocity with the RV measurements

fi2, ax2 = plt.subplots()
ax2.grid(which="major", lw=1, alpha=0.2)
ax2.grid(which="minor", lw=0.5, alpha=0.2)

nt = 10000
tmod = np.linspace(T0 - T/8, T0 + T/8, nt)
rv_mod = orb.getPrimaryRadialVelocity(tmod)
ax2.plot(Time(tmod, format="mjd").byear, rv_mod, color="b", label="Best-fit model")
ax2.scatter(Time(mjdRv, format="mjd").byear, 
            rv, color="r", label="Data (Miroshnichenko+ 2001)", marker=".")
ax2.set_ylim(-60,5)
ax2.legend()
ax2.set_xlim(2000, 2001.4)
ax2.set_xlabel("Time (yr)")
ax2.set_ylabel("Radial Velocity (km/s)")
plt.savefig(save_dir / "ExampleBinary_delsco_fitted_rv.png")
#%% Looking at Meilland et al. 2011 AMBER data measurements

data_dir = path /"data"/"RealData"/"AMBER"/"delSco"
fname = list(data_dir.glob("*.fits"))

data = oim.oimData(fname)

#binning the data are we are only interested here in the Low frequency signal
filt = oim.oimWavelengthBinningFilter(bin=25)
data.setFilter(filt)

figData = plt.figure(figsize=(15,6))#constrained_layout=True)
gs = figData.add_gridspec(2,5)
ax_uv = figData.add_subplot(gs[:, -2:],projection='oimAxes')
ax_v2 =  figData.add_subplot(gs[0, :-2],projection='oimAxes')
ax_cp =  figData.add_subplot(gs[1, :-2],projection='oimAxes')
data.uvplot(color="byBaseline",axe=ax_uv)
data.plot("SPAFREQ","VIS2DATA",color="byBaseline",xunit="cycle/mas",axe=ax_v2,legend=True,errorbar=True)
data.plot("SPAFREQ","T3PHI",color="byBaseline",xunit="cycle/mas",axe=ax_cp,legend=True,errorbar=True)
figData.tight_layout()
#%% First a "normal" two components model with the primary as a UD
ud = oim.oimUD()
pt = oim.oimPt()
model = oim.oimModel(ud,pt)
model.normalizeFlux()
ud.params["d"].set(min=0,max=5)
pt.params["x"].set(min=-150,max=150,free=True)
pt.params["y"].set(min=-150,max=150,free=True)
ud.params["f"].set(min=0.9,max=1)

#%% Fit on the model using a standard mcmc fitter 
fit2 = oim.oimFitterEmcee(data,model,nwalkers=20,dataTypes=["VIS2DATA","T3PHI"])
fit2.prepare()
fit2.run(nsteps=10000,progress=True)
#%% Print and plot the results
figWalker, axWalker = fit2.walkersPlot(chi2limfact=5)
figwCorner, axCorner = fit2.cornerPlot(chi2limfact=1.2,discard=8000)

fit2.printResults(chi2limfact=1.2,discard=8000)
figwFit, axFit = fit2.simulator.plotWithResiduals(["VIS2DATA","T3PHI"],xunit="cycle/arcsec")
axFit[2].set_ylim(-20,20)
x_fit, y_fit = pt.params["x"].value,pt.params["y"].value
x_fit_err, y_fit_err = pt.params["x"].error,pt.params["y"].error
#%% Grid explorer around the best-fit position
grid = oim.oimFitterRegularGrid(data,model,dataTypes=["VIS2DATA","T3PHI"])
grid.prepare(params=[pt.params["x"],pt.params["y"]],min=[-10,70],max=[50,130],steps=[0.2,0.2])
grid.run()
#%% Plot the grid of
figGrid,axGrid = grid.plotMap()
#%% Compare the position found fitting the AMBER data and the orbit from Meilland et al. 2011

fig,ax = plt.subplots()
ax.grid(which="major", lw=1, alpha=0.2)
ax.grid(which="minor", lw=0.5, alpha=0.2)
T0=Time(2000.6941,format="byear").mjd
T = 10.811*u.yr.to(u.day)

orb = oim.oimBinaryOrbit(e  = 0.94,            # Eccentricity
                         a  = 98.74,           # semi-major axis (mas)
                         T0 = 2000.6941*u.yr,  # Time Periastron passage (MJD by default or decimal year)
                         T  = 10.811*u.yr,     # Period (days by default or any compatible astropy unit if specified)
                         i  = 30.2,            # inclination angle (deg)
                         O  = 174,             # Longitude of ascending node (deg)
                         o  = 0.7,             # Argument of periastron
                         Ka = 23.9,            # Radial Velocity semi-amplitude of the first component (km/s by default)
                         V0 = -6.7             # Systemic velocity  (km/s by default)
                         )

morb = oim.oimModel(orb)

MJD = data.data[0]["OI_VIS2"].data["MJD"][0]
dateText=Time(MJD,format="mjd").iso.split()[0]

#plotting the orbit from Meiulland+ 2001
nt = 10000
tmod = np.linspace(T0 - T/2, T0 + T/2, nt)
x_mod0, y_mod0 = orb.getSeparation(tmod, mas=True)
ax.plot(x_mod0,y_mod0,color="lightgrey",label="Best-fit Orbit",lw=5)
ax.axis('equal')

#overplotting the Chi2r_min+1 contour from the grid exploration
im = grid.chi2rMap
X,Y = np.meshgrid(grid.grid[0],grid.grid[1])
ax.contour(X,Y,np.transpose(im),levels=[im.min()+1],colors="r")

# Plottinf the expected position at AMBER data MJD according to Meilland+ 2001 orbit
x0,y0 = orb.getSeparation(MJD,mas=True)
ax.scatter(x0,y0,color="grey",marker="o",
           label=f"Expected position\nfor {dateText}",zorder=10)

#plotting the best-fit position from the mcmc run on the AMBER data
ax.errorbar(x_fit,y_fit,xerr=x_fit_err,yerr=y_fit_err,color="b",marker=".",
            ls="",label="Best-fit position (oimFitterEmcee)",zorder=10)

p,h = ax.get_legend_handles_labels()
p.append(Line2D([0], [0], color='r', lw=1, label='Line'))
h.append("min$\\left(\\chi^2_r\\right)$+1 contour (oimFitterRegularGrid)")


legend = ax.legend(p[0:2],h[0:2],title="Meilland+ (2011)", 
                   title_fontsize=13,loc=1)
ax.legend(p[2:],h[2:],title=f"Fit of the {dateText} AMBER data", 
          title_fontsize=13,loc=4)
ax.add_artist(legend)

ax.set_xlabel("$\\alpha$ (mas)")
ax.set_ylabel("$\\delta$ (mas)")
ax.set_xlim(X.min(),X.max())
ax.set_ylim(Y.min(),Y.max())

#%%
orb.primary=oim.oimUD(d=2.3,f=0.94)
orb.primary.params["d"].set(min=0,max=5)
orb.secondary.params["f"]=oim.oimParamNorm(orb.primary.params["f"])
orb.primary.params["d"].set(min=0,max=5)
orb.primary.params["f"].set(min=0.9,max=1)
morb.getFreeParameters()
#%%
fit3=oim.oimFitterEmcee(data, morb, nwalkers=50)
fit3.simulator.cprior = sepRvPrior 

