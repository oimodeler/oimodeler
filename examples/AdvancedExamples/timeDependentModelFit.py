from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


path = Path(__file__).parent.parent.parent
data_dir = path / "examples" / "testData" / "ASPRO_SPICA_GROWING_UD"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

files = list(data_dir.glob("*.fits"))
data = oim.oimData(files)

# %%
mjd_ref = [60076, 60080]
ud = oim.oimUD(d=oim.oimInterp("time", mjd=mjd_ref,
               values=[0, 0], extrapolate=True))
model = oim.oimModel(ud)
ud.params['d'].params[0].set(min=0, max=10)
ud.params['d'].params[1].set(min=0, max=10)
ud.params['f'].free = False

# %%
fit = oim.oimFitterEmcee(files, model, nwalkers=20)
fit.prepare(init="random")

f1 = oim.oimRemoveArrayFilter(target="all", arr=["OI_VIS", "OI_T3"])
fit.data.setFilter(oim.oimDataFilter([f1]))

# %%
fit.run(nsteps=300, progress=True)

# %%
figWalkers, axeWalkers = fit.walkersPlot(cmap="plasma_r",
                                         savefig=save_dir / "complexModel_timeDependentWalkers.png")
figCorner, axeCorner = fit.cornerPlot(discard=100,
                                      savefig=save_dir / "complexModel_timeDependentCorner.png")

# %%
median, err_l, err_u, err = fit.getResults(mode='median', discard=100)

# %%
fig, ax = plt.subplots(1, 1, subplot_kw=dict(
    projection='oimAxes'), figsize=(8, 8))
timecol = ax.oiplot(fit.simulator.data, "SPAFREQ", "VIS2DATA",
                    xunit="cycles/mas", label="Data", cname="MJD", lw=4)
ax.oiplot(fit.simulator.simulatedData, "SPAFREQ", "VIS2DATA",
          xunit="cycles/mas", color="k", ls=":", label="Model")
fig.colorbar(timecol, ax=ax, label="MJD (days)")
ax.legend()
plt.savefig(save_dir / "complexModel_timeDependentFit.png")

# %%
fig, ax = plt.subplots(1, 1, figsize=(8, 4))

diam_ref = ud.params['d'](0, mjd_ref)

mjd_Obs = np.unique(fit.data.vect_mjd)
diam_Obs = ud.params['d'](0, mjd_Obs)

mjd = np.linspace(60076, 60080, num=100)
diam = ud.params['d'](0, mjd)
ax.plot(mjd, diam, color="k", lw=2, label="interpolated parameter")
ax.scatter(mjd_ref, diam_ref, marker="o", s=100, color="r", zorder=10,
           label="Ref values for model fitting")
ax.scatter(mjd_Obs, diam_Obs, marker="X", s=100, c=mjd_Obs, zorder=10,
           label="Observations")
ax.set_xlabel("Time (MJD)")
ax.set_ylabel("Diameter (mas)")
ax.legend()
plt.savefig(save_dir / "complexModel_timeDependentParameter.png")
