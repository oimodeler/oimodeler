from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


path = Path(__file__).parent.parent.parent

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)


# %%
g = oim.oimGauss(fwhm=oim.oimInterp("wl", wl=[3e-6, 4e-6], values=[2, 8]))
mg = oim.oimModel(g)
#%%
pprint(g.params['fwhm'](wl=[3e-6, 3.5e-6, 4e-6, 4.5e-6]))

# %%

figGim, axGim, im = mg.showModel(256, 0.1, wl=[3e-6, 3.5e-6, 4e-6, 4.5e-6],
                                 legend=True,colorbar=False,
                                 figsize=(3.2, 3.2),normalize=True)
figGim.tight_layout()
figGim.savefig(save_dir / "complexModel_chromaticGaussian.png")


#%%
nB = 500  # number of baselines
nwl = 100  # number of walvengths

# Create some spatial frequencies
wl = np.linspace(3e-6, 4e-6, num=nwl)
B = np.linspace(1, 400, num=nB)

spf = (B[np.newaxis,:]/wl[:,np.newaxis]).flatten()
wls = (np.ones((1,nB))*wl[:,np.newaxis]).flatten()



# %%
vis = np.abs(mg.getComplexCoherentFlux(
    spf, spf*0, wls)).reshape(len(wl), len(B))
vis /= np.outer(np.max(vis, axis=1), np.ones(nB))

figGv, axGv = plt.subplots(1, 1, figsize=(8, 4))
sc = axGv.scatter(spf, vis, c=wls*1e6, s=0.2, cmap="plasma")
figGv.colorbar(sc, ax=axGv, label="$\\lambda$ ($\\mu$m)")
axGv.set_xlabel("B/$\\lambda$ (cycles/rad)")
axGv.set_ylabel("Visiblity")
axGv.margins(0, 0)

plt.savefig(save_dir / "complexModel_chromaticGaussianVis.png")

# %%
ud = oim.oimUD(d=0.5, f=oim.oimInterp("wl", wl=[3e-6, 4e-6], values=[2, 0.2]))
m2 = oim.oimModel([ud, g])
fig2im, ax2im, im2 = m2.showModel(256, 0.1, wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6],
                                  normPow=0.2, figsize=(3.2, 3.2), 
                                  legend=True,colorbar=False)
fig2im.tight_layout()
fig2im.savefig(save_dir / "complexModel_UDAndGauss.png")



vis = np.abs(m2.getComplexCoherentFlux(
    spf, spf*0, wls)).reshape(len(wl), len(B))
vis /= np.outer(np.max(vis, axis=1), np.ones(nB))

fig2v, ax2v = plt.subplots(1, 1, figsize=(14, 8))
sc = ax2v.scatter(spf, vis, c=wls*1e6, s=0.2, cmap="plasma")
fig2v.colorbar(sc, ax=ax2v, label="$\\lambda$ ($\\mu$m)")
ax2v.set_xlabel("B/$\\lambda$ (cycles/rad)")
ax2v.set_ylabel("Visiblity")
ax2v.margins(0, 0)
ax2v.set_ylim(0, 1)
plt.savefig(save_dir / "complexModel_UDAndGaussVis.png")

# %%
eg = oim.oimEGauss(fwhm=oim.oimInterp(
    "wl", wl=[3e-6, 4e-6], values=[2, 8]), elong=2, pa=90)
el = oim.oimEllipse(d=0.5, f=oim.oimInterp(
    "wl", wl=[3e-6, 4e-6], values=[2, 0.2]), elong=2, pa=90)

m3 = oim.oimModel([el, eg])
fig3im, ax3im, im3 = m3.showModel(256, 0.1, wl=[3e-6,  4e-6],legend=True,
                                  figsize=(3.5, 3.2), normPow=0.2,colorbar=False)
fig3im.tight_layout()
fig3im.savefig(save_dir / "complexModel_Elong.png")


# %%
fig3v, ax3v = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)

# East-West
vis = np.abs(m3.getComplexCoherentFlux(
    spf, spf*0, wls)).reshape(len(wl), len(B))
vis /= np.outer(np.max(vis, axis=1), np.ones(nB))
ax3v[0].scatter(spf, vis, c=wls*1e6, s=0.2, cmap="plasma")
ax3v[0].set_title("East-West Baselines")
ax3v[0].margins(0, 0)
ax3v[0].set_ylim(0, 1)
ax3v[0].set_xlabel("B/$\\lambda$ (cycles/rad)")
ax3v[0].set_ylabel("Visiblity")

# North-South
vis = np.abs(m3.getComplexCoherentFlux(
    spf*0, spf, wls)).reshape(len(wl), len(B))
vis /= np.outer(np.max(vis, axis=1), np.ones(nB))
sc = ax3v[1].scatter(spf, vis, c=wls*1e6, s=0.2, cmap="plasma")
ax3v[1].set_title("North-South Baselines")
ax3v[1].set_xlabel("B/$\\lambda$ (cycles/rad)")
fig3v.colorbar(sc, ax=ax3v.ravel().tolist(), label="$\\lambda$ ($\\mu$m)")

plt.savefig(save_dir / "complexModel_ElongVis.png")

# %% Check the number of free parameters
pprint(m3.getFreeParameters())

# %% linkning some parameters
eg.params['elong'] = el.params['elong']
eg.params['pa'] = el.params['pa']
pprint(m3.getFreeParameters())


# %%
er = oim.oimERing()
er.params['elong'] = eg.params['elong']
er.params['pa'] = oim.oimParamLinker(eg.params["pa"], "+", 90)
er.params['din'] = oim.oimParamLinker(el.params["d"], "*", 2)
er.params['dout'] = oim.oimParamLinker(el.params["d"], "*", 4)

m4 = oim.oimModel([el, eg, er])

fig4im, ax4im, im4 = m4.showModel(256, 0.1, wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6],
                                  colorbar=False,figsize=(3.2, 3.2), normPow=0.2,
                                  legend=True)

fig4im.tight_layout()
fig4im.savefig(save_dir /  "complexModel_link.png")

pprint(m4.getFreeParameters())

#%%
el.params['d'].value = 4
el.params['pa'].value = 45

fig5im, ax5im, im = m4.showModel(256, 0.1, wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6],
                                 figsize=(3.2, 3.2), colorbar=False, normPow=0.2,
                                 legend=True)
                                 
fig5im.tight_layout()
fig5im.savefig(save_dir / "complexModel_linkRotScale.png")

#%%
star1 = oim.oimUD(f=0.8,d=1)
star2 = oim.oimPt(f=0.15,x=5,y=5)
star3 = oim.oimPt(x=15,y=12)
mtriple = oim.oimModel(star1,star2,star3)
star2.params["x"].free=True
star2.params["y"].free=True
star3.params["x"].free=True
star3.params["y"].free=True

pprint(mtriple.getFreeParameters())

#%%
#star3.params["f"] = oim.oimParamNorm([star1.params["f"],star2.params["f"]])
mtriple.normalizeFlux()

pprint(mtriple.getFreeParameters())

print(star3.params["f"]())
star1.params["f"].value = 0.5
print(star3.params["f"]())
#%%
gd1 = oim.oimGauss(fwhm=oim.oimInterp("time", mjd=[0, 1, 3], values=[1, 4, 1]))
ud1 = oim.oimUD(d=oim.oimInterp(
    "wl", wl=[1e-6, 3e-6], values=[0.5, 2]), x=-4, y=0, f=0.1)

m6 = oim.oimModel(gd1, ud1)

wls = np.array([1, 2, 3])*1e-6
times = [0, 1, 2, 3, 4]

fig6im, ax6im, im6 = m6.showModel(256, 0.04, wl=wls, t=times, legend=True,
                                  figsize=(2.5, 2), fromFT=True, normalize=True)

fig6im.savefig(save_dir / "complexModel_time.png")
