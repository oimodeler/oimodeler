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
# %%
pprint(g.fwhm(wl=[3e-6, 3.5e-6, 4e-6, 4.5e-6]))

# %%

figGim, axGim, im = mg.showModel(
    256,
    0.1,
    wl=[3e-6, 3.5e-6, 4e-6, 4.5e-6],
    legend=True,
    colorbar=False,
    figsize=(3.2, 3.2),
    normalize=True,
)
figGim.tight_layout()
figGim.savefig(save_dir / "complexModel_chromaticGaussian.png")


# %%
nB = 500  # number of baselines
nwl = 100  # number of walvengths

# NOTE: Create some spatial frequencies
wl = np.linspace(3e-6, 4e-6, num=nwl)
B = np.linspace(1, 400, num=nB)

figGv, axGv = plt.subplots(1, 1, figsize=(8, 4))
mg.plotVis(B,wl,axe=axGv)
axGv.margins(0, 0)
figGv.savefig(save_dir / "complexModel_chromaticGaussianVis.png")


# %%
ud = oim.oimUD(d=0.5, f=oim.oimInterp("wl", wl=[3e-6, 4e-6], values=[2, 0.2]))
m2 = oim.oimModel([ud, g])
fig2im, ax2im, im2 = m2.showModel(
    256,
    0.1,
    wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6],
    normPow=0.2,
    figsize=(3.2, 3.2),
    legend=True,
    colorbar=False,
)
fig2im.tight_layout()
fig2im.savefig(save_dir / "complexModel_UDAndGauss.png")



fig2v, ax2v = plt.subplots(1, 1, figsize=(14, 8))
m2.plotVis(B,wl,axe=ax2v)
ax2v.margins(0, 0)
ax2v.set_ylim(0, 1)
plt.savefig(save_dir / "complexModel_UDAndGaussVis.png")

# %%
eg = oim.oimEGauss(
    fwhm=oim.oimInterp("wl", wl=[3e-6, 4e-6], values=[2, 8]), elong=2, pa=90
)
el = oim.oimEllipse(
    d=0.5,
    f=oim.oimInterp("wl", wl=[3e-6, 4e-6], values=[2, 0.2]),
    elong=2,
    pa=90,
)

m3 = oim.oimModel([el, eg])
fig3im, ax3im, im3 = m3.showModel(
    256,
    0.1,
    wl=[3e-6, 4e-6],
    legend=True,
    figsize=(3.5, 3.2),
    normPow=0.2,
    colorbar=False,
)
fig3im.tight_layout()
fig3im.savefig(save_dir / "complexModel_Elong.png")


# %%

m3.plotVis(B,wl,PA=[0,90])
plt.savefig(save_dir / "complexModel_ElongVis.png")

#%%
# NOTE: Check the number of free parameters
pprint(m3.getFreeParameters())

# NOTE: Link some parameters
eg.elong = el.elong
eg.pa = el.pa
pprint(m3.getFreeParameters())

# %%
er = oim.oimERing()
er.elong = eg.elong
er.pa = oim.oimParamLinker(eg.pa, "+", 90)
er.din = oim.oimParamLinker(el.d, "*", 2)
er.dout = oim.oimParamLinker(el.d, "*", 4)

m4 = oim.oimModel([el, eg, er])

fig4im, ax4im, im4 = m4.showModel(
    256,
    0.1,
    wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6],
    colorbar=False,
    figsize=(3.2, 3.2),
    normPow=0.2,
    legend=True,
)

fig4im.tight_layout()
fig4im.savefig(save_dir / "complexModel_link.png")

pprint(m4.getFreeParameters())

# %%
el.d.value = 4
el.pa.value = 45

fig5im, ax5im, im = m4.showModel(
    256,
    0.1,
    wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6],
    figsize=(3.2, 3.2),
    colorbar=False,
    normPow=0.2,
    legend=True,
)

fig5im.tight_layout()
fig5im.savefig(save_dir / "complexModel_linkRotScale.png")

#%%

innerRim = oim.oimInnerRim(dim=128,d=20,incl=60,h=3,y=0,pa=-90+67,f=1)
expRing = oim.oimExpRing(dim=64, fwhm=5,f=10,elong=1)
mdisk = oim.oimModel(innerRim,expRing)

expRing.d  = innerRim.d

expRing.pa = oim.oimParamLinker(innerRim.pa,operator="add",fact=90)


def func(p):
    return 1./np.cos(np.deg2rad(p))

expRing.elong=oim.oimParamLinkerFunction(innerRim.incl,func)

pprint(mdisk.getFreeParameters())


fig, ax , _ = mdisk.showModel(256,0.3,fromFT=True,normPow=1)

fig.savefig(save_dir / "complexModel_linkingFunction0.png")


#%%

fig, ax = plt.subplots(1,4,figsize=(15,5))

d= [20,15,20,17]
pa = [0,20,90,-50]
incl = [20,40,60,45]

for i in range(4):
    innerRim.d.value=d[i]
    innerRim.pa.value=pa[i]
    innerRim.incl.value=incl[i]
        
    mdisk.showModel(256,0.3,fromFT=True,normPow=1,axe=ax[i],colorbar=False)
    if i!=0:
        ax[i].get_yaxis().set_visible(False)
fig.tight_layout()

fig.savefig(save_dir / "complexModel_linkingFunction.png")


# %%
star1 = oim.oimUD(f=0.8, d=1)
star2 = oim.oimPt(f=0.15, x=5, y=5)
star3 = oim.oimPt(x=15, y=12)
mtriple = oim.oimModel(star1, star2, star3)
star2.x.free = True
star2.y.free = True
star3.x.free = True
star3.y.free = True

pprint(mtriple.getFreeParameters())
mtriple.normalizeFlux()

pprint(mtriple.getFreeParameters())

print(star3.f())
star1.f.value = 0.5
print(star3.f())

# %%
gd1 = oim.oimGauss(fwhm=oim.oimInterp("time", mjd=[0, 1, 3], values=[1, 4, 1]))
ud1 = oim.oimUD(
    d=oim.oimInterp("wl", wl=[1e-6, 3e-6], values=[0.5, 2]), x=-4, y=0, f=0.1
)

m6 = oim.oimModel(gd1, ud1)

wls = np.array([1, 2, 3]) * 1e-6
times = [0, 1, 2, 3, 4]

fig6im, ax6im, im6 = m6.showModel(
    256,
    0.04,
    wl=wls,
    t=times,
    legend=True,
    figsize=(2.5, 2),
    fromFT=True,
    normalize=True,
)

fig6im.savefig(save_dir / "complexModel_time.png")
