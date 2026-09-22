from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np

import oimodeler as oim

path = Path(__file__).parent.parent.parent

# NOTE: Change these path if you want to save the products at another location
save_dir = path / "images"
product_dir = path / "data"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# %%
pt = oim.oimPt(f=0.1)
ud = oim.oimUD(d=10, f=0.5)
g = oim.oimGauss(fwhm=5, f=1)
r = oim.oimIRing(d=5, f=0.5)

# %%
pprint(ud)
#%%
pprint(ud.params["d"])
pprint(ud.d)

# NOTE: Build a few models from the components
mPt = oim.oimModel(pt)
mUD = oim.oimModel(ud)
mG = oim.oimModel(g)
mR = oim.oimModel(r)
mUDPt = oim.oimModel(ud, pt)

# %%
params = mUDPt.getParameters()
pprint(params)

freeParams = mUDPt.getFreeParameters()
pprint(freeParams)

# %%
im = mUDPt.getImage(512, 0.1)
plt.figure()
plt.imshow(im**0.2)
plt.savefig(save_dir / "basicModel_imshow.png")

# %%
im = mUDPt.getImage(256, 0.1, toFits=True)
pprint(im)
pprint(im.header)
pprint(im.data)
im = mUDPt.saveImage(product_dir / "modelImage.fits", 256, 0.1)


# %%
figImg, axImg, Img = mUDPt.showModel(
    512,
    0.1,
    normPow=0.2,
    figsize=(5, 4),
    savefig=save_dir / "basicModel_showModel.png",
)

# NOTE: Set some spatial frequencies (Baselines from 0 to 300m at 2.1 microns)
wl = 2.1e-6
B = np.linspace(0.0, 300, num=200)
spf = B / wl
spf0 = spf * 0
pprint(spf)

# %%
ccf = mUDPt.getComplexCoherentFlux(spf, spf * 0)  # East-West baselines

# %%
v = np.abs(ccf)
v = v / v.max()
plt.figure()
plt.plot(spf, v)
plt.xlabel("spatial frequency (cycles/rad)")
plt.ylabel("Visbility")
plt.savefig(save_dir / "basicModel_vis0.png")

#%%

fig, ax = mUD.plotVis(B,wl,xunit="cycle/mas")
fig.savefig(save_dir / "basicModel_vis1.png")

#%%
# NOTE: Some components
models = [mPt, mUD, mG, mR, mUDPt]
mNames = [
"Point Source (Pt)",
"Uniform Disk (UD",
"Gausian",
"Ring",
"UD + Pt",
]

nmodel = len(models)

fig, ax = plt.subplots(4, nmodel, figsize=(10,10/ nmodel*4))

for i, m in enumerate(models):
    m.showModel(512, 0.1, normPow=0.2, axe=ax[0, i], colorbar=False)
    spfmax = 0.69 # 11 cycle/mas
    m.showFourier(512, spfmax, axe=ax[1, i],colorbar=False,display="amp",unit="cycle/mas")
    m.showFourier(512, spfmax, axe=ax[2, i],colorbar=False,display="phase",unit="cycle/mas")
    m.plotVis(B,wl,xunit="cycle/mas",axe=ax[3,i])
    ax[0, i].set_title(mNames[i])
    ax[3, i].set_ylim(-0.05,1.05)

    if i!=0:
        for j in range(4):
            ax[j,i].get_yaxis().set_visible(False)
    for j in range(3):
        ax[j,i].get_xaxis().set_visible(False)            
   
ax[0,0].text(0.05,0.9,"IMAGE",      transform=ax[0,0].transAxes,color="w",ha="left",va="top")  
ax[1,0].text(0.05,0.9,"FT MODULUS", transform=ax[1,0].transAxes,color="w",ha="left",va="top")   
ax[2,0].text(0.05,0.9,"FT PHASE",   transform=ax[2,0].transAxes,color="w",ha="left",va="top")   
ax[3,0].text(0.05,0.9,"VISIBILITY", transform=ax[3,0].transAxes,color="k",ha="left",va="top")   


fig.tight_layout()
fig.savefig(save_dir / "basicModel_all.png")

