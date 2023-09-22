from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import oimodeler as oim


path = Path(__file__).parent.parent.parent

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

mas2rad = u.mas.to(u.rad)

# Some components
gd = oim.oimGauss(fwhm=oim.oimInterp("wl", wl=[3e-6, 3.5e-6, 4e-6], values=[1, 3, 4]),
                  f=oim.oimInterp("wl", wl=[3e-6, 4e-6], values=[1, 4]))
ud = oim.oimUD(d=0.7, f=1)
ud2 = oim.oimUD(x=0, y=5, d=1.5, f=oim.oimInterp(
    "wl", wl=[3e-6, 4e-6], values=[0, 1]))
ud3 = oim.oimUD(x=-5, y=2, d=0.5, f=oim.oimInterp("wl",
                wl=[3e-6, 4e-6], values=[0.1, 1]))
eg = oim.oimEGauss(fwhm=1, elong=1.5, pa=oim.oimInterp("wl", wl=[3e-6, 4e-6],
                                                       values=[20, 60]), f=oim.oimInterp("wl", wl=[3e-6, 4e-6], values=[1, 0.1]))
er = oim.oimERing(din=8, dout=oim.oimInterp(
    "wl", wl=[3e-6, 4e-6], values=[15, 20]), elong=2, pa=0)

# linking some parameters, actually replacing them
er.params["elong"] = eg.params["elong"]
er.params["pa"] = eg.params["pa"]

# Building a few models from the components
m1 = oim.oimModel([gd])
m2 = oim.oimModel([ud, gd])
m3 = oim.oimModel([ud, ud2, ud3])
m4 = oim.oimModel([eg, er])
models = [m1, m2, m3, m4]

nB = 50  # number of baselines
nwl = 100  # number of walvengths

# Create some spatial frequencies
wl = np.linspace(3e-6, 4e-6, num=nwl)
B = np.linspace(1, 400, num=nB)
Bs = np.tile(B, (nwl, 1)).flatten()
wls = np.transpose(np.tile(wl, (nB, 1))).flatten()
spf = Bs/wls
spf0 = spf*0


# %%
fig, ax = plt.subplots(4, 4, figsize=(8, 8))
for i, m in enumerate(models):

    # Compute visibilities for the models and spatial frequencies
    cf1 = np.abs(m.getComplexCoherentFlux(
        spf, spf0, wls)).reshape(len(wl), len(B))
    cf2 = np.abs(m.getComplexCoherentFlux(
        spf0, spf, wls)).reshape(len(wl), len(B))

    for iwl, wli in enumerate(wl):
        ax[2, i].plot(B/wli*mas2rad, cf1[iwl, :]/cf1[iwl, 0],
                      color=plt.cm.plasma(iwl/nwl), alpha=0.5)
        ax[3, i].plot(B/wli*mas2rad, cf2[iwl, :]/cf2[iwl, 0],
                      color=plt.cm.plasma(iwl/nwl), alpha=0.5)

    sc = ax[3, i].scatter(spf0*mas2rad, cf2, c=wls*1e6, s=0.1, cmap="plasma")
    ax[3, i].set_xlabel("B/$\\lambda$ (cycles/mas)")
    ax[2, i].margins(0, 0)
    ax[2, i].set_ylim(0, 1)
    ax[3, i].margins(0, 0)
    ax[3, i].set_ylim(0, 1)
    # ax[2,i].set_title(model_names[i],fontsize=8)
    if i != 0:
        ax[2, i].get_yaxis().set_visible(False)
        ax[3, i].get_yaxis().set_visible(False)
ax[2, 0].set_ylabel("Vis. (East-West)")
ax[3, 0].set_ylabel("Vis. (North-South)")

# %%
dim = 256  # Number of pixels for the image
pix = 0.1  # size of the pixel in the image in mas
wls2 = [3e-6, 4e-6]  # wavelengths array for the images

for i, m in enumerate(models):
    im = m.getImage(dim, pix, wls2, fromFT=True)
    for iwl, wli in enumerate(wls2):
        imi = im[iwl, :, :]
        imi /= np.max(imi)
        imshow = ax[iwl, i].imshow(imi, extent=[dim/2*pix, -dim/2*pix, -dim/2*pix, dim/2*pix],
                                   norm=colors.TwoSlopeNorm(vmin=0., vcenter=1e-1, vmax=1), cmap="plasma", interpolation=None)
        ax[iwl, i].text(0, dim*pix/2, "\n{:.1f}$\\mu$m".format(wli*1e6),
                        color="w", va="top", ha="center")
        if i != 0:
            ax[iwl, i].get_yaxis().set_visible(False)
        ax[iwl, i].set_xlabel("$\\alpha$(mas)")
        if i == 0:
            ax[iwl, i].set_ylabel("$\\delta$(mas)")

fig.tight_layout()
plt.subplots_adjust(right=1.05)
fig.colorbar(imshow, ax=ax[:2, :].ravel().tolist(),
             label="Normalized Intensity")
norm = colors.Normalize(vmin=np.min(wl), vmax=np.max(wl))
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sc, ax=ax[2:, :].ravel().tolist(), label="$\\lambda$ ($\\mu$m)")

file_name = save_dir / "createModelChromatic.png"
fig.savefig(file_name)

# %%
nB = 200  # number of baselines
nwl = 200  # number of walvengths

# Create some spatial frequencies
wl = np.linspace(3e-6, 4e-6, num=nwl)
B = np.linspace(1, 400, num=nB)
Bs = np.tile(B, (nwl, 1)).flatten()
wls = np.transpose(np.tile(wl, (nB, 1))).flatten()
spf = Bs/wls
spf0 = spf*0


fig, ax = plt.subplots(1, 1)
cf2 = np.abs(models[-2].getComplexCoherentFlux(spf0,
             spf, wls)).reshape(len(wl), len(B))

for iwl, wli in enumerate(wl):
    ax.plot(B/wli*mas2rad, cf2[iwl, :]/cf2[iwl, 0],
            color=plt.cm.plasma(iwl/nwl), alpha=0.5, lw=1)

ax.axis('off')
file_name = save_dir / "imageForLogo.png"
fig.savefig(file_name, dpi=240, transparent=True)
