# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 07:28:17 2022

@author: Ame
"""
from pathlib import Path

import matplotlib.pyplot as plt
import oimodeler as oim


path = Path(__file__).parent.parent.parent
data_dir = path / "examples" / "testData" / "ASPRO_MATISSE2"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# TODO: After pathlib change of all `oimodeler` modules, remove str casting.
files = list(map(str, data_dir.glob("*.fits")))
data = oim.oimData(files)

# %%
fig1 = plt.figure()
ax1 = plt.subplot(projection='oimAxes')
ax1.uvplot(data)
plt.savefig(save_dir / "ExampleOimPlot_uv.png")

# %%
fig2 = plt.figure()
ax2 = plt.subplot(projection='oimAxes')
lamcol = ax2.oiplot(data, "SPAFREQ", "VIS2DATA", xunit="cycle/mas", label="Data",
                    cname="EFF_WAVE",cunit="micron", errorbar=True)

#plt.colorbar(lamcol, ax=ax2, label="$\\lambda$ ($\mu$m)")
#ax2.legend()

plt.savefig(save_dir / "ExampleOimPlot_v2.png")

# %%
fig3 = plt.figure()
ax3 = plt.subplot(projection='oimAxes')
ax3.oiplot(data, "EFF_WAVE", "VIS2DATA",xunit="micron", color="byConfiguration",
           errorbar=True, kwargs_error={"alpha": 0.3})
ax3.legend()

plt.savefig(save_dir / "ExampleOimPlot_v2Wl.png")

# using the projection='oimAsxes' for all subpltos allow to use oimPlot custom plots
fig4, ax4 = plt.subplots(2, 2, subplot_kw=dict(
    projection='oimAxes'), figsize=(8, 8))

# First we plot the uv plan coverage
ax4[0, 0].uvplot(data)

# Then we plot the V² as a function of spa. freq. with a wavelength-colorscale
lamcol = ax4[0, 1].oiplot(data, "SPAFREQ", "VIS2DATA", xunit="cycle/mas", label="Data",
                          cname="EFF_WAVE", cunit="micron", ls=":", errorbar=True)

# Adding the corresponding colobar of the walvength colorscale
#fig4.colorbar(lamcol, ax=ax4[0, 1], label="$\\lambda$ ($\mu$m)")

# Legending also work with multicolor plots
ax4[0, 1].legend()

# Plotting the V² as a function of the wavlength with errorbars
# kwargs_error keyword allows to pass keywords to the errorbar plot
ax4[1, 0].oiplot(data, "EFF_WAVE", "VIS2DATA", xunit="nm",color="byBaseline",
                 errorbar=True, kwargs_error={"alpha": 0.1})
ax4[1, 0].legend(fontsize=6)


# Finally we plot the Closure Phase with a fewstyling options
ax4[1, 1].oiplot(data, "SPAFREQ", "T3PHI", xunit="cycle/rad", errorbar=True,
                 lw=2, ls=":", color="byFile")
ax4[1, 1].legend(fontsize=4)


# setting log scale as for normal matplotlib plots
ax4[0, 1].set_yscale('log')

# autolim allows to directly set ylim to (0,1) for vis and (-180,180) for phase
ax4[1, 0].autolim()
ax4[1, 1].autolim()

fig4.tight_layout()

plt.savefig(save_dir / "ExampleOimPlot_multi.png")
