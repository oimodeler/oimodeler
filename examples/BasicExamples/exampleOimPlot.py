# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 07:28:17 2022

@author: Ame
"""
from pathlib import Path

import matplotlib.pyplot as plt
import oimodeler as oim

path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "ASPRO_MATISSE2"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

files = list(data_dir.glob("*.fits"))
data = oim.oimData(files)

# %%

fig1 = plt.figure()
ax1 = plt.subplot(projection="oimAxes")
ax1.uvplot(data)
plt.savefig(save_dir / "ExampleOimPlot_uv.png")

# %%

fig2, ax2 = plt.subplots(
    3, 2, subplot_kw=dict(projection="oimAxes"), figsize=(8, 10)
)

ax2[0, 0].uvplot(
    data.data, color="byBaseline", marker=".", legendkwargs=dict(fontsize=6)
)
ax2[0, 1].uvplot(
    data.data, color="byFile", facecolor="w", legendkwargs=dict(fontsize=3.8)
)
ax2[1, 0].uvplot(data.data, color="byConfiguration", colorTab=["r", "g", "b"])
ax2[1, 1].uvplot(
    data.data, color="byArrname", colorTab=["r", "g", "b"], marker="+"
)
ax2[2, 0].uvplot(data.data, label="custom label", unit="km")
ax2[2, 1].uvplot(
    data.data, unit="cycle/rad", cunit="micron", lw=3, color="byConfiguration"
)

fig2.tight_layout()
plt.savefig(save_dir / "ExampleOimPlot_uv2.png")

# %%
fig3 = plt.figure()
ax3 = plt.subplot(projection="oimAxes")
ax3.uvplot(data.data, unit="cycle/mas", cunit="micron",
           label="cmap on wavelength", lw=2, cmap="plasma",
)
plt.savefig(save_dir / "ExampleOimPlot_uv3.png")

# %%
fig4 = plt.figure()
ax4 = plt.subplot(projection="oimAxes")
lamcol = ax4.oiplot(
    data, "SPAFREQ", "VIS2DATA",
    xunit="cycle/mas", label="Data",
    cname="EFF_WAVE", cunit="micron", errorbar=True)
ax4.legend()

plt.savefig(save_dir / "ExampleOimPlot_v2.png")

# %%
fig5 = plt.figure()
ax5 = plt.subplot(projection="oimAxes")
ax5.oiplot(
    data, "EFF_WAVE", "VIS2DATA",
    xunit="micron", color="byConfiguration",
    errorbar=True, kwargs_error={"alpha": 0.3})
ax5.legend()
plt.savefig(save_dir / "ExampleOimPlot_v2Wl.png")
# %%

# using the projection='oimAsxes' for all subpltos allow to use oimPlot custom plots
fig6, ax6 = plt.subplots(
    2, 2, subplot_kw=dict(projection="oimAxes"), figsize=(8, 8))

# Then we plot the V² as a function of spa. freq. with a wavelength-colorscale
ax6[0, 0].oiplot(
    data, "SPAFREQ", "VIS2DATA",
    xunit="cycle/mas", label="Data", cname="EFF_WAVE",
    cunit="micron", ls=":", errorbar=True)

# Legending also work with multicolor plots
ax6[0, 0].legend()

# setting log scale as for normal matplotlib plots
ax6[0, 0].set_yscale("log")


# Plotting the V² as a function of the wavlength with errorbars
# kwargs_error keyword allows to pass keywords to the errorbar plot
ax6[0, 1].oiplot(
    data, "EFF_WAVE", "VIS2DATA",
    xunit="nm", color="byBaseline",
    errorbar=True, kwargs_error={"alpha": 0.1})
ax6[0, 1].legend(fontsize=6)


# Finally we plot the Closure Phase with a fewstyling options
ax6[1, 0].oiplot(
    data, "SPAFREQ", "T3PHI",
    xunit="cycle/rad", errorbar=True,
    lw=2, ls=":", color="byFile")
ax6[1, 0].legend(fontsize=4)


ax6[1, 1].oiplot(
    data, "EFF_WAVE", "T3PHI",
    xunit="m", cname="LENGTH",
    errorbar=True, kwargs_error={"alpha": 0.1})

fig6.tight_layout()
plt.savefig(save_dir / "ExampleOimPlot_multi.png")

# %%

data_path = path / "examples" / "data" / "AMBER_AlphaCol"
files = [
    data_path / "ALPHACOL_2010-01-09T00_58.fits",
    data_path / "ALPHACOL_2010-01-20T10_36.fits",
]
data = oim.oimData(files)


# %%
fig7, ax7 = plt.subplots(
    3, 1, subplot_kw=dict(projection="oimAxes"), figsize=(12, 10)
)
ax7[0].oiplot(
    data, "EFF_WAVE", "VIS2DATA", xunit="Angstrom", color="byBaseline"
)
ax7[0].legend()
ax7[1].oiplot(data, "EFF_WAVE", "VISPHI", xunit="Angstrom", color="byBaseline")
ax7[1].legend()
ax7[2].oiplot(data, "EFF_WAVE", "T3PHI", xunit="Angstrom", color="byBaseline")
ax7[2].legend()
fig7.tight_layout()
plt.savefig(save_dir / "ExampleOimPlot_AlphaCol0.png")

# %%

fig = plt.figure(FigureClass=oim.oimWlTemplatePlots, figsize=(12, 7))
fig.autoShape(data.data, shape=[["VIS2DATA", None], ["VISPHI", "T3PHI"]])
fig.set_xunit("micron")
fig.plot(
    data.data, plotFunction=plt.Axes.errorbar,
    plotFunctionkwarg=dict(color="gray", alpha=0.3))
fig.plot(data.data, plotFunctionkwarg=dict(color="tab:blue", lw=0.5))
fig.set_ylim(["VISPHI", "T3PHI"], -25, 25)
fig.set_ylim(["VIS2DATA"], 0, 1.2)
fig.set_xlim(2.16, 2.172)
fig.set_legends(
    0.5, 0.1, "$BASELINE$ $LENGTH$m $PA$$^o$",
    ["VIS2DATA", "VISPHI"], fontsize=12, ha="center")
fig.set_legends(0.5, 0.1, "$BASELINE$", ["T3PHI"], fontsize=12, ha="center")
fig.tight_layout()
plt.savefig(save_dir / "ExampleOimPlot_AlphaCol1.png")
