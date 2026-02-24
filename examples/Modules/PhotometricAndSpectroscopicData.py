# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:03:43 2025

@author: ame
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

import oimodeler as oim

path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "RealData" / "MATISSE" / "FSCMa"
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# %%%
filenames = list(data_dir.glob("*.fits"))
data = oim.oimData(filenames)

iso_spectrum_path = path / "data" / "NonInterferometricData"
iso_spectrum_fname = iso_spectrum_path / "HD45677_ISO.txt"

isodata = ascii.read(iso_spectrum_fname)

wl = isodata.columns["col1"].data * 1e-6  # in m
dwl = wl - np.roll(wl, 1)
dwl[0] = dwl[1]

flx = isodata.columns["col2"].data  # in Jy
err_flx = isodata.columns["col3"].data  # in Jy

oitarget = data.data[0]["OI_TARGET"].copy()

isoFluxData = oim.oimFluxData(oitarget, wl, dwl, flx, err_flx)

data.addData(isoFluxData)

# %%
data.info()
figFlux, axFlux = plt.subplots(
    figsize=(15, 5), subplot_kw=dict(projection="oimAxes")
)
data.plot(
    "EFF_WAVE",
    "FLUXDATA",
    color="byFile",
    xunit="micron",
    errorbar=True,
    axe=axFlux,
)
axFlux.legend()
axFlux.set_xscale("log")
axFlux.set_xlim(1, 100)
figFlux.savefig(save_dir / "oimDataExample_plot_oimFluxData.png")

