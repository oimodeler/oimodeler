# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 09:08:46 2022

@author: Ame
"""
from pathlib import Path

import matplotlib.pyplot as plt
import oimodeler as oim


# Path to a fake MATISSE-L-band binary observation (3 oifits) created with ASPRO
path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "FSCMa_MATISSE"

# NOTE: Change this path if you want to save the products at another location
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

files = list(data_dir.glob("*.fits"))

# %%
data = oim.oimData(files)

f1 = oim.oimWavelengthRangeFilter(targets="all", wlRange=[3.0e-6, 4e-6])

filters = oim.oimDataFilter([f1])
data.setFilter(filters)

# data.setFilters() #removing the filter
data.useFilter = False

# %%
fig = plt.figure()
ax = plt.subplot(projection='oimAxes')

data.useFilter = False
ax.oiplot(data, "SPAFREQ", "VIS2DATA", color="tab:blue",
          lw=3, alpha=0.2, label="unfiltered")

data.useFilter = True
ax.oiplot(data, "SPAFREQ", "VIS2DATA", color="tab:blue", label="filtered")

ax.set_yscale('log')
ax.legend()
ax.autolim()

fig.tight_layout()
fig.savefig(save_dir / "ExampleFilter_wavelengthCut.png")

# %%
f2 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX"])
f3 = oim.oimDataTypeFilter(targets="all", dataType=["T3AMP", "T3PHI"])
data.setFilter(oim.oimDataFilter([f1, f2, f3]))
