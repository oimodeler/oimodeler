# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 13:35:11 2025

@author: ame
"""

from pathlib import Path

import oimodeler as oim

path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "RealData" / "MATISSE" / "FSCMa"
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# NOTE: Load some data
path = Path(__file__).parent.parent.parent
dir0 = path / "data" / "RealData" / "MATISSE" / "FSCMa"
filenames = list(dir0.glob("*.fits"))
data = oim.oimData(filenames)


# %%
Bname = oim.getBaselineName(data.data[0], length=True, angle=True)
CPname = oim.getBaselineName(data.data[0], hduname="OI_T3")
confname = oim.getConfigName(data.data[0])
B, PA = oim.getBaselineLengthAndPA(data.data[0])
u, v = oim.get2DSpaFreq(data.data[0])
wl = oim.getWlFromOifits(data.data[0])

print(f"This oifits files contains data taken with the {confname} array")
print(rf"The wavelength range is {wl.min()*1e6:.2f}-{wl.max()*1e6:.2f}$\mu$m")

for i in range(len(Bname)):
    print(f"{Bname[i]}")

target = oim.createOiTargetFromSimbad("Gamma Cas")
target.data[0]
