# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:46:20 2025

@author: ame
"""

from pathlib import Path

import oimodeler as oim

path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "FSCMa_MATISSE"

save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

# NOTE: Load data
files = list(data_dir.glob("*.fits"))
data = oim.oimData(files)

# NOTE: Print the content of the data as a list of hdulists
print(data.data)

# NOTE: Number of oifits files in the oimData object
print(len(data.data))

# NOTE: Summarise the oifits data in the oimData object
data.info()

# NOTE: Summarise the inforomation on the first oifits files of the oimData object
data.data[0].info()

# NOTE: Direct modification of some data.
# Here we multiply by two the errors on *= 2
# ata.data[2]["OI_VIS2"].data["VIS2ERR"]*= 2

# NOTE: Print the shape of optimized vectors of u,v, wl, and time coordinates
print(data.vect_u.shape)
print(data.vect_v.shape)
print(data.vect_wl.shape)
print(data.vect_mjd.shape)

# NOTE: Print the optimized structure of data used for simulation
print(data.struct_nB)
print(data.struct_nwl)
print(data.struct_arrType)
print(data.struct_dataType)

# NOTE: Plot some data
figuv, axuv = data.uvplot(color="byConfiguration")
figuv.savefig(save_dir / "oimDataExample_uvplot.png")
figdata, axdata = data.plot(
    "SPAFREQ",
    ["VIS2DATA", "T3PHI"],
    cname="EFF_WAVE",
    cunit="micron",
    errorbar=True,
    xunit="cycle/mas",
)
axdata[0].set_yscale("log")
figdata.savefig(save_dir / "oimDataExample_plot.png")

