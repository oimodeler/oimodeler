# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""
from pathlib import Path
from pprint import pprint

import oimodeler as oim
import matplotlib.pyplot as plt


path = Path(__file__).parent.parent.parent
data_dir = path / "data" / "FSCMa_MATISSE"

files = list(data_dir.glob("*.fits"))
data = oim.oimData(files)
pprint(data.data)

data.prepareData()
pprint(data.vect_u)
pprint(data.vect_v)
pprint(data.vect_wl)
pprint(data.vect_u.shape)


pprint(data.struct_wl)

#%%

data.setFilter(oim.oimWavelengthRangeFilter(wlRange=[3.1e-6,4e-6]))
fig,ax = data.plot("SPAFREQ","VIS2DATA",cname="EFF_WAVE",cunit="micron",yscale="log",xunit="cycle/mas",removeFilter=True,alpha=0.1,showColorbar=False)
data.plot("SPAFREQ","VIS2DATA",axe=ax,cname="EFF_WAVE",cunit="micron",errorbar=True,yscale="log",xunit="cycle/mas")
data.uvplot(unit="cycle/mas",cunit="micron",label="cmap on wavelength",lw=3,cmap="plasma")