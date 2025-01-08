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
data_dir = path / "data"  / "RealData" / "MATISSE"/ "FSCMa"
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)
    
files = list(data_dir.glob("*.fits"))
data = oim.oimData(files)
pprint(data.data)

data.prepareData()
pprint(data.vect_u)
pprint(data.vect_v)
pprint(data.vect_wl)
pprint(data.vect_u.shape)
#%%

figuv, axuv = data.uvplot(unit="cycle/mas",cunit="micron",lw=3,cmap="plasma",
                           showLegend=False,savefig=save_dir / "oimDataExample_uvplot.png")


figdata, axdata = data.plot("SPAFREQ",["VIS2DATA","T3PHI"],cname="EFF_WAVE",cunit="micron",
                             errorbar=True,xunit="cycle/mas")
axdata[0].set_yscale("log")
figdata.savefig(save_dir / "oimDataExample_plot.png")


#%%
filt_wl = oim.oimWavelengthRangeFilter(wlRange=[3.1e-6, 4e-6])
data.setFilter(filt_wl)


figcut,axcut = data.plot("SPAFREQ","VIS2DATA",yscale="log",xunit="cycle/mas",removeFilter=True,label="Original data",color="orange",lw=4)
data.plot("SPAFREQ","VIS2DATA",axe=axcut,xunit="cycle/mas",label="Filtered data",color="k",lw=2)
axcut.legend()
axcut.set_title("Data cut in with 3.1<$\lambda$<4 microns") 
figcut.savefig(save_dir / "oimDataExample_plot_wlcut.png")
#%%
path = Path(__file__).resolve().parents[2] / "data" / "RealData" / "GRAVITY" / "HD58647"
filenames = list(path.glob("*.fits"))
data = oim.oimData(filenames)

filt_bin=oim.oimWavelengthBinningFilter(bin=100,normalizeError=False)
data.setFilter(filt_bin)

figbin,axbin = data.plot("SPAFREQ","VIS2DATA",yscale="log",xunit="cycle/mas",removeFilter=True,label="Original data",color="orange",lw=4)
data.plot("SPAFREQ","VIS2DATA",axe=axbin,xunit="cycle/mas",label="Filtered data",color="k",lw=2,showFlagged=True)
axbin.legend()
axbin.set_title("Data binned by a factor of 100")
figbin.savefig(save_dir / "oimDataExample_plot_bin.png")

#%%
path = Path(__file__).parent.parent.parent
dir0 = path / "data"  / "RealData" / "MATISSE"/ "FSCMa"
filenames = list(dir0.glob("*.fits"))
data = oim.oimData(filenames)

filt_length=oim.oimFlagWithExpressionFilter(expr="LENGTH>50")
data.setFilter(filt_length)

figflag,axflag = data.plot("SPAFREQ","VIS2DATA",xunit="cycle/mas",removeFilter=True,color="orange",label="Original data",lw=5)
data.plot("SPAFREQ","VIS2DATA",axe=axflag,xunit="cycle/mas",label="Filtered data",color="k",lw=2)
axflag.legend()
axflag.set_title("Removing (=flaggging) data where  B>50m")
figflag.savefig(save_dir / "oimDataExample_plot_flag.png")

#%%%
data = oim.oimData(filenames)
from astropy.io import ascii
import numpy as np

iso_spectrum_path=path /  "data"  / "NonInterferometricData"
iso_spectrum_fname=iso_spectrum_path / "HD45677_ISO.txt"

isodata=ascii.read(iso_spectrum_fname)

wl  = isodata.columns['col1'].data*1e-6 # in m
dwl = wl-np.roll(wl,1)
dwl[0] = dwl[1]

flx = isodata.columns['col2'].data      # in Jy
err_flx = isodata.columns['col3'].data  # in Jy

oitarget=data.data[0]["OI_TARGET"].copy()

isoFluxData = oim.oimFluxData(oitarget,wl,dwl,flx,err_flx)

data.addData(isoFluxData)

#%%
data.info()
figFlux,axFlux =plt.subplots(figsize=(15,5),subplot_kw=dict(projection="oimAxes"))
data.plot("EFF_WAVE","FLUXDATA",color="byFile",xunit="micron",errorbar=True,axe=axFlux)
axFlux.legend()
axFlux.set_xscale("log")
axFlux.set_xlim(1,100)
figFlux.savefig(save_dir / "oimDataExample_plot_oimFluxData.png")
