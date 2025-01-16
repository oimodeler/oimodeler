# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:00:58 2025

@author: ame
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
    
#%% filtering in wavelength

path = Path(__file__).parent.parent.parent
dir0 = path / "data"  / "RealData" / "MATISSE"/ "FSCMa"
filenames = list(dir0.glob("*.fits"))
data = oim.oimData(filenames)

filt_wl = oim.oimWavelengthRangeFilter(wlRange=[3.1e-6, 4e-6])
data.setFilter(filt_wl)


figcut,axcut = data.plot("SPAFREQ","VIS2DATA",yscale="log",xunit="cycle/mas",removeFilter=True,label="Original data",color="orange",lw=4)
data.plot("SPAFREQ","VIS2DATA",axe=axcut,xunit="cycle/mas",label="Filtered data",color="k",lw=2)
axcut.legend()
axcut.set_title("Data cut in with 3.1<$\lambda$<4 microns") 
figcut.savefig(save_dir / "oimDataExample_plot_wlcut.png")
#%% binning in wavelength
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

#%% filtering with an expression: here baseline > 50m
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


#%% filterning by baseline name
path = Path(__file__).parent.parent.parent
dir0 = path / "data"  / "RealData" / "MATISSE"/ "FSCMa"
filenames = list(dir0.glob("*.fits"))
data = oim.oimData(filenames)

baselines=["A0-B2","A0-D0"]
filt_baselines=oim.oimKeepBaselinesFilter(baselines=baselines,arr="OI_VIS2")
data.setFilter(filt_baselines)
figflag,axflag = data.plot("SPAFREQ","VIS2DATA",xunit="cycle/mas",removeFilter=True,color="orange",label="Original data",lw=5)
data.plot("SPAFREQ","VIS2DATA",axe=axflag,xunit="cycle/mas",label="Filtered data",color="k",lw=2)
axflag.legend()
axflag.set_title(f"Keep only baselines {baselines}")
figflag.savefig(save_dir / "oimDataExample_plot_keepBaselines.png")

