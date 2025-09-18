# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 13:35:11 2025

@author: ame
"""

from pathlib import Path
from pprint import pprint
import numpy as np
import oimodeler as oim
import matplotlib.pyplot as plt


path = Path(__file__).parent.parent.parent
data_dir = path / "data"  / "RealData" / "MATISSE"/ "FSCMa"
save_dir = path / "images"
if not save_dir.exists():
    save_dir.mkdir(parents=True)
    
#%% loading some data


path = Path(__file__).parent.parent.parent
dir0 = path / "data"  / "RealData" / "MATISSE"/ "FSCMa"
filenames = list(dir0.glob("*.fits"))
data = oim.oimData(filenames)


#%%

Bname       = oim.getBaselineName(data.data[0],length=True,angle=True)
CPname      =  oim.getBaselineName(data.data[0],hduname="OI_T3")
confname    = oim.getConfigName(data.data[0])
B,PA        = oim.getBaselineLengthAndPA(data.data[0])
u,v       = oim.get2DSpaFreq(data.data[0])
wl          = oim.getWlFromOifits(data.data[0])


print(f"This oifits files contains data taken with the {confname} array")
print(f"The wavelength range is {wl.min()*1e6:.2f}-{wl.max()*1e6:.2f}\mum")

for i in range(len(Bname)):
    print(f"{Bname[i]}")
   
#%%
"""


vis2 = oim.createOiVis2(OI_REVN=OI_REVN,
                        DATE-OBS=DATE-OBS,
                        ARRNAME=ARRNAME,
                        INSNAME=INSNAME,
                        TARGET_ID=TARGET_ID,
                        TIME=TIME,
                        MJD=MJD,
                        INT_TIME=INT_TIME,
                        VIS2DATA=VIS2DATA,
                        VIS2ERR=VIS2ERR,
                        UCOORD=UCOORD,
                        VCOORD=VCOORD,
                        STA_INDEX=STA_INDEX,
                        FLAG=FLAG)
"""

#%%

target = oim.createOiTargetFromSimbad("Gamma Cas")
target.data[0]

#%%

