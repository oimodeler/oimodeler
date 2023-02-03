# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 09:08:46 2022

@author: Ame
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim
from astropy.io import fits

path = os.path.dirname(oim.__file__)

# Path to a fake MATISSE-L-band binary observation (3 oifits) created with ASPRO
pathData=os.path.join(path,os.pardir,"examples","testData","FSCMa_MATISSE")
files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

#%%

data=oim.oimData(files)


f1=oim.oimWavelengthRangeFilter(targets="all",wlRange=[3.0e-6,4e-6])

filters=oim.oimDataFilter([f1])
data.setFilter(filters)

#data.setFilters() #removing the filter
data.useFilter = False


#%%

fig=plt.figure()
ax = plt.subplot(projection='oimAxes')

data.useFilter = False
ax.oiplot(data,"SPAFREQ","VIS2DATA",color="tab:blue",lw=3,alpha=0.2,label="unfiltered")

data.useFilter = True
ax.oiplot(data,"SPAFREQ","VIS2DATA",color="tab:blue",label="filtered")

ax.set_yscale('log')
ax.legend()
ax.autolim()

fig.tight_layout()
fig.savefig(os.path.join(path,os.pardir,"images","ExampleFilter_wavelengthCut.png"))

#%%

f2=oim.oimRemoveArrayFilter(targets="all",arr=["OI_VIS","OI_FLUX"])         
f3=oim.oimDataTypeFilter(targets="all",dataType=["T3AMP","T3PHI"])
data.setFilter(oim.oimDataFilter([f1,f2,f3]))

