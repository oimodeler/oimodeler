# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 09:08:46 2022

@author: Ame
"""


import oimodeler as oim
import matplotlib.pyplot as plt
import os
import numpy as np
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

#%%
f1=oim.oimRemoveArrayFilter(targets="all",arr=["OI_VIS","OI_FLUX"])         
f2=oim.oimWavelengthRangeFilter(targets="all",wlRange=[[3.0e-6,4e-6],[4.55e-6,4.95e-6]])
f3=oim.oimDataTypeFilter(targets="all",dataType=["T3AMP","T3PHI"])
"""
data.setFilter(oim.oimDataFilter([f1,f2]))
print(len(data.vect_u))


data.setFilter(oim.oimDataFilter([f1,f2,f3]))
print(len(data.vect_u))
"""


f4=oim.oimWavelengthRangeFilter(targets=[0],wlRange=[3.0e-6,4e-6])
f5=oim.oimWavelengthRangeFilter(targets=[1],wlRange=[3.5e-6,3.8e-6])
f6=oim.oimWavelengthRangeFilter(targets=[2],wlRange=[3.2e-6,3.9e-6])
data.setFilter(oim.oimDataFilter([f4,f5,f6]))
"""
plt.figure()
ax = plt.subplot(projection='oimAxes')
ax.oiplot(data.data,"SPAFREQ","VIS2DATA",color="byFile",alpha=0.1,label="Non filtered ")



ax.oiplot(data.data,"SPAFREQ","VIS2DATA",color="byFile",label="filtered ")
ax.legend()

ax.set_yscale('log')
ax.autolim()


#%%

plt.figure()
ax = plt.subplot(projection='oimAxes')
cb=ax.oiplot(data.data,"SPAFREQ","VIS2DATA",cname="VIS2ERR",label="Color by VIS2ERR")
ax.legend()
plt.colorbar(cb,ax=ax,label="VIS2ERR")
ax.autolim()


plt.figure()
ax = plt.subplot(projection='oimAxes')
ax.oiplot(data.data,"SPAFREQ","VIS2DATA",color="byConfiguration",label="My star ")
ax.legend()
ax.autolim()

plt.figure()
ax = plt.subplot(projection='oimAxes')
ax.oiplot(data.data,"SPAFREQ","VIS2DATA",color="byBaseline",)
ax.legend()
ax.autolim()

plt.figure()
ax = plt.subplot(projection='oimAxes')
ax.oiplot(data.data,"SPAFREQ","VIS2DATA",color="k",label="My label ")
ax.legend()
ax.autolim()

plt.figure()
ax = plt.subplot(projection='oimAxes')
cb=ax.oiplot(data.data,"SPAFREQ","VIS2DATA",cunitmultiplier=1e6,cname="EFF_WAVE",label="Color by Walvength")
ax.legend()
plt.colorbar(cb,ax=ax,label="Wavelength ($\\mu$m)")
ax.autolim()


plt.figure()
ax = plt.subplot(projection='oimAxes')
cb=ax.oiplot(data.data,"SPAFREQ","VIS2DATA",cname="VIS2DATA",label="Color by VIS2DATA")
ax.legend()
plt.colorbar(cb,ax=ax,label="VIS2DATA")
ax.autolim()

plt.figure()
ax = plt.subplot(projection='oimAxes')
cb=ax.oiplot(data.data,"PA","VIS2DATA",cname="B",label="Color by B")
ax.legend()
plt.colorbar(cb,ax=ax,label="Baseline Length in (m)")
ax.autolim()

"""
plt.figure()
ax = plt.subplot(projection='oimAxes')
ax.oiplot(data.data,"SPAFREQ","VIS2DATA",cmap="plasma",cname="EFF_WAVE",errorbar=True,kwargs_error=dict(color="green",alpha=0.2))
ax.legend()
ax.autolim()