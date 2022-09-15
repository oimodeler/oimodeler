# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 07:28:17 2022

@author: Ame
"""
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import fits
import oimodeler as oim



path="D:\\Travail\\GitHub\\oimodeler\\examples\\testData\\ASPRO_MATISSE2\\"
files=[os.path.abspath(os.path.join(path,fi)) for fi in os.listdir(path) if ".fits" in fi]
data=[fits.open(fi,mode="update") for fi in files]


fig, ax = plt.subplots(1,4, subplot_kw=dict(projection='oimAxes'),figsize=(17,4))

ax[0].uvplot(data)
lamcol=ax[1].oiplot(data,"SPAFREQ","VIS2DATA" ,xunit="cycles/mas",color="byWavelength",lw=2,ls=":",errorbar=True)
fig.colorbar(lamcol, ax=ax[1])

ax[2].oiplot(data,"EFF_WAVE","VIS2DATA",errorbar=True,kwargs_error={"alpha":0.1})
ax[3].oiplot(data,"SPAFREQ","T3PHI",xunit="cycles/mas",errorbar=True,lw=2,ls=":",colorTab=["red","green","blue"])


ax[1].set_yscale('log')
#ax[1].autolim()
ax[2].autolim()
ax[3].autolim()
fig.tight_layout()     

