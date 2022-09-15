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



path="D:\\Travail\\GitHub\\oimodeler\\examples\\testData\\FSCMa_MATISSE\\"
files=[os.path.abspath(os.path.join(path,fi)) for fi in os.listdir(path) if ".fits" in fi]
data=[fits.open(fi,mode="update") for fi in files]


fig, ax = plt.subplots(2,2, subplot_kw=dict(projection='oimAxes'),figsize=(8,8))

ax[0,0].uvplot(data)
lamcol=ax[0,1].oiplot(data,"SPAFREQ","VIS2DATA" ,xunit="cycles/mas",
                    cname="EFF_WAVE",cunitmultiplier=1e6,ls=":",errorbar=True)
fig.colorbar(lamcol, ax=ax[0,1],label="$\\lambda$ ($\mu$m)")

ax[1,0].oiplot(data,"EFF_WAVE","VIS2DATA",xunitmultiplier=1e6,errorbar=True,kwargs_error={"alpha":0.1})
ax[1,1].oiplot(data,"SPAFREQ","T3PHI",xunit="cycles/mas",errorbar=True,lw=2,ls=":",colorTab=["red","green","blue"])


ax[0,1].set_yscale('log')
#ax[1].autolim()
ax[1,0].autolim()
ax[1,1].autolim()
fig.tight_layout()     

