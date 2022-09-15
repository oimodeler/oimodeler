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

path = os.path.dirname(oim.__file__)
pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")

files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]
data=[fits.open(fi,mode="update") for fi in files]


#%%
# using the projection='oimAsxes' for all subpltos allow to use oimPlot custom plots
fig, ax = plt.subplots(2,2, subplot_kw=dict(projection='oimAxes'),figsize=(8,8))

#First we plot the uv plan coverage 
ax[0,0].uvplot(data)

#Then we plot the V² as a function of spa. freq. with a wavelength-colorscale  
lamcol=ax[0,1].oiplot(data,"SPAFREQ","VIS2DATA" ,xunit="cycles/mas",label="Data",
                    cname="EFF_WAVE",cunitmultiplier=1e6,ls=":",errorbar=True)

#Adding the corresponding colobar of the walvength colorscale
fig.colorbar(lamcol, ax=ax[0,1],label="$\\lambda$ ($\mu$m)")

#Legending also work with multicolor plots
ax[0,1].legend()

# Plotting the V² as a function of the wavlength with errorbars 
#kwargs_error keyword allows to pass keywords to the errorbar plot
ax[1,0].oiplot(data,"EFF_WAVE","VIS2DATA",xunitmultiplier=1e6,
               errorbar=True,kwargs_error={"alpha":0.1})


# Finally we plot the Closure Phase with a fewstyling options
ax[1,1].oiplot(data,"SPAFREQ","T3PHI",xunit="cycles/mas",errorbar=True,
               lw=2,ls=":")

#setting log scale as for normal matplotlib plots
ax[0,1].set_yscale('log')

#autolim allows to directly set ylim to (0,1) for vis and (-180,180) for phase
ax[1,0].autolim()
ax[1,1].autolim()

fig.tight_layout()     


#%%
#Saving the plot
filename=os.path.join(path,os.pardir,"images","oimodel_test_plots.png")
plt.savefig(filename)