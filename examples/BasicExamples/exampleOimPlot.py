# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 07:28:17 2022

@author: Ame
"""
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import oimodeler as oim

path = os.path.dirname(oim.__file__)
pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")

files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]
data=oim.oimData(files)


#%%
fig1 = plt.figure()
ax1 = plt.subplot(projection='oimAxes')
ax1.uvplot(data)
plt.savefig(os.path.join(path,os.pardir,"images","ExampleOimPlot_uv.png"))

#%%
fig2 = plt.figure()
ax2 = plt.subplot(projection='oimAxes')
lamcol=ax2.oiplot(data,"SPAFREQ","VIS2DATA" ,xunit="cycles/mas",label="Data",
                cname="EFF_WAVE",cunitmultiplier=1e6,errorbar=True)
                
plt.colorbar(lamcol, ax=ax2,label="$\\lambda$ ($\mu$m)")
ax2.legend()

plt.savefig(os.path.join(path,os.pardir,"images","ExampleOimPlot_v2.png"))

#%%
fig3= plt.figure()
ax3 = plt.subplot(projection='oimAxes')
ax3.oiplot(data,"EFF_WAVE","VIS2DATA",xunitmultiplier=1e6,color="byConfiguration",
               errorbar=True,kwargs_error={"alpha":0.3})

plt.savefig(os.path.join(path,os.pardir,"images","ExampleOimPlot_v2Wl.png"))

#%%
# using the projection='oimAsxes' for all subpltos allow to use oimPlot custom plots
fig4, ax4 = plt.subplots(2,2, subplot_kw=dict(projection='oimAxes'),figsize=(8,8))

#First we plot the uv plan coverage 
ax4[0,0].uvplot(data)

#Then we plot the V² as a function of spa. freq. with a wavelength-colorscale  
lamcol=ax4[0,1].oiplot(data,"SPAFREQ","VIS2DATA" ,xunit="cycles/mas",label="Data",
                    cname="EFF_WAVE",cunitmultiplier=1e6,ls=":",errorbar=True)

#Adding the corresponding colobar of the walvength colorscale
fig4.colorbar(lamcol, ax=ax4[0,1],label="$\\lambda$ ($\mu$m)")

#Legending also work with multicolor plots
ax4[0,1].legend()

# Plotting the V² as a function of the wavlength with errorbars 
#kwargs_error keyword allows to pass keywords to the errorbar plot
ax4[1,0].oiplot(data,"EFF_WAVE","VIS2DATA",xunitmultiplier=1e6,
               errorbar=True,kwargs_error={"alpha":0.1})


# Finally we plot the Closure Phase with a fewstyling options
ax4[1,1].oiplot(data,"SPAFREQ","T3PHI",xunit="cycles/mas",errorbar=True,
               lw=2,ls=":")

#setting log scale as for normal matplotlib plots
ax4[0,1].set_yscale('log')

#autolim allows to directly set ylim to (0,1) for vis and (-180,180) for phase
ax4[1,0].autolim()
ax4[1,1].autolim()

fig4.tight_layout()     

filename=os.path.join(path,os.pardir,"images","ExampleOimPlot_multi.png")
plt.savefig(filename)