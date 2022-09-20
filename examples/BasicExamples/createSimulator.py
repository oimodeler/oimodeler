# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""

import oimodeler as oim
import matplotlib.pyplot as plt
import os
from datetime import datetime



path = os.path.dirname(oim.__file__)

# Path to a fake MATISSE-L-band binary observation (3 oifits) created with ASPRO
pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")
files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

#Building a oimodeler model with the same parameters
ud=oim.oimUD(d=5,f=1,y=10)
pt=oim.oimPt(f=1)
model=oim.oimModel([ud,pt])

# Creating the simulator with the filename list and the model
sim=oim.oimSimulator(data=files,model=model)

# Preparing data (building vectors of coordinates and structure of data types)
sim.data.prepareData()

#Computing the complex corr flux from the model at the data spatial freq
# with option to compute chi2 and final simulated data in oifits format
start_time = datetime.now()
sim.compute(computeChi2=True,computeSimulatedData=True)
end_time = datetime.now()
print('Simulation computation time = {:.3f}ms'.format((end_time - start_time)
                                                      .total_seconds()*1000 ))

#Printing data model chi2r
print("Chi2r = {}".format(sim.chi2r))


#%%
# plotting  data and simulated data

#list of data type to be plotted
arr=["VIS2DATA","VISAMP","VISPHI","T3AMP","T3PHI"]

#Set the projection to oimAxes for all subplots to use oimodeler custom plots
fig,ax=plt.subplots(len(arr),1,sharex=True,figsize=(8,6),
                    subplot_kw=dict(projection='oimAxes'))

plt.subplots_adjust(left=0.09,top=0.98,right=0.98,hspace=0.14)

# Ploting loop :  plotting data and simulated data for each data type in arr
for iax,axi in enumerate(ax):
    
    #plotting the data with wavelength colorscale + errorbars vs spatial frequencies
    scale=axi.oiplot(sim.data.data,"SPAFREQ",arr[iax] ,xunit="cycles/mas",
            cname="EFF_WAVE",cunitmultiplier=1e6,lw=2,cmap="coolwarm",
            errorbar=True,label="ASPRO")

    #over-plotting the simulated data as a dotted line  vs spatial frequencies
    axi.oiplot(sim.simulatedData.data,"SPAFREQ",arr[iax] ,xunit="cycles/mas",
            color="k",ls=":",lw=1,label="oimodeler")
    
    if axi!=ax[-1]: axi.get_xaxis().set_visible(False)
    if axi==ax[0]:axi.legend()
    
    #automatic ylim => 0-1 for visibilties, -180,180 for phases
    axi.autolim()
    
#Create a colorbar for the data plotted with wavelength colorscale option
fig.colorbar(scale, ax=ax.ravel().tolist(),label="$\\lambda$ ($\mu$m)")

#%%
#Save the plot
filename=os.path.join(path,os.pardir,"images","oimodel_Create_simulator_data.png")
plt.savefig(filename)

