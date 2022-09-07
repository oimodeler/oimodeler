# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:16:59 2022

@author: Ame
"""

import oimodeler as oim
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import astropy.units as u
import os
from datetime import datetime
import oifitstools
from astropy.io import fits
import numpy as np

path = os.path.dirname(oim.__file__)
#pathData=os.path.join(path,os.pardir,"examples","testData","FSCMA_MATISSE")
pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE")
files0=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

    
text=["Complex Corr Flux only","Complex Corr Flux + Chi2","Complex Corr Flux + Sim. Data","Full Computation "]
computeChi2=[False,True,False,True]
computeSimulatedData=[False,False,True,True]

ndata=30
dt=np.ndarray([4,ndata])

start_time0 = datetime.now()
for idata in range(ndata):
    
    ud=oim.oimUD(d=20,f=4)
    pt=oim.oimPt(f=6)
    model=oim.oimModel([ud,pt])
    
    files=files0*(idata+1)
    
    sim=oim.OImSimulator(data=files,model=model)    
    sim.data.prepareData()
    
    if idata==0:
        x0=np.size(sim.data.vect_u)
    
    
    navg=1 
    for itype in range(4):
        start_time = datetime.now()
        for i in range(navg): 
            sim.compute(computeChi2=computeChi2[itype],
                        computeSimulatedData=computeSimulatedData[itype])
        end_time = datetime.now()
        dt[itype,idata]=(end_time - start_time).total_seconds() * 1000/navg
   
end_time0 = datetime.now()
print('Full Computation time {:.3f}s'.format((end_time0 - start_time0).total_seconds() ))

#%%
x=np.linspace(1,ndata,ndata)*x0

col=plt.rcParams['axes.prop_cycle'].by_key()['color']
for itype in range(4):
    y=np.poly1d(np.polyfit(x,dt[itype,:],2))(x)
    plt.plot(x,dt[itype,:],label=text[itype],marker="o",ls="",color=col[itype])    
    plt.plot(x,y,color=col[itype])
plt.legend()
plt.xlabel("Number of data points")
plt.ylabel("Computation time (ms)")
plt.xlim(0,np.max(x))
plt.ylim(0,np.max(dt))
txt=""
for c in model.components:
    txt+=c.__str__().split("x")[0]
    txt+="+ "
txt=txt[:-3]
    
plt.title("Computation time for a {} model".format(txt))
filename=os.path.join(path,os.pardir,"images","oimodel_test_simulator_speed.png")
plt.savefig(filename)
