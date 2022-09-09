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
from astropy.io import fits

path = os.path.dirname(oim.__file__)

# Path to a fake MATISSE-L-band binary observation (3 oifits) created with ASPRO
pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")
files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

#Building a oimodeler model with the same parameters
ud=oim.oimUD(d=5,f=1,y=10)
pt=oim.oimPt(f=1)
model=oim.oimModel([ud,pt])

# Creating the simulator with the filename list and the model
sim=oim.OImSimulator(data=files,model=model)

# Preparing data (building vectors of corrdinates and structure of data types)
sim.data.prepareData()

#Computing the complex corr flux from the model at the data spatial freq
# with option to compute chi2 and final simulated data in oifits format
start_time = datetime.now()
sim.compute(computeChi2=True,computeSimulatedData=True)
end_time = datetime.now()
print('Simulation computation time = {:.3f}ms'.format((end_time - start_time).total_seconds()*1000 ))

#Printing data model chi2r
print("Chi2r = {}".format(sim.chi2r))


#%%
# plotting  data and simulated data
#TODO : to be remplaced by some internal functions


fig,ax=plt.subplots(6,1,sharex=True,figsize=(7,7))

nfiles=len(sim.data.data)
  
for ifile,d in enumerate(sim.data.data):

        dsim=sim.simulatedData.data[ifile]
        spFreqV=oim.getSpaFreq(d,unit="cycles/mas")   
        V2=d["OI_VIS2"].data["VIS2DATA"]
        V2sim=dsim["OI_VIS2"].data["VIS2DATA"]     
        nBV2=np.shape(V2)[0]    

        V=d["OI_VIS"].data["VISAMP"]
        Vsim=dsim["OI_VIS"].data["VISAMP"]     
        nBV=np.shape(V)[0]    
        
        phi=d["OI_VIS"].data["VISPHI"]
        phiSim=dsim["OI_VIS"].data["VISPHI"]     
        nBPhi=np.shape(phiSim)[0]    
 
        spFreqCP=oim.getSpaFreq(d,"OI_T3",unit="cycles/mas")
        CP=d["OI_T3"].data["T3PHI"]
        CPsim=dsim["OI_T3"].data["T3PHI"]        
        flagCP=d["OI_T3"].data["FLAG"]    
        nBCP=np.shape(CP)[0]      
        
        T3=d["OI_T3"].data["T3AMP"]
        T3sim=dsim["OI_T3"].data["T3AMP"]        
        flagCP=d["OI_T3"].data["FLAG"]    
        nBT3=np.shape(T3)[0]
        
        try:
            flx=d["OI_FLUX"].data["FLUXDATA"]
            flxSim=dsim["OI_FLUX"].data["FLUXDATA"]     
            nBFlx=np.shape(flx)[0]    

        except:
            pass

        lam=d["OI_WAVELENGTH"].data["EFF_WAVE"]
        nlam=np.size(lam)
           
        for iB in range(nBV2):
            ax[0].plot(spFreqV[iB,:],V2[iB,:],color="red")
            ax[0].plot(spFreqV[iB,:],V2sim[iB,:],color="blue")
                              
        for iB in range(nBV):
            ax[1].plot(spFreqV[iB,:],V[iB,:],color="red")
            ax[1].plot(spFreqV[iB,:],Vsim[iB,:],color="blue")    

        for iB in range(nBPhi):

            ax[2].plot(spFreqV[iB,:],phi[iB,:],color="red")
            ax[2].plot(spFreqV[iB,:],phiSim[iB,:],color="blue")      
                                
                
        for iB in range(nBT3):
            ax[3].plot(spFreqCP[iB,:],T3[iB,:],color="red")
            ax[3].plot(spFreqCP[iB,:],T3sim[iB,:],color="blue")        
            
                
        for iB in range(nBCP):
            ax[4].plot(spFreqCP[iB,:],CP[iB,:],color="red")
            ax[4].plot(spFreqCP[iB,:],CPsim[iB,:],color="blue")              
            
        try:
            ax[5].plot(spFreqV[iB,:],flx[iB,:],color="red")
            ax[5].plot(spFreqV[iB,:],flxSim[iB,:],color="blue")     
        except:
            pass
        
ax[0].set_ylim(0,1)
ax[1].set_ylim(0,1)
ax[2].set_ylim(-180,180)
ax[4].set_ylim(-180,180)

ax[0].set_ylabel("VIS2DATA")
ax[1].set_ylabel("VISAMP")
ax[2].set_ylabel("VISPHI (deg)")
ax[3].set_ylabel("T3AMP")
ax[4].set_ylabel("T3PHI (deg)")
ax[5].set_ylabel("FLUXADATA")

for i in range(len(ax)-1):
    ax[i].get_xaxis().set_visible(False)


ax[-1].set_xlabel("B/$\lambda$ (cycles/mas)")

fig.suptitle("Data Simulated with ASPRO (red) and oimodeler (blue)")
fig.tight_layout() 
#%%
#Saving the plot

filename=os.path.join(path,os.pardir,"images","oimodel_Create_simulator_data.png")
plt.savefig(filename)

