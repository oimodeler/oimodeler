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

path = os.path.dirname(oim.__file__)
#pathData=os.path.join(path,os.pardir,"examples","testData","FSCMA_MATISSE")
#pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE")
pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")
files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]


#ASPRO T3AMP set to 0
"""
for fi in files:
    d=fits.open(fi,mode='update')
    d["OI_T3"].data["T3AMP"]*=0
    d.flush()
"""

#eg=oim.oimEGauss(fwhm=1,elong=1.5,pa=oim.oimInterpWl([3e-6,4e-6],[20,60]),f=oim.oimInterpWl([3e-6,4e-6],[1,0.1]))
#er=oim.oimERing(din= 8,dout=oim.oimInterpWl([3e-6,4e-6],[15,20]),elong=2,pa=0)

#ud=oim.oimUD(d=20,f=4)
#pt=oim.oimPt(f=6)
#model=oim.oimModel([ud,pt])


ud=oim.oimUD(d=5,f=1,y=10)
pt=oim.oimPt(f=1)
model=oim.oimModel([ud,pt])


sim=oim.OImSimulator(data=files,model=model)

sim.data.prepareData()

n=1



start_time = datetime.now()
for k in range(n): 
    sim.compute(computeChi2=False,computeSimulatedData=False)
end_time = datetime.now()
print('complexCorrFlux only {:.3f}ms'.format((end_time - start_time).total_seconds() * 1000/n))


start_time = datetime.now()
for k in range(n): 
    sim.compute(computeChi2=True,computeSimulatedData=False)
end_time = datetime.now()
print('CCF+Chi2  {:.3f}ms'.format((end_time - start_time).total_seconds() * 1000/n))


start_time = datetime.now()
for k in range(n): 
    sim.compute(computeChi2=False,computeSimulatedData=True)
end_time = datetime.now()
print('CCF+SimData  {:.3f}ms'.format((end_time - start_time).total_seconds() * 1000/n))


start_time = datetime.now()
for k in range(n): 
    sim.compute(computeChi2=True,computeSimulatedData=True)
end_time = datetime.now()
print('Full Computation {:.3f}ms'.format((end_time - start_time).total_seconds() * 1000/n))



#%%

fig,ax=plt.subplots(6,1,sharex=True)

nfiles=len(sim.data.data)
  
for ifile,d in enumerate(sim.data.data):

        dsim=sim.simulatedData.data[ifile]
        spFreqV=oifitstools.getSpaFreq(d,unit="cycles/mas")   
        V2=d["OI_VIS2"].data["VIS2DATA"]
        V2sim=dsim["OI_VIS2"].data["VIS2DATA"]     
        nBV2=np.shape(V2)[0]    

        V=d["OI_VIS"].data["VISAMP"]
        Vsim=dsim["OI_VIS"].data["VISAMP"]     
        nBV=np.shape(V)[0]    
        
        phi=d["OI_VIS"].data["VISPHI"]
        phiSim=dsim["OI_VIS"].data["VISPHI"]     
        nBPhi=np.shape(phiSim)[0]    
 
        spFreqCP=oifitstools.getSpaFreq(d,"OI_T3",unit="cycles/mas")
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
            for ilam2 in range(nlam-2):
                ax[0].plot(spFreqV[iB,ilam2:ilam2+2],V2[iB,ilam2:ilam2+2],color="red")#col(lam[ilam2]))
                ax[0].plot(spFreqV[iB,ilam2:ilam2+2],V2sim[iB,ilam2:ilam2+2],color="blue")#col(lam[ilam2]))
                
                
                   
        for iB in range(nBV):
            for ilam2 in range(nlam-2):
                ax[1].plot(spFreqV[iB,ilam2:ilam2+2],V[iB,ilam2:ilam2+2],color="red")#col(lam[ilam2]))
                ax[1].plot(spFreqV[iB,ilam2:ilam2+2],Vsim[iB,ilam2:ilam2+2],color="blue")#col(lam[ilam2]))        

        for iB in range(nBPhi):
            for ilam2 in range(nlam-2):
                ax[2].plot(spFreqV[iB,ilam2:ilam2+2],phi[iB,ilam2:ilam2+2],color="red")#col(lam[ilam2]))
                ax[2].plot(spFreqV[iB,ilam2:ilam2+2],phiSim[iB,ilam2:ilam2+2],color="blue")#col(lam[ilam2]))        
                                
                
        for iB in range(nBT3):
            for ilam2 in range(nlam-2):
                ax[3].plot(spFreqCP[iB,ilam2:ilam2+2],T3[iB,ilam2:ilam2+2],color="red")#col(lam[ilam2]))
                ax[3].plot(spFreqCP[iB,ilam2:ilam2+2],T3sim[iB,ilam2:ilam2+2],color="blue")#col(lam[ilam2]))                
            
                
        for iB in range(nBCP):
            for ilam2 in range(nlam-2):
                ax[4].plot(spFreqCP[iB,ilam2:ilam2+2],CP[iB,ilam2:ilam2+2],color="red")#col(lam[ilam2]))
                ax[4].plot(spFreqCP[iB,ilam2:ilam2+2],CPsim[iB,ilam2:ilam2+2],color="blue")#col(lam[ilam2]))                
            
        try:
            for iB in range(nBFlx):
                for ila52 in range(nlam-2):
                    ax[5].plot(spFreqV[iB,ilam2:ilam2+2],flx[iB,ilam2:ilam2+2],color="red")#col(lam[ilam2]))
                    ax[5].plot(spFreqV[iB,ilam2:ilam2+2],flxSim[iB,ilam2:ilam2+2],color="blue")#col(lam[ilam2]))        
        except:
            pass
ax[0].legend()
ax[0].set_ylim(0,1)
ax[1].set_ylim(0,1)
ax[2].set_ylim(-180,180)

ax[4].set_ylim(-180,180)
#ax[4].set_ylim(0,20)


ax[0].set_ylabel("VIS2DATA")
ax[1].set_ylabel("VISAMP")
ax[2].set_ylabel("VISPHI")
ax[3].set_ylabel("T3AMP")
ax[4].set_ylabel("T3PHI")
ax[5].set_ylabel("FLUXADATA")

for i in range(len(ax)-1):
    ax[i].get_xaxis().set_visible(False)


ax[-1].set_xlabel("B/$\lambda$ (cycles/mas)")

fig.suptitle("Data Simulated with ASPRO (red) and oimmodeler (blue)")
filename=os.path.join(path,os.pardir,"images","oimodel_Create_simulator_data.png")
plt.savefig(filename)

print("Chi2r={}".format(sim.chi2r))
