import oimodeler as oim
import numpy as np
import matplotlib.pyplot as plt
import os

path = os.path.dirname(oim.__file__)



nB=50 #number of baselines 
nwl=100 #number of walvengths

#Create some spatial frequencies
wl=np.linspace(3e-6,4e-6,num=nwl)
B=np.linspace(1,400,num=nB)
Bs=np.tile(B,(nwl,1)).flatten()
wls=np.transpose(np.tile(wl,(nB,1))).flatten()
spf=Bs/wls
spf0=spf*0



#%%
g=oim.oimGauss(fwhm=oim.oimInterpWl([3e-6,4e-6],[1,4]))
    
print(g.params['fwhm']([3e-6,3.5e-6,4e-6,4.5e-6]))
