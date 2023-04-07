import os

import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim

path = os.path.dirname(oim.__file__)

#%%
pt = oim.oimPt(f=0.1)
ud = oim.oimUD(d=10,f=0.5)
g  = oim.oimGauss(fwhm=5,f=1)
r  = oim.oimIRing(d=5,f=0.5)


#%%
print(ud)
print(ud.params['d'])

#%%


#Building a few models from the components 
mPt   = oim.oimModel(pt)
mUD   = oim.oimModel(ud)
mG    = oim.oimModel(g)
mR    = oim.oimModel(r)
mUDPt = oim.oimModel(ud,pt)



#%%
params = mUDPt.getParameters()
print(params)

freeParams = mUDPt.getFreeParameters()
print(freeParams)   
 

#%%   

im = mUDPt.getImage(512,0.1)
plt.figure()
plt.imshow(im**0.2)
plt.savefig(os.path.join(path,os.pardir,"images","basicModel_imshow.png"))

#%%
im=mUDPt.getImage(256,0.1,toFits=True)
print(im)
print(im.header)
print(im.data)
im=mUDPt.saveImage("modelImage.fits",256,0.1)


#%%
figImg,axImg,Img=mUDPt.showModel(512,0.1,normPow=0.2, figsize=(5,4),
    savefig=os.path.join(path,os.pardir,"images","basicModel_showModel.png"))



#%%
#Create some spatial frequencies (Baselines from 0 to 300m at 2.1 microns)
wl=2.1e-6
B=np.linspace(0.0,300,num=200)
spf=B/wl
spf0=spf*0
print(spf)

#%%


ccf = mUDPt.getComplexCoherentFlux(spf,spf*0) # East-West baselines

#%%

v = np.abs(ccf)
v=v/v.max()
plt.figure()
plt.plot(spf , v)
plt.xlabel("spatial frequency (cycles/rad)")
plt.ylabel("Visbility")

plt.savefig(os.path.join(path,os.pardir,"images","basicModel_vis0.png"))

#%%
#Some components
models = [mPt,mUD,mG,mR,mUDPt]
mNames=["Point Source","Uniform Disk","Gausian","Ring",
              "Uniform Disk + Point Source"]


fig,ax=plt.subplots(2,len(models),figsize=(3*len(models),6),sharex='row',sharey='row')

for i, m in enumerate(models):
    m.showModel(512,0.1,normPow=0.2,axe=ax[0,i],colorbar=False)
    
    v = np.abs(m.getComplexCoherentFlux(spf,spf*0)) 
    v=v/v.max()
    ax[1,i].plot(spf , v)
    
    ax[0,i].set_title(mNames[i])
    ax[1,i].set_xlabel("sp. freq. (cycles/rad)")
    

fig.savefig(os.path.join(path,os.pardir,"images","basicModel_all.png"))   
