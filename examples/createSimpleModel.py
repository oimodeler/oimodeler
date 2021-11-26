import oimodeler as oim
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

mas2rad=u.mas.to(u.rad)

#Create some spatial frequencies (Baselines from 0 to 400m at 2.1 microns)
lam=2.1e-6
B=np.linspace(0.01,400,num=1000)
spf=B/lam
spf0=spf*0


#Some components
pt=oim.oimPt(f=0.1)
ud=oim.oimUD(d=3,f=0.5)
g=oim.oimGauss(fwhm=5,f=1)
e=oim.oimEllipse(d=10,f=0.5,pa=0,elong=2)

#Building a few models from the components 
m1=oim.oimModel([ud])
m2=oim.oimModel([ud,pt])
m3=oim.oimModel([g])
m4=oim.oimModel([pt,e])

models=[m1,m2,m3,m4]
names=["UD","UD+PT","G","PT+EL"]
cols= plt.rcParams['axes.prop_cycle'].by_key()['color']


fig,ax=plt.subplots(1,1,figsize=(10,7))
for i,m in enumerate(models):
    #Compute visibilities for the models and spatial frequencies
    v1=np.abs(m.getComplexCoherentFlux(spf,spf0)) # East-West baselines
    v2=np.abs(m.getComplexCoherentFlux(spf0,spf)) # North-South baselines
    
    Ftot=v1[0]
    plt.plot(B/lam*mas2rad,v1/Ftot,color=cols[i],label=names[i])
    plt.plot(B/lam*mas2rad,v2/Ftot,color=cols[i],ls=":")
plt.xlabel("B/$\\lambda$ (cycle/mas)")
plt.ylabel("Visibility")  
plt.legend()
plt.margins(0,0)



#%%
# Create direct images of the model
dim=128   #Number of pixel for the image
pix=0.1  #size of the pixel in the image in mas
fig,ax=plt.subplots(1,len(models),figsize=(5*len(models),5))
for i,m in enumerate(models):    
    im=m.getImage(dim,pix)
    im=im/np.sum(im)
    ax[i].imshow(im,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix])
    ax[i].set_xlabel("$\\alpha$(mas)")
    ax[i].set_ylabel("$\\delta$(mas)")