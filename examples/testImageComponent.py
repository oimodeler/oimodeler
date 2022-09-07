import oimodeler as oim
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import astropy.units as u

mas2rad=u.mas.to(u.rad)

# Create some spatial frequencies (Baselines from 0 to 400m at 2.1 microns)
npts=1000
wl=0.65e-6
B=np.linspace(1e-10,100,num=npts)
spf=B/wl
spf0=spf*0



# Creating a Uniform disk image 
dim=800
pixSize=0.1 #mas
D=10 #mas
xx=np.tile((np.arange(dim)-dim/2)*pixSize,(dim,1))
yy=np.transpose(xx)
data=((xx**2+yy**2)<(D/2)**2).astype(float)
data=data/np.sum(data) 


# A few quantities useful to know if the object is well sampled

spfMax=1/pixSize  #max spatial frequency that can Be probed in cycles/mas
Bmax=spfMax*wl/mas2rad # Corresponding maximum baseline for wl
print("Maximum baseline with wl={}um and pixSize={}mas : {:.1f}m".format(wl*1e6,pixSize,Bmax))

fov=pixSize*dim # Field of view in mas
dspf= 1/fov # Sampling in spatial frequency 
dB=dspf*wl/mas2rad # sampling in B (m)
print("Baseline sampling with fov={}mas: {:.1f}m".format(fov,dB))



ud0=oim.oimUD(f=0.2,d=2,x=10)

# Building the model with the direct image of the UD  
imUD=oim.oimComponentImage(pixSize=pixSize)
imUD.setImage(data)
modelDirect=oim.oimModel([imUD,ud0])

# Building the model with the FT analytical component
ud=oim.oimUD(d=10)#,x=10)
modelFT=oim.oimModel([ud,ud0])

vDirect=np.abs(modelDirect.getComplexCoherentFlux(spf,spf0))

vFT=np.abs(modelFT.getComplexCoherentFlux(spf,spf0)) 


fig,ax=plt.subplots(2,1,figsize=(10,7))
plt.margins(0,0)

ax[0].plot(B,vDirect/vDirect[0],label="Direct image")
ax[0].plot(B,vFT/vFT[0],label="FT based")

ax[0].set_ylabel("Visibility") 
ax[0].legend()
 

residuals=(vDirect/vDirect[0]-vFT/vFT[0])/(vFT/vFT[0])*100
mresiduals=np.median(np.abs(residuals))
print("Mean residual between direct and FT : {0}%".format(mresiduals))

ax[1].scatter(B,residuals)
ax[1].set_ylabel("Residuals (%)") 
ax[1].set_xlabel("B (m)")
ax[1].set_ylim(-5*mresiduals,5*mresiduals)


#%%
# Create images of the models  from their FT formulas

normPow=1
fig,ax=plt.subplots(1,2,figsize=(10,5))

names=["2 FT UDs","Image UD + FT UD"]
for i,mi in enumerate([modelDirect,modelFT]):

    im=mi.getImage(dim,pixSize,fromFT=True)
    im=im/np.sum(im)
    ax[i].imshow(im,extent=[dim/2*pixSize,-dim/2*pixSize,-dim/2*pixSize,
                            dim/2*pixSize],norm=colors.PowerNorm(gamma=normPow))
    ax[i].set_xlabel("$\\alpha$(mas)")    
    
    ax[i].text(0,0.4*dim*pixSize,names[i],color="w",va="top",
               ha="center",fontsize=15)
    
ax[0].set_ylabel("$\\delta$(mas)")
