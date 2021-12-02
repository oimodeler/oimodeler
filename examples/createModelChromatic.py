import oimodeler as oim
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import astropy.units as u
import os

mas2rad=u.mas.to(u.rad)

#Create some spatial frequencies (Baselines from 0 to 400m at 2.1 microns)
lam=2.1e-6
B=np.linspace(0.01,400,num=1000)
spf=B/lam
spf0=spf*0


#Some components
gd=oim.oimGauss(fwhm=oim.oimInterpWl(wl=[3e-6,3.5e-6,4e-6],value=[1,4,4]),
                 f=oim.oimInterpWl([3e-6,4e-6],[1,4]))
ud=oim.oimUD(d=0.7,f=1)
ud2=oim.oimUD(x=0,y=5,d=1.5,f=oim.oimInterpWl(wl=[3e-6,4e-6],value=[1,0.]))
ud3=oim.oimUD(x=-5,y=2,d=0.5,f=oim.oimInterpWl(wl=[3e-6,4e-6],value=[0.1,1]))
eg=oim.oimEGauss(fwhm=1,elong=1.5,pa=oim.oimInterpWl([3e-6,4e-6],[20,60]),f=oim.oimInterpWl([3e-6,4e-6],[1,0.1]))
er=oim.oimERing(din= 8,dout=oim.oimInterpWl([3e-6,4e-6],[15,20]),elong=2,pa=0)



#linking some parameters, actually replacing them
er.params["elong"]=eg.params["elong"]
er.params["pa"]=eg.params["pa"]

#Building a few models from the components 
m1=oim.oimModel([gd])
m2=oim.oimModel([ud,gd])
m3=oim.oimModel([ud,ud2,ud3])
m4=oim.oimModel([eg,er])
models=[m1,m2,m3,m4]

nB=50 #number of baselines 
nwl=100 #number of walvengths

#Create some spatial frequencies
wl=np.linspace(3e-6,4e-6,num=nwl)
B=np.linspace(1,400,num=nB)
Bs=np.tile(B,(nwl,1)).flatten()
wls=np.transpose(np.tile(wl,(nB,1))).flatten()
spf=Bs/wls
spf0=spf*0



fig,ax=plt.subplots(4,4,figsize=(8,8))

for i,m in enumerate(models):
    
    #Compute visibilities for the models and spatial frequencies
    cf1=np.abs(m.getComplexCoherentFlux(spf,spf0,wls)).reshape(len(wl),len(B))
    cf2=np.abs(m.getComplexCoherentFlux(spf0,spf,wls)).reshape(len(wl),len(B))
    
    
    for iwl,wli in enumerate(wl):
        ax[2,i].plot(B/wli*mas2rad,cf1[iwl,:]/cf1[iwl,0],color=plt.cm.plasma(iwl/nwl),alpha=0.5)
        ax[3,i].plot(B/wli*mas2rad,cf2[iwl,:]/cf2[iwl,0],color=plt.cm.plasma(iwl/nwl),alpha=0.5)        
   
    sc=ax[3,i].scatter(spf0*mas2rad,cf2,c=wls*1e6,s=0.1,cmap="plasma")
    ax[3,i].set_xlabel("B/$\\lambda$ (cycles/mas)")
    ax[2,i].margins(0,0)
    ax[2,i].set_ylim(0,1)
    ax[3,i].margins(0,0)
    ax[3,i].set_ylim(0,1)
    #ax[2,i].set_title(model_names[i],fontsize=8)
    if i!=0:
        ax[2,i].get_yaxis().set_visible(False)
        ax[3,i].get_yaxis().set_visible(False) 
ax[2,0].set_ylabel("Vis. (East-West)")         
ax[3,0].set_ylabel("Vis. (North-South)")  

dim=256   #Number of pixels for the image
pix=0.1  #size of the pixel in the image in mas
wls=[3e-6,4e-6] #wavelengths array for the images

for i,m in enumerate(models):    
    im=m.getImage(dim,pix,wls,fromFT=True)
    for iwl,wli in enumerate(wls):     
        imi=im[iwl,:,:]
        imi/=np.max(imi)
        imshow=ax[iwl,i].imshow(imi,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
                                norm=colors.TwoSlopeNorm(vmin=0., vcenter=1e-1, vmax=1)
                                ,cmap="plasma", interpolation=None)
        ax[iwl,i].text(0,dim*pix/2,"\n{:.1f}$\\mu$m".format(wli*1e6),
                       color="w",va="top",ha="center")
        if i!=0:
            ax[iwl,i].get_yaxis().set_visible(False)
        ax[iwl,i].set_xlabel("$\\alpha$(mas)")          
        if i==0:
            ax[iwl,i].set_ylabel("$\\delta$(mas)")
            
fig.tight_layout()            
fig.colorbar(imshow, ax=ax[:2,:].ravel().tolist(),label="Normalized Intensity")
norm = colors.Normalize(vmin=np.min(wl),vmax=np.max(wl))
sm = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
fig.colorbar(sc, ax=ax[2:,:].ravel().tolist(),label="$\\lambda$ ($\\mu$m)")


#fig.savefig("createModelChromatic.png")