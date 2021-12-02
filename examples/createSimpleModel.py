import oimodeler as oim
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
eg=oim.oimEGauss(fwhm=8,f=0.5,pa=0,elong=2)
r=oim.oimRing(d=5,f=0.5)
er=oim.oimRing(d=10,f=0.5,pa=30,elong=2)


#Building a few models from the components 
m1=oim.oimModel([ud])
m2=oim.oimModel([ud,pt])
m3=oim.oimModel([g])
m4=oim.oimModel([pt,e])
m5=oim.oimModel([eg])
m6=oim.oimModel([r])
m7=oim.oimModel([er])

models=[m1,m2,m3,m4,m5,m6,m7]
names=["UD","UD+PT","G","PT+EL","EG","R","ER"]
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
# Create images of the model from their formula in direct space or using 
# the FT formula

normPow=0.2
dim=256   #Number of pixel for the image
pix=0.1  #size of the pixel in the image in mas
fig,ax=plt.subplots(3,len(models),figsize=(3*len(models),9))
for i,m in enumerate(models):    
    

    #get the image for the direct space formula
    im=m.getImage(dim,pix)
    im=im/np.sum(im)
    ax[0,i].imshow(im,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
                   norm=colors.PowerNorm(gamma=normPow))
    ax[0,i].set_xlabel("$\\alpha$(mas)")
    ax[0,i].set_ylabel("$\\delta$(mas)")
    
    ax[0,i].text(0,-dim/2*pix,names[i],va="bottom",ha="center",color="w",fontsize=20)
    
    #get the image for the FT formula
    im2=m.getImage(dim,pix,fromFT=True)
    im2=im2/np.sum(im2)
    
    ax[1,i].imshow(im2,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
                   norm=colors.PowerNorm(gamma=normPow))
    
   
    
    ax[2,i].imshow(im2-im,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
                   norm=colors.PowerNorm(gamma=normPow))
    ax[2,i].set_xlabel("$\\alpha$(mas)")
    
    
    if i!=0:
        for j in range(3):
            ax[j,i].get_yaxis().set_visible(False)
    else:
        for j in range(3):
            ax[j,i].set_ylabel("$\\delta$(mas)") 
    
    ax[0,i].get_xaxis().set_visible(False)
    ax[1,i].get_xaxis().set_visible(False)
plt.figtext(0.5,0.9,"Images generated using the formula defined in _imageFunction",ha="center", fontsize=20)
plt.figtext(0.5,0.62,"Images generated using the inverse Fourier transform and the formula defined in _visFunction",ha="center", fontsize=20)
plt.figtext(0.5,0.35,"Difference",ha="center", fontsize=20)