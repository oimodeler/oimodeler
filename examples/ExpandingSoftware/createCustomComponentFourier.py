"""
In this example we show how to implement a new model component using a formula
in the Fourier plan. The component derives from the  oimComponentFourier class.
The Fourier formula should be implemented in  _visFunction and optionally the 
formula in the image plan can be implemented using  _imageFunction.

The component implemented here is a rectangular box.

After the implementation it is used to build  two models: one with a single box
and a more complex one consisting of three boxes and a uniform disk. 
Complex visibility for North-South and East-West baselines between 0 and 100m 
are computed for both models and plotted as well as images.

"""

import oimodeler as oim
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import astropy.units as u

##############################################################################
"""
Example of implementation of a new component defined by formula in Fourier Plan
The 2D-box function (very useful for squared stars!)
"""
class box(oim.oimComponentFourier):
    name="2D Box"
    shortname = "BOX"
    def __init__(self,**kwargs): 
         # First call the __init__ function from the parent class
         super().__init__(**kwargs)
         # Then define the component parameters  
         self.params["dx"]=oim.oimParam(name="dx", value=1,description="Size in x",unit=u.mas)
         self.params["dy"]=oim.oimParam(name="dy", value=1,description="Size in y",unit=u.mas)
         # Finally call the _eval function that allow the parameters to be processed         
         self._eval(**kwargs)
    
    """
    Implementation of the complex visibility function that will be called when
    using the getComplexCoherentFlux function. The component parameters should 
    be called with (wl,t) to allow parameter chromaticity and time dependence.
  
    The parameter also contains a unit keyword allowing conversion.
    
    Here a simple 2D-sinc function as the FT of our 2D box component.
    Note that the size arte converted from the given unit (usually mas) to rad    
    """
    def _visFunction(self,ucoord,vcoord,rho,wl,t):
        
        x=self.params["dx"](wl,t)*self.params["dx"].unit.to(u.rad)*ucoord
        y=self.params["dy"](wl,t)*self.params["dy"].unit.to(u.rad)*vcoord
        
        return np.sinc(x)*np.sinc(y) 
    """
    Implementation of the image function that will be called whenusing the 
    getImage function.
    If not implemented the model will use the Fourier based formula to compute 
    the image. It will also be the case if the keyword fromFT is set to True
    when the getImage is called.
    
    The box image function is simply implemented with logical operations.
    """
    def _imageFunction(self,xx,yy,wl,t):
        
        return ((np.abs(xx)<=self.params["dx"](wl,t)/2) &
                (np.abs(yy)<=self.params["dy"](wl,t)/2)).astype(float)



##############################################################################

# Create components with newly define box class
b1=box(dx=40,dy=10)
#Note that the x, y and f parameters are implemented for all components
b2=box(dx=2,dy=2,x=20,y=0,f=0.5)
# So is the pa for elliptic model but that is also valid for rotation
b3=box(dx=10,dy=20,x=-30,y=10,pa=50,f=10)


# Building the model from the 1st box component 
model1=oim.oimModel([b1])

#Building a model with all the boxes 
#Let's even  add a circle to this list of boxes
c=oim.oimUD(d=10,x=-30,y=-10)
model2=oim.oimModel([b1,b2,b3,c])



#%%
"""
Creating complex coherent flux and extracting visibility from the models
and plotting them as the function of the spatial frequency.
"""

# Create some spatial frequencies (Baselines from 0 to 100m at 2.1 microns)
lam=2.1e-6
B=np.linspace(0.01,100,num=1000)
spf=B/lam
spf0=spf*0


fig,ax=plt.subplots(1,1,figsize=(10,7))

#Compute visibilities for the models and spatial frequencies
v11=np.abs(model1.getComplexCoherentFlux(spf,spf0)) # East-West baselines
v21=np.abs(model1.getComplexCoherentFlux(spf0,spf)) # North-South baselines
Ftot1=v11[0] # total flux for normalization 


v12=np.abs(model2.getComplexCoherentFlux(spf,spf0)) # East-West baselines
v22=np.abs(model2.getComplexCoherentFlux(spf0,spf)) # North-South baselines
Ftot2=v12[0] # total flux for normalization 

plt.plot(B/lam*u.mas.to(u.rad),v11/Ftot1,color="tab:blue",
         label="First Box : East-West")
plt.plot(B/lam*u.mas.to(u.rad),v21/Ftot1,ls=":",color="tab:blue",
         label="First Box : North-South")
plt.plot(B/lam*u.mas.to(u.rad),v12/Ftot2,color="tab:orange",
         label="All Components : East-West")
plt.plot(B/lam*u.mas.to(u.rad),v22/Ftot2,ls=":",color="tab:orange",
         label="All Components: North-South")
plt.xlabel("B/$\\lambda$ (cycle/mas)")
plt.ylabel("Visibility")  
plt.legend()


#%%
"""
Creating images of the models from their formula in direct space or using the 
FT formula
"""

normPow=0.5  # exponent for normalization of the image 

dim=1024   #Number of pixel for the image
pix=0.2  #size of the pixel in the image in mas

fig,ax=plt.subplots(2,3,figsize=(18,12))
models=[model1,model2]
for i,model in enumerate(models):
    #get the image for the direct space formula
    im0=model.getImage(dim,pix)
    im0=im0/np.sum(im0)
    ax[i,0].imshow(im0,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
                   norm=colors.PowerNorm(gamma=normPow))
    ax[i,0].text(0,0.9*dim/2*pix,"Direct image",va="top",ha="center",
               color="w",fontsize=14)
    
    #get the image for the Fourier space formula
    im1=model.getImage(dim,pix,fromFT=True)
    im1=im1/np.sum(im1)
    ax[i,1].imshow(im1,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
                   norm=colors.PowerNorm(gamma=normPow))
    ax[i,1].text(0,0.9*dim/2*pix,"Image ocmputed for the Fourier formula",va="top",ha="center",
               color="w",fontsize=14)
    # Plot the difference between the two to check
    im2=im0-im1
    ax[i,2].imshow(im2,extent=[dim/2*pix,-dim/2*pix,-dim/2*pix,dim/2*pix],
                   norm=colors.PowerNorm(gamma=normPow))
    ax[i,2].text(0,0.9*dim/2*pix,"Difference",va="top",ha="center",
               color="w",fontsize=14)
    
    ax[i,0].set_ylabel("$\\delta$(mas)")
  
    if i==1:
        for j in range(3):
            ax[i,j].set_xlabel("$\\alpha$(mas)")
    