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
import os

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import oimodeler as oim


path = os.path.dirname(oim.__file__)

##############################################################################
"""
Example of implementation of a new component defined by formula in Fourier Plan
The 2D-box function (very useful for squared stars!)
"""
class oimBox(oim.oimComponentFourier):
    name="2D Box"
    shortname = "BOX"
    def __init__(self,**kwargs): 
         # First call the __init__ function from the parent class
         super().__init__(**kwargs)
         # Then define the component parameters  
         self.params["dx"]=oim.oimParam(name="dx", value=1,description="Size in x",unit=u.mas)
         self.params["dy"]=oim.oimParam(name="dy", value=1,description="Size in y",unit=u.mas)
         # Finally call the _eval method that allows the parameters to be processed         
         self._eval(**kwargs)
    
    """
    Implementation of the complex visibility function that will be called when
    using the getComplexCoherentFlux function. The component parameters should 
    be called with (wl,t) to allow parameter chromaticity and time dependence.
  
    The parameter also contains a unit keyword allowing conversion.
    
    Here a simple 2D-sinc function as the FT of our 2D box component.
    Note that the size are converted from the given unit (usually mas) to rad    
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

#%%

b1=oimBox(dx=40,dy=10)
m1=oim.oimModel([b1])

fig, ax = plt.subplots(1,2,figsize=(10,5))
m1.showModel(512,0.2,axe=ax[0],colorbar=False)
m1.showModel(512,0.2,axe=ax[1],fromFT=True,colorbar=False)
ax[0].set_title("Image with _imageFunction")
ax[1].set_title("Image with FFT of _visFunction")
fig.savefig(os.path.join(path,os.pardir,"images","customCompBox1Image.png"))

#%%

b2=oimBox(dx=2,dy=2,x=-20,y=0,f=0.5)
b3=oimBox(dx=10,dy=20,x=30,y=-10,pa=-40,f=10)
c=oim.oimUD(d=10,x=30,y=10)
m2=oim.oimModel([b1,b2,b3,c])
m2.showModel(512,0.2,colorbar=False,figsize=(5,5),
    savefig=os.path.join(path,os.pardir,"images","customCompBoxesImage.png"))




#%%
b4=oimBox(dx=oim.oimInterp("wl", wl=[2e-6,2.4e-6], values=[5,10]),dy=2,x=20,y=0,f=0.5)
b4.params['dy']=oim.oimParamLinker(b4.params['dx'],'mult',4)
    
m3=oim.oimModel([b4])

m3.showModel(512,0.2,wl=[2e-6,2.2e-6,2.4e-6],colorbar=False,swapAxes=True,
    savefig=os.path.join(path,os.pardir,"images","customCompChromBoxImages.png"))



#%%
"""
Creating complex coherent flux and extracting visibility from the models
and plotting them as the function of the spatial frequency.
"""


nB = 200  # number of baselines
nwl = 50  # number of walvengths

# Create some spatial frequencies
wl = np.linspace(2e-6, 2.5e-6, num=nwl)
B = np.linspace(1, 100, num=nB)
Bs = np.tile(B, (nwl, 1)).flatten()
wls = np.transpose(np.tile(wl, (nB, 1))).flatten()
spf = Bs/wls
spf0 = spf*0

fig,ax=plt.subplots(3,2,figsize=(10,7))


models=[m1,m2,m3]
names =["1 Box", "Multi Boxes","Chromatic box"]

for i,m in enumerate(models):
    
    visWest=np.abs(m.getComplexCoherentFlux(spf,spf0,wls)).reshape(nwl, nB)
    visWest /= np.outer(np.max(visWest, axis=1), np.ones(nB))
    visNorth=np.abs(m.getComplexCoherentFlux(spf0,spf,wls)).reshape(nwl, nB)
    visNorth /= np.outer(np.max(visNorth, axis=1), np.ones(nB))

    cb=ax[i,0].scatter(spf, visWest, c=wls*1e6, s=0.2, cmap="plasma")
    ax[i,1].scatter(spf, visNorth, c=wls*1e6, s=0.2, cmap="plasma")

    
    ax[i,0].set_ylabel("Vis. of {}".format(names[i]))
    
    if i!=2:
        ax[i,0].get_xaxis().set_visible(False)
        ax[i,1].get_xaxis().set_visible(False)
        
    ax[i,1].get_yaxis().set_visible(False)
        

fig.colorbar(cb, ax=ax.ravel().tolist(),label="$\\lambda$ ($\\mu$m)")

ax[2,0].set_xlabel("B/$\\lambda$ (cycles/rad)")
ax[2,1].set_xlabel("B/$\\lambda$ (cycles/rad)")
ax[0,0].set_title("East-West baselines")
ax[0,1].set_title("North-South baselines")

fig.savefig(os.path.join(path,os.pardir,"images","customCompMultiBoxesVis.png"))
