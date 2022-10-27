..  _expandingSoftware:
 
Expanding the Software
----------------------

In this section we present examples that show how to expand the functionalities of the oimodeler sofwate by creating customs objects : oimComponents, oimFilterComponents, oimFitters, and custom plotting function or utils.

Creating new Fourier Components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the `createCustomComponentFourier.py <https://github.com/oimodeler/oimodeler/blob/main/examples/ExpandingSoftware/createCustomComponentFourier.py>`_ example we show how to implement a new model component using a formula in the Fourier plan. The component will inherit from the  **oimComponentFourier** class. The Fourier formula should be implemented in  ``_visFunction`` and optionally the formula in the image plan can be implemented using  ``_imageFunction``. 


For this example we will show how to implement a basic rectangular box component. We start by importing oimodeler and some other useful packages.

.. code-block:: python

    import oimodeler as oim
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    import astropy.units as u


Our new component will be named **oimBox**, and it will have two parameters, ``dx`` and ``dy`` the size of the box in the x and y directions. Le'ts start to implement the oimBox class and its ``__init__`` method.


.. code-block:: python

    class box(oim.oimComponentFourier):
    name="2D Box"
    shortname = "BOX"
    def __init__(self,**kwargs): 
        
         super().__init__(**kwargs)
         self.params["dx"]=oim.oimParam(name="dx", value=1,description="Size in x",unit=u.mas)
         self.params["dy"]=oim.oimParam(name="dy", value=1,description="Size in y",unit=u.mas)       
         self._eval(**kwargs)
         

The class inherit from **oim.oimComponentFourier**. In the ``__init__`` method is called with the ``**kwargs`` arguments that allows to pass keyword arguments. To inherit from the parent class, we first call its  initialization method with ``super()__init__``. Then we define the two new parameters, dx and dy which are instances of the oimParam class. Finally we need to call the ``_eval`` method that allows the parameters to be processed.

Now that the new class is created, we need to implement the ``_visFunction`` method, with the Fourier transform formula of our component.  This method is called when using the getComplexCoherentFlux method of the oimComponent class. 

Note that the component parameters should be called with (wl,t) to allow parameter chromaticity and time dependence. The parameters have a unit and this should also be used to allow the use of other units when creating instances of the component.

In our case the complex visibilty of a rectangle is quite easy to write. It is a simple 2D-sinc function. Note that the x and y sizes are converted from the given unit (usually mas) to rad 

.. code-block:: python

    def _visFunction(self,ucoord,vcoord,rho,wl,t):
        
        x=self.params["dx"](wl,t)*self.params["dx"].unit.to(u.rad)*ucoord
        y=self.params["dy"](wl,t)*self.params["dy"].unit.to(u.rad)*vcoord
        
        return np.sinc(x)*np.sinc(y) 
    

We also need to implement the image method that will be called whenusing the getImage method. If not implemented the model will use the Fourier based formula to compute the image. It will also be the case if the keyword fromFT is set to True when the getImage is called. However it is always interesting to implement the image method, at least for debugging purpose, to check that the image compute for the image formula and using the fromFT option gives compatible results. We will check that later in that example.

For our box, we can implement the image method with logical operations

.. code-block:: python

    def _imageFunction(self,xx,yy,wl,t):
            
            return ((np.abs(xx)<=self.params["dx"](wl,t)/2) &
                    (np.abs(yy)<=self.params["dy"](wl,t)/2)).astype(float)


The full code of the oimBox component is quite short.

.. code-block:: python

    class oimBox(oim.oimComponentFourier):
    name="2D Box"
    shortname = "BOX"
    
    def __init__(self,**kwargs):       
         super().__init__(**kwargs)
         self.params["dx"]=oim.oimParam(name="dx", value=1,description="Size in x",unit=u.mas)
         self.params["dy"]=oim.oimParam(name="dy", value=1,description="Size in y",unit=u.mas)       
         self._eval(**kwargs)

    def _visFunction(self,ucoord,vcoord,rho,wl,t): 
        x=self.params["dx"](wl,t)*self.params["dx"].unit.to(u.rad)*ucoord
        y=self.params["dy"](wl,t)*self.params["dy"].unit.to(u.rad)*vcoord      
        return np.sinc(x)*np.sinc(y) 

    def _imageFunction(self,xx,yy,wl,t):            
            return ((np.abs(xx)<=self.params["dx"](wl,t)/2) &
                    (np.abs(yy)<=self.params["dy"](wl,t)/2)).astype(float)


We can now use it as any other oimodeler components. Let's build our first model with it.

.. code-block:: python
    
    b1=oimBox(dx=40,dy=10)
    m1=oim.oimModel([b1])
    
  
Now we can create images of our model: 

- with the _imageFunction
- with the FFT of the _visFunction

Both can be created with the ``showModel`` method of the oimComponent. To create the image from the FFT of the visibilty function, we just need to set the ``fromFT`` keyword to True.

.. code-block:: python

    fig, ax = plt.subplots(1,2,figsize=(10,5))
    m1.showModel(512,0.2,axe=ax[0],colorbar=False)
    m1.showModel(512,0.2,axe=ax[1],fromFT=True,colorbar=False)
    ax[0].set_title("Image with _imageFunction")
    ax[1].set_title("Image with FFT of _visFunction")


.. image:: ../../images/customCompBox1Image.png
  :alt: Alternative text   

Of course as our oimBox inherit from the oimComponent, it has three addtionnal parameters : the positions ``x`` and ``y`` and the flux ``f``. All oimComponent can also be rotated using the ``pa`` parameter. Note that if not set at the component creation the ``pa`` parameters (and the ``elong`` one) are not added to the model.

Let's create a complex model with boxes and uniform disk.

.. code-block:: python

    b2=oimBox(dx=2,dy=2,x=20,y=0,f=0.5)
    b3=oimBox(dx=10,dy=20,x=-30,y=10,pa=50,f=10)
    c=oim.oimUD(d=10,x=-30,y=-10)
    m2=oim.oimModel([b1,b2,b3,c])
    m2.showModel(512,0.2,colorbar=False,figsize=(5,5))


.. image:: ../../images/customCompBoxesImage.png
  :alt: Alternative text  
  
We could also create a chromatic box component using the oimInterpWl class or link parameters with 

.. code-block:: python

    b4=oimBox(dx=oim.oimInterpWl([2e-6,2.4e-6],[5,10]),dy=2,x=20,y=0,f=0.5)
    b4.params['dy']=oim.oimParamLinker(b4.params['dx'],'mult',4)
    
    m3=oim.oimModel([b4])

    m3.showModel(512,0.2,wl=[2e-6,2.2e-6,2.4e-6],colorbar=False,swapAxes=True)

.. image:: ../../images/customCompChromBoxImages.png
  :alt: Alternative text   
    

Let's finish this example by plotting the visibility of such models for a set of East-West and North-South baselines and wavelengths in the K band.



.. code-block:: python
     
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

        ax[i,0].scatter(spf, visWest, c=wls*1e6, s=0.2, cmap="plasma")
        ax[i,1].scatter(spf, visNorth, c=wls*1e6, s=0.2, cmap="plasma")
        ax[i,0].scatter(spf, visWest, c=wls*1e6, s=0.2, cmap="plasma")
        ax[i,1].scatter(spf, visNorth, c=wls*1e6, s=0.2, cmap="plasma")
        
        ax[i,0].set_ylabel("Vis. of {}".format(names[i]))
        
        if i!=2:
            ax[i,0].get_xaxis().set_visible(False)
            ax[i,1].get_xaxis().set_visible(False)
            
        ax[i,1].get_yaxis().set_visible(False)
            
    ax[2,0].set_xlabel("B/$\\lambda$ (cycles/rad)")
    ax[2,1].set_xlabel("B/$\\lambda$ (cycles/rad)")
    ax[0,0].set_title("East-West baselines")
    ax[0,1].set_title("North-South baselines")
                  

.. image:: ../../images/customCompMultiBoxesVis.png
  :alt: Alternative text   
    
Of course, only the third model is chromatic.

Creating new Image Components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::
    Example will be added when te oimComponentImage will be fully implemented

Creating new Radial profile Components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. warning::
    Example will be added when te oimComponentRadialProfile will be fully implemented