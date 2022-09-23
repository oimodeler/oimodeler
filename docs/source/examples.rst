..  _examples:

Examples
========

All the following examples can be found in the examples subdirectories of the oimodeler github repository.
Before looking at these examples, you might want to check the :ref:`getting_started` page


Basic Examples
--------------

In this section we presents script we presents showing the basic functionalities of the oimodeler software.


Loading oifits data
^^^^^^^^^^^^^^^^^^^

The ``exampleOimData.py`` script show how to create a oimData object from a list of oifits files and how the data in organized in the oimData instance.


.. code-block:: python

    import oimodeler as oim
    import os

    path = os.path.dirname(oim.__file__)
    pathData=os.path.join(path,os.pardir,"examples","testData","FSCMa_MATISSE")
    files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData)]

    data=oim.oimData(files)


The oifits data, stored in the ``astropy.io.fits.hdulist`` format, can be accessed using the ``oimData.data`` variable

.. code-block:: python

    print(data.data)
    
.. code-block:: python

    >>[[<astropy.io.fits.hdu.image.PrimaryHDU object at 0x000002657CBD7CA0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E546AF0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E3EA970>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E3EAAC0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E406520>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E402EE0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E406FD0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E4600D0>],
    [<astropy.io.fits.hdu.image.PrimaryHDU object at 0x000002657E458F70>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x0000026500769BE0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002650080EA60>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x00000265007EA430>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x00000265007EAAF0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002650080EC40>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E4DC820>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E4ECFD0>],
    [<astropy.io.fits.hdu.image.PrimaryHDU object at 0x000002657E4DCCA0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x0000026500B7EB50>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E9F79D0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E5913A0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E591A60>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E591B20>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E5B7790>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E5BAEB0>]]
    
    
To learn more on the astropy.i OIFITS2

To be used in the oimSimulator and oiFitter data need to be optimized in a simpler vectorial/structure. Tis step is done automatically when using the simulator or fitter but can be done manually using the following command:
    
.. code-block:: python
    
    data.prepareData()
    
For instance this create single vectors fgor the data coordinate : ``data.vect_u``, ``data.vect_v``, ``data.vect_wl``

.. code-block:: python

    print(data.vect_u)
    print(data.vect_v)   
    print(data.vect_wl)  
    print(data.vect_u.shape)
    
.. code-block:: python
    
    [0. 0. 0. ... 0. 0. 0.]
    [0. 0. 0. ... 0. 0. 0.]
    [4.20059359e-06 4.18150239e-06 4.16233070e-06 ... 2.75303296e-06
     2.72063039e-06 2.68776785e-06]
    (5376,)
    
    
Basic models
^^^^^^^^^^^^

The ``basicModel.py`` script demonstrate the basic functionalities of the oimModel and oimComponents object.


First we import the oimodeler and numpy.


.. code-block:: python

    import oimodeler as oim
    import numpy as np
    
    
A model is a collection of components. All components derived from the oimComponent class. The components may be described in the image plan by their intensity distribution or directly in the Fourier plan for components with known analytical Fourier transforms. In these example we will only focus on this later type of component which all derived from the oimFourierComponent class. In the table below is a list of the currently implemented oimFourierComponents:

+---------------+--------------------------------+------------------------------+
| class         | description                    | parameters                   |
+===============+================================+==============================+
| oimPt         | Point source                   | x,y,f                        |
+---------------+--------------------------------+------------------------------+
| oimBackground | Background                     | x,y,f                        |
+---------------+--------------------------------+------------------------------+
| oimUD         | Uniform Disk                   | x,y,f,d                      |
+---------------+--------------------------------+------------------------------+
| oimEllipse    | Uniform Ellipse                | x,y,f,d,pa,elong             |
+---------------+--------------------------------+------------------------------+
| oimGauss      | Gaussian Disk                  | x,y,f,fwhm                   |
+---------------+--------------------------------+------------------------------+
| oimEGauss     | Point source                   | x,y,f,fwhm,pa,elong          |
+---------------+--------------------------------+------------------------------+
| oimIRing      | Infinitesimal Ring             | x,y,f,d                      |
+---------------+--------------------------------+------------------------------+
| oimEIRing     | Ellitical infinitesimal ring   | x,y,f,d,pa,elong             |
+---------------+--------------------------------+------------------------------+
| oimRing       | Ring                           | x,y,f,din,dout               |
+---------------+--------------------------------+------------------------------+
| oimERing      | Ellitical  ring                | x,y,f,din,dout,pa,elong      |
+---------------+--------------------------------+------------------------------+
| ESKRing       | Skewed Ellitical ring          | x,y,f,d,skw,skwPa,pa,elong   |
+---------------+--------------------------------+------------------------------+


To create models we must first create the components. Let's create a few simple components.


.. code-block:: python

    pt = oim.oimPt(f=0.1)
    ud = oim.oimUD(d=10,f=0.5)
    g  = oim.oimGauss(fwhm=5,f=1)
    r  = oim.oimIRing(d=5,f=0.5)

    
Here we have create a point source components, a 10 mas uniform disk, a Gaussian distribution with a 5 mas fwhm and a 5 mas infinitesimal ring. 

Note that the model parameters which are not set explicitly during the components creation are set to their default values (i.e., f=1 x=y=0).

We can print the description of the component easily


.. code-block:: python

    print(ud)

.. code-block::
    
    >>Uniform Disk x=0.00 y=0.00 f=0.50 d=10.00

Or you want to print the details of a parameter:

.. code-block:: python

    print(ud.params['d'])
 
.. code-block:: 
    
    >>oimParam d = 10 ± 0 mas range=[-inf,inf] free

Note that the components parameters are instances of the oimParam class which hold not only the parameter value stored in oimParam.value but also : 

- oimParam.error : the parameters uncertainties (for model fitting)
- oimParam.unit : the unit as a astropy.unit object
- oimParam.min : minimum possible value (for model fitting)
- oimParam.max : minimum possible value (for model fitting)
- oimParam.free : True=free parameter and False=fixed parameter (for model fitting)
- oimParam.description : A string that describes the model parameter

We can now create our first models uinsg the oimModel class.


.. code-block:: python

    mPt   = oim.oimModel([pt])
    mUD   = oim.oimModel([ud])
    mG    = oim.oimModel([g])
    mR    = oim.oimModel([r])
    mUDPt = oim.oimModel([ud,pt])
    
    

we now have 4 one-component models and 1 2-components models.

We can get the parameters of our models using the getParameter method of the oimModel class. 

.. code-block:: python
    
    params=mUDPt.getParameters()
    print(params)
        

.. code-block::

    {'c1_UD_x': oimParam at 0x23de5c62fa0 : x=0 ± 0 mas range=[-inf,inf] free=False ,
    'c1_UD_y': oimParam at 0x23de5c62580 : y=0 ± 0 mas range=[-inf,inf] free=False , 
    'c1_UD_f': oimParam at 0x23de5c62400 : f=0.5 ± 0  range=[-inf,inf] free=True ,
    'c1_UD_d': oimParam at 0x23debc1abb0 : d=10 ± 0 mas range=[-inf,inf] free=True , 
    'c2_Pt_x': oimParam at 0x23debc1a8b0 : x=0 ± 0 mas range=[-inf,inf] free=False , 
    'c2_Pt_y': oimParam at 0x23debc1ab80 : y=0 ± 0 mas range=[-inf,inf] free=False , 
    'c2_Pt_f': oimParam at 0x23debc1ac10 : f=0.1 ± 0  range=[-inf,inf] free=True }

getParameters returns a dict of all parameters of the components of the model. The keys are defined as x{num of component}_{short Name of component}_{param name}.

Alternatively we can get the free parameters using the getFreeParameters method:

.. code-block:: python
    
    freeParams=mUDPt.getParameters()
    print(freeParams)
        
.. code-block::

    {'c1_UD_f': oimParam at 0x23de5c62400 : f=0.5 ± 0  range=[-inf,inf] free=True ,
    'c1_UD_d': oimParam at 0x23debc1abb0 : d=10 ± 0 mas range=[-inf,inf] free=True ,
    'c2_Pt_f': oimParam at 0x23debc1ac10 : f=0.1 ± 0  range=[-inf,inf] free=True }


The oiModel can return an image of the model using the getImage method. It takes two arguments, the image size in pixels and the pixel size in mas.

.. code-block:: python
    
    im=mUDPt.getImage(512,1)
    plt.imshow(im**0.2)

.. image:: ../../images/basicModel_imshow.png
  :alt: Alternative text   
  

We plot the image with a 0.2 power-law to make the uniform disk components visible: both components have the same total flux but the UD is spread on much more pixels.

Alternatively we can use the method showModel which take the same argument as the getImage, but directly create a plot with proper axes and colorbar.

.. code-block:: python

    figImg,axImg=mUDPt.showModel(512,0.2,normPow=0.1


.. image:: ../../images/basicModel_showModel.png
  :alt: Alternative text  


In other examples, we use  oimModel and oimData objects within a oimSimulator to simulate interferometric quantities from the model at the spatial frequencies from the data.  Without the oimSimulator the oimModel can only produce complex coherent flux (i.e. non normalized complex visibility) for a vector of spatial frequecies and wavelengths. 

.. code-block:: python

    wl=2.1e-6
    B=np.linspace(0.0,300,num=200)
    spf=B/wl

Here we have create a vector of 200 spatial frequencies for baselines ranging from 0 to 300 m  and for an observing wavelength of 2.1 microns.

    We can now use this vector to get the complex coherent flux (CCF) from our model. 
    

.. code-block:: python

    ccf = mUDPt.getComplexCoherentFlux(spf,spf*0) 
    
The getComplexCoherentFlux take three parameters : the spatial frequencies along the east-west axis, the spatial frequencies along the North-South axis, and optionally, the wavelength. Here we are dealing with grey models so we don't need to specify the wavelength. And, as our models are circular, we don't care about the baseline orientation and a set the North-South component of the spatial frequencies to zero.


We can now plot the visibility from the CCF as the function of the spatial frequencies:

.. code-block:: python

    
    v = np.abs(ccf)
    v=v/v.max()
    plt.plot(spf , v)
    plt.xlabel("spatial frequency (cycles/rad)")
    plt.ylabel("Visbility")

.. image:: ../../images/basicModel_vis0.png
  :alt: Alternative text  


Let's finish this example by creating a figure with the image and visibility for all the previously created models.

.. code-block:: python

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
        

.. image:: ../../images/basicModel_all.png
  :alt: Alternative text 

.. _createModelChromatic:

Complex models
^^^^^^^^^^^^^^

In this examples we create and play with more complex Fourier-based models with includes:

- flatenning of some components
- linked parameters between components
- Chromaticity of some parameters

First we import the useful packages and create a set of spatial frequencies and wavelengths to be used to generate visibilities.

.. code-block:: python

    import oimodeler as oim
    import numpy as np
    import matplotlib.pyplot as plt
    
    
    
    nB=100 #number of baselines 
    nwl=100 #number of walvengths

    #Create some spatial frequencies
    wl=np.linspace(3e-6,4e-6,num=nwl)
    B=np.linspace(0,150,num=nB)
    Bs=np.tile(B,(nwl,1)).flatten()
    
    wls=np.transpose(np.tile(wl,(nB,1))).flatten()
    spf=Bs/wls
    
Unlike in the previous example with the grey data, we create a 2D-array for the spatial frequencies of ``nB`` baselines by ``nwl`` wavelengths. The wavlength vector is tiled itself to have the same length as the spatial frequency vector.

Let's create our first chromatic components. Chromaticity can added to grey Fourier-based model by using the oimInterpWl when creating a new component.

.. code-block:: python

    g=oim.oimGauss(fwhm=oim.oimInterpWl([3e-6,4e-6],[2,8]))
    
We have created a Gaussian component with a fwhm growing from 2 mas at 3 microns to 8 mas at 4 microns.
We can access to the interpolated value of the parameters using the call operator ().


.. code-block:: python

    print(g.params['fwhm']([3e-6,3.5e-6,4e-6,4.5e-6]))

.. code-block:: python
    
    >>[2. 5. 8. 8.]
    
The values are interpolated within the wavelength range [3e-6,4e-6] and fixed beyond these range.

Let's build a simple model with this component and plot the images at few wavelengths and the visibilities for the baselines we created before.

.. code-block:: python

    vis=np.abs(mg.getComplexCoherentFlux(spf,spf*0,wls)).reshape(len(wl),len(B))
    vis/=np.outer(np.max(vis,axis=1),np.ones(nB))

    figGv,axGv=plt.subplots(1,1,figsize=(14,8))
    sc=axGv.scatter(spf,vis,c=wls*1e6,s=0.2,cmap="plasma")
    figGv.colorbar(sc, ax=axGv,label="$\\lambda$ ($\\mu$m)")
    axGv.set_xlabel("B/$\\lambda$ (cycles/rad)")
    axGv.set_ylabel("Visiblity")
    axGv.margins(0,0)
    

.. image:: ../../images/complexModel_chromaticGaussian.png
  :alt: Alternative text 

.. image:: ../../images/complexModel_chromaticGaussianVis.png
  :alt: Alternative text 

Now let's add a second component: a uniform disk with a chromatic flux.

.. code-block:: python
    
    ud=oim.oimUD(d=0.5,f=oim.oimInterpWl([3e-6,4e-6],[2,0.2]))
    m2=oim.oimModel([ud,g])

    fig2im,ax2im = m2.showModel(256,0.1,wl=[3e-6,3.25e-6,3.5e-6,4e-6],figsize=(14,2.5))
    vis=np.abs(m2.getComplexCoherentFlux(spf,spf*0,wls)).reshape(len(wl),len(B))
    vis/=np.outer(np.max(vis,axis=1),np.ones(nB))

    fig2v,ax2v=plt.subplots(1,1,figsize=(14,8))
    sc=ax2v.scatter(spf,vis,c=wls*1e6,s=0.2,cmap="plasma")
    fig2v.colorbar(sc, ax=ax2v,label="$\\lambda$ ($\\mu$m)")
    ax2v.set_xlabel("B/$\\lambda$ (cycles/rad)")
    ax2v.set_ylabel("Visiblity")
    ax2v.margins(0,0)
    ax2v.set_ylim(0,1)


.. image:: ../../images/complexModel_UDAndGauss.png
  :alt: Alternative text 

.. image:: ../../images/complexModel_UDAndGaussVis.png
  :alt: Alternative text 
    


Now let's create a similar model but with elongated components. We will replace the uniform disk by an ellipse and the Gaussian by an elongated Gaussian.

.. code-block:: python

    eg=oim.oimEGauss(fwhm=oim.oimInterpWl([3e-6,4e-6],[2,8]),elong=2,pa=90)
    el=oim.oimEllipse(d=0.5,f=oim.oimInterpWl([3e-6,4e-6],[2,0.1]),elong=2, pa=90)

    m3=oim.oimModel([el,eg])
    fig3im,ax3im = m3.showModel(256,0.1,wl=[3e-6,3.25e-6,3.5e-6,4e-6],figsize=(14,2.5),normPow=0.2)

.. image:: ../../images/complexModel_Elong.png
  :alt: Alternative text

Now that our model is no more circular, we need to take care of the baselines orientations. Let's plot both North-South and East-West baselines.

.. code-block:: python

    fig3v,ax3v=plt.subplots(1,2,figsize=(14,5),sharex=True,sharey=True)

    # East-West
    vis = np.abs(m3.getComplexCoherentFlux(spf, spf*0, wls)).reshape(len(wl), len(B))
    vis /= np.outer(np.max(vis, axis=1), np.ones(nB))
    ax3v[0].scatter(spf, vis, c=wls*1e6, s=0.2, cmap="plasma")
    ax3v[0].set_title("East-West Baselines")
    ax3v[0].margins(0, 0)
    ax3v[0].set_ylim(0, 1)
    ax3v[0].set_xlabel("B/$\\lambda$ (cycles/rad)")
    ax3v[0].set_ylabel("Visiblity")

    # North-South
    vis = np.abs(m3.getComplexCoherentFlux(spf*0, spf, wls)).reshape(len(wl), len(B))
    vis /= np.outer(np.max(vis, axis=1), np.ones(nB))
    sc = ax3v[1].scatter(spf, vis, c=wls*1e6, s=0.2, cmap="plasma")
    ax3v[1].set_title("North-South Baselines")
    ax3v[1].set_xlabel("B/$\\lambda$ (cycles/rad)")
    fig3v.colorbar(sc, ax=ax3v.ravel().tolist(), label="$\\lambda$ ($\\mu$m)")
    
.. image:: ../../images/complexModel_ElongVis.png
  :alt: Alternative text
  
  
Let's have a look at our last model free parameters.

.. code-block:: python

    print(m3.getFreeParameters())
    
   
.. code-block::   
  
    >>{'c1_eUD_f_interp1': oimParam at 0x23d9e7194f0 : f=2 ± 0  range=[-inf,inf] free=True ,
    'c1_eUD_f_interp2': oimParam at 0x23d9e719520 : f=0.2 ± 0  range=[-inf,inf] free=True ,
    'c1_eUD_elong': oimParam at 0x23d9e7192e0 : elong=2 ± 0  range=[-inf,inf] free=True ,
    'c1_eUD_pa': oimParam at 0x23d9e719490 : pa=90 ± 0 deg range=[-inf,inf] free=True ,
    'c1_eUD_d': oimParam at 0x23d9e7193a0 : d=0.5 ± 0 mas range=[-inf,inf] free=True ,
    'c2_EG_f': oimParam at 0x23d9e7191c0 : f=1 ± 0  range=[-inf,inf] free=True ,
    'c2_EG_elong': oimParam at 0x23d9e7191f0 : elong=2 ± 0  range=[-inf,inf] free=True ,
    'c2_EG_pa': oimParam at 0x23d9e719220 : pa=90 ± 0 deg range=[-inf,inf] free=True ,
    'c2_EG_fwhm_interp1': oimParam at 0x23d9e7192b0 : fwhm=2 ± 0 mas range=[-inf,inf] free=True ,
    'c2_EG_fwhm_interp2': oimParam at 0x23d9e719340 : fwhm=8 ± 0 mas range=[-inf,inf] free=True }
  
We see here that for the Ellipse (C1_eUD) the f parameter has been replaced by two independent parameters called ``c1_eUD_f_interp1`` and ``c1_eUD_f_interp2``. They represent the value of the flux at 3 and 4 microns. We could have added more reference wavelengths in our model and would have ended with more parameters. The same happens for the elongated Gaussian (C2_EG) fwhm.

Currently our model has 10 free parameters. In certain cases we might want to link or share two or more parameters. In our case, we might consider that the two components have the same ``pa`` and ``elong``. The can be done easily. To share a parameter you can just replace one parameter by another.

.. code-block:: python
   
    eg.params['elong']=el.params['elong']
    eg.params['pa']=el.params['pa']
    
    print(m3.getFreeParameters())
    
.. code-block::  

    {'c1_eUD_f_interp1': oimParam at 0x23d9e7194f0 : f=2 ± 0  range=[-inf,inf] free=True ,
    'c1_eUD_f_interp2': oimParam at 0x23d9e719520 : f=0.2 ± 0  range=[-inf,inf] free=True ,
    'c1_eUD_elong': oimParam at 0x23d9e7192e0 : elong=2 ± 0  range=[-inf,inf] free=True ,
    'c1_eUD_pa': oimParam at 0x23d9e719490 : pa=90 ± 0 deg range=[-inf,inf] free=True ,
    'c1_eUD_d': oimParam at 0x23d9e7193a0 : d=0.5 ± 0 mas range=[-inf,inf] free=True ,
    'c2_EG_f': oimParam at 0x23d9e7191c0 : f=1 ± 0  range=[-inf,inf] free=True ,
    'c2_EG_fwhm_interp1': oimParam at 0x23d9e7192b0 : fwhm=2 ± 0 mas range=[-inf,inf] free=True ,
    'c2_EG_fwhm_interp2': oimParam at 0x23d9e719340 : fwhm=8 ± 0 mas range=[-inf,inf] free=True }
    
    
That way we have reduced our number of free parameters to 8. If you change the eg.params['elong'] or el.params['elong'] values it will change both parameters are they are actually the same instance of the oimParam class.

Let's create a new model which include a elongated ring perpendicular to the Gaussian and Ellipse pa and with a inner and outer radii equals to 2 and 4 times the ellipse diameter, respectively.

.. code-block:: python

    er = oim.oimERing()

    er.params['elong']=eg.params['elong']
    er.params['pa']=oim.oimParamLinker(eg.params["pa"],"add",90)
    er.params['din']=oim.oimParamLinker(el.params["d"],"mult",2)
    er.params['dout']=oim.oimParamLinker(el.params["d"],"mult",4)

    m4= oim.oimModel([el, eg,er])

    m4.showModel(256, 0.1, wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6], figsize=(14, 2.5), normPow=0.2)
       
 
.. image:: ../../images/complexModel_link.png
  :alt: Alternative text 
    
Although quite complex this models only have 9 free parameters. If we change the ellipse diameter and its position angle, the components will scale (except the Gaussian that fwhm is independent) and rotate.

.. code-block:: python

    el.params['d'].value=4
    el.params['pa'].value=45
        
    m4.showModel(256, 0.1, wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6], figsize=(14, 2.5), normPow=0.2)    
    
    
.. image:: ../../images/complexModel_linkRotScale.png
  :alt: Alternative text  


    
.. _createSimulator:

Create a Simulator
^^^^^^^^^^^^^^^^^^




Plotting data from oifits files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Expanding the Software
----------------------

In this section we present examples that show how to expand the functionalities of the oimodeler sofwate by crezating customs objects : oimComponents, oimFilterComponents, oimFitters, and custom plotting function or utils.

Performance Tests
-----------------

Scripts concerning performance tests are presented in this section.

Data for tests
--------------

