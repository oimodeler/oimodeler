..  _examples:

Examples
========

All the following examples can be found in the examples subdirectories of the oimodeler github repository.
Before looking at these examples, you might want to check the :ref:`getting_started` page


Basic Examples
--------------

In this section we presents script we presents showing the basic functionalities of the oimodeler software.

..  _exampleOimData:

Loading oifits data
^^^^^^^^^^^^^^^^^^^

The `exampleOimData.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimData.py>`_ script show how to create a oimData object from a list of oifits files and how the data in organized in the oimData instance.


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

The `basicModel.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/basicModel.py>`_ script demonstrate the basic functionalities of the oimModel and oimComponents object.


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

In the example `complexModel.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/complexModels.py>`_ we create and play with more complex Fourier-based models with includes:

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

    el.params['d'].value = 4
    el.params['pa'].value = 45
        
    m4.showModel(256, 0.1, wl=[3e-6, 3.25e-6, 3.5e-6, 4e-6], figsize=(14, 2.5), normPow=0.2)    
      
.. image:: ../../images/complexModel_linkRotScale.png
  :alt: Alternative text  




.. _createSimulator:

Comparing data and model with the oimSimulator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the `exampleOimSimulator.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOiSimulator.py>`_ script, we use the oimSimulator class to compare some oifits data with a model. We will compute the reduced chi2 and plot the comparison between the data an the simulated data from the model.

Let's start by importing the needed modules and setting ``files`` to the list of the same oifits files as in the :ref:`exampleOimData` example. 

.. code-block:: python

    import oimodeler as oim
    import matplotlib.pyplot as plt
    import os
    
    path = os.path.dirname(oim.__file__)
    pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")
    files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

These oifits were simulated with ASPRO as a MATISSE observation of a partly resolved binary star. 

We set a model a binary star with one component resolved. It consists in two components : a uniform disk and a point source.

.. code-block:: python

    ud=oim.oimUD(d=3,f=1,x=10,y=20)
    pt=oim.oimPt(f=0.5)
    model=oim.oimModel([ud,pt])

We now create a oimSimulator with the oimModel and the data. The data can either be :

- an oimData instance previously created
- a list of previously opened astropy.io.fits.hdulist
- a list of filenames to the oifits files (list of string)

.. code-block:: python

    sim=oim.oimSimulator(data=files,model=model)
    
Before using the simulator we need to prepare the data using the `prepare` method. This call the `prepare` method of the created oimData instance within the oimSimulator instance. The function is used to create vectorized coordinates for the data (spatial frequencies in x and y and wavelengths) to be passed to the oimModel instance to compute the complex Coherent Flux (ccf) using the oimModel.getComplexCoherentFlux method, and some structures to go back from the ccf to the measured interferometric quantities contained in the oifits files: VIS2DATA, VISAMP, VISPHI, T3AMP, T3PHI, and FLUXDATA.

.. code-block:: python

    sim.data.prepareData()

Once the data is prepared we can call the compute method to compute the chi2 and the simulatedData.

.. code-block:: python

    sim.compute(computeChi2=True,computeSimulatedData=True)
    print("Chi2r = {}".format(sim.chi2r))

.. code-block:: python

    Chi2r = 5674.502111807307


Our model isn't fitting well the data. Let's plot the data model comparison for all interferometric quantities contained in the oifits files.

.. code-block:: python

    fig0,ax0= sim.plot(["VIS2DATA","VISAMP","VISPHI","T3AMP","T3PHI"])
  
  
.. image:: ../../images/ExampleOimSimulator_model0.png
  :alt: Alternative text  


You can try to fit the model to the data "by hand", or go to the next example where we use a oimFitter subclass to automatically find the good parameters.


Running a mcmc fit
^^^^^^^^^^^^^^^^^^

In the `exampleOimFitterEmcee.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimFitterEmcee>`_ script, we perform a complete emcee run to determine the values of the parameters of the same binary as in the :ref:`createSimulator` example.

We start by setting up the script with imports, data list and a binary model. We don't need to specify values for the biary parameters as they will be fitted.

.. code-block:: python

    import oimodeler as oim
    import os

    path = os.path.dirname(oim.__file__)

    pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")
    files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

    ud=oim.oimUD()
    pt=oim.oimPt()
    model=oim.oimModel([ud,pt])


Before starting the run we need to specify which parameters are free and what are there range. By dafault all parameters are free but the components coordinates x and y. For a binary we need to set them to free for one of the components. As we only deal with relative fluxes, we can set the flux of one of the component to be fixed to one.

.. code-block:: python

    ud.params['d'].set(min=0.01,max=20)
    ud.params['x'].set(min=-50,max=50,free=True)
    ud.params['y'].set(min=-50,max=50,free=True)
    ud.params['f'].set(min=0.,max=10.)
    pt.params['f'].free=False
    
    print(model.getFreeParameters())
    

.. code-block::

    {'c1_UD_x': oimParam at 0x23d940e4850 : x=0 ± 0 mas range=[-50,50] free=True , 
    'c1_UD_y': oimParam at 0x23d940e4970 : y=0 ± 0 mas range=[-50,50] free=True ,
    'c1_UD_f': oimParam at 0x23d940e4940 : f=0.5 ± 0  range=[0.0,10.0] free=True ,
    'c1_UD_d': oimParam at 0x23d940e4910 : d=3 ± 0 mas range=[0.01,20] free=True }

We have 4 free-parameters, the position (x,y) flux and diameters of the uniform disk component.

Now we can create a fitter with our model and our filenames list of oifits files. We use the emcee fitter that have only one parameter, the number of walkers that will explore the parameters space. If you are not confident with emcee, you should have a look at the documentation `here <https://emcee.readthedocs.io/en/stable/>`_

.. code-block:: python
    
    fit=oim.oimFitterEmcee(files,model,nwalkers=32)
    

We need to initialize the fitter using its prepare method. The an emcee run that mainly mean setting the initial values of the walkers. The default method is to set them to random values within the parameters space.

.. code-block:: python
    
    fit.prepare(init="random")
    print(fit.initialParams)
    
.. code-block::  
 
    >>[[-37.71319618 -49.22761731   9.3299391   15.51294277]
       [-12.92392301  17.49431852   7.76169304   9.23732472]
       [-31.62470824 -11.05986877   8.71817772   0.34509237]
       [-36.38546264  33.856871     0.81935324   9.04534926]
       [ 45.30227534 -38.50625408   4.89978551  14.93004   ]
       [-38.01416866  -6.24738348   5.26662714  13.16349304]
       [-21.34600438 -14.98116997   1.20948714   8.15527356]
       [-17.14913499  10.40965493   0.37541088  18.81733973]
       [ -9.61039318 -12.02424002   6.81771974  16.22898422]
       [ 49.07320952 -34.48933488   1.75258006  19.96859116]]
       
 
We can now run the fit. We choose to run 2000 as a start and show interactively the progress as a progress bar. The fit should take a minutes on a standard computer to compute 64000 models (``nwalkers`` x ``nsteps``).

.. code-block:: python

    fit.run(nsteps=2000,progress=True)
 
The oimFitterEmcee instance store the emcee sampler as a member variable oimFitterEmcee.sampler. you can, for example, acces the chain of walkers and the log of probability directly.  

.. code-block:: python

    sampler = fit.sampler
    chain   = fit.sampler.chain
    lnprob  = fit.sampler.lnprobability
    
We can manipulate yourself these data. But the oimFitterEmcee implements varoius methods to retrieve and plot the results of the mcmc run.

The walkers position as the function of the steps can be plotted using the walkersPlot method.

.. code-block:: python

    figWalkers,axeWalkers=fit.walkersPlot(cmap="plasma_r")


.. image:: ../../images/exampleOimFitterEmceeWalkers.png
  :alt: Alternative text  


After a few hundred steps most walkers converge to a position with a good reduced chi2. However, from that figure will clearly see that:

- not all walkers have converge after 2000 steps
- some walkers converge to a solution that gives significantly worst chi2

In optical interferometry there are often local minimas in the chi2 and it seems that some of our walkers are locked there. In our case, this minimum is due to the fact that object is close be symmetrical if not for the fact than one of the component is resolved. Neverless, the chi2 of the local minimum is about 20 times worst the one of the global minimum.

We can plot the famous corner plot with the 1D and 2D density distribution. oimodel use the `corner.py <https://corner.readthedocs.io/en/latest/>`_ library for that purpose. We will discard the 1000 first steps as most of the walkers have converge after that. By default, the corner plot remove also the data with a chi2 greater than 20 times those of the best model. This option can be changed using the keyword ``chi2limfact`` 

.. code-block:: python

    figCorner,axeCorner=fit.cornerPlot(discard=1000)
    
   
.. image:: ../../images/exampleOimFitterEmceeCorner.png
  :alt: Alternative text  
  
  
We now can get the result of our fit. The oimFitterEmcee fitter can either return the ``best``, the ``mean`` or the ``median`` model. It return uncertainties estimated from the density distribution (see emcee doc for more details. 

.. code-block:: python
    
    median,err_l,err_u,err=fit.getResults(mode='median',discard=1000)

To compute the median and mean model we have to remove, as in the corner plot, the walkers that didn't converge with the ``chi2limitfact`` keyword (default in 20) and remove the steps of the bruning phase with the ``discard`` option.

When asking for the results, the simulatedData with these value are also produced in the fitter internal simulator. We can plot again the data/model and compute the final reduced chi2:

.. code-block:: python 
    
    figSim,axSim=fit.simulator.plot(["VIS2DATA","VISAMP","VISPHI","T3AMP","T3PHI"])
    print("Chi2r = {}".format(fit.simulator.chi2r))
    
.. image:: ../../images/ExampleOimFitterEmcee_fittedData.png
  :alt: Alternative text 

Filtering data
^^^^^^^^^^^^^^

Filtering can be applied to the oimData using the oimDataFilter class. The oimDataFilter is basically a stack of filters derived from the oimDataFilterComponent abstract class. The example presented here comes from the `exampleOimDataFilter.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimDataFilter>`_ script.

As for other example we will start by importing oimodeler and other useful packages and create a list of oifits files.  

.. code-block:: python 
    
    import oimodeler as oim
    import matplotlib.pyplot as plt
    import os

    path = os.path.dirname(oim.__file__)
    pathData=os.path.join(path,os.pardir,"examples","testData","FSCMa_MATISSE")
    files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]

We create an oimData object which will contain the oifits data. 

.. code-block:: python 
    
    data=oim.oimData(files)

We now create a simple filter to cut data to a specific wavelength range in the ``oimWavelengthRangeFilter`` class. 

.. code-block:: python 
    
    f1=oim.oimWavelengthRangeFilter(targets="all",wlRange=[3.0e-6,4e-6])
    
The ``oimWavelengthRangeFilter`` has two parameters :

- ``targets`` : which is common to all filter components : it specify the targeted files within the data structure to which the filter apply. Possible values are : "all" for all files, a single file specify by its index, or a list of indexes. Here we specify that we want to apply our filter to all data files.

- ``wlRange`` : the wavelength range to cut as a two elements list (min wavelength and max wavelength), or a list of multiple two elements list if you want to cut multiple wavelengths ranges simultaneously. In our example you have selected wavelength between 3 and 4 microns. Wavelengths outside this range will be removed from the data.
    
Now we can create a filter stack with this single filter and apply it to our data.

.. code-block:: python 

    filters=oim.oimDataFilter([f1])
    data.setFilter(filters)
    

By default the filter will be automatically activated as soon as a filter is set using the ``setFilter`` method of the oimData class. This means that the call to oimData.data will return the filtered data, and that if using the oimData class within a oimSimulator or a oimFitter, the filtered data will be used instead of the unfiltered data. 

.. note::
    The unfiltered data can always be accessed using oimData._data and the filtered data, that may be None if no filter have been set, using oimData._filteredData
   
To switch off a filter we can either call the setFilter without parameters (this will remove the filter completely) or set the useFilter variable to False.

.. code-block:: python 

    #data.setFilters() #removing the filter
    data.useFilter = False
    
Let's plot the unfiltered and filtered data using the oimPlot method.

.. code-block:: python 

    fig=plt.figure()
    ax = plt.subplot(projection='oimAxes')

    data.useFilter = False
    ax.oiplot(data,"SPAFREQ","VIS2DATA",color="tab:blue",lw=3,alpha=0.2,label="unfiltered")

    data.useFilter = True
    ax.oiplot(data,"SPAFREQ","VIS2DATA",color="tab:blue",label="filtered")

    ax.set_yscale('log')
    ax.legend()
    ax.autolim()
    

.. image:: ../../images/ExampleFilter_wavelengthCut.png
  :alt: Alternative text 
  
The other simple filters for data selection are :

- ``oimRemoveArrayFilter`` : removing array (such as OI_VIS, OI_T3...) from the data. 
- ``oimDataTypeFilter`` : removing data type (such as VISAMP, VISPHI, T3AMP...) from the data.

.. note::
    Actually oimDataTypeFilter doesn't remove the columns with the data type from any array as these column are complusory in the the oifits format definition. Instead it is setting all the values of the column to zero which oimodeler will recognize as emplty for data simulation and model fitting. 

.. code-block:: python 

    f2=oim.oimRemoveArrayFilter(targets="all",arr=["OI_VIS","OI_FLUX"])         
    f3=oim.oimDataTypeFilter(targets="all",dataType=["T3AMP","T3PHI"])
    data.setFilter(oim.oimDataFilter([f1,f2,f3]))

Here we create a new filter stack with the previous wavelength filter (f1), a filter (f2) removing the array OI_VIS and OI_FLUX from the data, and a filter (f3) removing the columns T3AMP and T3PHI. Basically, we only have VIS2DATA left in our oifits structure.

.. note::
    Removing T3AMP and T3PHI from the OI_T3 is equivalent for model-fitting to remove the array OI_T3 for model-fitting. 


Plotting data from oifits files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Beyond the specific plots shown in the previous example the oimPlot module allow to plot most of the oifits data in a very simple way. The example presented here comes from the `exampleOimPlot.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimPlot.py>`_ script.

Let's start by setting up the project with imports, path, and some data.

.. code-block:: python 

    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from astropy.io import fits
    import oimodeler as oim

    path = os.path.dirname(oim.__file__)
    pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_MATISSE2")

    files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]
    data=[fits.open(fi,mode="update") for fi in files]
    
oimodeler comes with the oimAxes class that subclass the standard matplotlib.pytplotAxes class (base class for all matplotlib plots). To use it you simply need to specify it as a projection (actually it calls the subclass) when creating the axe or axes.

.. code-block:: python 

    fig, ax = plt.subplots(subplot_kw=dict(projection='oimAxes'))
   
First we can plot the classic uv coverage using the uvplot method by passing the oifits data.

.. code-block:: python 

    ax[0,0].uvplot(data)
    
.. image:: ../../images/ExampleOimPlot_uv.png
  :alt: Alternative text     
    
We can use the oiplot method of the oimAxes to plot any quantity inside an oifits file as a function of another one. For instance let's plot the squared visibilities as a function of the spatial frequencies with the wavelength as a colorscale

.. code-block:: python
   
    ax = plt.subplot(projection='oimAxes')
    lamcol=ax.oiplot(data,"SPAFREQ","VIS2DATA" ,xunit="cycles/mas",label="Data",
                    cname="EFF_WAVE",cunitmultiplier=1e6,errorbar=True)
                    
    plt.colorbar(lamcol, ax=ax,label="$\\lambda$ ($\mu$m)")
    ax.legend()
    
.. image:: ../../images/ExampleOimPlot_v2.png
  :alt: Alternative text     
  
  
We can also plot the square visibility as the function of the wavelength.

.. code-block:: python

    ax.oiplot(data,"EFF_WAVE","VIS2DATA",xunitmultiplier=1e6,
               errorbar=True,kwargs_error={"alpha":0.3})
  
.. image:: ../../images/ExampleOimPlot_v2Wl.png
  :alt: Alternative text       
  
Finally, we can create a 2x2 figure with multiple plots. The projection keyword have to be set for all Axes using the subplot_kw keyword in the subplots method.

.. code-block:: python

    fig, ax = plt.subplots(2,2, subplot_kw=dict(projection='oimAxes'),figsize=(8,8))
   
    ax[0,0].uvplot(data)

    lamcol=ax[0,1].oiplot(data,"SPAFREQ","VIS2DATA" ,xunit="cycles/mas",label="Data",
                        cname="EFF_WAVE",cunitmultiplier=1e6,ls=":",errorbar=True)
    fig.colorbar(lamcol, ax=ax[0,1],label="$\\lambda$ ($\mu$m)")
    ax[0,1].legend()
    ax[0,1].set_yscale('log')   

    ax[1,0].oiplot(data,"EFF_WAVE","VIS2DATA",xunitmultiplier=1e6,
                   errorbar=True,kwargs_error={"alpha":0.3})
    ax[1,0].autolim()

    ax[1,1].oiplot(data,"SPAFREQ","T3PHI",xunit="cycles/mas",errorbar=True,
                   lw=2,ls=":")
    ax[1,1].autolim()
    
.. image:: ../../images/ExampleOimPlot_multi.png
  :alt: Alternative text   
    

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
    m2.showModel(512,0.2,colorbar=False)


.. image:: ../../images/customCompBoxesImage.png
  :alt: Alternative text  
  
We could also create a chromatic box component using the oimInterpWl class or link parameters with 

.. code-block:: python

    b4=oimBox(dx=oim.oimInterpWl([2e-6,2.4e-6],[5,10]),dy=2,x=20,y=0,f=0.5)
    b4.params['dy']=oim.oimParamLinker(b4.params['dx'],'mult',4)
    
    m3=oim.oimModel([b4])

    m3.showModel(512,0.2,wl=[2e-6,2.2e-6,2.4e-6],colorbar=False)

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
    Example will be added when te oimComponentImage will be fully implemented

Performance Tests
-----------------



Scripts concerning performance tests are presented in this section.

Data for tests
--------------

