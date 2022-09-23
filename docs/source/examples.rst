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
    
Unlike in the previous example with the grey data, we create a 2D-array for the spatial frequencies of ``nB`` baselines by ``nwl`` wavelengths. The wavlength vector is tiled itself to have the same




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

