:tocdepth: 2

..  _models:

Building and using models
=========================

.. _basics of models:

The basics of models
--------------------

This complete code corresponding to this section is available in `TheBasicsOfModels.py <https://github.com/oimodeler/oimodeler/blob/main/examples/Modules/TheBasicsOfModels.py>`_ 

Models, Components and Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the **oimodeler** framework, a model is and instance of the :func:`oimModel <oimodeler.oimModel.oimModel>` class. 
It contains a collection of components, which all derived from the :func:`oimComponent <oimodeler.oimComponent.oimComponent>` 
semi-abstract class. The components may be described in the image plane, by their 1D or 2D intensity distribution,
or directly in the Fourier plane, for the most simple components with known analytical Fourier transforms. 
Each components is described by a set of parameters which are instances of the :func:`oimParam <oimodeler.oimParam.oimParam>` class.

Thus, building models in **oimodeler** relies on three classes:

- :func:`oimModel <oimodeler.oimModel.oimModel>`: the model class 

- :func:`oimComponent <oimodeler.oimComponent.oimComponent>`: the abstract class from which all components derive

- :func:`oimParam <oimodeler.oimParam.oimParam>`: the parameter class


To create models we must first create some components.
Let's create a few simple components.

.. code-block:: ipython3

    pt = oim.oimPt(f=0.1)
    ud = oim.oimUD(d=10, f=0.5)
    g  = oim.oimGauss(fwhm=5, f=1)
    r  = oim.oimIRing(d=5, f=0.5)

Here, we have create a point source, a 10 mas uniform disk, a Gaussian distribution 
with a 5 mas fwhm and a 5 mas infinitesimal ring. 
The comprehensive list of components available is **oimodeler** is given in the next section. 

The model parameters which are not set explicitly during the components creation
are set to their default values (i.e., f=1, x=y=0).

We can print the description of the component easily:

.. code:: ipython3

    print(ud.params)

.. parsed-literal::
    
    Uniform Disk x=0.00 y=0.00 f=0.50 d=10.00

Or if you want to print the details of a parameter:

.. code-block:: ipython3

    print(ud.params['d'])

 
.. parsed-literal::
    
    oimParam d = 10 ± 0 mas range=[-inf,inf] free


Note that the components parameters are instances of the
:func:`oimParam <oimodeler.oimParam.oimParam>` class which hold not only the
parameter value stored in the ``oimParam.value`` attribute, but in addition to it
the following attributes: 

- ``oimParam.error``: the parameters uncertainties (for model fitting).
- ``oimParam.unit``: the unit as a ``astropy.units`` object.
- ``oimParam.min``: minimum possible value (for model fitting).
- ``oimParam.max``: minimum possible value (for model fitting).
- ``oimParam.free``: Describes a free parameter for ``True``
  and a fixed parameter for ``False`` (for model fitting).
- ``oimParam.description``: A string that describes the model parameter.


Building Models
~~~~~~~~~~~~~~~

We can now create our first models using the
:func:`oimModel <oimodeler.oimModel.oimModel>` class.

.. code-block:: ipython3

    mPt   = oim.oimModel(pt)
    mUD   = oim.oimModel(ud)
    mG    = oim.oimModel(g)
    mR    = oim.oimModel(r)
    mUDPt = oim.oimModel(ud, pt)
    

Now, we have four one-component models and one two-component model.

We can get the parameters of our models using the 
:func:`oimModel.getParameters <oimodeler.oimModel.oimModel.getParameters>`
method.


.. code-block:: ipython3
    
    params = mUDPt.getParameters()
    print(params)
        

.. parsed-literal::

    {'c1_UD_x': oimParam at 0x23de5c62fa0 : x=0 ± 0 mas range=[-inf,inf] free=False ,
     'c1_UD_y': oimParam at 0x23de5c62580 : y=0 ± 0 mas range=[-inf,inf] free=False , 
     'c1_UD_f': oimParam at 0x23de5c62400 : f=0.5 ± 0  range=[-inf,inf] free=True ,
     'c1_UD_d': oimParam at 0x23debc1abb0 : d=10 ± 0 mas range=[-inf,inf] free=True , 
     'c2_Pt_x': oimParam at 0x23debc1a8b0 : x=0 ± 0 mas range=[-inf,inf] free=False , 
     'c2_Pt_y': oimParam at 0x23debc1ab80 : y=0 ± 0 mas range=[-inf,inf] free=False , 
     'c2_Pt_f': oimParam at 0x23debc1ac10 : f=0.1 ± 0  range=[-inf,inf] free=True }

The method returns a dict of all parameters of the model components.
The keys are defined as 

    ``x{num of component}_{short Name of component}_{param name}``.

Alternatively, we can get the free parameters using the
:func:`getFreeParameters <oimodeler.oimModel.oimModel.getFreeParameters>` method:

.. code-block:: ipython3
    
    freeParams = mUDPt.getParameters()
    print(freeParams)
    
.. parsed-literal::

    {'c1_UD_f': oimParam at 0x23de5c62400 : f=0.5 ± 0  range=[-inf,inf] free=True ,
     'c1_UD_d': oimParam at 0x23debc1abb0 : d=10 ± 0 mas range=[-inf,inf] free=True ,
     'c2_Pt_f': oimParam at 0x23debc1ac10 : f=0.1 ± 0  range=[-inf,inf] free=True }

The two main methods of an :func:`oimModel <oimodeler.oimModel.oimModel>` object are:

- :func:`getImage <oimodeler.oimModel.oimModel.getImage>`: which returns an image of the model 
- :func:`oimModel.getComplexCoherentFlux <oimodeler.oimModel.oimModel.getComplexCoherentFlux>` which returns the complex Coherent Flux of the model 

Althought the :func:`getImage <oimodeler.oimModel.oimModel.getImage>`  is only used to vizualize the model intensity 
distribution and is not used for  model-fitting, :func:`getComplexCoherentFlux <oimodeler.oimModel.oimModel.getComplexCoherentFlux>` is
 at the base of the computation of all interferometric observables and thus of the data-model comparison.


Getting the model image
~~~~~~~~~~~~~~~~~~~~~~~

Let's first have a look at the :func:`oimModel.getImage <oimodeler.oimModel.oimModel.getImage>` method.

It takes two arguments, the image's size in pixels and the pixel size in mas.

.. code-block:: ipython3
    
    im = mUDPt.getImage(512, 0.1)
    plt.figure()
    plt.imshow(im**0.2)

.. image:: ../../images/basicModel_imshow.png
  :alt: Alternative text   
  
We plot the image with a 0.2 power-law to make the uniform disk components visible:
Both components have the same total flux but the uniform disk is spread on many more
pixels.

The image can also be returned as an ``astropy hdu`` object (instead of a ``numpy array``)
setting the ``toFits`` keyword to ``True``.
The image will then contained a header with the proper fits image keywords
(NAXIS, CDELT, CRVAL, etc.).

.. code-block:: ipython3
    
    im = mUDPt.getImage(256, 0.1, toFits=True)
    print(im)
    print(im.header)
    print(im.data.shape)


.. parsed-literal::
  
    ... <astropy.io.fits.hdu.image.PrimaryHDU object at 0x000002610B8C22E0>
    
    SIMPLE  =                    T / conforms to FITS standard                      
    BITPIX  =                  -64 / array data type                                
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                  256                                                  
    NAXIS2  =                  256                                                  
    EXTEND  =                    T                                                  
    CDELT1  = 4.84813681109536E-10                                                  
    CDELT2  = 4.84813681109536E-10                                                  
    CRVAL1  =                    0                                                  
    CRVAL2  =                    0                                                  
    CRPIX1  =                128.0                                                  
    CRPIX2  =                128.0                                                  
    CUNIT1  = 'rad     '                                                            
    CUNIT2  = 'rad     '                                                            
    CROTA1  =                    0                                                  
    CROTA2  =                    0                                                 
    
    (256, 256)
    

.. note::

    Currently only **regular** grids in wavelength and time are allowed when exporting
    to fits-image format. If specified, the **wl** and **t** vectors need to be regularily
    sampled. The easiest way is to use the 
    `numpy.linspace <https://numpy.org/doc/stable/reference/generated/numpy.linspace.html>`_
    function.

    If their sampling is irregular an error will be raised.


    
Using the :func:`oimModel.saveImage <oimodeler.oimModel.oimModel.saveImage>` method
will also return an image in the fits format and save it to the specified fits file. 

.. code-block:: ipython3
   
    im = mUDPt.saveImage("modelImage.fits", 256, 0.1)


.. note::

    The returned image in fits format will be 2D, if  time and wavelength are not
    specified, or if they are numbers, 3D if one of them is an array, and 4D if both
    are arrays.


Alternatively, we can use the :func:`oimModel.showModel <oimodeler.oimModel.oimModel.showModel>`
method which take the same argument as the getImage, but directly create a plot with
proper axes and colorbar.

.. code-block:: ipython3

    figImg, axImg = mUDPt.showModel(512, 0.1, normPow=0.2)


.. image:: ../../images/basicModel_showModel.png
  :alt: Alternative text  

Getting the model Complex Coherent Flux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In most of the cases the user won't use directly the :func:`oimModel.getComplexCoherentFlux <oimodeler.oimModel.oimModel.getComplexCoherentFlux>` 
method to retrieve the model complex coherent flux for a set of coordinates but will create  :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>`
or a  :func:`oimSimulator <oimodeler.oimFitter.oimFitter>` that will contain the instance of :func:`oimModel <oimodeler.oimModel.oimModel>`
and some interferometric data in an :func:`oimData <oimodeler.oimData.oimData>` to simulate interferometric quantities from the model at the 
spatial frequenciesfrom our data.  This will be covered in the XXXXXXXXXXX section.

Nevertheless, in some cases and for explanatory purposes we will directly use this methods in the following example.
Without the :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` class, the :func:`oimModel <oimodeler.oimModel.oimModel>`
can only produce complex coherent flux (i.e., non normalized complex visibility) for a vector of spatial frequecies and wavelengths. 

.. code-block:: ipython3

    wl = 2.1e-6
    B = np.linspace(0.0, 300, num=200)
    spf = B/wl


Here, we have created a vector of 200 spatial frequencies, for baselines ranging from 0 to 300 m at an observing wavelength of 2.1 microns.

We can now use this vector to get the complex coherent flux (CCF) from our model. 
    

.. code-block:: ipython3

    ccf = mUDPt.getComplexCoherentFlux(spf, spf*0) 

    
The :func:`oimModel.getComplexCoherentFlux <oimodeler.oimModel.oimModel.getComplexCoherentFlux>`
method takes four parameters: 

- the spatial frequencies along the East-West axis (u coordinates in cycles/rad), 
- the spatial frequencies along the North-South axis (v coordinates in cycles/rad), 

and optionally,

- the wavelength (in meters)
- time (mjd)

Here, we are dealing with grey and time-independent models so we don't need to specify the wavelength. 
Additionnally, as our models are circular, we don't care about the baseline orientation.
That why we set the North-South component of the spatial frequencies to zero.

We can now plot the visibility from the CCF as the function of the spatial frequencies:

.. code-block:: ipython3

    v = np.abs(ccf)
    v = v/v.max()
    plt.figure()
    plt.plot(spf, v)
    plt.xlabel("spatial frequency (cycles/rad)")
    plt.ylabel("Visbility")


.. image:: ../../images/basicModel_vis0.png
  :alt: Alternative text  


Let's finish this example by creating a figure with the image and visibility
for all the previously created models.

.. code-block:: ipython3

    models = [mPt, mUD, mG, mR, mUDPt]
    mNames = ["Point Source", "Uniform Disk", "Gausian", "Ring",
              "Uniform Disk + Point Source"]

    fig, ax = plt.subplots(2, len(models), figsize=(
        3*len(models), 6), sharex='row', sharey='row')

    for i, m in enumerate(models):
        m.showModel(512, 0.1, normPow=0.2, axe=ax[0, i], colorbar=False)
        v = np.abs(m.getComplexCoherentFlux(spf,  spf*0))
        v = v/v.max()
        ax[1, i].plot(spf, v)
        ax[0, i].set_title(mNames[i])
        ax[1, i].set_xlabel("sp. freq. (cycles/rad)")


.. image:: ../../images/basicModel_all.png
  :alt: Alternative text 



Types of components
-------------------

The code corresponding to this section is available in `TypesOfComponents.py <https://github.com/oimodeler/oimodeler/blob/main/examples/Modules/TypesOfComponents.py>`_

**oimodeler** components are of three different types:

| 1. the components defined in the Fourier space by an analytical formula.
 | They inherit from the  :func:`oimComponentFourier <oimodeler.oimcomponent.oimComponentFourier>` class
|2. the components defined by their 2D intensity map in the image space.
 | They inherit from the  :func:`oimComponentImage <oimodeler.oimcomponent.oimComponentImage>` class
|3. the components defined by their 1D intensity profile in the image space.
 | They inherit from the  :func:`oimComponentRadialProfile <oimodeler.oimcomponent.oimComponentRadialProfile>` class


Basic Fourier components
~~~~~~~~~~~~~~~~~~~~~~~~

In the table below is a list of the current Fourier-based components, which all derived from
the :func:`oimComponentFourier <oimodeler.oimComponent.oimComponentFourier>` semi-abstract class.

.. csv-table:: Available Fourier based components
   :file: table_components_fourier.csv
   :header-rows: 1  
   :delim: |
   :widths: auto

To print the comprehensive list of Fourier-based compnents you can type:

.. code-block:: ipython3

    print(oim.listComponents(componentType="fourier"))

.. parsed-literal::

    ['oimComponentFourier', 'oimPt', 'oimBackground', 'oimUD', 'oimEllipse', 'oimGauss', 'oimEGauss', 'oimIRing',
     'oimEIRing', 'oimRing', 'oimRing2', 'oimERing', 'oimERing2', 'oimESKIRing', 'oimESKGRing', 'oimESKRing', 'oimLorentz',
     'oimELorentz', 'oimLinearLDD', 'oimQuadLDD', 'oimPowerLawLDD', 'oimSqrtLDD', 'oimAEIRing', 'oimAERing', 'oimBox',
     'oimGaussLorentz', 'oimStarHaloGaussLorentz', 'oimStarHaloIRing']
     
.. note:: If you want to have more information on a component (for instance, on its paramaters) you can use python **help** function.

Although simple, these components can allow to build complex models, For instance, Chromaticity and/or time-dependency 
can be added to any parameters of these components to build more complex models.
We will see this in details in the :ref:`Advanced parameters` section.

.. note:: 
    Models using Fourier-based components are usually faster to run as they use a simple function to compute the 
    complex Coherent Flux whereas imaged-based used FFT or Hankel-Transform (for radial profile) 


Image components
~~~~~~~~~~~~~~~~

**oimodeler** allows to use components described in the image space. This can be done by subclassing the semi-abstract
:func:`oimComponentImage <oimodeler.oimcomponent.oimComponentImage>` class.

In the table below is a list of the current image-plan components:

.. csv-table:: Available Image plane components
   :file: table_components_image.csv
   :header-rows: 1  
   :delim: |
   :widths: auto

To print the comprehensive list of image-based compnents you can type:

.. code-block:: ipython3

    print(oim.listComponents(componentType="image"))


Describing an object by its intensity distribution instead of its Fourier transform can be useful in three cases:

1. the component cannot be described using an analytical formula in the Fourier space but can be described by an analytical formula in the image plan
2. the component cannot be described by a simple analytical formula even in the image space but an image can easily be computed, for instance with a iterative code
3. the user want to use external code such as images from a radiative transfert model


Here are three examples of these three kind of image components implemented in **oimodeler**.

The first one is a spiral implemented as :func:`oimSpiral <oimodeler..oimCustomComponents.oimSpiral.oimSpiral>`.
Its implementation is descrbided in details in the ex

.. code-block:: ipython3

    spiral = oim.oimSpiral(dim=256, fwhm=20, P=0.1, width=0.2, pa=30, elong=2)
    mspiral = oim.oimModel(spiral)


The second one is a simple simulation of a fast-rotator using the Roche model and including the gravity-darkening.
It is implemented as :func:`oimFastRotator <oimodeler..oimCustomComponents.oimFastRotator.oimFastRotator>` and a full description is given in XXXXXXXXXXXXXXXXXXXX.

.. code-block:: ipython3

    frot = oim.oimFastRotator(dpole=5, dim=128, incl=-50,rot=0.99, Tpole=20000, beta=0.25,pa=20)
    mfrot = oim.oimModel(frot)

Finally, the last one is an output from the radiative transfer code  `RADMC3D <https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/>`_ simulating the inner part of a dusty disk around the B[e] star FS CMa.
The simulation was made from 1.5 to 13μm. and the output was saved as a chromatic image-cube in the fits  format with proper axes descruibed in the header (size of pixel in x, y and wavelength). We use the  :func:`oimComponentFitsImage <oimodeler..oimComponents.oimComponentFitsImage>` class to load the image as a image-components. 

.. code-block:: ipython3

    radmc3D_fname = product_dir / "radmc3D_model.fits"
    radmc3D = oim.oimComponentFitsImage(radmc3D_fname,pa=180)
    mradmc3D = oim.oimModel(radmc3D)

For more information on using fits images in oimodeler, read the next section.

.. image:: ../../images/componentImages_images.png
  :alt: Alternative text 

Unlike when using



Fits images component
~~~~~~~~~~~~~~~~~~~~~

Radial-Profile components
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Available radial profile components
   :file: table_components_radial.csv
   :header-rows: 1  
   :delim: |
   :widths: auto


.. _Advanced parameters:

Advanced parameters
-------------------

Linking parameters
~~~~~~~~~~~~~~~~~~

chromatic & time-dependent Interpolator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

