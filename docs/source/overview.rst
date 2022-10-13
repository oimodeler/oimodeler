Overview
========

oimodeler is a project aiming at developping a modular and easily expandable python-based modelling software for optical interferometry. The project started end of 2021, and the software is currently at a very early stage of development. 


Software modularity
-------------------

As shown in the diagram below, the oimodeler software is composed of 7 modules. Models can be created with the **oimModel** module and various **oimComponents**.  Interferometric data can be loaded into a **oimData** object from standard oifits files. This module also allows to load flux/spectra in various format.**oimData** can optionallly be filtered using the **oimDataFilter** class with various **oimFilterComponents**.**oimData** and **oimModel** are plugged into a **oimSimulator** to simulate data from the model at the same spatial/spectral coordinates than the data. The module also allow to compute the model/data chi2.The **oimSimulator** ( or directly **oimData** and **oimModel**) can be plugged into a **oimFitter** to perform model fitting. Currently only a simple emcee fitter is implemented. **oimPlots** contains plotting functions for oifits data and oimodeler objects.. **oimUtils** contains various functions to manipulate oifits data

.. image:: _static/diagram.png
  :alt: Alternative text



Modules 
-------

oiModel
^^^^^^^

The oimModel module is dedicated to the creation of models for optical interferometry. Models are modular and composed of one or many oimComponents.
The oimComponents can be analytical models defined in the Fourier plan (e.g. Uniform Disks, 2D-Gaussian distribution...) or image plan (useful if no analytical fromula in the Fourier plan is available). In the image plan, the code allow to build component in 1D (i.e. radial profile), 2D (images, or chromatic radial profile), or 3D (chromatic image-cubes). oimComponent can be easily subclassed to create new customs components.


.. warning::
    Only basic analytical models in the Fourier plan are fully implemented. The current implementation of the image-plan model may evolved in upcoming versions of the code.

oimData
^^^^^^^

The oimDate module allows to encapsulate interferometric (but also photometric and spectroscopic) data. The oimData allow to retrieve the original data as a list of astropy.io.fits.hdulist (oimData.data) but also provide optimization of the data as vector/structure for faster model fitting. OimDataFilter can be plugged to the oimData to filter the data.

oimSimulator
^^^^^^^^^^^^

The oimSimulator contains an oimModel object (oimSimulator.model) and oimData object (oimSimulator.data). It allows to simulate a new dataset (stored in oimSimulator.simulatedData) with the same quantities and spatial/spectral coordinated of an oimData and using the oimModel. It also allows to compute chi2 for data and model comparison.

oimFitter
^^^^^^^^^

The oimFitter module is dedicated to model fitting. The parent class oimFitter is an abstract class that can be subclassed to implement various fitter. Currently only a simple emcee-based fitter is implemented. 

oimPlots
^^^^^^^^

oimPlots contains various plotting tools for oifits data and oimodeler objects. The oimAxes is a subclass of matplotlib.pyplot.Axes class with methods dedicated to produce plot data oifits data in the astropy.io.fits.hdulist format.

oimDataFilter
^^^^^^^^^^^^^
oimDataFilter is dedicated to filtering and modifying data. It allows data selection (truncating, removing arrays, and flagging) based on various criteria (wavelengths, SNR ...), and other data manipulation such as smoothing and binning of the data.

oimUtils
^^^^^^^^

oimUtils containscontains various functions to manipulate oifits data such as : retrieving baselines names, length, orientation, spatial frequencies, and creating new oifits array in astropy.io.fits.hudlist format.

Expandability
-------------

oimodeler is written to allow easy implementation of new model component, fitters, data filters or loader (e.g. for non-oifits format data), and plots. Feel free to contact us if you develop custom features and want them to be included in the oimodeler distribution.

