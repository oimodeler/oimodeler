:tocdepth: 1

Modularity & Expandability
==========================


Description oimodeler modules
-----------------------------

As described below and shown in the diagram, **oimodeler** is a modular software:

- Models are created with the :mod:`oimModel <oimodeler.oimModel>` module and base components in 
  :mod:`oimComponent <oimodeler.oimComponent>`, which include model parameters from 
  :mod:`oimParam <oimodeler.oimParam>`.

- Interferometric data can be loaded from standard OIFITS files using the :mod:`oimData <oimodeler.oimData>` 
  module. This also supports loading flux/spectra in various formats with the 
  :func:`oimData.oimFluxData <oimodeler.oimFluxData>` class, which can be optionally filtered by the 
  :func:`oimData.oimDataFilter <oimodeler.oimData.oimDataFilter>` class, using filters from the 
  :mod:`oimDataFilter <oimodeler.oimDataFilter>` module.

- Data simulation and calculation are handled by the :mod:`oimSimulator <oimodeler.oimSimulator>` module, 
  which takes :func:`oimData.oimData <oimodeler.oimData.oimData>` and 
  :func:`oimModel.oimModel <oimodeler.oimModel.oimModel>` objects as input to simulate data at the same 
  spatial and spectral coordinates as the observations. It also computes model/data chi2.

- Model fitting is performed by fitters in the :mod:`oimFitter <oimodeler.oimFitter>` module, which also 
  takes :func:`oimData.oimData <oimodeler.oimData.oimData>` and 
  :func:`oimModel.oimModel <oimodeler.oimModel.oimModel>` classes as input.

- The :mod:`oimPlots <oimodeler.oimPlots>` module provides plotting functions for OIFITS data and 
  **oimodeler** objects.

- The :mod:`oimUtils <oimodeler.oimUtils>` module contains various utility functions for manipulating 
  OIFITS data.
.. image:: _static/diagram.png
  :alt: Alternative text


oimModel
^^^^^^^^

The :mod:`oimModel <oimodeler.oimModel>` module focuses on creating models for optical interferometry. 
Models are modular and consist of one or more 
:func:`oimComponent.oimComponent <oimodeler.oimComponent.oimComponent>` objects. They can generate 
complex coherent fluxes and images, which can then be integrated into an 
:func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` object and/or any fitter in the 
:mod:`oimFitter <oimodeler.oimFitter>` module for data analysis and modeling. See the :ref:`model` 
section for more details.


oimComponent
^^^^^^^^^^^^

The :mod:`oimComponent <oimodeler.oimComponent>` module manages model components that can be defined 
analytically in the Fourier plan (e.g., Uniform Disks, 2D Gaussian distributions) or in the image plane 
(useful when no analytical Fourier formula exists). An :mod:`oimComponent <oimodeler.oimComponent>` can 
also wrap external code, such as functions computing images, radial profiles, or hyperspectral cubes, 
or image-fit files (e.g., from radiative transfer models). Additionally, components can be easily 
inherited to create new custom components.

oimParam
^^^^^^^^

The :mod:`oimParam <oimodeler.oimParam>` module contains basic model parameters. Its 
:func:`oimParam.oimParam <oimodeler.oimParam.oimParam>` class defines component parameters 
(based on base classes from the :mod:`oimComponent <oimodeler.oimComponent>` module). It also 
includes parameter linkers (:func:`oimParam.oimParamLinker <oimodeler.oimParam.oimParamLinker>`), 
normalizers (:func:`oimParam.oimParamNormalize <oimodeler.oimParam.oimParamNormalize>`), and advanced 
interpolators (:func:`oimParam.oimParamInterpolator <oimodeler.oimParam.oimParamInterpolator>`) 
that enable building chromatic and time-dependent models.


oimData
^^^^^^^

The :mod:`oimData<oimodeler.oimData>` module allows to encapsulate interferometric
(also photometric and spectroscopic) data. The :func:`oimData.oimData <oimodeler.oimData.oimData>`
class contains the original OIFITS data as a list of 
`astropy.io.fits.hdulist <https://docs.astropy.org/en/stable/io/fits/api/hdulists.html>`_ 
but also provide optimization of the data as vector/structure for faster model fitting. 


oimFluxData
^^^^^^^^^^^

The :mod:`oimData <oimodeler.oimData>` module encapsulates interferometric, photometric, and spectroscopic data. 
The :func:`oimData.oimData <oimodeler.oimData.oimData>` class holds the original OIFITS data as a list of 
`astropy.io.fits.hdulist <https://docs.astropy.org/en/stable/io/fits/api/hdulists.html>`_ objects and 
also provides optimized data structures (vectors) for faster model fitting.


oimDataFilter
^^^^^^^^^^^^^

The :mod:`oimDataFilter <oimodeler.oimDataFilter>` module handles filtering and modifying data within 
:func:`oimData.oimData <oimodeler.oimData.oimData>` classes. It supports data selection (truncation, 
array removal, flagging) based on criteria like wavelength or SNR, as well as data manipulation 
methods such as smoothing, binning, and error computation.


oimSimulator
^^^^^^^^^^^^

The :mod:`oimSimulator <oimodeler.oimSimulator>` module is the core for comparing models 
(:func:`oimModel.oimModel <oimodeler.oimModel.oimModel>`) and data 
(:func:`oimData.oimData <oimodeler.oimData.oimData>`). It simulates datasets with matching 
quantities and spatial/spectral coordinates from both data and model. It also computes :math:`\chi^2_r`
for comparison. The :func:`oimSimulator.oimSimulator <oimodeler.oimSimulator.oimSimulator>` class 
is fully compatible with OIFITS2 and can simulate any data type from an OIFITS file 
(e.g., VIS2DATA, VISAMP in absolute, differential, and correlated flux).

oimFitter

^^^^^^^^^
The :mod:`oimFitter <oimodeler.oimFitter>` module is dedicated to model fitting. The parent class 
:func:`oimFitter.oimFitter <oimodeler.oimFitter.oimFitter>` is an abstract class to be inherited 
for implementing various fitters. Current fitters include an MCMC sampler based on the popular 
emcee library, a simple grid search, and a minimizer using the scipy minimize function.

oimPlots
^^^^^^^^

The :mod:`oimPlots <oimodeler.oimPlots>` module provides plotting tools for OIFITS data and 
**oimodeler** objects. The :func:`oimPlots.oimAxes <oimodeler.oimPlots.oimAxes>` subclass extends 
`matplotlib.pyplot.Axes <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axes.html>`_ 
with methods tailored to plot OIFITS data from the 
`astropy.io.fits.hdulist <https://docs.astropy.org/en/stable/io/fits/api/hdulists.html>`_ format.

oimUtils
^^^^^^^^

The :mod:`oimUtils <oimodeler.oimUtils>` module contains functions to manipulate OIFITS data, such as:

- Retrieving baseline names, lengths, orientations, and spatial frequencies
- Creating new OIFITS arrays in 
  `astropy.io.fits.hdulist <https://docs.astropy.org/en/stable/io/fits/api/hdulists.html>`_ format
and more.


Expandability
-------------

**oimodeler** is designed for easy implementation of new model components, fitters, data filters, 
parameter interpolators, data importers (e.g., for non-OIFITS formats), and plotting tools.

Feel free to contact `Anthony Meilland <mailto://ame@oca.eu>`_ if you develop custom features and 
would like them included in the **oimodeler** distribution. Alternatively, submit a pull request 
on the `GitHub repository <https://github.com/oimodeler/oimodeler>`_ and become an **oimodeler** 
contributor.

For bug reports and feature requests, please use the 
`GitHub issue tracker <https://github.com/oimodeler/oimodeler/issues>`_.

