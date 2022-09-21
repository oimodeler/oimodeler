.. oimodeler documentation master file, created by
   sphinx-quickstart on Wed Nov 24 16:00:55 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

oimodeler 
=========

oimodeler is a project aiming at developping a modular and easily expandable python-based modelling software for optical interferometry. The project started end of 2021, and the software is currently at a very early stage of development. 


It will allow to manipulate data in the oifits format, build complex models from various components, simulate data from the model at the spatial frequencies of your observations, computed chi2, perform model fitting (using mcmc or other fitters), and plot results easily. It will allow many model components including gray or chromatic analytical models defined in Fourier or image plan, imported radial profile, images, image-cubes (with chromaticity), or the use of precomputed grids of models providing that the output is in fits image format. As the software is modular and object oriented, it is easy to expand it by creating new components by deriving abstract classes. 


Modules
--------

- **oimModel** : Create models with various components as bricks 
- **oimData** :  Handle interferometric, spectroscopic and photometric data
- **oimDataFilter** : Filtering and modifying data (wavelength-range cut, smoothing, removing flags...)  
- **oimSimulator** : Main class holding evertyhing together and producing final results :plots, tables...
- **oimFitter** : Define and perform model-fitting   
- **oimPlot** : Plotting tools
- **oimUtils** : Various utility for optical-interferometry


.. warning::

    The software is in very early development. Partially implemented modules are the following:  
    
    - oimModel: Working with gray and chromatic models defined in Fourier plan, early implementation of Image-plan based models
    - oimData class: No filtering, no optimization of data      
    - oimSimulator: Simulated data and chi2 computation (no filtering yet) 
    - oimFitter : Implementation of a basic emcee-based fitter
    - oimPlot: Basics plots of oifits data and uv-plan plot    
    - oimUtils: Spatial frequencies, baseline name, length and PA, create oifits arrays    
    
    No module is complete and have been fully verified up to now

A Few examples
--------------

Here are some plots for the :ref:`createModelChromatic.py <createModelChromatic>` example showing various chromatic-geometric models and the corresponding simulated Visibilities.

.. image:: ../../images/createModelChromatic.png
  :alt: Alternative text


Here is an example from the :ref:`createSimulator.py <createSimulator>` script showing high-end plots of some MATISSE LM-band data and a model create with oimodeler . In that case the data were simulated using the `APSRO <https://www.jmmc.fr/english/tools/proposal-preparation/aspro/>`_ software from `JMMC <http://www.jmmc.fr/>`_.

.. image:: ../../images/oimodel_Create_simulator_data.png
  :alt: Alternative text



Here is an example from the :ref:`simpleFitEmcee.py <createSimulator>` script showing :

- the parameters values (and corresponding chi2 as a colorscale) for walkers from a emcee run on a binary model
- the famous corners plots for the 4 free parameters: x and y positions of the binary, diameter and flux

 Again, the data is simulated with `APSRO <https://www.jmmc.fr/english/tools/proposal-preparation/aspro/>`.

.. image:: ../../images/SimpleFitWalkers.png
  :alt: Alternative text
  
.. image:: ../../images/SimpleFitCorner.png
  :alt: Alternative text

.. toctree::
    overview
    installation
    getting_started
    examples
    api
    :maxdepth: 2
   




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
