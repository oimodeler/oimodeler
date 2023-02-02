.. oimodeler documentation master file, created by
   sphinx-quickstart on Wed Nov 24 16:00:55 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: ../../images/oimodelerlogo.png
  :alt: Alternative text

oimodeler
=========

The oimodeler project aims at developping a modular and easily expandable python-based modelling software for optical interferometry. The project started end of 2021, and the software is currently at an early stage of development. 

It allows to manipulate data in the oifits format, build complex models from various components, simulate data from the model at the spatial frequencies of your observations, compute chi2, perform model fitting (using mcmc or other fitters), and plot results easily. Components can be defined in the image or Fourier plan using analytical formulas or precomputed images. Components or model parameters can be chromatic and/or time dependent. The software is modular and object oriented, in order to make it easy to expand it by creating new model components or other features from abstract classes. 

.. warning::

    The software is in early development:  
        - Models : in Fourier and Image plans with chromaticity and time dependence
        - Data : interferometric data only. No photometric or spectroscopic data.
        - Data Filters : Filtering wavelength range, and data type (VIS2DATA, VISAMP...)    
        - Fitters : Implementation of a basic emcee-based fitter with plots for results
        - Plots : Basics plots of oifits data and uv-plan plot
        - Utils : miscs utilities for oifits data (creating and modifying array, getting info..)   

    No module is complete and have been fully verified up to now!

A Few examples
==============

Here are some plots for the `createModelChromatic.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/createModelChromatic.py>`_ example showing various chromatic-geometric models and the corresponding simulated Visibilities.

.. image:: ../../images/createModelChromatic.png
  :alt: Alternative text


Here is an example from the :ref:`createSimulator.py <createSimulator>` script showing high-end plots of some MATISSE LM-band data and a model create with oimodeler . In that case the data were simulated using the `APSRO <https://www.jmmc.fr/english/tools/proposal-preparation/aspro/>`_ software from `JMMC <http://www.jmmc.fr/>`_.

.. image:: ../../images/oimodel_Create_simulator_data.png
  :alt: Alternative text

Here is an example from the :ref:`simpleFitEmcee.py <createSimulator>` script showing :

- the parameters values (and corresponding chi2 as a colorscale) for walkers from a emcee run on a binary model
- the famous corners plots for the 4 free parameters: x and y positions of the binary, diameter and flux

 Again, the data is simulated with `APSRO <https://www.jmmc.fr/english/tools/proposal-preparation/aspro/>`_.

.. image:: ../../images/SimpleFitWalkers.png
  :alt: Alternative text
  
.. image:: ../../images/SimpleFitCorner.png
  :alt: Alternative text


Here is a plot showing a model consisting of a fast rotating star plus a uniform disk. Chromatic images of fast rotator are computed with an external function encapsulated into a oimodeler component. The uniform disk is a simple Fourier-based component. The code is in the createCustomComponentImageFastRotator.py

.. image:: ../../images/customCompImageFastRotatorImageAndVis.png
  :alt: Alternative text


License 
=======

Copyright 2021-2022 Anthony Meilland and `contributors <https://github.com/oimodeler/oimodeler/graphs/contributors>`_.

oimodeler is a free software distributed under GNU General Public License. 

For details see the `LICENSE <https://github.com/oimodeler/oimodeler/blob/main/LICENSE>`_.

Acknowledgment
==============

The oimodeler package is developed with support from the `VLTI/MATISSE <https://www.matisse.oca.eu/fr/accueil-matisse>`_ consortium and the ANR project `MASSIF <https://www.anr-massif.fr>`_ .


Contact
=======

If you have any question or comments contact `Anthony Meilland <mailto://ame@oca.eu>`_.

You can use `Github <https://github.com/oimodeler/oimodeler/issues>`_ interface for bug reports.

Table of Contents
=================

.. toctree::
    overview
    installation
    getting_started
    examples
    api
    :maxdepth: 4


* :ref:`genindex`
* :ref:`modindex`

