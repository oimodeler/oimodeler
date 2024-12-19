.. image:: ../../images/oimodeler_title.png
  :alt: Alternative text

.. rst-class:: hide-me

    Introduction
    ============

The oimodeler project aims at developing a modular and easily expandable python-based
modelling software for optical interferometry. The project started end of 2021, and the software is now in a fully functional beta version.

It allows to manipulate data in the OIFITS format, build complex models from various
components, simulate data from the model at the spatial frequencies of your observations,
compute chi2, perform model fitting (using mcmc or other fitters), and plot results easily. Components can be defined in the image or Fourier plan using analytical formulas or precomputed images. Components or model parameters can be chromatic and/or time dependent. The software is modular and object oriented, in order to make it easy to expand it by creating new model components or other features from abstract classes. 

.. image:: _static/oimComponents.png
  :alt: Alternative text


.. toctree::
   :hidden:
   :caption: Introduction
    
   afewexamples
   modularity
   installation
   getting_started
    
.. toctree:: 
   :hidden:
   :caption: Modules Description
   
   data
   models
   simulator
   fitter
   plot
   utils
   
.. toctree::
   :hidden: 
   :caption: Examples   
   
   notebooks
   
.. toctree:: 
   :hidden: 
   :caption: References
   
   api   
   news
   ackn

