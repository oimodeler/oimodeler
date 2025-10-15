.. image:: ../../images/oimodeler_title.png
  :alt: Alternative text

.. raw:: html

   <div style="height: 0; visibility: hidden;">

Index
=====

.. raw:: html

   </div>

.. rst-class:: hide-me

    Introduction
    ============

The **oimodeler** project aims to develop modular, easily expandable Python software for optical
interferometry modeling. Initiated in late 2021, it now offers a fully functional beta version.

It supports manipulation of OIFITS format data, building complex models from diverse components,
simulating data at observed spatial frequencies, computing chi2, performing fits (MCMC or others),
and easy result plotting.

Components can be defined in image or Fourier space via analytical formulas or precomputed images.
Model parameters or components may be chromatic and/or time-dependent.

The software is modular and object-oriented, facilitating expansion by creating new components or
features from abstract classes.

**oimodeler** was presented at the 2024 SPIE conference in Yokohama (paper on
(`HAL <https://cnrs.hal.science/hal-04797236>`_, `ADS <https://ui.adsabs.harvard.edu/abs/2024SPIE13095E..2WM/abstract>`_,
`BibTeX <https://ui.adsabs.harvard.edu/abs/2024SPIE13095E..2WM/exportcitation>`_).

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
   :caption: Examples   
   
   notebooks
    
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
   :caption: Expanding oimodeler

   expanding


   
.. toctree:: 
   :hidden: 
   :caption: References
   
   api   
   news
   ackn

