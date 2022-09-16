.. oimodeler documentation master file, created by
   sphinx-quickstart on Wed Nov 24 16:00:55 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

oimodeler documentation
=======================

OiModeler is a modular modeling tool for optical interferometric. It allows to load data in oifits format, build complex models from various components, simulate data from the model at the spatial frequencies of your observations, computed chi2, perform model fitting, and plot results easily. The model components includes gray or chromatic analytical models defined in Fourier or image plan, imported radial profile, image, image-cube (with chromaticity), or precomputed grids of models providing that the outpu is in fits image format. As the software is modular and object oriented, it is easy to expand it by creating new components by deriving abstract classes. 


.. image:: ../../images/createModelChromatic.png
  :width: 4600
  :alt: Alternative text


.. toctree::
   installation
   getting_started
   examples
   api
   :maxdepth: 2
   :caption: Contents:

Check out the :doc:`usage` section for further information, including how to
:ref:`install <installation>` the project.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
