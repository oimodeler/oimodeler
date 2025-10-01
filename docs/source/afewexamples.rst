A Few examples
==============

Here are some plots from the 
`createModelChromatic.py <https://github.com/oimodeler/oimodeler/blob/main/examples/Other/createModelChromatic.py>`_
example, showcasing various chromatic-geometric models and their simulated visibilities.

.. image:: ../../images/createModelChromatic.png
  :alt: Alternative text

Here is an example from the `exampleOimSimulator.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimSimulator.py>`_ 
script showing high-quality plots of some MATISSE LM-band data and a model created with **oimodeler**. 
In this case, the data were simulated using the `APSRO <https://www.jmmc.fr/english/tools/proposal-preparation/aspro/>`_ 
software from the `JMMC <http://www.jmmc.fr/>`_.


.. image:: ../../images/oimodel_Create_simulator_data.png
  :alt: Alternative text


Here is an example from the `exampleOimFitterEmcee.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimFitterEmcee.py>`_ script showing:

- Parameter values (with chi2 color scale) for walkers from an emcee run on a binary model
- Corner plots for the 4 free parameters: binary’s x and y positions, diameter, and flux

The data is again simulated using `APSRO <https://www.jmmc.fr/english/tools/proposal-preparation/aspro/>`_.


.. image:: ../../images/SimpleFitWalkers.png
  :alt: Alternative text

  
.. image:: ../../images/SimpleFitCorner.png
  :alt: Alternative text


Here is a plot of a model combining a fast-rotating star and a uniform disk. Chromatic images of 
the fast rotator are computed using an external function wrapped in an **oimodeler** component. 
The uniform disk is a simple Fourier-based component. The code is in 
`createCustomComponentImageFastRotator.py <https://github.com/oimodeler/oimodeler/blob/main/examples/ExpandingSoftware/createCustomComponentImageFastRotator.py>`_.

.. image:: ../../images/customCompImageFastRotatorImageAndVis.png
  :alt: Alternative text
