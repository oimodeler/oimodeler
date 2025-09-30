A Few examples
==============

Here are some plots for the `createModelChromatic.py <https://github.com/oimodeler/oimodeler/blob/main/examples/Other/createModelChromatic.py>`_
example showing various chromatic-geometric models and the corresponding simulated
visibilities.

.. image:: ../../images/createModelChromatic.png
  :alt: Alternative text


Here is an example from the `exampleOimSimulator.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimSimulator.py>`_ script showing
high-end plots of some MATISSE LM-band data and a model create with **oimodeler**.
In that case the data were simulated using the `APSRO <https://www.jmmc.fr/english/tools/proposal-preparation/aspro/>`_
software from `JMMC <http://www.jmmc.fr/>`_.

.. image:: ../../images/oimodel_Create_simulator_data.png
  :alt: Alternative text


Here is an example from the `exampleOimFitterEmcee.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimFitterEmcee.py>`_ script showing:

- The parameters values (and corresponding chi2 as a colorscale) for walkers from a emcee run on a binary model
- The famous corners plots for the 4 free parameters: x and y positions of the binary, diameter and flux

Again, the data is simulated with `APSRO <https://www.jmmc.fr/english/tools/proposal-preparation/aspro/>`_.

.. image:: ../../images/SimpleFitWalkers.png
  :alt: Alternative text

  
.. image:: ../../images/SimpleFitCorner.png
  :alt: Alternative text


Here is a plot showing a model consisting of a fast rotating star plus a uniform disk.
Chromatic images of fast rotator are computed with an external function encapsulated
into a **oimodeler** component. The uniform disk is a simple Fourier-based component.
The code is in the `createCustomComponentImageFastRotator.py <https://github.com/oimodeler/oimodeler/blob/main/examples/ExpandingSoftware/createCustomComponentImageFastRotator.py>`_

.. image:: ../../images/customCompImageFastRotatorImageAndVis.png
  :alt: Alternative text
