:tocdepth: 1

..  _nimodeler:

nimodeler
=========


**nimodeler** is an experimental Python package designed to simulate nulling interferometry data compatible with the NIFITS format.

It is built on top of the following libraries:

- **NIFITS**: https://github.com/rlaugier/nifits
- **oimodeler**: https://github.com/oimodeler/oimodeler/

Project Status
--------------

This package is currently at a very early stage of development.

Current Features
----------------

At its current stage, **nimodeler** allows users to:

- Load a single NIFITS file
- Plot instrument responses for various channels: Photometric, Additive or Nulling
- Extract data from these channels
- Extract the differential nulling channel, when available
- Compute simulated data using an **oimodeler** model with the same instrument response
- Perform data–model comparisons using a simple chi-squared method

Some plots from  nimodeler and some VLTI/NOTT simulated data
------------------------------------------------------------


.. figure:: https://github.com/oimodeler/nimodeler/tree/main/images/test_nobackground_v3_channels_responses.png

    Channels response for a exoplanet simulation

.. figure::https://github.com/oimodeler/nimodeler/tree/main/images/flux_allchannels_nobackground_v3.png

    Data/Model comparison for an exoplanet simulation


.. figure:: https://github.com/oimodeler/nimodeler/tree/main/images/diff_channel_flux_nobackground_v3.png

    Data/Model comparison of differential null for an exoplanet simulation observed at 20 different hour angle position
    (full night of observation)

.. figure:: https://github.com/oimodeler/nimodeler/tree/main/imagess/position_exploration_nobackground_v3.png

    Grid exploration of the position of an exoplanet
