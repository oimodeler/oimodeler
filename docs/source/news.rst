.. _news:

News and Changelog
==================

2025-10-01
----------

v0.9.0:  2nd Beta Version
^^^^^^^^^^^^^^^^^^^^^^^^^

New Features:
:::::::::::::

- added :math:`\chi^2_r` grid exploration: :func:`oimFitterGrid <oimodeler.oimFitter.oimFitterGrid>`

- added simple :math:`\chi^2_r` minimizer based on Minimize scipy function:  :func:`oimFitterMinimize <oimodeler.oimFitter.oimFitterMinimize>`

- added Dynamic Nested Sampling fitter based on  `dynesty <https://dynesty.readthedocs.io/>`_ package :  :func:`oimFitterDynesty <oimodeler.oimFitter.oimFitterDynesty>`

- New models: updated fast rotator, temperature gradient disk,

- new residual plots for simulator: :func:`plotWithResidual <oimodeler.oimSimulator.oimSimulator.plotWithResidual>`

- new documentation including examples from the 2024 VLTI School

2023-07-12
----------

v0.8.0:  First Beta Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Basic features of the oimodeler package are almost fully implemented although not fully tested and verified. Features still missing before v1.0:

- Components & models : Exploration of grid of models, binary orbit, rotating star, temperature gradient disk, blackbody-based flux interpolator, saving of model and parameters

- Fitter : Levenberg-Marquardt fitter with various global search patterns, more options for the emcee-base fitter, chain of fitters, saving fitter & results


New Features
::::::::::::

- add support to pathlib Path for :func:`oimData <oimodeler.oimData.oimData>` and :func:`oimData <oimodeler.oimPlot>` 

- new class :func:`oimFluxData <oimodeler.oimFluxData.oimFluxData>`: importing photometric and spectroscopic measurement into oimodeler

- New plot options for :func:`oimSimulator.plot <oimodeler.oimSimulator.oimSimulator.plot>`

- New method :func:`plotWlTemplate <oimodeler.oimSimulator.oimSimulator.plotWlTemplate>` for generating multi-plots (one per baseline) for interferometric as function of the wavelength for the :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>`  class

- added option `dataType` in :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` to allow to choose which data type (VIS2DATA, VISPHI...) to include in the chi2 computation without the use of data filtering

- added kwargs options from the corner.py module in the :func:`oimFitterEmcee.cornerPlot <oimodeler.oimFitter.oimFitterEmcee.cornerPlot>` method

- added many data filters/modifiers : 
    - :func:`oimKeepDataType <oimodeler.oimDataFilter.oimKeepDataType>` : specifying which data type to keep : VIS2DATA, VISPHI, ...
    - :func:`oimWavelengthShiftFilter <oimodeler.oimDataFilter.oimWavelengthShiftFilter>` : shifting wavelengths 
    - :func:`oimWavelengthSmoothingFilter <oimodeler.oimDataFilter.oimWavelengthSmoothingFilter>` : "Smoothing" data by convolution on x walvength-pixels
    - :func:`oimWavelengthBinningFilter  <oimodeler.oimDataFilter.oimWavelengthBinningFilter>`: binning data in wavelength 
    - :func:`oimFlagWithExpressionFilter  <oimodeler.oimDataFilter.oimFlagWithExpressionFilter>`: filtering out or in with an expression including columns from the oifits table (ex "VIS2ERR/VIS2DATA>0.05" will remove VIS2DATA with more than 5% error)
    - :func:`oimDiffErrFilter  <oimodeler.oimDataFilter.oimDiffErrFilter>` : computing differential errors from the rms in a specified spectral band
    - :func:`oimSetMinErrFilter  <oimodeler.oimDataFilter.oimSetMinErrFilter>` : set a minimum error on the data in % (for visibilities) or deg (for phase)
    
- rewritting of :func:`uvPlot <oimodeler.oimPlots.uvPlot>` function: now allow colorscale for baselines, configuration, array, file, and allow plotting as a function of spatial frequency instead of length.

- added new :func:`oimWlTemplatePlots <oimodeler.oimPlots.oimWlTemplatePlots>` class to produce figure with multiple wavelength-plots per baseline of all interferometric quantities.

- added functions to create all oifits extensions: :func:`createOiTarget <oimodeler.oimUtils.createOiTarget>` :func:`createOiArray <oimodeler.oimUtils.createOiArray>` :func:`createOiWavelength <oimodeler.oimUtils.createOiWavelength>`  :func:`createOiVis2 <oimodeler.oimUtils.createOiVis2>` :func:`createOiVis <oimodeler.oimUtils.createOiVis>`  :func:`createOiT3 <oimodeler.oimUtils.createOiT3>` :func:`createOiFlux <oimodeler.oimUtils.createOiFlux>` :func:`createOiTargetFromSimbad <oimodeler.oimUtils.createOiTargetFromSimbad>` 

Bugs fixed
::::::::::

- rewritting of :func:`oimFitterEmcee.walkerPlot <oimodeler.oimFitter.oimFitterEmcee.walkersPlot>` method to speed it up and adding ncolors option for number of color in colorscale (generation time of the plot is proportionnal to the number of colors)

- corrected many bugs in the :func:`oimPlot <oimodeler.oimPlots.oimPlot>` function : possibility to plot data (or set coloscale) as function of PA, baseline LENGTH ...





