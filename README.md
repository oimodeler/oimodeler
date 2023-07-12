<h1 align="center">
<img src="images/oimodelerlogo.png" width="600">
</h1><br>

[![License: GNU](https://img.shields.io/badge/License-GNU-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
![Lifecycle:
Early Development](https://img.shields.io/badge/lifecycle-EarlyDevelopment-orange.svg)

Oimodeler is a modular modelling tool for optical interferometry


* Installation: https://oimodeler.readthedocs.io/en/latest/installation.html
* Usage: https://oimodeler.readthedocs.io/en/latest/getting_started.html
* Documentation https://oimodeler.readthedocs.io/en/latest/
* Contributing: https://oimodeler.readthedocs.io/en/latest/expandingSoftware.html
* Source code: https://github.com/oimodeler/oimodeler
* Bug reports: https://github.com/oimodeler/oimodeler/issues

## Features
### Modules:
* **oimModel** : Create models with various components as bricks
* **oimComponent** : component for building models
* **oimParam** : Component parameters and parameters interpolators
* **oimData** :  Handle interferometric, spectroscopic and photometric data
* **oimDataFilter** : Filtering and modifying data (wavlengths range cut, smoothing, removing flags...)
* **oimSimulator** : Main class holding evertyhing together and producing final results :plots, tables...
* **oimFitter** : Define and perform model-fitting
* **oimPlot** : Plotting tools
* **oimUtils** : Various utility for optical-interferometry

### :warning: Under Development :warning:
The following modules have only been partially implemented
* oimModel:  Basic functionnalities for building model from all kind of components 
* oimComponent: Many basics Fourier-based components, possiblity to build radial profile and images based components, import of image at fits format and a few advanced components (kinematic disk, fast-rotator)
* oimParam: Possibility to link parameters, normalize them and various interpolators in wavelength and time
* oimData: Interferometric data and basic wrapper for photometric/spectroscopic data
* oimDataFilter: Wavelength range cut and shift, data smoothing and binning, datatype selection, flagging based on criteria
* oimSimulator: Simulated data and chi2 computation with filtering
* oimFitter:  Basic MCMC fitter based on emcee module.
* oimPlot: Various plots for oifits data (see examples below)
* oimUtils: Helpers for oifits format: retrieve baselines info, create new tables, manipulated data...
)

## Examples

For more examples see: https://oimodeler.readthedocs.io/en/latest/examples.html

> Various example scripts are available in the `examples` directory.<br>

This plot has been created with the `createModelChromatic.py` script:

![boo](./images/createModelChromatic.png)

Here is the resulting plot from the `oimodel_Create_simulator_data.py` script.
It plots data of a partly resolved binary created with:
- the [ASPRO](https://www.jmmc.fr/english/tools/proposal-preparation/aspro/) software from JMMC (including realistic noise)
- oimodeler using a shifted uniform disk + unresolved component

![boo](./images/oimodel_Create_simulator_data.png)

Here is a plot showing a model consisting of a fast rotating star plus a uniform disk.
Chromatic images of fast rotator are computed with an external function encapsulated into a oimodeler component.
The uniform disk is a simple Fourier-based component. The code is in the `createCustomComponentImageFastRotator.py`

![boo](./images/customCompImageFastRotatorImageAndVis.png)

## Contact
[Anthony Meilland](https://github.com/AnthonyMeilland)
