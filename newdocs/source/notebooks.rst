:tocdepth: 2


.. _notebooks:
    
The basics of oimodeler 
=======================

Here are a few very python scripts from **oimodeler** 0.8 documentation demonstrating the basics of the package. 
Theses examples uses simulated dataset computed with the `APSRO <https://www.jmmc.fr/english/tools/proposal-preparation/aspro/>`_
software from the `JMMC <http://www.jmmc.fr/>`_. 

Here is the list of examples:

- `Loading oifits data <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimData.py>`_
- `Basic models <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/basicModels.py>`_
- `Precomputed fits-formated image <https://github.com/oimodeler/oimodeler/tree/main/examples/BasicExamples/FitsImageModel.py>`_
- `Data/model comparison <https://github.com/oimodeler/oimodeler/tree/main/examples/BasicExamples/exampleOimSimulator.py>`_
- `Running a mcmc fit <https://github.com/oimodeler/oimodeler/tree/main/examples/BasicExamples/exampleOimFitterEmcee.py>`_
- `Plotting oifits data <https://github.com/oimodeler/oimodeler/tree/main/examples/BasicExamples/exampleOimPlot.py>`_
- `Filtering data <https://github.com/oimodeler/oimodeler/tree/main/examples/BasicExamples/exampleOimDataFilter.py>`_
- `Building Complex models <https://github.com/oimodeler/oimodeler/tree/main/examples/AdvancedExamples/complexModels.py>`_
- `Precomputed chromatic image-cubes <https://github.com/oimodeler/oimodeler/tree/main/examples/AdvancedExamples/FitsImageCubeModels.py>`_
- `Parameters Interpolators <https://github.com/oimodeler/oimodeler/tree/main/examples/AdvancedExamples/paramInterpolators.py>`_
- `Fitting a chromatic model <https://github.com/oimodeler/oimodeler/tree/main/examples/AdvancedExamples/ChromaticModelFit.py>`_

Modelling real VLTI data
========================

These tutorials were initially developed for the `12th VLTI School of Interferometry <https://vltischool2024.sciencesconf.org/>`_.
Using real datasets from VLTI instruments (PIONIER, MATISSE, and GRAVITY) they demonstrate how to use 
**oimodeler** to perform model fitting, including chromaticity and kinematics, for stars and 
circumstellar environments. 

The tutorials are provided as notebooks, containing Python code, explanations, and questions related to 
model-fitting issues.



.. list-table:: 
   :class: borderless
   :widths: 1 1 

   * - .. figure:: _static/notebooks/notebook_Ex1_thumbnail.png

            `PIONIER observation of Canopus 
            <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex1_canopus.ipynb>`_ 
            
            In this first example we will demonstrate how to use oimodeler to create a simple model and perform 
            model-fitting. We will also introduce the concepts of :math:`\chi^2_r` minimization and of local 
            and global minima.
            
            We will use VLTI/PIONIER data obtained on the giant star Canopus as an example.

     - .. figure:: _static/notebooks/notebook_Ex2_thumbnail.png

            `The Binary star 94 Aqr MATISSE observation 
            <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex2_94_Aqr.ipynb>`_

            In this second exercise we will use VLTI/MATISSE data obtained on the binary star 94 Aqr. We will 
            demonstrate probblematic around binary star fitting including local minima and degeneracy.
            
            We will use MCMC and grid fitters. 



   * - .. figure:: _static/notebooks/notebook_Ex3_thumbnail.png

            `Building chromatic model for the YSO HD 179218 
            <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex3_HD179278.ipynb>`_ 

            With the advent of routine spectro-interferometry, the use of simple achromatic models to reproduce
            interferometric observations often appears insufficient. **oimodeler** provides different tools to build 
            chromatic models. 
            
            One method is to use parameter interpolators. As we will see, such feature can be very useful when 
            modelling highly chromatic objects. To illustrate the concept, let us consider the case of L-band MATISSE 
            data obtained on the Herbig star HD17921.


     - .. figure:: _static/notebooks/notebook_Ex4_thumbnail.png

            `MWC297: continuum vs line emission 
            <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex4_MWC297.ipynb>`_ 

            Here, we illustrate the use of a chromatic interpolator in the case of high resolution MATISSE data 
            obtained on the massive and luminous YSO MWC 297. We will focus on the gaseous Br:math:`_\gamma` 
            emission around 4.05:math:`_\mu`m, and its adjacent continuum (produced by the dust).

            The aim is to determine the size of those two emitting components, thus informing us about the
            physical process at the origin of the gas emission.


   * - .. figure:: _static/notebooks/notebook_Ex5_thumbnail.png

            `Kinematics of the Gaseous disk around HD58647 
            <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex5_HD58647.ipynb>`_ 

            Spectro-interferometry in atomic lines allows to constrain dynamics thanks to the Doppler effect. 
            In this section we will focus on rotating disk around stars. Such disk are found around YSO, some 
            evolved stars and of course classical Be stars.

            The idea of this section is not to provide a comprehensive guide of how to model spectro-interferometric 
            measurements of rotating disks but to introduce the subject and show how to use the oimKinematicDisk 
            component implemented in **oimodeler**.


     - .. figure:: _static/notebooks/notebook_Ex6_thumbnail.png

            `Radiative Transfer model of FS CMa 
            <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex6_FSCMa.ipynb>`_ 

            **oimodeler** contains a specific class, oimComponentFitsImage, to import precomputed images or 
            hyperspectral image-cubes such as outputs from radiative-transfer models in a fits-image format.
            The loaded image can be shifted, rotated and scaled if necessary.

            In this example we will load one image-cube created with the radiative transfer code RADMC3D and compare 
            it to the MATISSE observations of the B[e] star FS CMa 