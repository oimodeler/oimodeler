:tocdepth: 2


.. _notebooks:

The basics of oimodeler
=======================

Here are a few Python scripts from **oimodeler** 0.8 documentation demonstrating the basics of the package.  
These examples use simulated datasets computed with the `ASPRO <https://www.jmmc.fr/english/tools/proposal-preparation/aspro/>`_  
software from the `JMMC <http://www.jmmc.fr/>`_.

Below is the list of basic examples, with direct links to their source code:

- `Loading oifits data <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimData.py>`_
- `Basic models <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/basicModels.py>`_
- `Precomputed fits-formatted image <https://github.com/oimodeler/oimodeler/tree/main/examples/BasicExamples/FitsImageModel.py>`_
- `Data/model comparison <https://github.com/oimodeler/oimodeler/tree/main/examples/BasicExamples/exampleOimSimulator.py>`_
- `Running a MCMC fit <https://github.com/oimodeler/oimodeler/tree/main/examples/BasicExamples/exampleOimFitterEmcee.py>`_
- `Plotting oifits data <https://github.com/oimodeler/oimodeler/tree/main/examples/BasicExamples/exampleOimPlot.py>`_
- `Filtering data <https://github.com/oimodeler/oimodeler/tree/main/examples/BasicExamples/exampleOimDataFilter.py>`_
- `Building complex models <https://github.com/oimodeler/oimodeler/tree/main/examples/AdvancedExamples/complexModels.py>`_
- `Precomputed chromatic image cubes <https://github.com/oimodeler/oimodeler/tree/main/examples/AdvancedExamples/FitsImageCubeModels.py>`_
- `Parameters interpolators <https://github.com/oimodeler/oimodeler/tree/main/examples/AdvancedExamples/paramInterpolators.py>`_
- `Fitting a chromatic model <https://github.com/oimodeler/oimodeler/tree/main/examples/AdvancedExamples/ChromaticModelFit.py>`_

*We recommend starting with the basic examples before moving on to the notebooks using real data.*

Modelling real VLTI data
========================

These tutorials were initially developed for the `12th VLTI School of Interferometry <https://vltischool2024.sciencesconf.org/>`_.  
Using real datasets from VLTI instruments (PIONIER, MATISSE, and GRAVITY), they demonstrate how to use **oimodeler** to perform model fitting,  
including chromaticity and kinematics, for stars and circumstellar environments.

The tutorials are provided as notebooks, containing Python code, explanations, and questions related to model-fitting issues.

.. list-table:: 
   :class: borderless
   :widths: 1 1 

   * - .. figure:: _static/notebooks/notebook_Ex1_thumbnail.png
            :target: https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex1_canopus.ipynb
            :alt: PIONIER observation of Canopus

            `PIONIER observation of Canopus 
            <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex1_canopus.ipynb>`_ 

            In this first example, we demonstrate how to use oimodeler to create a simple model and perform model fitting.  
            We introduce the concepts of reduced chi-square ( :math:`\chi^2_r` ) minimization and of local and global minima.

            We use VLTI/PIONIER data obtained on the giant star Canopus.

     - .. figure:: _static/notebooks/notebook_Ex2_thumbnail.png
              :target: https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex2_94_Aqr.ipynb
              :alt: The Binary star 94 Aqr MATISSE observation     

              `The binary star 94 Aqr MATISSE observation <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex2_94_Aqr.ipynb>`_

              In this exercise, we use VLTI/MATISSE data obtained on the binary star 94 Aqr.  
              We demonstrate issues around binary star fitting including local minima and degeneracy.

              We use MCMC and grid fitters.

   * - .. figure:: _static/notebooks/notebook_Ex3_thumbnail.png
              :target: https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex3_HD179278.ipynb
              :alt: Chromatic model for the YSO HD 179218

              `Chromatic model for the YSO HD 179218 <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex3_HD179278.ipynb>`_

              With the advent of routine spectro-interferometry, simple achromatic models are often insufficient to reproduce observations.  
              **oimodeler** provides tools to build chromatic models, including parameter interpolators.

              This example uses L-band MATISSE data on the Herbig star HD 179218 to illustrate the concept.

     - .. figure:: _static/notebooks/notebook_Ex4_thumbnail.png
              :target: https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex4_MWC297.ipynb
              :alt: MWC297: continuum vs line emission      

              `MWC297: continuum vs line emission <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex4_MWC297.ipynb>`_

              This notebook illustrates the use of a chromatic interpolator with high-resolution MATISSE data of the massive YSO MWC 297.  
              Focus is on the gaseous Br:math:`_\gamma` emission at 4.05 :math:`\mu`m and its adjacent dust-produced continuum.

              The goal is to determine the sizes of these two components, informing about the physical processes producing the gas emission.

   * - .. figure:: _static/notebooks/notebook_Ex5_thumbnail.png
              :target: https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex5_HD58647.ipynb
              :alt: Kinematics of the gaseous disk around HD58647

              `Kinematics of the gaseous disk around HD58647 <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex5_HD58647.ipynb>`_

              Spectro-interferometry in atomic lines can constrain dynamics via the Doppler effect.  
              Here, we focus on rotating disks around stars, found around YSOs, some evolved stars, and classical Be stars.

              This section introduces the **oimKinematicDisk** component implemented in **oimodeler**.

     - .. figure:: _static/notebooks/notebook_Ex6_thumbnail.png
              :target: https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex6_FSCMa.ipynb
              :alt: Radiative Transfer model of FS CMa
              
              `Radiative Transfer model of FS CMa <https://github.com/oimodeler/oimodeler/tree/main/examples/notebooks/oimodeler_Ex6_FSCMa.ipynb>`_

              **oimodeler** includes a class, **oimComponentFitsImage**, to import precomputed images or hyperspectral image cubes  
              (e.g., outputs from radiative-transfer models in FITS format). The loaded images can be shifted, rotated, and scaled.

              This example loads an image cube created with the radiative transfer code RADMC3D and compares it to MATISSE observations  
              of the B[e] star FS CMa.
