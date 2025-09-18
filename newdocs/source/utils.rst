..  _utils:

Oifits Helpers & Utils
======================

Many functions and classes developed for **oimodeler** to filter data, perform fitting, or produce plots can 
be accessed directly by the user. Most of these helper functions are implemented in the 
:func:`oimUtils <oimodeler.oimUtils>` module. A complete description is available in the **API** section.

Here, we only describe some of the features that might be of interest to regular **oifits** and/or **oimodeler** users.

Getting information on an oifits file
--------------------------------------

The functions listed below where developped to retrieve some information from the oifits tables:

- :func:`getBaselineName <oimodeler.oimUtils.getBaselineName>` : list of all baselines for an oifits array 
- :func:`getConfigName <oimodeler.oimUtils.getConfigName>`  : array configuration name (i.e., ATs position, UTs, CHARA telescopes) 
- :func:`getBaselineLengthAndPA <oimodeler.oimUtils.getBaselineLengthAndPA>` :  length (m) and PA (deg) of the baselines
- :func:`getWlFromOifits <oimodeler.oimUtils.getWlFromOifits>` the wavelength vector (m)
- :func:`get2DSpaFreq <oimodeler.oimUtils.get2DSpaFreq>` :   spatial frequencies of the baselines (with a format nB,n:math:`\lambda` ) 

ll of these functions operate similarly and return information for a single array (OI_VIS, OI_VIS2, or OI_T3) 
from one oifits file.

When working with an :func:`oimData <oimodeler.oimData.oimData>` object, the user needs to access
a single oifits file contained in the ``data`` variable.

Here are some examples of using these function on the first oifits file contained in an 
:func:`oimData <oimodeler.oimData.oimData>` instance.

.. code-block:: ipython3

    Bname       = oim.getBaselineName(data.data[0], length=True, angle=True)
    CPname      =  oim.getBaselineName(data.data[0],hduname="OI_T3")
    confname    = oim.getConfigName(data.data[0])
    B,PA        = oim.getBaselineLengthAndPA(data.data[0])
    u,v       = oim.get2DSpaFreq(data.data[0])
    wl          = oim.getWlFromOifits(data.data[0])


    print(f"This oifits files contains data taken with the {confname} array")
    print(f"The wavelength range is {wl.min()*1e6:.2f}-{wl.max()*1e6:.2f}\mum")

    for i in range(len(Bname)):
        print(f"{Bname[i]}")
    

.. parsed-literal::

    This oifits files contains data taken with the A0-B2-C1-D0 array
    The wavelength range is 2.69-4.20\mum
    D0-C1 23m -149$^o$
    A0-B2 24m 147$^o$
    B2-D0 34m 31$^o$
    B2-C1 11m 31$^o$
    A0-D0 32m 74$^o$
    A0-C1 22m 119$^o$

Note that, by default the functions return information on the first OI_VIS2 extension. Other extension 
can be accessed 
- using the ``hduname`` option to access extensions with other name such as  OI_T3, OI_VIS
- using the ``extver`` to specify a specific extension by it EXTVER number

The :func:`getBaselineName <oimodeler.oimUtils.getBaselineName>` can return the baselines length and orientation
 as seen in the above example.


Modifying oifits arrays
-----------------------

The :func:`oimUtils <oimodeler.oimUtils>` module also contains functions which are at the core of the data
filtering available using the :func:`oimDataFilter <oimodeler.oimDataFilter.oimDataFilter>` class.

Instead of using these feature through the :func:`oimDataFilter <oimodeler.oimDataFilter.oimDataFilter>` 
interface one can directly modify oifits files opened with astrpy.io.fits module. 

The main functions are the following :

- :func:`shiftWavelength <oimodeler.oimUtils.shiftWavelength>`: Shift the wavelength of an oifits file
- :func:`spectralSmoothing <oimodeler.oimUtils.spectralSmoothing>`: Smooth the spectral data of an oifits file
- :func:`binWavelength <oimodeler.oimUtils.binWavelength>`: Bin the wavelength of an oifits file
- :func:`oifitsFlagWithExpression <oimodeler.oimUtils.oifitsFlagWithExpression>` : Flag the data with an expression
- :func:`computeDifferentialError <oimodeler.oimUtils.computeDifferentialError>` : Compute the differential error
- :func:`setMinimumError <oimodeler.oimUtils.setMinimumError>` : Set the minimum error of a given data type to a given value


Creating oifits arrays
----------------------

the :func:`oimUtils <oimodeler.oimUtils>` module also contains helpers function to create fits tables compatibles 
with the  OIFITS2 (Optical Interferometry FITS) standard defined  in `Duvert et al. (2017) <https://www.aanda.org/articles/aa/pdf/2017/01/aa26405-15.pdf>`_. 

- :func:`createOiTarget <oimodeler.oimUtils.createOiTarget>`
- :func:`createOiArray <oimodeler.oimUtils.createOiArray>`
- :func:`createOiWavelength <oimodeler.oimUtils.createOiWavelength>`
- :func:`createOiVis <oimodeler.oimUtils.createOiVis>`
- :func:`createOiVis2 <oimodeler.oimUtils.createOiVis2>`
- :func:`createOiT3 <oimodeler.oimUtils.createOiT3>`
- :func:`createOiFlux <oimodeler.oimUtils.createOiFlux>`

For instance, to create a OI_VIS2 table, ones should provide the complete 
list required KEYWORDS and COLUMNS as defined in the OFITS2 standard:


.. code-block:: ipython3

    vis2 = oim.createOiVis2(OI_REVN=OI_REVN,
                            DATE-OBS=DATE-OBS,
                            ARRNAME=ARRNAME,
                            INSNAME=INSNAME,
                            TARGET_ID=TARGET_ID,
                            TIME=TIME,
                            MJD=MJD,
                            INT_TIME=INT_TIME,
                            VIS2DATA=VIS2DATA,
                            VIS2ERR=VIS2ERR,
                            UCOORD=UCOORD,
                            VCOORD=VCOORD,
                            STA_INDEX=STA_INDEX,
                            FLAG=FLAG)

Note that the data provided as numpy array should also have the homogeneous dimension in term of 
number of baselines and number of wavelengths.


A particular case is the :func:`createOiTargetFromSimbad <oimodeler.oimUtils.createOiTargetFromSimbad>` function,
which can generate an ``OI_TARGET`` table filled with values retrieved from the SIMBAD database when given 
a target name.

.. code-block:: ipython3

    target = oim.createOiTargetFromSimbad("Gamma Cas")
    print(target.data[0])

.. parsed-literal::

    (1, 'Gamma Cas', 14.1772125, 60.71674, 2000.0, 1.9444444276928152e-08, 
    1.9444444276928152e-08, 0.0, 'UNKNOWN', 'OPTICAL', 6.991666666666668e-06,
     -1.088888888888889e-06, 2.2222222284540294e-08, 2.2222222284540294e-08, 
     1.65e-06, 3.3333333e-08, 'B0.5IVpe')

Listing oimodeler features
--------------------------

the :func:`oimUtils <oimodeler.oimUtils>` module also contains helpers to list all **oimodeler** 
classes deriving from a parent class.

.. csv-table:: oimodeler class listing functions
   :file: table_listingFunctions.csv
   :header-rows: 1  
   :delim: |
   :widths: auto


These functions can be used to test which components are available in your version of **oimodeler**
