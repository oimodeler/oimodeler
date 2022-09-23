..  _examples:

Examples
========

All the following examples can be found in the examples subdirectories of the oimodeler github repository.


Basic Examples
--------------

In this section we presents script we presents showing the basic functionalities of the oimodeler software.


Loading Data
^^^^^^^^^^^^

The ``loadData.py`` script show how to create a oimData object from a list of oifits files.


.. code-block:: python

    import oimodeler as oim
    import os

    path = os.path.dirname(oim.__file__)
    pathData=os.path.join(path,os.pardir,"examples","testData","FSCMa_MATISSE")
    files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData)]

    data=oim.oimData(files)


The oifits data, stored in the ``astropy.io.fits.hdulist`` format, can be accessed using the ``oimData.data`` variable

.. code-block:: python

    print(data.data)
    
.. code-block:: python

    >>[[<astropy.io.fits.hdu.image.PrimaryHDU object at 0x000002657CBD7CA0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E546AF0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E3EA970>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E3EAAC0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E406520>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E402EE0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E406FD0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E4600D0>],
    [<astropy.io.fits.hdu.image.PrimaryHDU object at 0x000002657E458F70>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x0000026500769BE0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002650080EA60>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x00000265007EA430>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x00000265007EAAF0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002650080EC40>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E4DC820>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E4ECFD0>],
    [<astropy.io.fits.hdu.image.PrimaryHDU object at 0x000002657E4DCCA0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x0000026500B7EB50>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E9F79D0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E5913A0>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E591A60>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E591B20>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E5B7790>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x000002657E5BAEB0>]]
    
    
To learn more on the astropy.i OIFITS2

To be used in the oimSimulator and oiFitter data need to be optimized in a simpler vectorial/structure. Tis step is done automatically when using the simulator or fitter but can be done manually using the following command:
    
.. code-block:: python:
    
    data.prepareData()
    
For instance this create single vectors fgor the data coordinate : ``data.vect_u``, ``data.vect_v``, ``data.vect_wl``

.. code-block:: python:

    print(data.vect_u)
    print(data.vect_v)   
    print(data.vect_wl)  
    print(data.vect_u.shape)
    
.. code-block:: python:
    
    [0. 0. 0. ... 0. 0. 0.]
    [0. 0. 0. ... 0. 0. 0.]
    [4.20059359e-06 4.18150239e-06 4.16233070e-06 ... 2.75303296e-06
     2.72063039e-06 2.68776785e-06]
    (5376,)
    
    
Create a simple model
^^^^^^^^^^^^^^^^^^^^^




.. _createModelChromatic:

Create a chromatic model
^^^^^^^^^^^^^^^^^^^^^^^^



.. _createSimulator:

Create a Simulator
^^^^^^^^^^^^^^^^^^




Plotting data from oifits files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Expanding the Software
----------------------

In this section we present examples that show how to expand the functionalities of the oimodeler sofwate by crezating customs objects : oimComponents, oimFilterComponents, oimFitters, and custom plotting function or utils.

Performance Tests
-----------------

Scripts concerning performance tests are presented in this section.

Data for tests
--------------

