..  _manipulatingOifitsData:

Manipulating oifits data
------------------------

In this section we introduce tools aiming at simplifying the manipulation of oifits data


Plotting oifits data
^^^^^^^^^^^^^^^^^^^^

Beyond the specific plots shown in the previous examples, the
:func:`oimPlot <oimodeler.oimPlot.oimPlot>` module allows to plot most of the
OIFITS data in a very simple way. The example presented here comes from the
`exampleOimPlot.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimPlot.py>`_
script.

Let's start by setting up the project with imports, path, and some data.

.. code-block:: python 

    from pathlib import Path

    import matplotlib.pyplot as plt
    import oimodeler as oim


    path = Path(__file__).parent.parent.parent
    data_dir = path / "examples" / "testData" / "ASPRO_MATISSE2"
    save_dir = path / "images"
    if not save_dir.exists():
        save_dir.mkdir(parents=True)

    files = list(map(str, data_dir.glob("*.fits")))
    
The ``oimodeler`` package comes with the :func:`oimAxes <oimodeler.oimPlot.oimAxes>`
class that is a subclass of the standard `matplotlib.pytplot.Axes <https://matplotlib.org/stable/api/axes_api.html>`_
class (the base class for all matplotlib plots). To use it, you simply need to
specify it as a projection (actually this calls the subclass) when creating
an axe or axes.

.. code-block:: python 

    fig1 = plt.figure()
    ax1 = plt.subplot(projection='oimAxes')

   
First, we can plot the classic uv coverage using the
:func:`oimAxes.uvplot <oimodeler.oimPlot.oimAxes.uvplot>` method by passing the
list of OIFITS files (filename or opened) or an instance of a :func:`oimData <oimodeler.oimData.oimData>`
class.

.. code-block:: python 

    ax1.uvplot(data)

    
.. image:: ../../images/ExampleOimPlot_uv.png
  :alt: Alternative text     

    
We can use the :func:`oiplot <oimodeler.oimPlot.oimAxes.oiplot>` method of
the :func:`oimAxes <oimodeler.oimPlot.oimAxes>` to plot any quantity inside
an OIFITS file as a function of another one. For instance, let's plot the
squared visibilities as a function of the spatial frequencies with the wavelength
(in microns) as a colorscale.

.. code-block:: python
   
    fig2 = plt.figure()
    ax2 = plt.subplot(projection='oimAxes')
    lamcol=ax2.oiplot(data, "SPAFREQ", "VIS2DATA", xunit="cycles/mas", label="Data",
                      cname="EFF_WAVE", cunitmultiplier=1e6, errorbar=True)
                
    plt.colorbar(lamcol, ax=ax2, label="$\\lambda$ ($\mu$m)")
    ax2.legend()

    
.. image:: ../../images/ExampleOimPlot_v2.png
  :alt: Alternative text     
  
  
We can also plot the square visibility as the function of the wavelength while 
colouring the curves by the interferometer configurations (i.e., the list of all
baselines within one file). Note that we can pass parameters to the error plots
with the ``kwargs_error`` keyword.

.. code-block:: python

    fig3 = plt.figure()
    ax3 = plt.subplot(projection='oimAxes')
    ax3.oiplot(data, "EFF_WAVE", "VIS2DATA", xunitmultiplier=1e6, color="byConfiguration",
               errorbar=True, kwargs_error={"alpha": 0.3})
    ax3.legend()

  
.. image:: ../../images/ExampleOimPlot_v2Wl.png
  :alt: Alternative text       


.. note::
    Special values of the color option are ``"byFile"``, ``"byConfiguration"``,
    ``"byArrname"``, or ``"byBaseline"``. Other values will be interpreted as a
    standard `matplotlib colorname <https://matplotlib.org/stable/gallery/color/named_colors.html>`_.
    When using one of these values, the corresponding labels are added to the plots.
    Using the :func:`oimAxes.legend <oimodeler.oimPlot.oimAxes.legend>` method
    will automatically add the proper names.

  
Finally, we can create a ``2x2`` figure with multiple plots. The projection keyword
has to be set for all :func:`oimAxes <oimodeler.oimPlot.oimAxes>`
using the ``subplot_kw`` keyword in the
`matplotlib.pyplot.subplots <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html>`_
method.

.. code-block:: python

    fig4, ax4 = plt.subplots(2, 2, subplot_kw=dict(
        projection='oimAxes'), figsize=(8, 8))

    ax4[0, 0].uvplot(data)

    lamcol = ax4[0, 1].oiplot(data, "SPAFREQ", "VIS2DATA", xunit="cycles/mas", label="Data",
                              cname="EFF_WAVE", cunitmultiplier=1e6, ls=":", errorbar=True)

    fig4.colorbar(lamcol, ax=ax4[0, 1], label="$\\lambda$ ($\mu$m)")
    ax4[0, 1].legend()
    ax4[1, 0].oiplot(data, "EFF_WAVE", "VIS2DATA", xunitmultiplier=1e6, color="byBaseline",
                     errorbar=True, kwargs_error={"alpha": 0.1})

    ax4[1, 0].legend(fontsize=6)
    ax4[1, 1].oiplot(data, "SPAFREQ", "T3PHI", xunit="cycles/mas", errorbar=True,
                     lw=2, ls=":", color="byFile")

    ax4[1, 1].legend(fontsize=4)
    ax4[0, 1].set_yscale('log')
    ax4[1, 0].autolim()
    ax4[1, 1].autolim()

    
.. image:: ../../images/ExampleOimPlot_multi.png
  :alt: Alternative text   
    


Filtering data
^^^^^^^^^^^^^^

Filtering can be applied to the :func:`oimData <oimodeler.oimData.oimData>` class
using the :func:`oimDataFilter <oimodeler.oimDataFilter.oimDataFilter>` class.
It is basically a stack of filters derived from the 
:func:`oimDataFilterComponent <oimodeler.oimDataFilter.oimDataFilterComponent>`
abstract class. The example presented here comes from the
`exampleOimDataFilter.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimDataFilter>`_
script.

As done before the required packages and create a list of the OIFITS files.  

.. code-block:: python 

    from pathlib import Path

    import matplotlib.pyplot as plt
    import oimodeler as oim

    path = Path(__file__).parent.parent.parent
    data_dir = path / "examples" / "testData" / "FSCMa_MATISSE"
    save_dir = path / "images"
    if not save_dir.exists():
        save_dir.mkdir(parents=True)

    files = list(map(str, data_dir.glob("*.fits")))

We create an :func:`oimData <oimodeler.oimData.oimData>` object which will contain
the OIFITS data. 

.. code-block:: python 
    
    data = oim.oimData(files)


We now create a simple filter to cut the data to a specific wavelength range with
the :func:`oimWavelengthRangeFilter <oimodeler.oimDataFilter.oimWavelengthRangeFilter>`
class. 

.. code-block:: python 
    
    f1 = oim.oimWavelengthRangeFilter(targets="all", wlRange=[3.0e-6, 4e-6])

    
The :func:`oimWavelengthRangeFilter <oimodeler.oimDataFilter.oimWavelengthRangeFilter>`
has two keywords:

- ``targets``: Which is common to all filter components: It specifies the targeted
  files within the data structure to which the filter applies.

  - Possible values are: ``"all"`` for all files (which we use in this example).
  - A single file specify by its index.
  - Or a list of indexes.

- ``wlRange``: The wavelength range to cut as a two elements list
  (min wavelength and max wavelength), or a list of multiple two-elements lists
  if you want to cut multiple wavelengths ranges simultaneously. In our example
  you have selected wavelength between 3 and 4 microns. Wavelengths outside this
  range will be removed from the data.
    
Now we can create a filter stack with this single filter and apply it to our data.

.. code-block:: python 

    filters = oim.oimDataFilter([f1])
    data.setFilter(filters)
    

By default the filter will be automatically activated as soon as a filter is set using
the :func:`oimData.setFilter <oimodeler.oimData.oimData.setFilter>` method
of the :func:`oimData <oimodeler.oimData.oimData>` class.
This means that querying the ``oimData.data`` attribute will return the filtered data,
and that when using the :func:`oimData <oimodeler.oimData.oimData>` class within
an :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` or an
:func:`oimFitter <oimodeler.oimFitter.oimFitter>`, the filtered data will be used
instead of the unfiltered data. 

.. note::

    The unfiltered data can always be accessed using the ``oimData._data`` attribute
    and, in a similar manner, also the filtered data (that may be ``None`` if no filters
    have been applied) using the private attribute ``oimData._filteredData``.

   
To switch off a filter we can either call the :func:`oimData.setFilter <oimodeler.oimData.oimData.setFilter>`
method without any arguments (this will remove the filter completely),

.. code-block:: python 

    data.setFilters()


or set the ``useFilter`` attribute to ``False``.

.. code-block:: python 

    data.useFilter = False

    
Let's plot the unfiltered and filtered data using the :func:`oimPlot <oimodeler.oimPlot.oimPlot>`
method.

.. code-block:: python 

    fig = plt.figure()
    ax = plt.subplot(projection='oimAxes')

    data.useFilter = False
    ax.oiplot(data, "SPAFREQ", "VIS2DATA", color="tab:blue", lw=3, alpha=0.2, label="unfiltered")

    data.useFilter = True
    ax.oiplot(data, "SPAFREQ", "VIS2DATA", color="tab:blue", label="filtered")

    ax.set_yscale('log')
    ax.legend()
    ax.autolim()
    

.. image:: ../../images/ExampleFilter_wavelengthCut.png
  :alt: Alternative text 

  
Other filters for data selection are:

- ``oimRemoveArrayFilter``: Removes array(s) (such as OI_VIS, OI_T3, etc.) from the data. 
- ``oimDataTypeFilter``: Removes data type(s) (such as VISAMP, VISPHI, T3AMP, etc.)
  from the data.

.. note::
    Actually, :func:`oimDataTypeFilter <oimodeler.oimDataFilter.oimDataTypeFilter>`
    doesn't remove the columns with the data type from
    any array as these columns are complusory in the the OIFITS format definition.
    Instead, it is setting all the values of the column to zero which ``oimodeler``
    will recognize as empty for data simulation and model fitting. 


.. code-block:: python 

    f2 = oim.oimRemoveArrayFilter(targets="all", arr=["OI_VIS", "OI_FLUX"])         
    f3 = oim.oimDataTypeFilter(targets="all", dataType=["T3AMP"," T3PHI"])
    data.setFilter(oim.oimDataFilter([f1, f2, f3]))


Here, we create a new filter stack with the previous wavelength filter `f1`,
a filter `f2` for removing the array OI_VIS and OI_FLUX from the data, and a filter
`f3` removing the columns T3AMP and T3PHI. Basically, we only have the VIS2DATA left
in our OIFITS structure.

.. note::
    Removing T3AMP and T3PHI from the OI_T3 is equivalent for model-fitting to removing
    the array OI_T3. 