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
    data_dir = path / "examples" / "data" / "ASPRO_MATISSE2"
    save_dir = path / "images"
    if not save_dir.exists():
        save_dir.mkdir(parents=True)

    files = list(map(str, data_dir.glob("*.fits")))
    
.. note:
    The data consist in three OIFITS simulated with ASPRO as a MATISSE observation of a partly resolved binary star for the VLTI SMALL, MEDIUM and LARGE VLTI configurations of Auxiliary Telescopes (ATs).

   
The ``oimodeler`` package comes with the :func:`oimAxes <oimodeler.oimPlot.oimAxes>`
class that is a subclass of the standard `matplotlib.pytplot.Axes <https://matplotlib.org/stable/api/axes_api.html>`_
class (the base class for all matplotlib plots). To use it, you simply need to
specify it as a projection (actually this calls the subclass) when creating
an axe or axes.

.. code-block:: python 

    fig1 = plt.figure()
    ax1 = plt.subplot(projection='oimAxes')

   
First, we can plot the simple uv coverage using the
:func:`oimAxes.uvplot <oimodeler.oimPlot.oimAxes.uvplot>` method by passing the
list of OIFITS files (filename or opened) or an instance of a :func:`oimData <oimodeler.oimData.oimData>`
class.

.. code-block:: python 

    ax1.uvplot(data)

    
.. image:: ../../images/ExampleOimPlot_uv.png
  :alt: Alternative text     

There are several colouring options for the uvplot:

- ``color="byBaseline"`` : color the plots by baseline name 
- ``color="byFile"`` : color the plots by oifits file name 
- ``color="byConfiguration"`` : color the plots by Array configuration ("A0-B2-C1-D0" ...)
- ``color="byArray"`` : color the plots by Array name (VLTI, CHARA...) 

The data can be plotted in unit of length (by default) or in unit of spatial frequency using the ``cunit`` keyword.

Here are a few examples of uv plots with different options. To use the :func:`oimAxes <oimodeler.oimPlot.oimAxes>` 
in a multiple plot created by `matplotlib.pyplot.subplots <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html>`_
method, one need to set the  projection  using the ``subplot_kw`` keyword in when creating the axes.

.. code-block:: python 
 
    fig2, ax2 = plt.subplots(2, 3, subplot_kw=dict(projection='oimAxes'), figsize=(18, 9))

    ax2[0,0].uvplot(data.data,color="byBaseline",marker=".",legendkwargs=dict(fontsize=8))
    ax2[0,1].uvplot(data.data,color="byFile",facecolor="w",legendkwargs=dict(fontsize=5.))
    ax2[0,2].uvplot(data.data,color="byConfiguration",colorTab=["r","g","b"])
    ax2[1,0].uvplot(data.data,color="byArrname",colorTab=["r","g","b"],marker="+")
    ax2[1,1].uvplot(data.data,label="custom label",unit="km")
    ax2[1,2].uvplot(data.data,unit="cycle/rad",cunit="micron",lw=2,color="byConfiguration")
    fig2.tight_layout()
    
.. image:: ../../images/ExampleOimPlot_uv2.png
  :alt: Alternative text     
  
If the data is plotted in unit of spatial frequency without a color code specified, a colormap based on the wavelength will be used and a scalorscale will be added to the figure.

.. code-block:: python 

    fig3 = plt.figure()
    ax3 = plt.subplot(projection='oimAxes')
    ax3.uvplot(data.data,unit="cycle/mas",cunit="micron",
               label="cmap on wavelength",lw=3,cmap="plasma")
 
.. image:: ../../images/ExampleOimPlot_uv3.png
  :alt: Alternative text  

We can use the :func:`oiplot <oimodeler.oimPlot.oimAxes.oiplot>` method of
the :func:`oimAxes <oimodeler.oimPlot.oimAxes>` to produce plots of the following quantities:

.. csv-table:: Quantities plottable with oiplot method
   :file: table_plotData.csv
   :header-rows: 1  
   :delim: ;
   :widths: auto
   


For instance, let's plot the square visibilities (and corresponding errors) as a
function of the spatial frequency with the wavelength (converted in microns)
as a colorscale. 

.. code-block:: python
   
    fig4 = plt.figure()
    ax4 = plt.subplot(projection='oimAxes')
    ax4.oiplot(data, "SPAFREQ", "VIS2DATA", xunit="cycle/mas", label="Data",
                    cname="EFF_WAVE",cunit="micron", errorbar=True)
    ax4.legend()

   
.. image:: ../../images/ExampleOimPlot_v2.png
  :alt: Alternative text     
  

As for ``uvplot``, the color code can alternatively set using the ``color`` keyword.
Here we plot the square visibility as the function of the wavelength while 
colouring it by interferometer configurations (i.e., the list of all
telescopes). Note that here,  we are passing parameters to the error plot function
using the ``kwargs_error`` keyword.

.. code-block:: python

    fig5 = plt.figure()
    ax5 = plt.subplot(projection='oimAxes')
    ax5.oiplot(data, "EFF_WAVE", "VIS2DATA", xunit="micron",color="byConfiguration",
               errorbar=True,kwargs_error={"alpha": 0.3})
    ax5.legend()


.. image:: ../../images/ExampleOimPlot_v2Wl.png
  :alt: Alternative text       


.. note::
    Special values of the color option are ``"byFile"``, ``"byConfiguration"``,
    ``"byArrname"``, or ``"byBaseline"``. Other values will be interpreted as a
    standard `matplotlib colorname <https://matplotlib.org/stable/gallery/color/named_colors.html>`_.
    When using one of these values, the corresponding labels are added to the plots.
    Using the :func:`oimAxes.legend <oimodeler.oimPlot.oimAxes.legend>` method
    will automatically add the proper names.

  
Let's create a figure with multiple oiplots. As for uvplot, the projection keyword
has to be set for all :func:`oimAxes <oimodeler.oimPlot.oimAxes>`
using the ``subplot_kw`` keyword in the
`matplotlib.pyplot.subplots <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html>`_
method.

.. code-block:: python

    fig6, ax6 = plt.subplots(2, 2, subplot_kw=dict(
    projection='oimAxes'), figsize=(8, 8))

    
    ax6[0, 0].oiplot(data, "SPAFREQ", "VIS2DATA", xunit="cycle/mas", label="Data",
                              cname="EFF_WAVE", cunit="micron", ls=":", errorbar=True)
    ax6[0, 0].legend()
    ax6[0, 0].set_yscale('log')

    ax6[0, 1].oiplot(data, "EFF_WAVE", "VIS2DATA", xunit="nm",color="byBaseline",
                     errorbar=True, kwargs_error={"alpha": 0.1})
    ax6[0, 1].legend(fontsize=6)


    ax6[1, 0].oiplot(data, "SPAFREQ", "T3PHI", xunit="cycle/rad", errorbar=True,
                     lw=2, ls=":", color="byFile")
    ax6[1, 0].legend(fontsize=4)
 

    ax6[1, 1].oiplot(data, "EFF_WAVE", "T3PHI", xunit="m",cname="LENGTH",
                     errorbar=True, kwargs_error={"alpha": 0.1})



.. image:: ../../images/ExampleOimPlot_multi.png
  :alt: Alternative text   
    

Let's have a look at another data set: two VLTI/AMBER observations of the classical Be star Alpha Col.
Observation were centered on the BrGamma Emission line.


.. code-block:: python
    
    data_path = path / "examples" / "data" / "AMBER_AlphaCol"
    files = [data_path / "ALPHACOL_2010-01-09T00_58.fits",
              data_path / "ALPHACOL_2010-01-20T10_36.fits"]
    data=oim.oimData(files)
    

We can plot VIS2DATA, VISPHI, T3PHI as a function of the wavelength throught the emission line.

.. code-block:: python

    fig7, ax7 = plt.subplots(3, 1, subplot_kw=dict(projection='oimAxes'), figsize=(12, 10))
    ax7[0].oiplot(data, "EFF_WAVE", "VIS2DATA", xunit="Angstrom",color="byBaseline")
    ax7[0].legend()
    ax7[1].oiplot(data, "EFF_WAVE", "VISPHI", xunit="Angstrom",color="byBaseline")
    ax7[1].legend()
    ax7[2].oiplot(data, "EFF_WAVE", "T3PHI", xunit="Angstrom",color="byBaseline")
    ax7[2].legend()

We clearly see some interesting signal in the emission line but it is hard to disantangle 
signal from each baseline.

.. image:: ../../images/ExampleOimPlot_AlphaCol0.png
  :alt: Alternative text   
    

We have included in **oimodeler** a new template to produce easily per-baseline plots: :func:`oimWlTemplatePlots <oimodeler.oimPlot.oimWlTemplatePlots>`. It derives from the matplotlib.figure. Figure class and can be used by specifiying 
``FigureClass = oim.oimWlTemplatePlot`` in the figure creation.

.. code-block:: python

    fig=plt.figure(FigureClass=oim.oimWlTemplatePlots, figsize=(12, 7))
    
First we need to define what we want to plot by passing the oimData (or list of oifits files) to the :func:`autoshape <oimodeler.oimPlot.oimWlTemplatePlots.autoshape>` method of the newly created figure. The function also require a shape with a list 
of what data types we want to include in the figure. For instance, to create a figure with the VIS2DATA on the first row and the VISPHI and T3PHI on the second one :

.. code-block:: python

    fig.autoShape(data.data, shape=[["VIS2DATA",None],["VISPHI","T3PHI"]])
    fig.set_xunit("micron")

Here we have also specified that we sant the plots x-axis to be in microns. We can now plot the data using the basic plot function from matplotlib or custom one. We can pass keyword to the plotting function using the  ``plotFunctionkwarg`` dictionary. Here we use the standard errorbar and plot functions of matplotlib.

.. code-block:: python

    fig.plot(data.data, plotFunction=plt.Axes.errorbar, 
             plotFunctionkwarg=dict(color="gray", alpha=0.3))
    fig.plot(data.data, plotFunctionkwarg=dict(color="tab:blue",lw=0.5))

Finally we can set the plot limits and legends. Legends text can include per baseline information such as $BASELINE$, $LENGTH$ or $PA$ which respectively return the baseline name, length and position angle.
    
.. code-block:: python

    fig.set_ylim(["VISPHI","T3PHI"],-25,25)
    fig.set_ylim(["VIS2DATA"],0,1.2)
    fig.set_xlim(2.16,2.172)
    fig.set_legends(0.5,0.1, "$BASELINE$ $LENGTH$m $PA$$^o$", ["VIS2DATA","VISPHI"],
                    fontsize=12, ha="center")
    fig.set_legends(0.5,0.1, "$BASELINE$", ["T3PHI"], fontsize=12, ha="center")


.. image:: ../../images/ExampleOimPlot_AlphaCol1.png
  :alt: Alternative text   
    

Note that the :func:`oimodel <oimodeler.oimModel.oimModel>`, :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` and  :func:`oimFitter <oimodeler.oimFitter.oimFitter>` classes also contain plotting methods of their own that are described in their respective section of this documentation. 


Filtering data
^^^^^^^^^^^^^^

Filters can be applied to the :func:`oimData <oimodeler.oimData.oimData>` class
using the :func:`oimDataFilter <oimodeler.oimDataFilter.oimDataFilter>` class.
This object is basically a stack of filters derived from the 
:func:`oimDataFilterComponent <oimodeler.oimDataFilter.oimDataFilterComponent>`
abstract class. The example presented here comes from the
`exampleOimDataFilter.py <https://github.com/oimodeler/oimodeler/blob/main/examples/BasicExamples/exampleOimDataFilter>`_
script.


.. csv-table:: Available filter components
   :file: table_dataFilter.csv
   :header-rows: 1  
   :delim: ;
   :widths: auto
   

Let's first include the required packages and create a list of the OIFITS files from **oimodeler** data directory. 

.. code-block:: python 

    from pathlib import Path

    import matplotlib.pyplot as plt
    import oimodeler as oim

    path = Path(__file__).parent.parent.parent
    data_dir = path / "examples" / "data" / "FSCMa_MATISSE"
    save_dir = path / "images"
    if not save_dir.exists():
        save_dir.mkdir(parents=True)

    files = list(map(str, data_dir.glob("*.fits")))

We create an :func:`oimData <oimodeler.oimData.oimData>` object which will contain
the OIFITS data. 

.. code-block:: python 
    
    data = oim.oimData(files)


Here we create a simple filter to cut the data to a specific wavelength range with
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

    filters = oim.oimDataFilter(f1)
    data.setFilter(filters)
    

Alternatively, as we have only one filter component in our stack we can directly assigned it 
to the :func:`setFilter <oimodeler.oimData.oimData.setFilter>` method.

.. code-block:: python 

    data.setFilter(f1)

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

    
Let's plot the unfiltered and filtered data using the :func:`oiplot <oimodeler.oimPlot.oimAxe.oiplot>`
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
