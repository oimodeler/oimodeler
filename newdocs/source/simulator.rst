..  _simulator:

Data/Model comparison
=====================

In **oimodeler** the :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` class is the main class to do data/model
comparison. In these section we will present:

- the basics use and functionalities of this class
- some details on how the simluated interferometric data are computed
- some details on the :math:`\chi^2_r` computation

The code for this section is in
`SimulatingData.py <https://github.com/oimodeler/oimodeler/tree/main/examples/Modules/SimulatingData.py>`_


The basics of the simulator
---------------------------

Let's start by importing the needed modules and setting the variable ``files`` to the
list of the same OIFITS files as in the :ref:`exampleOimData` example.

.. code-block:: ipython3

    from pathlib import Path
    from pprint import pprint

    import oimodeler as oim

    path = Path(__file__).parent.parent.parent
    data_dir = path / "examples" / "data" / "ASPRO_MATISSE2"
    save_dir = path / "images"
    if not save_dir.exists():
        save_dir.mkdir(parents=True)

    files = list(map(str, data_dir.glob("*.fits")))


These OIFITS were simulated with ASPRO as a MATISSE observation of a partly resolved
binary star.

We make a model of a binary star where one of the components is resolved.
It consists of two components: A uniform disk and a point source.

.. code-block:: python

    ud = oim.oimUD(d=3, f=1, x=10, y=20)
    pt = oim.oimPt(f=0.5)
    model = oim.oimModel([ud, pt])


We now create an
:func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` object and feed it
with the data and our model.

The data can either be:

- A previously created :func:`oimData <oimodeler.oimData.oimData>`.
- A list of previously opened `astropy.io.fits.hdulist <https://docs.astropy.org/en/stable/io/fits/api/hdulists.html#astropy.io.fits.HDUList>`_.
- A list of paths to the OIFITS files (list of strings).

.. code-block:: python

    sim = oim.oimSimulator(data=files, model=model)


When creating the simulator, it automatically calls the :func:`oimData.prepareData <oimodeler.oimData.oimData.prepareData>`
method of the created :func:`oimData <oimodeler.oimData.oimData>` instance within
the simulator instance. This calls the
:func:`oimData.prepare <oimodeler.oimData.oimData.prepare>` method of :func:`oimData <oimodeler.oimData.oimData>`
The function is used to create vectorized coordinates for the data (spatial frequencies
in x and y and wavelengths) to be passed to the :func:`oimModel <oimodeler.oimModel.oimModel>`
instance to compute the Complex Coherent Flux (CCF) using the
:func:`oimModel.getComplexCoherentFlux <oimodeler.oimModel.oimModel.getComplexCoherentFlux>`
method, and some structures to go back from the ``ccf`` to the measured interferometric
quantities contained in the OIFITS files: VIS2DATA, VISAMP, VISPHI, T3AMP, T3PHI,
and FLUXDATA.

Once the data is prepared, we can call the :func:`oimSimulator.compute <oimodeler.oimSimulator.oimSimulator.compute>`
method to compute the :math:`\chi^2` and the simulated data.

.. code-block:: python

    sim.compute(computeChi2=True, computeSimulatedData=True)
    pprint("Chi2r = {}".format(sim.chi2r))


.. code-block::

    ... Chi2r = 11356.162973124885


Our model isn't fitting the data well. Let's take a closer look and
plot the data-model comparison for all interferometric quantities contained
in the OIFITS files.

.. code-block:: python

    fig0, ax0= sim.plot(["VIS2DATA", "VISAMP", "VISPHI", "T3AMP", "T3PHI"])


.. image:: ../../images/ExampleOimSimulator_model0.png
  :alt: Alternative text

.. warning::

    By default the simulator uses all data types to compute the chi2. In the case of our ASPRO simulated data, this is OK as all
    datatypes are computed. But for most real interferometric instruments, some data type should be ignore. It is often the case
    of the closure-ampltiude (T3AMP). For some instruments like MATISSE, one should choose between using VISAMP or VIS2DATA.

We can force the chi2r computation to only a subset of datatype using the dataTypes option of :func:`oimSimulator.compute
<oimodeler.oimSimulator.oimSimulator.compute>`  method. For instance, in the following we only compute the chi2r
on the square visibliity and closure-phase.

.. code-block::

    sim.compute(computeChi2=True, dataTypes=["VIS2DATA","T3PHI"])
    pprint(f"Chi2r = {sim.chi2r}")

.. code-block::

    ... Chi2r = 24393.17539459703

We could now try to fit the model "by hand", or by making a loop on some parameters
 and looking at the chi2r. But oimodeler implement various fitter class to perform automatic
 model fitting as shown in the next example
where we use a fitter from the :mod:`oimFitter <oimodeler.oimFitter>`
module to automatically find a good fit (and thus well fitting parameters).


Simulating data
---------------



Computing Chi2
--------------
