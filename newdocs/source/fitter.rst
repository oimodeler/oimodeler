..  _fitter:

Model fitting
=============

In the **oimodeler** framework, model-fitting is performed by classes deriving from the abstract class
:func:`oimFitter <oimodeler.oimFitter.oimFitter>`. Here we describe the common functionalities and methods 
of these fitting classes. 

Theses classes encapsulate an :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` instance, which determines
the :math:`\chi^2` value between the data and the model as seen in the :ref:`simulator` section. As for the simulator 
class, :math:`\chi^2` can be computed only on a subset of datatypes setting the ``dataTypes`` kewyord. 

When creating a :func:`oimFitter <oimodeler.oimFitter.oimFitter>` the user should either pass an instance 
of :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` 

.. code:: ipython3

      fit = oimFitter(sim)

or an instance of :func:`oimData <oimodeler.oimData.oimData>` and  :func:`oimModel <oimodeler.oimModel.oimModel>`  

.. code:: ipython3

         fit = oimFitter(data, model)

The fitting classes have four common methods:

   - :func:`prepare <oimodeler.oimFitter.oimFitter.prepare>` : prepare the fitter: defining the parameter space...
   - :func:`run <oimodeler.oimFitter.oimFitter.run>` : launch the model-fitting run
   - :func:`getResults <oimodeler.oimFitter.oimFitter.getResults>` : return the best-fit parameters and uncertainties 
   - :func:`printResults <oimodeler.oimFitter.oimFitter.printResults>` : print the result of the model-fitting

.. warning::
   After calling the :func:`getResults <oimodeler.oimFitter.oimFitter.getResults>` or  
   :func:`printResults <oimodeler.oimFitter.oimFitter.printResults>`  the parameters of the 
   :func:`oimModel <oimodeler.oimModel.oimModel>` is set to their best-fit values.

The various classes of model-fitting also includes specific plotting functions that will be described 
in details in their recpective sub-sections.

In the current version, **oimodeler** includes four model-fitting algorithms.

.. csv-table:: Available Model-Fitters
   :file: table_fitters.csv
   :header-rows: 1
   :delim: |
   :widths: auto

These various algorithms allow the user to find the best-fit values of all free parameters of the model
(minimum of :math:`\chi^2`) and, depending on their nature, to evaluate their statistic. 

For instance, uncertainties can be estimated using :
   - the posterior probability function in the case of MCMC or DNS
   - the covariant matrix for the gradient-descent methods such as Minimize one

.. warning::
   It should be noted that no model fitting algorithm can guarantee convergence to the global minimum 
   of the chi-squared statistic.

Emcee fitter
------------

**A few words about the MCMC algorithm**

The Markov Chain Monte Carlo (MCMC) algorithm is a method used to sample from complex probability distributions
when direct sampling is difficult. It works by constructing a Markov chain, a sequence of random variables
where each variable depends only on the previous one.

Unlike optimization algorithms that seek the global minimum, MCMC methods do not directly aim to find it.
Instead, they converge toward a probability distribution, which may be concentrated near the global minimum
if that region has high probability. 

The burn-in phase refers to the initial iterations of the MCMC process, during which the chain has not yet reached 
the target distribution and the samples may not be representative
of the true posterior.

After a burning phase, the chain "wanders" over time through the sample space in such a way that the frequency 
of visits to each region reflects the target distribution, allowing for approximate estimation of expectations
and probabilities. 


**Description of the oimFitterEmcee class**

To implement a MCMC sampler, oimodeler use the python library **emcee**. 
If you are not confident with this package, you should have a look at the documentation 
`here <https://emcee.readthedocs.io/en/stable/>`_.

The **emcee** sampler is encapsulated into the :func:`oimFitterEmcee <oimodeler.oimFitter.oimFitterEmcee>` class.

At the creation of the fitter the number of desired walker exploring the paramter space can be specified using the 
keyword ``nwalkers``. The default number is 20.

.. code:: ipython3

   fit = oim.oimFitterEmcee(data,model, nwalkers=10, dataTypes=["VIS2DATA","T3PHI"])

The :func:`prepare <oimodeler.oimFitter.oimFitterEmcee.prepare>` method should then be called to set the initial 
walkers positions. Two options are available depending on the value of the keyword ``init``:
   - **random** : (default) Uniformly random positions within the parameter space limited by the values of :func:`oimParam.min <oimodeler.oimParam.oimParam.min>` and :func:`oimParam.max <oimodeler.oimParam.oimParam.max>`
   - **gaussian** : random positions from a normal (Gaussian) distribution around the current position defined by the :func:`oimParam.value <oimodeler.oimParam.oimParam.value>`  with standard deviation of :func:`oimParam.error <oimodeler.oimParam.oimParam.error>`

The :func:`oimFitterEmcee <oimodeler.oimFitter.oimFitterEmcee>` class offers two additionnal functionalities 
of the **emcee** package. 

The first one is the possiblity of changing the walker 
`moves <https://emcee.readthedocs.io/en/stable/user/moves/>`_ to optimize the parameter space exploration.

The user can also load and save the sampler using the `HDF5 backend <https://emcee.readthedocs.io/en/stable/user/backends/>`_.
This can be done by specifying an hdf5 file with the ``samplerFile`` keyword. If the file exists, it is 
loaded into the sampler. Its probability distribution can be explored and the best-fit parameters determined.
If the file doesn't exist, it will be created and the results of the MCMC run will be saved in this file.

.. code:: ipython3

   fit.prepare(init="gaussian", moves = moves.StretchMove, samplerFile=mySampler.h5)

After initializing the walkers, the MCMC run can be performed using the :func:`run <oimodeler.oimFitter.oimFitterEmcee.run>`
method. The number of iterations of the run is set by the ``nsteps`` keyword. The ``progress`` keyword can be used to show a 
progress bar.

.. code:: ipython3

   fit.run(nsteps=5000, progress=True)

After the MCMC run, the results can be plotted with three methods:



**Example on MIRCX data of a binary star**

To demonstrate the use of the :func:`oimFitterEmcee <oimodeler.oimFitter.oimFitterEmcee>` class we will use a single
`MIRCX <http://www.astro.ex.ac.uk/people/kraus/mircx.html>`_ observation of the binary star :math:`\beta` Ari. 
The code for this section is in `emceeFitting.py <https://github.com/oimodeler/oimodeler/tree/main/examples/Modules/emceeFitting.py>`_




Chi2 Minimizer
--------------

Regular Grid exploration
------------------------

Dynesty fitter
--------------

