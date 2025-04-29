..  _fitter:

Model fitting
=============

In the **oimodeler** framework, model-fitting is performed by classes deriving from the abstract class
:func:`oimFitter <oimodeler.oimFitter.oimFitter>`.


This class encapsulates an :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` instance, which determines
the :math:`\chi^2` value between the data and the model as seen in the :ref:`simulator` section.

In the current version, **oimodeler** includes four model-fitting algorithms.

.. csv-table:: Available Model-Fitters
   :file: table_fitters.csv
   :header-rows: 1
   :delim: |
   :widths: auto


These various algorithms allow the user to find the best-fit values of all free parameters of the model
(minimum of :math:`\chi^2`) and, depending on their nature, to evaluate their statistic. For instance, uncertainties can
be estimated using the posterior probability function in the case of MCMC or dynamic nested samplers,or using the
covariant matrix for the Minimize one.



Emcee fitter
------------

Chi2 Minimizer
--------------

Regular Grid exploration
------------------------

Dynesty fitter
--------------

