:tocdepth: 1

Installation
============

The ``oimodeler`` package can be installed directly using the pip install command:

.. code-block:: console

    $ pip install git+https://github.com/oimodeler/oimodeler


.. warning::

    The examples including  dedicated data won't be installed when using pip in this way.
    They are available on the `Github  <https://github.com/oimodeler/oimodeler/tree/main/examples/>`_
    archive for manual download.


Alternatively, you can install the complete ``oimodeler`` package (including examples
and data) by cloning the `Github repository <https://github.com/oimodeler/oimodeler>`_:

.. code-block:: console

    $ git clone https://github.com/oimodeler/oimodeler
    $ cd oimodeler/

and installing it. As a non-editable install:

.. code-block:: console

    $ pip install .


Or as an editable (development) install:

.. code-block:: console

    $ pip install -e .
    
    
Dependancies
------------

The ``oimodeler`` package requires the following packages:

- `numpy <https://numpy.org/>`_
- `scipy <https://scipy.org/>`_
- `matplotlib <https://matplotlib.org/>`_
- `astropy <https://www.astropy.org/>`_
- `astroquery <https://astroquery.readthedocs.io/en/latest/>`_
- `emcee <https://emcee.readthedocs.io/en/stable/>`_
- `corner <https://corner.readthedocs.io/en/latest/>`_
- `tqdm <https://tqdm.github.io/>`_
- `pyFFTW <https://pypi.org/project/pyFFTW/>`_ (optional)

And for development some optional packages might be useful:

- `pytest <https://docs.pytest.org/en/7.3.x/>`_
- `pytest-cov <https://pytest-cov.readthedocs.io/en/latest/index.html>`_
- `sphinx <https://www.sphinx-doc.org/>`_
- `sphinx-autobuild <https://github.com/executablebooks/sphinx-autobuild>`_
- `sphinx-autodoc-typehints <https://github.com/tox-dev/sphinx-autodoc-typehints>`_
- `sphinx_rtd_theme <https://sphinx-rtd-theme.readthedocs.io/en/stable/index.html>`_
- `numpydoc <https://numpydoc.readthedocs.io/en/latest/>`_

These packages (except the development dependencies) are automatically installed if missing
when installing ``oimodeler``.


Checking installation
---------------------

The check if ``oimodeler`` is properly installed you can run the :ref:`getting_started`
or other :ref:`examples` scripts.
