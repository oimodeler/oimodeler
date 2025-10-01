:tocdepth: 1

Installation
============

The **oimodeler** package can be installed directly using the pip install command:

.. code-block:: console

    $ pip install git+https://github.com/oimodeler/oimodeler


.. warning::

    The examples and their dedicated data are not installed via pip. They are available for manual 
    download on the `GitHub <https://github.com/oimodeler/oimodeler/tree/main/examples/>`_ archive.


Alternatively, we recommend installing the full **oimodeler** package (including examples and data) 
by cloning the `GitHub repository <https://github.com/oimodeler/oimodeler>`_:

.. code-block:: console

    $ git clone https://github.com/oimodeler/oimodeler
    $ cd oimodeler/

and installing it as a non-editable package:

.. code-block:: console

    $ pip install .

Or as an editable (development) install:

.. code-block:: console

    $ pip install -e .


Dependencies
------------

The **oimodeler** package requires the following dependencies:

- `numpy <https://numpy.org/>`_
- `scipy <https://scipy.org/>`_
- `matplotlib <https://matplotlib.org/>`_
- `astropy <https://www.astropy.org/>`_
- `astroquery <https://astroquery.readthedocs.io/en/latest/>`_
- `emcee <https://emcee.readthedocs.io/en/stable/>`_
- `corner <https://corner.readthedocs.io/en/latest/>`_
- `tqdm <https://tqdm.github.io/>`_
- `pyFFTW <https://pypi.org/project/pyFFTW/>`_ (optional)
- `dynesty <https://dynesty.readthedocs.io/>`_ (optional)

For development, the following optional packages might be useful:

- `pytest <https://docs.pytest.org/en/7.3.x/>`_
- `pytest-cov <https://pytest-cov.readthedocs.io/en/latest/index.html>`_
- `sphinx <https://www.sphinx-doc.org/>`_
- `sphinx-autobuild <https://github.com/executablebooks/sphinx-autobuild>`_
- `sphinx-autodoc-typehints <https://github.com/tox-dev/sphinx-autodoc-typehints>`_
- `sphinx_rtd_theme <https://sphinx-rtd-theme.readthedocs.io/en/stable/index.html>`_
- `numpydoc <https://numpydoc.readthedocs.io/en/latest/>`_

Except for development dependencies, these packages are automatically installed if missing when 
installing **oimodeler** via pip.


Checking Installation
---------------------

To check if **oimodeler** is properly installed, run the :ref:`getting_started` guide or try some 
examples from the :ref:`notebooks` section.

