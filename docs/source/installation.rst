Installation
============


oimodeler can be installed directly using the pip install command.

.. code-block:: console

    $ pip install git+https://github.com/oimodeler/oimodeler


.. note::
    The examples including  dedicated data won't be instaled when using pip. They are available on the `Github  <https://github.com/oimodeler/oimodeler/tree/main/examples/>`_ archive.

Alternatively, you can install the complete oimodeler package including examples and data by cloning the git repository.

.. code-block:: console

    $ git clone https://github.com/oimodeler/oimodeler
    
    
Dependancies
------------

oimodeler python library depends only on the following packages:

- `numpy <https://numpy.org/>`_
- `scipy <https://scipy.org/>`_
- `matplotlib <https://matplotlib.org/>`_
- `astropy <https://www.astropy.org/>`_
- `astroquery <https://astroquery.readthedocs.io/en/latest/>`_
- `emcee <https://emcee.readthedocs.io/en/stable/>`_
- `corner <https://corner.readthedocs.io/en/latest/>`_
- `tqdm <https://tqdm.github.io/>`_
- `pyFFTW <https://pypi.org/project/pyFFTW/>`_ (optional)

    
These packages (except the optional ones) are automatically installed if missing when installing oimodeler with the pip command but manual installation are required when cloning the repository.


Checking installation
---------------------

The check that oimodeler is properly installed you can run :ref:`getting_started` or other :ref:`examples` scripts