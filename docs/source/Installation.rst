Installation
============

Using pip
---------

FULMAR can be installed conveniently using ``pip``:

.. code-block:: bash

    pip install fulmar-astro

If you have multiple versions of Python and pip on your machine, make sure to use ``pip3``. Try:

.. code-block:: bash

    pip3 install fulmar-astro

Requirements
------------

FULMAR has the following requirements:

* `Arviz <https://arviz-devs.github.io/arviz/>`_
* `Astropy <https://www.astropy.org/>`_ >4.1, <5
* `celerite2 <https://celerite2.readthedocs.io/en/latest/>`_
* `corner <https://github.com/dfm/corner.py>`_
* `exoplanet <https://docs.exoplanet.codes/en/latest/>`_
* `Lightkurve <https://docs.lightkurve.org/>`_ >= 2
* `Matplotlib <https://matplotlib.org/>`_
* `NumPy <https://www.numpy.org/>`_
* `pymc3-ext <https://github.com/exoplanet-dev/pymc3-ext>`_
* `TransitLeastSquares <https://github.com/hippke/tls>`_ >= 1.0.31

These should be installed automatically if you use ``pip`` to install FULMAR.

From Source
-----------
The source code from FULMAR can be pulled from `github <https://github.com/astrojose9/fulmar>` :

.. code-block:: bash

    git clone https://github.com/astrojose9/Fulmar.git
    cd fulmar
    python setup.py install

If the command ``python`` does not point to Python 3 on your machine, you can try to replace the last line with ``python3 setup.py install``. If you don't have ``git`` on your machine, you can find installation instructions `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.

