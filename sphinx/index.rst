.. pyGrating documentation master file, created by
   sphinx-quickstart on Sat Jan 28 22:54:37 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyGrating's documentation!
=====================================

``pyGrating`` is a python module to reduce Ultraviolet Imaging Telescope (UVIT/AstroSat) Grating spectra. An image reduced using `CCDLAB`_ is required as input. The ``GratingImage`` module can identify all sources in the image and ``GratingSpectum`` module converts the 2D spectrum to 1D calibrated (wavelength and flux) spectrum.

Installation
------------
Add the `pygrating.py`_ file and the `data`_ folder to your working directory .


Requirements
~~~~~~~~~~~~
I recommend that you install the `Anaconda Distribution`_.

``pyGrating`` has following requirements:

* ``numpy`` >= 1.21.5
* ``scipy`` >= 1.7.3
* ``astropy`` >= 5.0.4
* ``specutils`` >= 1.7.0
* ``matplotlib`` >= 3.5.1

Caution!
~~~~~~~~
The software is under active development.

Citation
--------

.. image:: https://zenodo.org/badge/594533919.svg
   :target: https://zenodo.org/badge/latestdoi/594533919

Jadhav 2023, pyGrating, DOI: 10.5281/zenodo.7579951 [1]_

Tandon, et al. 2020, AJ, 159, 158 [2]_

Dewangan, 2021, JApA, 42, 49 [3]_


.. [1] https://zenodo.org/badge/latestdoi/594533919
.. [2] https://ui.adsabs.harvard.edu/abs/2020AJ....159..158T/abstract
.. [3] https://ui.adsabs.harvard.edu/abs/2021JApA...42...49D/abstract
.. _pygrating.py: https://github.com/jikrant3/pyGrating/blob/main/pygrating/pygrating.py
.. _data: https://github.com/jikrant3/pyGrating/tree/main/data
.. _CCDLAB: https://github.com/user29A/CCDLAB
.. _Anaconda Distribution: https://www.anaconda.com/products/individual


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Home <self>
   getting_started
   modules
   changlog


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
