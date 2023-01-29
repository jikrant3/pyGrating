# pyGrating
``pyGrating`` is a python module to reduce Ultraviolet Imaging Telescope (UVIT/AstroSat) Grating spectra. An image reduced using [CCDLAB](https://github.com/user29A/CCDLAB) is required as input. The ``GratingImage`` module can identify all sources in the image and ``GratingSpectum`` module converts the 2D spectrum to 1D calibrated (wavelength and flux) spectrum.

## pgGrating documentation is at [jikrant3.github.io/pyGrating/](https://jikrant3.github.io/pyGrating)

## Installation
Add the ``pygrating.py`` and ``data`` folder to your working directory ([GitHub](https://github.com/jikrant3/pyGrating/blob/main/pygrating/pygrating.py)).


### Requirements
I recommend that you install the [Anaconda Distribution](https://www.anaconda.com/products/individual).

``pyGrating`` has following requirements:

* ``numpy`` >= 1.21.5
* ``scipy`` >= 1.7.3
* ``astropy`` >= 5.0.4
* ``specutils`` >= 1.7.0
* ``matplotlib`` >= 3.5.1

### Caution!
The software is under active development.

## Citation
Jadhav 2023, pyGrating, DOI: 10.5281/zenodo.7579951 [![DOI](https://zenodo.org/badge/594533919.svg)](https://zenodo.org/badge/latestdoi/594533919)

Tandon, et al. [2020, AJ, 159, 158](https://ui.adsabs.harvard.edu/abs/2020AJ....159..158T/abstract)

Dewangan, [2021, JApA, 42, 49](https://ui.adsabs.harvard.edu/abs/2021JApA...42...49D/abstract)

