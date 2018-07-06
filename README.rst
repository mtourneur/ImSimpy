=======
ImSimpy
=======


.. image:: https://img.shields.io/pypi/v/ImSimpy.svg
        :target: https://pypi.python.org/pypi/ImSimpy

.. image:: https://img.shields.io/travis/dcorre/ImSimpy.svg
        :target: https://travis-ci.org/dcorre/ImSimpy

.. image:: https://readthedocs.org/projects/ImSimpy/badge/?version=latest
        :target: https://ImSimpy.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Image Simulator for optical/NIR telescope


* Free software: MIT license
* Documentation: https://ImSimpy.readthedocs.io.


Features
--------

ImSimpy is an Image Simulator for optical/NIR telescope developed in Python3.
The telescope characteristics are given as input using a hjson file. The telescope performance is first computed with `pyETC`_.
An other hjson file contains the condition of observations: seeing, filter band, exposure time, sky rightness, target elevation for instance but also parameters to compute or load the PSF, to download or load the catalogue, and to control the effects to include in the simulated image (Readout Noise, vignetting,hot and dead pixels, cosmic rays,...).
The simulated image is aved in a FITS file.

.. _pyETC: https://github.com/dcorre/pyETC

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
