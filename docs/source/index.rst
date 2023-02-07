Discorpy's documentation
========================

Discorpy is an open-source Python package for correcting radial distortion
with sub-pixel accuracy as required by tomography detector systems :cite:`Vo:2015`.
It is used to calculate parameters of a correction model, which are the center
of distortion and the polynomial coefficients, using a grid pattern image. From
version 1.4, perspective distortion correction and methods for processing
line-pattern and chessboard (checkerboard) images were added to the package.

A key feature of Discorpy is that radial distortion, the center of distortion,
and perspective distortion are determined and corrected independently using a
single calibration image. Discorpy was developed for calibrating lens-coupled
detectors of tomography systems but it also can be used for commercial cameras.

**Showcases**: https://discorpy.readthedocs.io/en/latest/usage.html#demonstrations

.. image:: img/index/showcase.png
  :width: 100 %
  :align: center

**Source code:** https://github.com/DiamondLightSource/discorpy

**Author:** Nghia T. Vo - NSLS-II, Brookhaven National Laboratory, US; Diamond Light Source, UK.

**Keywords:** Camera calibration, radial lens distortion, perspective distortion,
distortion correction, tomography.

Contents
========


.. toctree::
   :maxdepth: 3
   :numbered:

   install
   tutorials
   usage
   api
   credit
