Methods for correcting distortions
==================================

Introduction
------------

For correcting radial and/or perspective distortion, one needs to know a model to
map between the distorted space/plane and the undistorted space/plane. Mapping
from the undistorted space to the distorted space is the forward mapping (Fig. 1).
The reverse process is the backward mapping or inverse mapping (Fig. 2).

.. figure:: figs/methods/fig1.jpg
  :figwidth: 95 %
  :align: center
  :figclass: align-center

  Figure 1. Forward mapping.

.. figure:: figs/methods/fig2.jpg
  :figwidth: 95 %
  :align: center
  :figclass: align-center

  Figure 2. Backward mapping.

There are many models which can be chosen from literature :cite:`Clarke_and_Fryer:1998,
Ricolfe-Viala:2010, Criminisi:1999` such as polynomial, logarithmic, field-of-view,
or matrix-based models to describe the relationship between the undistorted space
and distorted space. Some models were proposed for only one type of distortion while
others are for both distortion types including the location of the optical center.
From a selected model, one can find a practical approach to calculate the parameters
of this model.

To calculate parameters of a distortion model, one needs to be able to determine
the coordinates of reference points in the distorted space and their positions
in the undistorted space, correspondingly. Reference points can be extracted
using an acquired image of a `calibration object <https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=7501>`_
giving a line or dot-pattern image (Figs. 1,2), which is distorted. Using conditions
that lines of these points must be straight, equidistant, parallel, or perpendicular
one can estimate the locations of these reference-points in the undistorted
space with high-accuracy.

*Discorpy* is the Python implementation of radial distortion correction methods
presented in :cite:`Vo:2015`. These methods employ polynomial models and use a
calibration image for calculating coefficients of the models where the optical
center is determined independently. The reason of using these models and a
calibration image is to achieve sub-pixel accuracy as strictly required by
parallel-beam tomography systems. The methods were developed and used internally
at the beamline I12, Diamond Light Source-UK, as Mathematica codes. In 2018, they
were converted to Python codes and packaged as open-source software :cite:`Vo:2018`
under the name Vounwarp. The name was changed to Discorpy in 2021. From version 1.4,
methods for correcting perspective distortion :cite:`Criminisi:1999` and extracting
reference points from a line-pattern image were added to the software. A key
feature of methods developed and implemented in Discorpy is that radial distortion,
center-of-distortion (optical center), and perspective distortion are determined
independently using a single calibration image. The following sections explain
methods implemented in Discorpy

Extracting reference-points from a calibration image
----------------------------------------------------

The purpose of a calibration-image (Fig. 3(a,b,c)) is to provide reference-points
(Fig. 3(d)) which can be extracted from the image using some image processing
techniques. As shown in Fig. 3, there are a few calibration-images can be used
in practice. A dot-pattern image (Fig. 3(a)) is the easiest one to process because one just
needs to segment the dots and calculate the center-of-mass of each dot. For a
line-pattern image (Fig. 3(b)), a line-detection technique is needed. Points on the detected
lines or the crossing points between these lines can be used as the reference-points.
For a chessboard image (Fig. 3(c)), one can employ some corner-detection techniques or
apply a gradient filter to the image and use a line-detection technique.

In practice, acquired calibration images do not always look nice as shown in Fig. 3.
Some are very challenging to get reference-points. The following sub-sections present
practical approaches to process calibration images in such cases:

.. toctree::

    methods/dot_pattern
    methods/line_pattern


.. figure:: figs/methods/fig3.jpg
  :figwidth: 95 %
  :align: center
  :figclass: align-center

  Figure 3. (a) Dot-pattern image. (b) Line-pattern image. (c) Chessboard image. (d)
  Extracted reference-points from the image (a),(b), and (c).

Grouping reference-points into horizontal lines and vertical lines
------------------------------------------------------------------

Different techniques of calculating parameters of a distortion-model use
reference-points differently. The techniques :cite:`Bailey:2002, Vo:2015`
implemented in Discorpy group reference-points into horizontal lines and
vertical lines, represent them by the coefficients of parabolic fits, and
use these coefficients for calculating distortion-parameters.

.. figure:: figs/methods/fig4.png
  :figwidth: 95 %
  :align: center
  :figclass: align-center

  Figure 4. (a) Points are grouped into horizontal lines. (b) Points are grouped
  into vertical lines.
