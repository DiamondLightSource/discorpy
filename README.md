# Vounwarp

## Distortion correction package


**Vounwarp** is an open-source Python package for calculating distortion coefficients
from a grid pattern image. It is the implementation of published methods, Nghia T. Vo et al.
"Radial lens distortion correction with sub-pixel accuracy for X-ray micro-tomography"
Optics Express 23, 32859-32868 (2015). https://doi.org/10.1364/OE.23.032859

Features
========
- Pre-processing methods for finding coordinates of dot-centroids, grouping them
 into lines, removing non-dot objects or misplaced dots.
- Procesing methods for calculating distortion coefficients, which are the center of distortion
  and polynomial coefficients, of the backward model, the forward model, and the
  backward-from-forward model.
- Post-processing methods for: unwarping lines of points, images, or slices of a 3D dataset; evaluating the accuracy of the correction model.
- Some methods may be useful for other applications:
  * Correct non-uniform background using a FFT-based filter and a median filter.
  * Select binary objects in a certain range of values.
  * Unwarp slices of a 3D dataset. 

Install
=======
- vounwarp is available on the conda-forge channel. To install:

  $ conda install -c conda-forge vounwarp
