# Discorpy
(**Dis**)tortion (**Cor**)rection (**Py**)thon-package

*Distortion calibration and correction software for lens-based detector systems*

**Discorpy** is an open-source Python package for distortion correction 
with sub-pixel accuracy as required by tomography detector systems. 
It calculates parameters of correction models using a grid pattern image. 
The package implements methods published by Nghia T. Vo et al.
"Radial lens distortion correction with sub-pixel accuracy for X-ray micro-
tomography," Optics Express, 23, 32859-32868 (2015), https://doi.org/10.1364/OE.23.032859. 
From version 1.4, correction for perspective distortion has been added to the package.

**Author:** Nghia Vo, *Diamond Light Source, UK.*

Update notes
============
- 25/02/2021: the package name is changed from Vounwarp to Discorpy. The old-name 
package is still available at https://github.com/nghia-vo/vounwarp 
and installable using conda-forge channel: https://anaconda.org/conda-forge/vounwarp
- 21/11/2021: correction for perspective distortion, pre-processing methods for
line-pattern images and chessboard images was added to the package. Version 1.4 
was released.

Features
========
- Pre-processing methods for: extracting reference-points from a dot-pattern image,
line-pattern image, and chessboard image; grouping these points line-by-line. 
- Processing methods for calculating the optical center, coefficients of polynomial 
models for correcting radial distortion, and parameters of a model for correcting 
perspective distortion.
- Post-processing methods for: undistorting lines of points, images, or slices of 
a 3D dataset; and evaluating the accuracy of the correction results.
- Some methods may be useful for other applications:
  * Correct non-uniform background using a FFT-based filter and a median filter.
  * Select binary objects in a certain range of values.
  * Unwarp slices of a 3D dataset.

Documentation
=============

- https://discorpy.readthedocs.io/en/latest/

Installation
============
- https://discorpy.readthedocs.io/en/latest/install.html

How to use
==========
- https://discorpy.readthedocs.io/en/latest/usage.html

Demonstrations
==============
- Apply to a visible dot-target collected at [Beamline I12](https://www.diamond.ac.uk/Instruments/Imaging-and-Microscopy/I12/Detectors-at-I12.html),
Diamond Light Source, UK:

    ![I12_before_after1](https://github.com/DiamondLightSource/discorpy/raw/master/data/demo/i12_data_1.jpg)

    ![I12_before_after2](https://github.com/DiamondLightSource/discorpy/raw/master/data/demo/i12_data_2.jpg)

- Apply to an X-ray dot-target collected at [Beamline I13](https://www.diamond.ac.uk/Instruments/Imaging-and-Microscopy/I13/Diamond-Manchester_Imaging_Branchline/Facilities_and_equipment_Imaging.html),
Diamond Light Source, UK:

    ![I13_before_after1](https://github.com/DiamondLightSource/discorpy/raw/master/data/demo/i13_data_1.jpg)

    ![I13_before_after2](https://github.com/DiamondLightSource/discorpy/raw/master/data/demo/i13_data_2.jpg)

- Improvement of a tomographic reconstructed image after distortion correction.
  + Before the correction:
    
    ![tomo_before](https://github.com/DiamondLightSource/discorpy/raw/master/data/demo/recon_before.jpg)
    
  + After the correction:
    
    ![tomo_before](https://github.com/DiamondLightSource/discorpy/raw/master/data/demo/recon_after.jpg)


- Apply to a hazard camera of the [Mars Perseverance Rover](https://mars.nasa.gov/mars2020/multimedia/raw-images/).
Details of how to estimate distortion coefficients of that camera without using
a calibration target are shown [here](https://github.com/DiamondLightSource/discorpy/blob/master/examples/Perseverance_distortion_correction/Distortion_correction_for_Perseverance_camera.md)  

    ![Mars_rover](https://github.com/DiamondLightSource/discorpy/raw/master/data/demo/Mars_Rover_camera.jpg)
