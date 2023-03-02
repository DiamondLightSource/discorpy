# Discorpy
(**Dis**)tortion (**Cor**)rection (**Py**)thon-package

*Camera calibration and distortion correction software for lens-based detector systems*

![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/DiamondLightSource/discorpy/discorpy_ga.yml) 
[![Downloads](https://static.pepy.tech/personalized-badge/discorpy?period=total&units=international_system&left_color=grey&right_color=blue&left_text=Pypi-downloads)](https://pepy.tech/project/discorpy) 
[![former_vounwarp_downloads](https://anaconda.org/conda-forge/vounwarp/badges/downloads.svg)](https://anaconda.org/conda-forge/vounwarp) 
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/discorpy/badges/downloads.svg)](https://anaconda.org/conda-forge/discorpy) 
[![Documentation Status](https://readthedocs.org/projects/discorpy/badge/?version=latest)](https://discorpy.readthedocs.io/en/latest/?badge=latest) 
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/discorpy/badges/platforms.svg)](https://anaconda.org/conda-forge/discorpy) 
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/DiamondLightSource/discorpy) 
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/discorpy/badges/license.svg)](https://anaconda.org/conda-forge/discorpy)
![Coverage](https://github.com/DiamondLightSource/discorpy/raw/master/docs/coverage_report/coverage.svg)


**Discorpy** is an open-source Python package for camera calibration and distortion 
correction with sub-pixel accuracy. It calculates parameters of correction models 
using a grid pattern image. The package mainly implements methods published in 
[Optics Express](https://doi.org/10.1364/OE.23.032859). It provides methods in 
a full pipeline of data processing. From version 1.4, perspective distortion 
correction was added to the package.

**Author and maintainer:** Nghia Vo, *NSLS-II, Brookhaven National Laboratory, US; Diamond Light Source, UK*

Major updates
=============
- 25/02/2021: the package name was changed from Vounwarp to Discorpy. The old-name 
package is still available at https://github.com/nghia-vo/vounwarp 
and installable using conda-forge channel: https://anaconda.org/conda-forge/vounwarp
- 21/11/2021: Version 1.4 was released with new features: perspective distortion 
correction, pre-processing methods for line-pattern images and chessboard images. 

Features
========
- Pre-processing methods for: extracting reference-points from a dot-pattern image,
line-pattern image, and chessboard (checkerboard) image; grouping these points line-by-line. 
- Processing methods for calculating the optical center, coefficients of polynomial 
models for correcting radial distortion, and parameters of a model for correcting 
perspective distortion.
- Post-processing methods for: unwarping lines of points, images, or slices of 
a 3D dataset; and evaluating the accuracy of the correction results.
- Some methods may be useful for other applications:
  * Correct non-uniform background of an image.
  * Select binary objects in a certain range of values.
  * Unwarp slices of a 3D dataset.
- Summarized by an AI chatbot: "It is a Python library for camera distortion 
  correction that is designed to be easy to use and accessible to both computer 
  vision experts and novice users. The library provides a simple API for 
  correcting lens distortion in images, and it is based on the popular NumPy 
  library for numerical computing. Discorpy also has good documentation and 
  examples to help users get started with the library".


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

- https://discorpy.readthedocs.io/en/latest/usage.html#demonstrations 

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
a calibration target are shown [here](https://github.com/DiamondLightSource/discorpy/blob/master/examples/Perseverance_distortion_correction/Distortion_correction_for_Perseverance_camera.md).  

    ![Mars_rover](https://github.com/DiamondLightSource/discorpy/raw/master/data/demo/Mars_Rover_camera.jpg)
