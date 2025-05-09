.. _usage:

=====
Usage
=====

Resources
---------

-   A peer-reviewed paper describing the methods developed and implemented for Discorpy is available on
    `Journal of Synchrotron Radiation <https://doi.org/10.1107/S1600577525002267>`__.
-   Examples of how to use the package are in the "/examples" folder at the `github page <https://github.com/DiamondLightSource/discorpy/tree/master/examples>`__.
-   Coefficients determined by the package can be used by other tomographic
    software such as `Tomopy <https://tomopy.readthedocs.io/en/latest/api/tomopy.prep.alignment.html>`__,
    `Savu <https://github.com/DiamondLightSource/Savu/blob/master/savu/plugins/corrections/distortion_correction.py>`__,
    or `Algotom <https://github.com/algotom/algotom>`__ for correction.

Notes related to Python programming
-----------------------------------

In the :ref:`workflow <methods>` of calibrating a distortion target, some functions
work with Python lists, some with Numpy arrays. They are different objects and
can cause confusion or give rise to errors if users don't use these objects
properly. There are many tutorials online to explain the difference between them.
Here we only show some examples related to how they may be used with Discorpy.

-   Initialization:

    .. code-block:: python

        import numpy as np

        # Declare a Python list
        points_n1 = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        # Declare a Numpy array
        points_n2 = np.asarray(points_n1, dtype=np.float32)
        # or
        points_n2 = np.zeros((3,2), dtype=np.float32)

-   Indexing:

    .. code-block:: python

        # For a Python list
        element1 = points_n1[0][1]
        # For a Numpy array
        element2 = points_n2[0,1]

        # For a Python list
        sub_set1 = points_n1[0][0:2]
        # For a Numpy array
        sub_set2 = points_n2[0,0:2]

-   Math operations:

    .. code-block:: python

        # For a Python list
        points_m1 = points_n1 + 1.0 # >> Raise an error
        # For a Numpy array
        points_m2 = points_n2 + 1.0 # >> Add 1.0 to every element of the array.

-   Storing:

    .. code-block:: python

        # A Python list can store different objects with different size.
        points_n1 =([[[0.0, 1.0], [1.0, 2.0]], np.ones((3,2))])
        # A Numpy array can only store Numpy objects with the same size.
        points_n2 = np.asarray(([[[0.0, 1.0], [1.0, 2.0]], np.ones((2,2))]))

In Discorpy, Python lists are mainly used to store Numpy arrays with different
size. For example, in the :ref:`grouping <methods>` step, the number of extracted
reference-points on each line are not the same.

.. _axes_origin:

It is important to be aware that different methods may use different coordinate systems.
In the image coordinate system, the origin is at the top-left corner of the image.
For example, the `center-of-mass <https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.center_of_mass.html>`__
of an object refers to this origin. In the `plotting <https://matplotlib.org/stable/tutorials/introductory/pyplot.html>`__
coordinate system, the origin is at the bottom-left. When overlaying reference points on a calibration image using
Matplotlib, be aware of this difference.

The coordinates of reference points extracted from a calibration image refer to the image coordinate system
(top-left origin). The values for the center of distortion are also referenced from this origin.

Demonstrations
--------------

.. toctree::

    usage/demo_01
    usage/demo_02
    usage/demo_03
    usage/demo_04
    usage/demo_05
    usage/demo_06
    usage/demo_07
    usage/demo_08
    usage/tips
