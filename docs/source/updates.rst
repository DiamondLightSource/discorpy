Update notes
============

-   11/04/2025:

    +   Added methods for calibrating strong distortion (fisheye images):

        *   Grouping points on strongly curved lines.
        *   Calculating the center of distortion using vanishing point approaches.
        *   Correcting perspective distortion.

    +   Updated loader-saver module to support Windows file paths.
    +   Added functions for saving and loading .json and .pickle files.
    +   Improved methods for detecting points in line-pattern and chessboard images.
    +   Added masking methods for removing unwanted points from calibration images.
    +   Release version 1.7.

-   20/11/2023:

    +   Add module of utility methods:

        *   Generate grid-pattern images.
        *   Find corresponding points between the distorted and undistorted space.
        *   Unwarp a color image with an option to keep the original field of view.
        *   Unwarping an image or video using the 'remap' function in Opencv for fast performance.

    +   Release version 1.6.

-   02/03/2023:

    +   Improve methods for perspective distortion correction.
    +   Add missing options to enable/disable optimization.
    +   Improve/add tests.
    +   Realise version 1.5

-   21/11/2021:

    +   Add methods for correcting perspective distortion.
    +   Add pre-processing methods for a line-pattern image and a chessboard image.
    +   Add Readthedocs documentation page.
    +   Realise version 1.4.

- 01/06/2021:

    +   Change name to Discorpy.
    +   Release version 1.3.