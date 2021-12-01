.. _demo_01:

Process a high-quality calibration image
========================================

The following workflow shows how to use Discorpy to process a high-quality
calibration `image <https://github.com/DiamondLightSource/discorpy/blob/master/data/dot_pattern_01.jpg>`_
acquired at Beamline I12, Diamond Light Source.

- Load the image:

    .. code-block:: python

        import numpy as np
        import discorpy.losa.loadersaver as io
        import discorpy.prep.preprocessing as prep
        import discorpy.proc.processing as proc
        import discorpy.post.postprocessing as post

        # Initial parameters
        file_path = "../data/dot_pattern_01.jpg"
        output_base = "E:/correction/"
        num_coef = 5  # Number of polynomial coefficients
        mat0 = io.load_image(file_path) # Load image
        (height, width) = mat0.shape

    .. figure:: figs/demo_01/fig1.jpg
        :name: fig_23
        :figwidth: 75 %
        :align: center
        :figclass: align-center

        Visible dot-target.

- Extract reference-points:

    .. code-block:: python
        :emphasize-lines: 14,15

        # Segment dots
        mat1 = prep.binarization(mat0)
        # Calculate the median dot size and distance between them.
        (dot_size, dot_dist) = prep.calc_size_distance(mat1)
        # Remove non-dot objects
        mat1 = prep.select_dots_based_size(mat1, dot_size)
        # Remove non-elliptical objects
        mat1 = prep.select_dots_based_ratio(mat1)
        io.save_image(output_base + "/segmented_dots.jpg", mat1) # Save image for checking
        # Calculate the slopes of horizontal lines and vertical lines.
        hor_slope = prep.calc_hor_slope(mat1)
        ver_slope = prep.calc_ver_slope(mat1)
        print("Horizontal slope: {0}. Vertical slope {1}".format(hor_slope, ver_slope))

        #>> Horizontal slope: 0.01124473800478091. Vertical slope -0.011342266682773354

    .. figure:: figs/demo_01/fig2.jpg
        :name: fig_24
        :figwidth: 75 %
        :align: center
        :figclass: align-center

        Segmented dots.

- As can be seen from the highlighted output above, the slopes of horizontal lines and
  vertical lines are nearly the same (note that their signs are opposite due to
  the use of different origin-axis). This indicates that perspective distortion
  is negligible. The next step is to group points into lines. After the grouping
  step, perspective effect; :ref:`if there is <Correcting perspective effect>`;
  can be corrected simply by adding a single line of command.

    .. code-block:: python

        # Group points to horizontal lines
        list_hor_lines = prep.group_dots_hor_lines(mat1, hor_slope, dot_dist)
        # Group points to vertical lines
        list_ver_lines = prep.group_dots_ver_lines(mat1, ver_slope, dot_dist)
        # Optional: remove horizontal outliners
        list_hor_lines = prep.remove_residual_dots_hor(list_hor_lines, hor_slope)
        # Optional: remove vertical outliners
        list_ver_lines = prep.remove_residual_dots_ver(list_ver_lines, ver_slope)
        # Save output for checking
        io.save_plot_image(output_base + "/horizontal_lines.png", list_hor_lines, height, width)
        io.save_plot_image(output_base + "/vertical_lines.png", list_ver_lines, height, width)

        # Optional: correct perspective effect. Only available from Discorpy 1.4
        # list_hor_lines, list_ver_lines = proc.regenerate_grid_points_parabola(
        #                      list_hor_lines, list_ver_lines, perspective=True)


    .. figure:: figs/demo_01/fig3.png
        :name: fig_25
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        Reference-points (center-of-mass of dots) after grouped into
        horizontal lines (a) and vertical lines (b).

- We can check the straightness of lines of points by calculating the distances
  of the grouped points to their fitted straight lines. This helps to assess if
  the distortion is significant or not. Note that the origin of the coordinate
  system before distortion correction is at the top-left corner of an image.

    .. code-block:: python

        list_hor_data = post.calc_residual_hor(list_hor_lines, 0.0, 0.0)
        list_ver_data = post.calc_residual_ver(list_ver_lines, 0.0, 0.0)
        io.save_residual_plot(output_base + "/hor_residual_before_correction.png",
                              list_hor_data, height, width)
        io.save_residual_plot(output_base + "/ver_residual_before_correction.png",
                              list_ver_data, height, width)

    .. figure:: figs/demo_01/fig4.png
        :name: fig_26
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        Plot of the distances of the dot-centroids from their fitted
        straight line against their distances from the axes origin. (a) For
        horizontal lines. (b) For vertical lines.

- As shown in :numref:`fig_26`, the residual is more than 12 pixels which means that
  distortion is significant and needs to be corrected. The next step is to
  calculate the center of distortion (COD) and the coefficients of the backward
  mapping for :ref:`a radial distortion model <methods>`.

    .. code-block:: python
        :emphasize-lines: 12-14

        # Calculate the center of distortion
        (xcenter, ycenter) = proc.find_cod_coarse(list_hor_lines, list_ver_lines)
        # Calculate coefficients of the correction model
        list_fact = proc.calc_coef_backward(list_hor_lines, list_ver_lines,
                                            xcenter, ycenter, num_coef)
        # Save the results for later use.
        io.save_metadata_txt(output_base + "/coefficients_radial_distortion.txt",
                             xcenter, ycenter, list_fact)
        print("X-center: {0}. Y-center: {1}".format(xcenter, ycenter))
        print("Coefficients: {0}".format(list_fact))
        """
        >> X-center: 1252.1528590042283. Y-center: 1008.9088499595639
        >> Coefficients: [1.00027631e+00, -1.25730878e-06, -1.43170401e-08,
                          -1.65727563e-12, 7.89109870e-16]
        """

- Using the determined parameters of the correction model, we can unwarp the
  lines of points and check the correction results.

    .. code-block:: python

        # Apply correction to the lines of points
        list_uhor_lines = post.unwarp_line_backward(list_hor_lines, xcenter, ycenter,
                                                    list_fact)
        list_uver_lines = post.unwarp_line_backward(list_ver_lines, xcenter, ycenter,
                                                    list_fact)
        # Save the results for checking
        io.save_plot_image(output_base + "/unwarpped_horizontal_lines.png", list_uhor_lines,
                           height, width)
        io.save_plot_image(output_base + "/unwarpped_vertical_lines.png", list_uver_lines,
                           height, width)
        # Calculate the residual of the unwarpped points.
        list_hor_data = post.calc_residual_hor(list_uhor_lines, xcenter, ycenter)
        list_ver_data = post.calc_residual_ver(list_uver_lines, xcenter, ycenter)
        # Save the results for checking
        io.save_residual_plot(output_base + "/hor_residual_after_correction.png",
                              list_hor_data, height, width)
        io.save_residual_plot(output_base + "/ver_residual_after_correction.png",
                              list_ver_data, height, width)


    .. figure:: figs/demo_01/fig5.png
        :name: fig_27
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        . (a) Unwarpped horizontal lines. (b) Unwarpped vertical lines.


    .. figure:: figs/demo_01/fig6.png
        :name: fig_28
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        Residual of the unwarpped points. Note that the origin of the
        coordinate system is at the center of distortion. (a) For horizontal lines.
        (b) For vertical lines.

- As can be seen from :numref:`fig_28` the accuracy of the correction results is sub-pixel.
  The last step of the workflow is to use the determined model for correcting images.

    .. code-block:: python

        # Load coefficients from previous calculation if need to
        # (xcenter, ycenter, list_fact) = io.load_metadata_txt(
        #     output_base + "/coefficients_radial_distortion.txt")
        # Correct the image
        corrected_mat = post.unwarp_image_backward(mat0, xcenter, ycenter, list_fact)
        # Save results. Note that the output is 32-bit numpy array. Convert to lower-bit if need to.
        io.save_image(output_base + "/corrected_image.tif", corrected_mat)
        io.save_image(output_base + "/difference.tif", corrected_mat - mat0)

    .. figure:: figs/demo_01/fig7.jpg
      :name: fig_29
      :figwidth: 75 %
      :align: center
      :figclass: align-center

      Corrected image.

    .. figure:: figs/demo_01/fig8.jpg
      :name: fig_30
      :figwidth: 75 %
      :align: center
      :figclass: align-center

      Difference between images before (:numref:`fig_23`) and after (:numref:`fig_29`)
      the correction.

Click :download:`here <./codes/demo_01.py>` to download the Python codes.