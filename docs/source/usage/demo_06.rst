Calibrate a camera using a chessboard image
===========================================

The following workflow shows how to calibrate a commercial camera using a chessboard
image. First of all, the chessboard pattern was created using this
`website <https://markhedleyjones.com/projects/calibration-checkerboard-collection>`__.
The generated pattern was printed and stuck to a flat wall (:numref:`fig_65`).
Then, its image was taken using the front-facing camera of a Microsoft Surface Pro
Laptop. The pincushion distortion is visible in the image.

- Load the `chessboard image <https://github.com/DiamondLightSource/discorpy/tree/master/data/laptop_camera>`__,
  convert to a line-pattern image.

    .. code-block:: python

        import numpy as np
        import discorpy.losa.loadersaver as io
        import discorpy.prep.preprocessing as prep
        import discorpy.prep.linepattern as lprep
        import discorpy.proc.processing as proc
        import discorpy.post.postprocessing as post

        # Initial parameters
        file_path = "C:/data/laptop_camera/chessboard.jpg"
        output_base = "./for_demo_06/"
        num_coef = 5  # Number of polynomial coefficients
        mat0 = io.load_image(file_path) # Load image
        (height, width) = mat0.shape

        # Convert the chessboard image to a line-pattern image
        mat1 = lprep.convert_chessboard_to_linepattern(mat0)
        io.save_image(output_base + "/line_pattern_converted.jpg", mat1)

        # Calculate slope and distance between lines
        slope_hor, dist_hor = lprep.calc_slope_distance_hor_lines(mat1, radius=15, sensitive=0.5)
        slope_ver, dist_ver = lprep.calc_slope_distance_ver_lines(mat1, radius=15, sensitive=0.5)
        print("Horizontal slope: ", slope_hor, " Distance: ", dist_hor)
        print("Vertical slope: ", slope_ver, " Distance: ", dist_ver)

    .. figure:: figs/demo_06/fig1.jpg
        :name: fig_65
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        . (a) Chessboard image. (b) Line-pattern image generated from image (a).

- Extract reference-points (adjust parameters: *radius* and/or *sensitive* if need to).
  Note that users can use other methods available in `Scikit-image <https://scikit-image.org/docs/dev/auto_examples/features_detection/plot_corner.html>`__
  to perform this step directly on the chessboard image. There are many unwanted points
  (:numref:`fig_66`) but they can be removed by the grouping method (:numref:`fig_67`)
  in Discorpy.

    .. code-block:: python

        # Extract reference-points
        list_points_hor_lines = lprep.get_cross_points_hor_lines(mat1, slope_ver, dist_ver,
                                                                 ratio=0.3, norm=True, offset=450,
                                                                 bgr="bright", radius=15,
                                                                 sensitive=1.0, denoise=True,
                                                                 subpixel=True)
        list_points_ver_lines = lprep.get_cross_points_ver_lines(mat1, slope_hor, dist_hor,
                                                                 ratio=0.3, norm=True, offset=150,
                                                                 bgr="bright", radius=15,
                                                                 sensitive=1.0, denoise=True,
                                                                 subpixel=True)
        io.save_plot_points(output_base + "/ref_points_horizontal.png", list_points_hor_lines,
                            height, width, color="red")
        io.save_plot_points(output_base + "/ref_points_vertical.png", list_points_ver_lines,
                            height, width, color="blue")

        # Group points into lines
        list_hor_lines = prep.group_dots_hor_lines(list_points_hor_lines, slope_hor, dist_hor,
                                                   ratio=0.1, num_dot_miss=2, accepted_ratio=0.8)
        list_ver_lines = prep.group_dots_ver_lines(list_points_ver_lines, slope_ver, dist_ver,
                                                   ratio=0.1, num_dot_miss=2, accepted_ratio=0.8)

        # Remove residual dots
        list_hor_lines = prep.remove_residual_dots_hor(list_hor_lines, slope_hor, 2.0)
        list_ver_lines = prep.remove_residual_dots_ver(list_ver_lines, slope_ver, 2.0)

        # Save output for checking
        io.save_plot_image(output_base + "/horizontal_lines.png", list_hor_lines, height, width)
        io.save_plot_image(output_base + "/vertical_lines.png", list_ver_lines, height, width)
        list_hor_data = post.calc_residual_hor(list_hor_lines, 0.0, 0.0)
        list_ver_data = post.calc_residual_ver(list_ver_lines, 0.0, 0.0)
        io.save_residual_plot(output_base + "/hor_residual_before_correction.png",
                              list_hor_data, height, width)
        io.save_residual_plot(output_base + "/ver_residual_before_correction.png",
                              list_ver_data, height, width)

    .. figure:: figs/demo_06/fig2.png
        :name: fig_66
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        Extracted reference points from :numref:`fig_65` (b). (a) For horizontal
        lines. (b) For vertical lines.

    .. figure:: figs/demo_06/fig3.png
        :name: fig_67
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        Grouped points. (a) Horizontal lines. (b) Vertical lines.

    .. figure:: figs/demo_06/fig4.png
        :name: fig_68
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        Residual of distorted points. (a) Horizontal lines. (b) Vertical lines.

- Next steps are straightforward. Coefficients of the radial-distortion model are calculated
  where the perspective effect is :ref:`corrected <Correcting perspective effect>` before that.
  The results of applying the model for unwarping lines and images can be seen in
  :numref:`fig_69` and :numref:`fig_70`

    .. code-block:: python

        # Regenerate grid points after correcting the perspective effect.
        list_hor_lines, list_ver_lines = proc.regenerate_grid_points_parabola(
            list_hor_lines, list_ver_lines, perspective=True)

        # Calculate parameters of the radial correction model
        (xcenter, ycenter) = proc.find_cod_coarse(list_hor_lines, list_ver_lines)
        list_fact = proc.calc_coef_backward(list_hor_lines, list_ver_lines,
                                            xcenter, ycenter, num_coef)
        io.save_metadata_txt(output_base + "/coefficients_radial_distortion.txt",
                             xcenter, ycenter, list_fact)
        print("X-center: {0}. Y-center: {1}".format(xcenter, ycenter))
        print("Coefficients: {0}".format(list_fact))

        # Check the correction results:
        # Apply correction to the lines of points
        list_uhor_lines = post.unwarp_line_backward(list_hor_lines, xcenter, ycenter,
                                                    list_fact)
        list_uver_lines = post.unwarp_line_backward(list_ver_lines, xcenter, ycenter,
                                                    list_fact)
        # Calculate the residual of the unwarpped points.
        list_hor_data = post.calc_residual_hor(list_uhor_lines, xcenter, ycenter)
        list_ver_data = post.calc_residual_ver(list_uver_lines, xcenter, ycenter)
        # Save the results for checking
        io.save_plot_image(output_base + "/unwarpped_horizontal_lines.png",
                           list_uhor_lines, height, width)
        io.save_plot_image(output_base + "/unwarpped_vertical_lines.png",
                           list_uver_lines, height, width)
        io.save_residual_plot(output_base + "/hor_residual_after_correction.png",
                              list_hor_data, height, width)
        io.save_residual_plot(output_base + "/ver_residual_after_correction.png",
                              list_ver_data, height, width)

    .. figure:: figs/demo_06/fig5.png
        :name: fig_69
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        Residual of unwarped points. (a) Horizontal lines. (b) Vertical lines.

    .. figure:: figs/demo_06/fig6.jpg
        :name: fig_70
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        . (a) Unwarped image of :numref:`fig_65` (a). (b) Difference between the images
        before and after unwarping.

- Calculated coefficients of the correction model can be used to unwarp
  `another image <https://github.com/DiamondLightSource/discorpy/tree/master/data/laptop_camera>`__
  taken by the same camera as demonstrated in :numref:`fig_71`. For a color image, we have
  to correct each channel of the image.


    .. code-block:: python

        # Load coefficients from previous calculation
        (xcenter, ycenter, list_fact) = io.load_metadata_txt(
            output_base + "/coefficients_radial_distortion.txt")

        # Load an image and correct it.
        img = io.load_image("../../../data/laptop_camera/test_image.jpg", average=False)
        img_corrected = np.copy(img)
        for i in range(img.shape[-1]):
            img_corrected[:, :, i] = post.unwarp_image_backward(img[:, :, i], xcenter,
                                                                ycenter, list_fact)
        io.save_image(output_base + "/test_image_corrected.jpg", img_corrected)

    .. figure:: figs/demo_06/fig7.jpg
        :name: fig_71
        :figwidth: 100 %
        :align: center
        :figclass: align-center

        . (a) Test image taken from the same camera. (b) Unwarped image. The red
        straight line is added for reference.

Click :download:`here <./codes/demo_06.py>` to download the Python codes.
