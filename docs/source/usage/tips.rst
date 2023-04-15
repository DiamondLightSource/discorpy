Useful tips
===========

1)  Apply distortion correction to a color image. Given distortion coefficients
    from the calibration process shown in previous demos, we have to apply unwarping
    to each channel of the image as follows:

    .. code-block:: python

        import numpy as np
        import discorpy.losa.loadersaver as io
        import discorpy.proc.processing as proc
        import discorpy.post.postprocessing as post

        file_path = "../../data/percy_cam/F_R_hazcam.png"
        output_base = "E:/output/unwarping/"

        img = io.load_image(file_path, average=False)
        # Define a function for convenient use.
        def unwarp_color_image(img, xcenter, ycenter, list_bfact, mode='reflect'):
            if len(img.shape) != 3:
                raise ValueError("Check if input is a color image!!!")
            img_corrected = []
            for i in range(img.shape[-1]):
                img_corrected.append(post.unwarp_image_backward(img[:, :, i], xcenter,
                                                                ycenter, list_bfact, mode=mode))
            img_corrected = np.moveaxis(np.asarray(img_corrected), 0, 2)
            return img_corrected

        # Suppose that we get xcenter, ycenter, list_bfact (backward mode) to this point
        img_corrected = unwarp_color_image(img, xcenter, ycenter, list_bfact)
        io.save_image(output_base + "/corrected_F_R_hazcam.png", img_corrected)

2)  Find the coordinates of a point in one space given its corresponding
    coordinates in another space. Suppose that we get distortion coefficients of
    the backward model; given a point in a distorted image, we want to find its
    position in the undistorted image. This can be done as the following:

    .. code-block:: python

        import numpy as np
        import discorpy.losa.loadersaver as io
        import discorpy.proc.processing as proc
        import discorpy.post.postprocessing as post

        img = io.load_image("E:/data/image.jpg", average=True)
        # Suppose that we get coefficients of the backward model (xcenter, yxenter, list_bfact).
        # We need to find the forward transformation using the given backward model.
        (height, width) = img.shape
        ref_points = [[i - ycenter, j - xcenter] for i in np.linspace(0, height, 40) for j in
                      np.linspace(0, height, 40)]
        list_ffact = proc.transform_coef_backward_and_forward(list_bfact, ref_points=ref_points)

        # Define the function to calculate corresponding points between distorted and undistorted space
        # This function can be used both ways:
        # In: distorted point, forward model :-> Out: undistorted point
        # In: undistorted point, backward model :-> Out: distorted point
        def find_point_to_point(points, xcenter, ycenter, list_fact):
            """
            points : (row_index, column_index) of the point.
            """
            xi, yi = points[1] - xcenter, points[0] - ycenter
            ri = np.sqrt(xi * xi + yi * yi)
            factor = np.float64(np.sum(list_fact * np.power(ri, np.arange(len(list_fact)))))
            xo = xcenter + factor * xi
            yo = ycenter + factor * yi
            return xo, yo

        # Find top-left point in the undistorted space given top-left point in the distorted space.
        xu_top_left, yu_top_left = find_point_to_point((0, 0), xcenter, ycenter, list_ffact)
        # Find bottom-right point in the undistorted space given bottom-right point in the distorted space.
        xu_bot_right, yu_bot_right = find_point_to_point((height - 1, width - 1), xcenter, ycenter,
                                                         list_ffact)

3)  Maintain the original field of view of an image. In correction for the barrel distortion,
    some parts of the image edges will be cropped. However, there are two possible methods to
    incorporate these areas back into the undistorted image:

        +   Apply padding to the original image before unwarping. The appropriate
            amount of padding can be determined using the results from the example mentioned above.
            Note that we have to update the center of distortion corresponding to the padding.

            .. code-block:: python

                # Calculate padding width for each side.
                pad_top = int(np.abs(yu_top_left))
                pad_bot = int(yu_bot_right - height)
                pad_left = int(np.abs(xu_top_left))
                pad_right = int(xu_bot_right - width)

                img_pad = np.pad(img, ((pad_top, pad_bot), (pad_left, pad_right), (0, 0)), mode="constant")
                img_corrected = unwarp_color_image(img_pad, xcenter + pad_left, ycenter + pad_top,
                                                   list_bfact, mode='constant')
                io.save_image(output_base + "/F_R_hazcam_unwarped_padding.jpg", img_corrected)

            .. image:: figs/tips/img1.jpg
                :width: 100 %
                :align: center

        +   Rescale the distortion coefficients and resize the original image correspondingly.

            .. code-block:: python

                import scipy.ndimage as ndi

                zoom = 2.0
                list_bfact1 = zoom * list_bfact
                xcenter1 = xcenter * zoom
                ycenter1 = ycenter * zoom
                img_corrected = []
                for i in range(img.shape[-1]):
                    img_tmp = ndi.zoom(img[:, :, i], zoom)
                    img_tmp = post.unwarp_image_backward(img_tmp, xcenter1, ycenter1, list_bfact1)
                    img_corrected.append(img_tmp)
                img_corrected = np.moveaxis(np.asarray(img_corrected), 0, 2)
                io.save_image(output_base + "/F_R_hazcam_unwarped_zoom.jpg", img_corrected)

            .. image:: figs/tips/img2.jpg
                :width: 100 %
                :align: center
