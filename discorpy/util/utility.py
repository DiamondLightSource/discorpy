# ============================================================================
# Copyright (c) 2023 Diamond Light Source Ltd. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ============================================================================
# Author: Nghia T. Vo
# E-mail:
# Description: Module for utility methods
# Publication date: 10 July 2023
# ============================================================================
# Contributors:
# ============================================================================


"""
Module of utility methods:

-   Generate a dot-pattern, line-pattern, and chessboard image.
-   Find corresponding points between the distorted and undistorted space.
-   Unwarp a color image with an option to keep the original field of view.
-   Unwarping an image or video using the 'remap' function in Opencv for fast
    performance.
"""

import numpy as np
import scipy.ndimage as ndi
from scipy.ndimage import map_coordinates
import discorpy.proc.processing as proc


def make_circle_mask(width, ratio):
    """
    Create a circle mask.

    Parameters
    -----------
    width : int
        Width of a square array.
    ratio : float
        Ratio between the diameter of the mask and the width of the array.

    Returns
    ------
    array_like
         Square array.
    """
    mask = np.zeros((width, width), dtype=np.float32)
    center = width // 2
    radius = ratio * center
    y, x = np.ogrid[-center:width - center, -center:width - center]
    mask_check = x * x + y * y <= radius * radius
    mask[mask_check] = 1.0
    return mask


def make_dot_pattern(height=1800, width=2000, dot_distance=90,
                     dot_size=15, margin=150):
    """
    Generate a dot-pattern image.

    Parameters
    ----------
    height : int
        Height of the image.
    width : int
        Width of the image.
    dot_distance : int
        Distance between two dots.
    dot_size : int
        Size of each dot.
    margin : int
        Blank area between the dots and the edges.

    Returns
    -------
    array_like
        Dot-pattern image.
    """
    dot_size = np.clip(dot_size, 1, min(height, width) // 8)
    if dot_distance < dot_size:
        raise ValueError("Dot size must be smaller than the dot-distance!!!")
    mat = np.zeros((height, width), dtype=np.float32)
    if isinstance(margin, tuple) or isinstance(margin, list):
        marg_ver, marg_hor = margin[0:2]
    else:
        marg_ver = marg_hor = margin
    half_dot = dot_size // 2 + 1
    mask = make_circle_mask(dot_size, 1.0)
    mat[marg_ver + half_dot: height - marg_ver - half_dot:dot_distance,
        marg_hor + half_dot: width - marg_hor - half_dot:dot_distance] = 1
    dot_pattern = np.float32(ndi.binary_dilation(mat,
                                                 iterations=1, structure=mask))
    return 1 - dot_pattern


def make_line_pattern(height=1800, width=2000, line_distance=90,
                      line_size=7, margin=100):
    """
    Generate a dot-pattern image.

    Parameters
    ----------
    height : int
        Height of the image.
    width : int
        Width of the image.
    line_distance : int
        Distance between two dots.
    line_size : int
        Size of each dot.
    margin : int
        Blank area between the lines and the edges.

    Returns
    -------
    array_like
        Dot-pattern image.
    """
    line_size = np.clip(line_size, 1, min(height, width) // 8)
    mat = np.zeros((height, width), dtype=np.float32)
    if isinstance(margin, tuple) or isinstance(margin, list):
        marg_ver, marg_hor = margin[0:2]
    else:
        marg_ver = marg_hor = margin
    list_i = np.arange(marg_ver, height - marg_ver - line_size,
                       line_distance)
    list_j = np.arange(marg_hor, width - marg_hor - line_size,
                       line_distance)
    for i in list_i:
        mat[i:i + line_size, list_j[0]:list_j[-1] + line_size] = 1
    for j in list_j:
        mat[list_i[0]:list_i[-1] + line_size, j:j + line_size] = 1
    return 1 - mat


def make_chessboard(height=1800, width=2000, size=100, margin=100,
                    margin_grayscale=0.95):
    """
    Generate a chessboard image.

    Parameters
    ----------
    height : int
        Height of the image.
    width : int
        Width of the image.
    size : int
        Size of each cell.
    margin : int
        Blank area between the chessboard and the edges.
    margin_grayscale : float
        Gray scale of margin area (0: black 1: white).

    Returns
    -------
    array_like
        Chessboard image.
    """
    mat = margin_grayscale * np.ones((height, width), dtype=np.float32)
    num = 0
    for i in range(size + margin, height - margin - size, size):
        if num % 2 == 0:
            num1 = 0
            for j in range(size + margin, width - margin - size, size):
                if num1 % 2 == 0:
                    mat[i: i + size, j: j + size] = 1.0
                else:
                    mat[i: i + size, j: j + size] = 0.0
                num1 = num1 + 1
        else:
            num1 = 0
            for j in range(size + margin, width - margin - size, size):
                if num1 % 2 == 0:
                    mat[i: i + size, j: j + size] = 0.0
                else:
                    mat[i: i + size, j: j + size] = 1.0
                num1 = num1 + 1
        num = num + 1
    return mat


def find_point_to_point(points, xcenter, ycenter, list_fact,
                        output_order="xy"):
    """
    Calculate corresponding points between distorted and undistorted space.
    This function can be used both ways:
    - Given the input is a distorted point and a forward model, the output
    is the undistorted point.
    - Given the input is an undistorted point and a backward model, the output
    is the distorted point.

    Parameters
    ----------
    points : tuple of point indexes.
        Tuple of (row_index, column_index) of the point. Note that the format
        is ij-index not xy-coordinate.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_fact : list of float
        Polynomial coefficients of the backward/forward model.
    output_order : {"xy", "yx"}
        To select the order of the output. "yx" <-> "ij".

    Returns
    -------
    tuple of float
        x- and y-coordinate of the point in another space.
    """
    xi, yi = points[1] - xcenter, points[0] - ycenter
    ri = np.sqrt(xi * xi + yi * yi)
    factor = np.float64(np.sum(
        list_fact * np.power(ri, np.arange(len(list_fact)))))
    xo = xcenter + factor * xi
    yo = ycenter + factor * yi
    if output_order == "xy":
        return xo, yo
    else:
        return yo, xo


def _calc_pad(pad, height, width, xcenter, ycenter, list_fact):
    """
    Supplementary method for calculating pad-width.
    """
    t_pad, b_pad, l_pad, r_pad = 0, 0, 0, 0
    if isinstance(pad, bool):
        if pad is True:
            ref_points = [[i - ycenter, j - xcenter] for i in
                          np.linspace(0, height, 40) for j in
                          np.linspace(0, width, 40)]
            list_tfact = proc.transform_coef_backward_and_forward(
                                            list_fact, ref_points=ref_points)
            xu_tl, yu_tl = find_point_to_point((0, 0), xcenter, ycenter,
                                               list_tfact)
            xu_tr, yu_tr = find_point_to_point((0, width - 1), xcenter,
                                               ycenter, list_tfact)
            xu_br, yu_br = find_point_to_point((height - 1, width - 1),
                                               xcenter, ycenter, list_tfact)
            xu_bl, yu_bl = find_point_to_point((height - 1, 0),
                                               xcenter, ycenter, list_tfact)
            l_val = min(xu_tl, xu_bl)
            if l_val < 0:
                l_pad = int(-l_val)
            r_val = max(xu_tr, xu_br)
            if r_val > width:
                r_pad = int(r_val - width)
            t_val = min(yu_tl, yu_tr)
            if t_val < 0:
                t_pad = int(-t_val)
            b_val = max(yu_bl, yu_br)
            if b_val > height:
                b_pad = int(b_val - height)
    elif isinstance(pad, int):
        t_pad = b_pad = l_pad = r_pad = pad
    elif isinstance(pad, tuple) or isinstance(pad, list):
        if len(pad) != 4:
            raise ValueError("Incorrect format!!! Please use a tuple/list of "
                             "(top_pad, bottom_pad, left_pad, right_pad)")
        else:
            t_pad, b_pad, l_pad, r_pad = pad
    else:
        raise ValueError("Invalid format of the 'pad' parameter!!!")
    return t_pad, b_pad, l_pad, r_pad


def unwarp_color_image_backward(mat, xcenter, ycenter, list_fact, order=1,
                                mode="reflect", pad=False,
                                pad_mode='constant'):
    """
    Unwarp a color image using a backward model.

    Parameters
    ----------
    mat : array_like
        2D/3D array.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_fact : list of float
        Polynomial coefficients of the backward model.
    order : int, optional.
        The order of the spline interpolation.
    mode : {'reflect', 'grid-mirror', 'constant', 'grid-constant', 'nearest',
           'mirror', 'grid-wrap', 'wrap'}, optional
        To determine how to handle image boundaries.
    pad : bool, int, or tuple of int.
        Use to keep the original view. If pad is True, the pad width is
        calculated automatically.
    pad_mode : str
        To select a method for padding: 'reflect', 'edge', 'mean',
        'linear_ramp', 'symmetric',...

    Returns
    -------
    array_like
        2D/3D array. Distortion-corrected image.
    """
    (height, width) = mat.shape[:2]
    t_pad, b_pad, l_pad, r_pad = _calc_pad(pad, height, width, xcenter,
                                           ycenter, list_fact)
    num_dim = len(mat.shape)
    if num_dim == 2:
        pad_width = [(t_pad, b_pad), (l_pad, r_pad)]
    else:
        pad_width = [(t_pad, b_pad), (l_pad, r_pad), (0, 0)]
    mat_pad = np.pad(mat, pad_width, mode=pad_mode)
    (height, width) = mat_pad.shape[:2]
    xcenter = xcenter + l_pad
    ycenter = ycenter + t_pad
    xu_list = np.arange(width) - xcenter
    yu_list = np.arange(height) - ycenter
    xu_mat, yu_mat = np.meshgrid(xu_list, yu_list)
    ru_mat = np.sqrt(xu_mat ** 2 + yu_mat ** 2)
    fact_mat = np.sum(np.asarray(
        [factor * ru_mat ** i for i, factor in enumerate(list_fact)]), axis=0)
    xd_mat = np.float32(np.clip(xcenter + fact_mat * xu_mat, 0, width - 1))
    yd_mat = np.float32(np.clip(ycenter + fact_mat * yu_mat, 0, height - 1))
    indices = np.reshape(yd_mat, (-1, 1)), np.reshape(xd_mat, (-1, 1))
    if num_dim == 2:
        mat_corr = np.reshape(map_coordinates(
            mat_pad, indices, order=order, mode=mode), (height, width))
    else:
        mat_corr = []
        for i in range(mat_pad.shape[-1]):
            mat_corr.append(np.reshape(map_coordinates(
                mat_pad[:, :, i], indices, order=order, mode=mode),
                (height, width)))
        mat_corr = np.moveaxis(np.asarray(mat_corr), 0, 2)
    return mat_corr


def mapping_cv2(mat, xmat, ymat, method=None, border=None):  # pragma: no cover
    """
    Apply a geometric transformation to a 2D array using Opencv.

    Parameters
    ----------
    mat : array_like.
        2D/3D array.
    xmat : array_like
        2D array of the x-coordinates. Origin is at the left of the image.
    ymat : array_like
        2D array of the y-coordinates. Origin is at the top of the image.
    method : opencv-object
        To select interpolation method. Note to use the prefix: cv2.<method>
        https://tinyurl.com/3afmv6jc
    border : opencv-object
        To select method for boundary handling. Note to use the prefix:
        cv2.<method> https://tinyurl.com/52xzkwt2

    Returns
    -------
    array_like
    """
    try:
        import cv2
    except ImportError:
        raise ValueError("You must install Opencv before using this function!")
    if method is None:
        method = cv2.INTER_LINEAR
    if border is None:
        border = cv2.BORDER_CONSTANT
    mat = cv2.remap(mat, xmat, ymat, interpolation=method, borderMode=border)
    return mat


def unwarp_image_backward_cv2(mat, xcenter, ycenter, list_fact, method=None,
                              border=None, pad=False,
                              pad_mode='constant'):  # pragma: no cover
    """
    Unwarp an image using the backward model with the Opencv 'remap' function
    for fast performance.

    Parameters
    ----------
    mat : array_like
        Input image. Can be a color image.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_fact : list of float
        Polynomial coefficients of the backward model.
    method : opencv-object
        To select interpolation method. Note to use the prefix: cv2.<method>
        https://tinyurl.com/3afmv6jc
    border : opencv-object
        To select method for boundary handling. Note to use the prefix:
        cv2.<method> https://tinyurl.com/52xzkwt2
    pad : bool, int, or tuple of int.
        Use to keep the original view. If pad is True, the pad width is
        calculated automatically.
    pad_mode : str
        To select a method for padding: 'reflect', 'edge', 'mean',
        'linear_ramp', 'symmetric',...

    Returns
    -------
    array_like
        Distortion-corrected image.
    """
    (height, width) = mat.shape[:2]
    t_pad, b_pad, l_pad, r_pad = _calc_pad(pad, height, width, xcenter,
                                           ycenter, list_fact)
    num_dim = len(mat.shape)
    if num_dim == 2:
        pad_width = [(t_pad, b_pad), (l_pad, r_pad)]
    else:
        pad_width = [(t_pad, b_pad), (l_pad, r_pad), (0, 0)]
    mat_pad = np.pad(mat, pad_width, mode=pad_mode)
    (height, width) = mat_pad.shape[:2]
    xu_list = np.arange(width) - (xcenter + l_pad)
    yu_list = np.arange(height) - (ycenter + t_pad)
    xu_mat, yu_mat = np.meshgrid(xu_list, yu_list)
    ru_mat = np.sqrt(xu_mat ** 2 + yu_mat ** 2)
    fact_mat = np.sum(np.asarray(
        [factor * ru_mat ** i for i, factor in enumerate(list_fact)]), axis=0)
    xd_mat = np.float32(
        np.clip(xcenter + l_pad + fact_mat * xu_mat, 0, width - 1))
    yd_mat = np.float32(
        np.clip(ycenter + t_pad + fact_mat * yu_mat, 0, height - 1))
    mat = mapping_cv2(mat_pad, xd_mat, yd_mat, method=method, border=border)
    return mat


def unwarp_video_cv2(cam_obj, xcenter, ycenter, list_fact, method=None,
                     border=None, pad=True,
                     pad_mode='constant'):  # pragma: no cover
    """
    Unwarp frames from Opencv video object using the backward model.

    Parameters
    ----------
    cam_obj : obj
        Opencv camera object. e.g. cv2.VideoCapture(0)
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_fact : list of float
        Polynomial coefficients of the backward model.
    method : opencv-object
        To select interpolation method. Note to use the prefix: cv2.<method>
        https://tinyurl.com/3afmv6jc
    border : opencv-object
        To select method for boundary handling. Note to use the prefix:
        cv2.<method> https://tinyurl.com/52xzkwt2
    pad : bool, int, or tuple of int.
        Use to keep the original view. If pad is True, the pad width is
        calculated automatically.
    pad_mode : str
        To select a method for padding: 'reflect', 'edge', 'mean',
        'linear_ramp', 'symmetric',...

    Returns
    -------
    Generator
    """
    try:
        import cv2
    except ImportError:
        raise ValueError("You must install Opencv before using this function!")
    width = int(cam_obj.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cam_obj.get(cv2.CAP_PROP_FRAME_HEIGHT))
    t_pad, b_pad, l_pad, r_pad = _calc_pad(pad, height, width, xcenter,
                                           ycenter, list_fact)
    width = width + l_pad + r_pad
    height = height + t_pad + b_pad
    xcenter = xcenter + l_pad
    ycenter = ycenter + t_pad
    xu_list = np.arange(width) - xcenter
    yu_list = np.arange(height) - ycenter
    xu_mat, yu_mat = np.meshgrid(xu_list, yu_list)
    ru_mat = np.sqrt(xu_mat ** 2 + yu_mat ** 2)
    fact_mat = np.sum(np.asarray(
        [factor * ru_mat ** i for i, factor in enumerate(list_fact)]), axis=0)
    xd_mat = np.float32(np.clip(xcenter + fact_mat * xu_mat, 0, width - 1))
    yd_mat = np.float32(np.clip(ycenter + fact_mat * yu_mat, 0, height - 1))
    check, frame = cam_obj.read()
    if check:
        num_dim = len(frame.shape)
        if num_dim == 2:
            pad_width = [(t_pad, b_pad), (l_pad, r_pad)]
        else:
            pad_width = [(t_pad, b_pad), (l_pad, r_pad), (0, 0)]
    while True:
        check, frame = cam_obj.read()
        if check:
            if (pad is False) or (pad == 0):
                uframe = mapping_cv2(frame, xd_mat, yd_mat, method=method,
                                     border=border)
            else:
                frame_pad = np.pad(frame, pad_width, mode=pad_mode)
                uframe = mapping_cv2(frame_pad, xd_mat, yd_mat, method=method,
                                     border=border)
            cv2.imshow("Press ESC to exit", uframe)
            c = cv2.waitKey(1)
            if c == 27:
                break
    cam_obj.release()
    cv2.destroyAllWindows()
