# ============================================================================
# Copyright (c) 2018 Diamond Light Source Ltd. All rights reserved.
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
# E-mail: nghia.vo@diamond.ac.uk
# Description: Python implementation of the author's methods of
# distortion correction, Nghia T. Vo et al "Radial lens distortion
# correction with sub-pixel accuracy for X-ray micro-tomography"
# Optics Express 23, 32859-32868 (2015), https://doi.org/10.1364/OE.23.032859
# Publication date: 10th July 2018
# ============================================================================

"""
Module of post-processing methods:
- Unwarp a line of dots, an image.
- Generate unwarped slices of a 3D dataset.
- Calculate the residual of unwarped dots.
"""
import numpy as np
from scipy import optimize
from scipy.ndimage import map_coordinates


def unwarp_line_forward(list_lines, xcenter, ycenter, list_fact):
    """
    Unwarp lines of dot-centroids using a forward model.

    Parameters
    ----------
    list_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each line.
    list_fact : list of floats
        Polynomial coefficients of the forward model.

    Returns
    -------
    list_ulines : list of 2D arrays
        List of the unwarped coordinates of dot-centroids on each line.
    """
    list_ulines = []
    list_expo = np.arange(len(list_fact), dtype=np.int16)
    for _, line in enumerate(list_lines):
        uline = np.zeros_like(line)
        for j, dot in enumerate(line):
            xd = dot[1] - xcenter
            yd = dot[0] - ycenter
            rd = np.sqrt(xd * xd + yd * yd)
            factor = np.sum(list_fact * np.power(rd, list_expo))
            uline[j, 1] = xcenter + factor * xd
            uline[j, 0] = ycenter + factor * yd
        list_ulines.append(uline)
    return list_ulines


def _func_diff(ru, rd, *list_fact):
    return (rd - ru * np.sum(np.asarray([fact * ru ** i for i,
                                                            fact in
                                         enumerate(list_fact)]))) ** 2


def unwarp_line_backward(list_lines, xcenter, ycenter, list_fact):
    """
    Unwarp lines of dot-centroids using a backward model. The method finds the
    coordinates of undistorted points from the coordinates of distorted
    points using numerical optimzation.

    Parameters
    ----------
    list_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each line.
    list_fact : list of floats
        Polynomial coefficients of the backward model.

    Returns
    -------
    list_ulines : list of 2D arrays
        List of the unwarped coordinates of dot-centroids on each line.
    """
    list_ulines = []
    for _, line in enumerate(list_lines):
        uline = np.zeros_like(line)
        for j, dot in enumerate(line):
            xd = dot[1] - xcenter
            yd = dot[0] - ycenter
            rd = np.sqrt(xd * xd + yd * yd)
            list_arg = [rd]
            list_arg.extend(list_fact)
            minimum = optimize.minimize(_func_diff, rd, args=tuple(list_arg))
            ru = minimum.x[0]
            factor = ru / rd
            uline[j, 1] = xcenter + factor * xd
            uline[j, 0] = ycenter + factor * yd
        list_ulines.append(uline)
    return list_ulines


def unwarp_image_backward(mat, xcenter, ycenter, list_fact):
    """
    Unwarp an image using a backward model.

    Parameters
    ----------
    mat : array_like
        2D array.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_fact : list of float
        Polynomial coefficients of the backward model.

    Returns
    -------
    array_like
        2D array. Distortion-corrected image.
    """
    (height, width) = mat.shape
    xu_list = np.arange(width) - xcenter
    yu_list = np.arange(height) - ycenter
    xu_mat, yu_mat = np.meshgrid(xu_list, yu_list)
    ru_mat = np.sqrt(xu_mat ** 2 + yu_mat ** 2)
    fact_mat = np.sum(
        np.asarray([factor * ru_mat ** i for i,
                                             factor in enumerate(list_fact)]),
        axis=0)
    xd_mat = np.float32(np.clip(xcenter + fact_mat * xu_mat, 0, width - 1))
    yd_mat = np.float32(np.clip(ycenter + fact_mat * yu_mat, 0, height - 1))
    indices = np.reshape(yd_mat, (-1, 1)), np.reshape(xd_mat, (-1, 1))
    mat = map_coordinates(mat, indices, order=1, mode='reflect')
    return mat.reshape((height, width))


def unwarp_image_forward(mat, xcenter, ycenter, list_fact):
    """
    Unwarp an image using a forward model. Should be used only for assessment
    due to the problem of vacant pixels.

    Parameters
    ----------
    mat : float
        2D array.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_fact : list of floats
        Polynomial coefficients of the forward model.

    Returns
    -------
    array_like
        2D array. Distortion-corrected image.
    """
    (height, width) = mat.shape
    xd_list = np.arange(width) - xcenter
    yd_list = np.arange(height) - ycenter
    xd_mat, yd_mat = np.meshgrid(xd_list, yd_list)
    rd_mat = np.sqrt(xd_mat ** 2 + yd_mat ** 2)
    fact_mat = np.sum(
        np.asarray([factor * rd_mat ** i for i,
                                             factor in enumerate(list_fact)]),
        axis=0)
    xu_mat = np.intp(
        np.round(np.clip(xcenter + fact_mat * xd_mat, 0, width - 1)))
    yu_mat = np.intp(
        np.round(np.clip(ycenter + fact_mat * yd_mat, 0, height - 1)))
    mat_unw = np.zeros_like(mat)
    mat_unw[yu_mat, xu_mat] = mat
    return mat_unw


def unwarp_slice_backward(mat3D, xcenter, ycenter, list_fact, index):
    """
    Generate an unwarped slice [:,index.:] of a 3D dataset, i.e.
    one unwarped sinogram of a 3D tomographic data.

    Parameters
    ----------
    mat3D : array_like
        3D array. Correction is applied along axis 1.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_fact : list of floats
        Polynomial coefficients of a backward model.
    index : int
        Index of the slice

    Returns
    -------
    array_like
        2D array. Distortion-corrected slice.
    """
    if len(mat3D.shape) < 3:
        raise ValueError("Input must be a 3D data")
    (depth, height, width) = mat3D.shape
    xu_list = np.arange(0, width) - xcenter
    yu = index - ycenter
    ru_list = np.sqrt(xu_list ** 2 + yu ** 2)
    flist = np.sum(
        np.asarray([factor * ru_list ** i for i,
                                              factor in enumerate(list_fact)]),
        axis=0)
    xd_list = np.clip(xcenter + flist * xu_list, 0, width - 1)
    yd_list = np.clip(ycenter + flist * yu, 0, height - 1)
    yd_min = np.int16(np.floor(np.amin(yd_list)))
    yd_max = np.int16(np.ceil(np.amax(yd_list))) + 1
    yd_list = yd_list - yd_min
    sino = np.zeros((depth, width), dtype=np.float32)
    indices = yd_list, xd_list
    for i in range(depth):
        sino[i] = map_coordinates(
            mat3D[i, yd_min:yd_max, :], indices, order=1, mode='reflect')
    return sino


def _mapping(mat, xmat, ymat):
    """
    Apply a geometric transformation to a 2D array

    Parameters
    ----------
    mat : array_like
        2D array.
    xmat : array_like
        2D array of x-coordinates.
    ymat : array_like
        2D array of y-coordinates.

    Returns
    -------
    array_like
        2D array.
    """
    coord = np.vstack((np.ndarray.flatten(ymat), np.ndarray.flatten(xmat)))
    mat = map_coordinates(mat, coord, order=1, mode='reflect')
    return mat.reshape(xmat.shape)


def unwarp_chunk_slices_backward(mat3D, xcenter, ycenter, list_fact,
                                 start_index, stop_index):
    """
    Generate a chunk  of unwarped slices [:,start_index: stop_index, :] used
    for tomographic data.

    Parameters
    ----------
    mat3D : array_like
        3D array. Correction is applied along axis 1.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_fact : list of floats
        Polynomial coefficients of a backward model.
    start_index : int
        Starting index of slices.
    stop_index : int
        Stopping index of slices.

    Returns
    -------
    array_like
        3D array. Distortion-corrected slices.
    """
    if (len(mat3D.shape) < 3):
        raise ValueError("Input must be a 3D data")
    (depth, height, width) = mat3D.shape
    index_list = np.arange(height, dtype=np.int16)
    if stop_index == -1:
        stop_index = height
    if (start_index not in index_list) or (stop_index not in index_list):
        raise ValueError("Selected index is out of the range")
    xu_list = np.arange(0, width) - xcenter
    yu1 = start_index - ycenter
    ru_list = np.sqrt(xu_list ** 2 + yu1 ** 2)
    flist = np.sum(
        np.asarray([factor * ru_list ** i for i,
                                              factor in enumerate(list_fact)]),
        axis=0)
    yd_list1 = np.clip(ycenter + flist * yu1, 0, height - 1)
    yu2 = stop_index - ycenter
    ru_list = np.sqrt(xu_list ** 2 + yu2 ** 2)
    flist = np.sum(
        np.asarray([factor * ru_list ** i for i,
                                              factor in enumerate(list_fact)]),
        axis=0)
    yd_list2 = np.clip(ycenter + flist * yu2, 0, height - 1)
    yd_min = np.int16(np.floor(np.amin(yd_list1)))
    yd_max = np.int16(np.ceil(np.amax(yd_list2))) + 1
    yu_list = np.arange(start_index, stop_index + 1) - ycenter
    xu_mat, yu_mat = np.meshgrid(xu_list, yu_list)
    ru_mat = np.sqrt(xu_mat ** 2 + yu_mat ** 2)
    fact_mat = np.sum(
        np.asarray([factor * ru_mat ** i for i,
                                             factor in enumerate(list_fact)]),
        axis=0)
    xd_mat = np.float32(np.clip(xcenter + fact_mat * xu_mat, 0, width - 1))
    yd_mat = np.float32(
        np.clip(ycenter + fact_mat * yu_mat, 0, height - 1)) - yd_min
    sino_chunk = np.asarray(
        [_mapping(mat3D[i, yd_min: yd_max, :],
                  xd_mat, yd_mat) for i in range(depth)])
    return sino_chunk


def calc_residual_hor(list_ulines, xcenter, ycenter):
    """
    Calculate the distances of unwarped dots (on each horizontal line) to
    each fitted straight line which is used to assess the straightness of
    unwarped lines.

    Parameters
    ----------
    list_ulines : list of 2D arrays
        List of the coordinates of dot-centroids on each unwarped horizontal
        line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.

    Returns
    -------
    array_like
        2D array. Each element has two values: 1) Distance of a dot to the
        center of distortion; 2) Distance of this dot to the nearest fitted
        straight line.
    """
    list_data = []
    for i, line in enumerate(list_ulines):
        y_list = line[:, 0] - ycenter
        x_list = line[:, 1] - xcenter
        (a_fact, b_fact) = np.polyfit(x_list, y_list, 1)
        dist_list = np.abs(
            a_fact * x_list - y_list + b_fact) / np.sqrt(a_fact ** 2 + 1)
        radi_list = np.sqrt(x_list ** 2 + y_list ** 2)
        list_tmp = np.asarray(list(zip(radi_list, dist_list)))
        list_data.extend(list_tmp)
    list_data = np.asarray(list_data)
    return list_data[list_data[:, 0].argsort()]


def calc_residual_ver(list_ulines, xcenter, ycenter):
    """
    Calculate the distances of unwarped dots (on each vertical line) to each
    fitted straight line which is used to assess the straightness of unwarped
    lines.

    Parameters
    ----------
    list_ulines : list of 2D arrays
        List of the coordinates of dot-centroids on each unwarped vertical line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.

    Returns
    -------
    array_like
        2D array. Each element has two values: 1) Distance of a dot to the
        center of distortion; 2) Distance of this dot to the nearest fitted
        straight line.
    """
    list_data = []
    for i, line in enumerate(list_ulines):
        y_list = line[:, 0] - ycenter
        x_list = line[:, 1] - xcenter
        (a_fact, b_fact) = np.polyfit(y_list, x_list, 1)
        dist_list = np.abs(
            a_fact * y_list - x_list + b_fact) / np.sqrt(a_fact ** 2 + 1)
        radi_list = np.sqrt(x_list ** 2 + y_list ** 2)
        list_tmp = np.asarray(list(zip(radi_list, dist_list)))
        list_data.extend(list_tmp)
    list_data = np.asarray(list_data)
    return list_data[list_data[:, 0].argsort()]


def check_distortion(list_data):
    """
    Check if the distortion is significant or not. If the number of dots
    having the residual greater than 1 pixel is greater than 15% of the total
    number of dots, there's distortion.

    Parameters
    ----------
    list_data : array_like
        List of [radius, residual] of each dot.

    Returns
    -------
    bool
    """
    check = False
    res_list = list_data[:, 1]
    perc_err = (
            1.0 * len(res_list[res_list > 1.0]) / len(res_list))
    if perc_err > 0.15:
        check = True
    return check
