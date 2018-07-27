#============================================================================
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
#============================================================================
# Author: Nghia T. Vo
# E-mail: nghia.vo@diamond.ac.uk
# Description: Python implementation (2.7) of the author's methods of
# distortion correction, Nghia T. Vo et al "Radial lens distortion
# correction with sub-pixel accuracy for X-ray micro-tomography"
# Optics Express 23, 32859-32868 (2015), https://doi.org/10.1364/OE.23.032859
# Publication date: 10th July 2018
#============================================================================

"""
Module of post-processing methods:
- Unwarp a line of dots, an image.
- Generate unwarped slices of a 3D dataset.
- Calculate the residual of the corrected dots.
"""
import numpy as np
import cv2
from scipy import interpolate
from scipy import optimize


def unwarp_line_forward(list_lines, xcenter, ycenter, list_fact):
    """
    Unwarp lines of dot-centroids using the forward model.
    ---------
    Parameters: - list_lines: List of the coordinates of dot-centroids
                            on the lines.
                - list_fact: Polynomial coefficients of the forward model.
    ---------
    Return:     - list_ulines: List of the corrected coordinates of
                             dot-centroids on the lines.
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
    return (rd - ru * np.sum(np.asarray([fact * ru**i for i,
                                         fact in enumerate(list_fact)])))**2


def unwarp_line_backward(list_lines, xcenter, ycenter, list_fact):
    """
    Unwarp lines of dot-centroids using the forward model.
    ---------
    Parameters: - list_lines: List of the coordinates of dot-centroids
                            on the lines.
                - list_fact: Polynomial coefficients of the backward model.
    ---------
    Return:     - list_ulines: List of the corrected coordinates of
                             dot-centroids on the lines.
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


def _mapping(mat, xmat, ymat):
    """
    Apply a geometric transformation to a 2D array using Opencv
    ---------
    Parameters: - mat: 2D array.
                - xmat: 2D array of the x-coordinates.
                - ymat: 2D array of the y-coordinates.
    ---------
    Return:     - 2D array. 
    """
    mat = cv2.remap(mat, xmat, ymat, interpolation=cv2.INTER_LINEAR)
    return mat


def unwarp_image_backward_cv(mat, xcenter, ycenter, list_fact):
    """
    Unwarp a 2D array using the backward model. Use Opencv library
    for fast performance.
    ---------
    Parameters: - mat: 2D array.
                - xcenter: Center of distortion in x-direction.
                - ycenter: Center of distortion in y-direction.
                - list_fact: Polynomial coefficients of the backward model.
    ---------
    Return:     - 2D array, distortion corrected.
    """
    (height, width) = mat.shape
    xu_list = np.arange(width) - xcenter
    yu_list = np.arange(height) - ycenter
    xu_mat, yu_mat = np.meshgrid(xu_list, yu_list)
    ru_mat = np.sqrt(xu_mat**2 + yu_mat**2)
    fact_mat = np.sum(
        np.asarray([factor * ru_mat**i for i,
                    factor in enumerate(list_fact)]), axis=0)
    xd_mat = np.float32(np.clip(xcenter + fact_mat * xu_mat, 0, width - 1))
    yd_mat = np.float32(np.clip(ycenter + fact_mat * yu_mat, 0, height - 1))
    mat = _mapping(mat, xd_mat, yd_mat)
    return mat


def unwarp_image_backward(mat, xcenter, ycenter, list_fact):
    """
    Unwarp a 2D array using the backward model.
    ---------
    Parameters: - mat: 2D array.
                - xcenter: Center of distortion in x-direction.
                - ycenter: Center of distortion in y-direction.
                - list_fact: Polynomial coefficients of the backward model.
    ---------
    Return:     - 2D array, distortion corrected.
    """
    (height, width) = mat.shape
    xu_list = np.arange(width) - xcenter
    yu_list = np.arange(height) - ycenter
    xu_mat, yu_mat = np.meshgrid(xu_list, yu_list)
    ru_mat = np.sqrt(xu_mat**2 + yu_mat**2)
    fact_mat = np.sum(
        np.asarray([factor * ru_mat**i for i,
                    factor in enumerate(list_fact)]), axis=0)
    xd_mat = np.float32(np.clip(xcenter + fact_mat * xu_mat, 0, width - 1))
    yd_mat = np.float32(np.clip(ycenter + fact_mat * yu_mat, 0, height - 1))
    finter = interpolate.RectBivariateSpline(
        yu_list + ycenter, xu_list + xcenter, mat, kx=1, ky=1)
    mat = finter.ev(yd_mat, xd_mat)
    return mat


def unwarp_image_forward(mat, xcenter, ycenter, list_fact):
    """
    Unwarp a 2D array using the forward model.
    Should be used only for testing due to the problem of vacant pixels.
    ---------
    Parameters: - mat: 2D array.
                - xcenter: Center of distortion in x-direction.
                - ycenter: Center of distortion in y-direction.
                - list_fact: Polynomial coefficients of the forward model.
    ---------
    Return:     - 2D array, distortion corrected.
    """
    (height, width) = mat.shape
    xd_list = np.arange(width) - xcenter
    yd_list = np.arange(height) - ycenter
    xd_mat, yd_mat = np.meshgrid(xd_list, yd_list)
    rd_mat = np.sqrt(xd_mat**2 + yd_mat**2)
    fact_mat = np.sum(
        np.asarray([factor * rd_mat**i for i,
                    factor in enumerate(list_fact)]), axis=0)
    xu_mat = np.int16(
        np.round(np.clip(xcenter + fact_mat * xd_mat, 0, width - 1)))
    yu_mat = np.int16(
        np.round(np.clip(ycenter + fact_mat * yd_mat, 0, height - 1)))
    mat_unw = np.zeros_like(mat)
    mat_norm = np.zeros_like(mat)
    for i in range(height):
        for j in range(width):
            posi = yu_mat[i, j]
            posj = xu_mat[i, j]
            mat_unw[posi, posj] += mat[i, j]
            mat_norm[posi, posj] += 1
    mat_norm[mat_norm == 0.0] = 1
    return mat_unw / mat_norm


def unwarp_slice_backward(mat3D, xcenter, ycenter, list_fact, index):
    """
    Generate an unwarped slice [:,index.:] of a 3D dataset, i.e
    one unwarped sinogram of a 3D tomographic data.
    ---------
    Parameters: - mat3D: 3D array.
                - xcenter: Center of distortion in x-direction.
                - ycenter: Center of distortion in y-direction.
                - list_fact: Polynomial coefficients of the backward model.
                - index: Index of the slice
    ---------
    Return:     - 2D array, distortion corrected.
    """
    if (len(mat3D.shape) < 3):
        raise ValueError("Input must be a 3D data")
    (depth, height, width) = mat3D.shape
    xu_list = np.arange(0, width) - xcenter
    yu = index - ycenter
    ru_list = np.sqrt(xu_list**2 + yu**2)
    flist = np.sum(
        np.asarray([factor * ru_list**i for i,
                    factor in enumerate(list_fact)]), axis=0)
    xd_list = np.clip(xcenter + flist * xu_list, 0, width - 1)
    yd_list = np.clip(ycenter + flist * yu, 0, height - 1)
    yd_min = np.int16(np.floor(np.amin(yd_list)))
    yd_max = np.int16(np.ceil(np.amax(yd_list))) + 1
    xlist = np.arange(0, width)
    ylist = np.arange(yd_min, yd_max)
    sino = np.zeros((depth, width), dtype=np.float32)
    for i in range(depth):
        finter = interpolate.RectBivariateSpline(
            ylist, xlist, mat3D[i, yd_min:yd_max, :], kx=1, ky=1)
        sino[i] = finter.ev(yd_list, xd_list)
    mean_val = np.mean(sino)
    sino[sino == 0.0] = mean_val
    return sino


def unwarp_chunk_slices_backward(mat3D, xcenter, ycenter, list_fact,
                                 start_index, stop_index):
    """
    Generate a chunk  of unwarped slices [:,start_index: stop_index, :].
    Useful for correcting 3D tomographic data.
    ---------
    Parameters: - mat3D: 3D array.
                - xcenter: Center of distortion in x-direction.
                - ycenter: Center of distortion in y-direction.
                - list_fact: Polynomial coefficients of the backward model.
                - start_index: Starting index of the slice.
                - stop_index: Stopping index of the slice.
    ---------
    Returns:    - 3D array, distortion corrected.
    """
    if (len(mat3D.shape) < 3):
        raise ValueError("Input must be a 3D data")
    (depth, height, width) = mat3D.shape
    index_list = np.arange(depth, dtype=np.int16)
    if stop_index == -1:
        stop_index = depth
    if (start_index not in index_list) or (stop_index not in index_list):
        raise ValueError("Selected index is out of the range")
    xu_list = np.arange(0, width) - xcenter
    yu1 = start_index - ycenter
    ru_list = np.sqrt(xu_list**2 + yu1**2)
    flist = np.sum(
        np.asarray([factor * ru_list**i for i,
                    factor in enumerate(list_fact)]), axis=0)
    yd_list1 = np.clip(ycenter + flist * yu1, 0, height - 1)
    yu2 = stop_index - ycenter
    ru_list = np.sqrt(xu_list**2 + yu2**2)
    flist = np.sum(
        np.asarray([factor * ru_list**i for i,
                    factor in enumerate(list_fact)]), axis=0)
    yd_list2 = np.clip(ycenter + flist * yu2, 0, height - 1)
    yd_min = np.int16(np.floor(np.amin(yd_list1)))
    yd_max = np.int16(np.ceil(np.amax(yd_list2))) + 1
    yu_list = np.arange(start_index, stop_index + 1) - ycenter
    xu_mat, yu_mat = np.meshgrid(xu_list, yu_list)
    ru_mat = np.sqrt(xu_mat**2 + yu_mat**2)
    fact_mat = np.sum(
        np.asarray([factor * ru_mat**i for i,
                    factor in enumerate(list_fact)]), axis=0)
    xd_mat = np.float32(np.clip(xcenter + fact_mat * xu_mat, 0, width - 1))
    yd_mat = np.float32(
        np.clip(ycenter + fact_mat * yu_mat, 0, height - 1)) - yd_min
    sino_chunk = np.asarray(
        [_mapping(mat3D[i, yd_min: yd_max, :],
                  xd_mat, yd_mat) for i in range(depth)])
    return sino_chunk


def calc_residual_hor(list_ulines, xcenter, ycenter):
    """
    Calculate the distances of corrected dots (on the horizontal lines)
    to fitted straight lines.
    Useful to check the straightness of the unwarped lines.     
    ---------
    Parameters: - list_ulines: List of the unwarped horizontal lines.
                - xcenter: Center of distortion in x-direction.
                - ycenter: Center of distortion in y-direction.
    ---------
    Returns:    - 2D array: each element has two values: 1) Distance
                of a dot to the center of distortion; 2) Distance
                of a dot to the fitted straight line. 
    """
    list_data = []
    for i, line in enumerate(list_ulines):
        y_list = line[:, 0] - ycenter
        x_list = line[:, 1] - xcenter
        (a_fact, b_fact) = np.polyfit(x_list,  y_list, 1)
        dist_list = np.abs(
            a_fact * x_list - y_list + b_fact) / np.sqrt(a_fact**2 + 1)
        radi_list = np.sqrt(x_list**2 + y_list**2)
        list_tmp = np.asarray(zip(radi_list, dist_list))
        list_data.extend(list_tmp)
    list_data = np.asarray(list_data)
    return list_data[list_data[:, 0].argsort()]


def calc_residual_ver(list_ulines, xcenter, ycenter):
    """
    Calculate the distances of corrected dots (on the vertical lines)
    to fitted straight lines.
    Useful to check the straightness of the unwarped lines.     
    ---------
    Parameters: - list_ulines: List of the unwarped vertical lines.
                - xcenter: Center of distortion in x-direction.
                - ycenter: Center of distortion in y-direction.
    ---------
    Returns:    - 2D array: each element has two values: 1) Distance
                of a dot to the center of distortion; 2) Distance
                of a dot to the fitted straight line. 
    """
    list_data = []
    for i, line in enumerate(list_ulines):
        y_list = line[:, 0] - ycenter
        x_list = line[:, 1] - xcenter
        (a_fact, b_fact) = np.polyfit(y_list,  x_list, 1)
        dist_list = np.abs(
            a_fact * y_list - x_list + b_fact) / np.sqrt(a_fact**2 + 1)
        radi_list = np.sqrt(x_list**2 + y_list**2)
        list_tmp = np.asarray(zip(radi_list, dist_list))
        list_data.extend(list_tmp)
    list_data = np.asarray(list_data)
    return list_data[list_data[:, 0].argsort()]


def check_distortion(list_data):
    """
    Check if the distortion is significant or not.
    If the number of dots having the residual greater than 1 pixel
    is greater than 15% of the total number of dots, there's distortion.     
    ---------
    Parameters: - list_data: List of [radius, residual] of the dots.             
    ---------
    Returns:    - Boolean.
    """
    check = False
    res_list = list_data[:, 1]
    perc_err = (
        1.0 * len(res_list[res_list > 1.0]) / len(res_list))
    if perc_err > 0.15:
        check = True
    return check
