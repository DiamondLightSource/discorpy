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
# Contributors:
# ============================================================================

"""
Module of pre-processing methods:
- Normalize, binarize an image.
- Determine the median dot-size, median distance between two nearest dots,
  and the slopes of grid-lines of a dot-pattern image.
- Remove non-dot objects or misplaced dots.
- Group dot-centroids into horizontal lines and vertical lines.

"""

import numpy as np
from scipy import ndimage as ndi
from scipy.fftpack import fft2, ifft2
from scipy.ndimage.measurements import center_of_mass
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
import skimage.morphology as morph
import skimage.measure as meas
from skimage.transform import radon


def normalization(mat, size=51):
    """
    Correct a non-uniform background of an image using the median filter.

    Parameters
    ----------
    mat : array_like
        2D array.
    size : int
        Size of the median filter.

    Returns
    -------
    array_like
        2D array. Corrected background.
    """
    mat_bck = ndi.median_filter(mat, size, mode="reflect")
    mean_val = np.mean(mat_bck)
    try:
        mat_cor = mean_val * mat / mat_bck
    except ZeroDivisionError:
        mat_bck[mat_bck == 0.0] = mean_val
        mat_cor = mean_val * mat / mat_bck
    return mat_cor


def _make_window(height, width, sigma=10):
    """
    Create a 2D Gaussian window.

    Parameters
    ----------
    height : int
        Height of the window.
    width : int
        Width of the window.
    sigma : int
        Sigma of the Gaussian window.

    Returns
    -------
    array_like
        2D array.
    """
    xcenter = (width - 1.0) / 2.0
    ycenter = (height - 1.0) / 2.0
    y, x = np.ogrid[-ycenter:height - ycenter, -xcenter:width - xcenter]
    num = 2.0 * sigma * sigma
    window = np.exp(-(x * x / num + y * y / num))
    return window


def _apply_fft_filter(mat, sigma, pad):
    """
    Apply a Fourier Gaussian filter.

    Parameters
    ----------
    mat : array_like
        2D array.
    sigma : int
        Sigma of the Gaussian filter.
    pad : int
        Pad width.

    Returns
    -------
    array_like
        2D array. Filtered image.
    """
    mat = np.pad(mat, ((pad, pad), (pad, pad)), mode='symmetric')
    (height, width) = mat.shape
    window = _make_window(height, width, sigma)
    xlist = np.arange(0, width)
    ylist = np.arange(0, height)
    x, y = np.meshgrid(xlist, ylist)
    matsign = np.power(-1.0, x + y)
    mat = np.real(ifft2(fft2(mat * matsign) * window) * matsign)
    return mat[pad:height - pad, pad:width - pad]


def normalization_fft(mat, sigma=10, pad=30):
    """
    Correct a non-uniform background image using a Fourier Gaussian filter.

    Parameters
    ----------
    mat : array_like
        2D array.
    sigma : int
        Sigma of the Gaussian.
    pad : int
        Pad width.

    Returns
    -------
    array_like
        2D array. Corrected background image.
    """
    mat_bck = _apply_fft_filter(mat, sigma, pad)
    mean_val = np.mean(mat_bck)
    try:
        mat_cor = mean_val * mat / mat_bck
    except ZeroDivisionError:
        mat_bck[mat_bck == 0.0] = mean_val
        mat_cor = mean_val * mat / mat_bck
    return mat_cor


def _select_roi(mat, ratio, square=False):
    """
    Select ROI around the middle of an image.

    Parameters
    ----------
    mat : array_like
        2D array.
    ratio : float
        Ratio between the ROI size and the image size.
    square : bool, optional
        To get a square area or not.

    Returns
    -------
    array_like
        2D array.
    """
    (height, width) = mat.shape
    ratio = np.clip(ratio, 0.05, 1.0)
    if square is True:
        c_hei = height // 2
        c_wid = width // 2
        radi = int(ratio * min(height, width)) // 2
        mat_roi = mat[c_hei - radi:c_hei + radi, c_wid - radi: c_wid + radi]
    else:
        depad_hei = int((height - ratio * height) / 2)
        depad_wid = int((width - ratio * width) / 2)
        mat_roi = mat[depad_hei:height - depad_hei, depad_wid:width - depad_wid]
    return mat_roi


def _invert_dots_contrast(mat):
    """
    Invert the contrast of a 2D binary array to make sure that dots are white.

    Parameters
    ----------
    mat : array_like
        2D binary array.

    Returns
    -------
    array_like
        2D array.
    """
    (height, width) = mat.shape
    ratio = np.sum(mat) / (height * width)
    max_val = np.max(mat)
    if ratio > 0.5:
        mat = max_val - mat
    return mat


def binarization(mat, ratio=0.3, thres=None, denoise=True):
    """
    Apply a list of operations: binarizing an 2D array; inverting the contrast
    of dots if needs to; removing border components; cleaning salty noise; and
    filling holes.

    Parameters
    ----------
    mat : array_like
        2D array.
    ratio : float
        Used to select the ROI around the middle of the image for calculating
        threshold.
    thres : float, optional
        Threshold for binarizing. Automatically calculated if None.
    denoise : bool, optional
        Apply denoising to the image if True.

    Returns
    -------
    array_like
        2D binary array.
    """
    if denoise:
        mat = ndi.median_filter(np.abs(mat), (2, 2))
    if thres is None:
        thres = threshold_otsu(_select_roi(mat, ratio), nbins=512)
    mat = np.asarray(mat > thres, dtype=np.float32)
    mat = _invert_dots_contrast(mat)
    mat = clear_border(mat)
    mat = morph.opening(mat, morph.disk(1))
    mat = np.int16(ndi.binary_fill_holes(mat))
    return mat


def check_num_dots(mat):
    """
    Check if the number of dots is not enough for parabolic fit.

    Parameters
    ----------
    mat : array_like
        2D binary array.

    Returns
    -------
    bool
        True means not enough.
    """
    check = False
    _, num_dots = ndi.label(mat)
    if num_dots < 5 * 5:
        print("WARNING!!! Number of detected dots: {}".format(num_dots))
        print("is not enough for the algorithm to work!")
        check = True
    return check


def calc_size_distance(mat, ratio=0.3):
    """
    Find the median size of dots and the median distance between two nearest
    dots.

    Parameters
    ----------
    mat : array_like
        2D binary array.
    ratio : float
        Used to select the ROI around the middle of an image.

    Returns
    -------
    dot_size : float
        Median size of the dots.
    dot_dist : float
        Median distance between two nearest dots.
    """
    mat = _select_roi(mat, ratio)
    mat = np.int16(clear_border(mat))
    mat_label, num_dots = ndi.label(mat)
    list_index = np.arange(1, num_dots + 1)
    list_sum = ndi.measurements.sum(mat, labels=mat_label, index=list_index)
    dot_size = np.median(list_sum)
    list_cent = np.asarray(
        center_of_mass(mat, labels=mat_label, index=list_index))
    list_dist = [np.sort(np.sqrt((dot[0] - list_cent[:, 0]) ** 2
                                 + (dot[1] - list_cent[:, 1]) ** 2))[1]
                 for dot in list_cent]
    dot_dist = np.median(list_dist)
    return dot_size, dot_dist


def _check_dot_size(mat, min_size, max_size):
    """
    Check if the size of a dot is in a range.

    Parameters
    ----------
    mat : array_like
        2D binary array.
    min_size : float
        Minimum size.
    max_size : float
        Maximum size.

    Returns
    -------
    bool
    """
    check = False
    dot_size = mat.sum()
    if (dot_size >= min_size) and (dot_size <= max_size):
        check = True
    return check


def select_dots_based_size(mat, dot_size, ratio=0.3):
    """
    Select dots having a certain size.

    Parameters
    ----------
    mat : array_like
        2D binary array.
    dot_size : float
        Size of the standard dot.
    ratio : float
        Used to calculate the acceptable range.
        [dot_size - ratio*dot_size; dot_size + ratio*dot_size]

    Returns
    -------
    array_like
        2D array. Selected dots.
    """
    min_size = np.clip(dot_size - ratio * dot_size, 0, None)
    max_size = dot_size + ratio * dot_size
    mat_label, _ = ndi.label(np.int16(mat))
    list_dots = ndi.find_objects(mat_label)
    dots_selected = [dot for dot in list_dots
                     if _check_dot_size(mat[dot], min_size, max_size)]
    mat1 = np.zeros_like(mat, dtype=np.int16)
    for _, j in enumerate(dots_selected):
        mat1[j] = mat[j]
    return mat1


def _check_axes_ratio(mat, ratio):
    """
    Check if the ratio of the axes length of a fitted ellipse is smaller than a
    threshold.

    Parameters
    ----------
    mat : array_like
        2D binary array.
    ratio : float
        Threshold value.

    Returns
    -------
    bool
    """
    check = False
    (height, width) = mat.shape
    if (height < 2) or (width < 2):
        return check
    minor_axis = meas.regionprops(mat)[0].minor_axis_length
    major_axis = meas.regionprops(mat)[0].major_axis_length
    if (height > 1) and (width > 1):
        val = np.abs(major_axis / minor_axis - 1.0)
        if val < ratio:
            check = True
    else:
        check = False
    return check


def select_dots_based_ratio(mat, ratio=0.3):
    """
    Select dots having the ratio between the axes length of the fitted ellipse
    smaller than a threshold.

    Parameters
    ----------
    mat : array_like
        2D binary array.
    ratio : float
        Threshold value.

    Returns
    -------
    array_like
        2D array. Selected dots.
    """
    mat = np.int16(mat)
    mat_label, _ = ndi.label(mat)
    list_dots = ndi.find_objects(mat_label)
    dots_selected = [dot for dot in list_dots
                     if _check_axes_ratio(mat[dot], ratio)]
    mat1 = np.zeros_like(mat)
    for _, j in enumerate(dots_selected):
        mat1[j] = mat[j]
    return mat1


def select_dots_based_distance(mat, dot_dist, ratio=0.3):
    """
    Select dots having a certain range of distance to theirs neighbouring dots.

    Parameters
    ----------
    mat : array_like
        2D array.
    dot_dist : float
        Median distance of two nearest dots.
    ratio : float
        Used to calculate acceptable range.

    Returns
    -------
    array_like
        2D array. Selected dots.
    """
    mat = np.int16(mat)
    mat_label, num_dots = ndi.label(mat)
    list_dots = ndi.find_objects(mat_label)
    list_index = np.arange(1, num_dots + 1)
    list_cent = np.asarray(
        center_of_mass(mat, labels=mat_label, index=list_index))
    list_dist = np.asarray([np.sort(np.sqrt(
        (dot[0] - list_cent[:, 0]) ** 2 + (dot[1] - list_cent[:, 1]) ** 2))[1:4]
                            for dot in list_cent])
    mat1 = np.zeros_like(mat)
    for i, j in enumerate(list_dots):
        dist = list_dist[i]
        num = dist // dot_dist
        dist_error = (dist - num * dot_dist) / dot_dist
        if any(dist_error < ratio):
            mat1[j] = mat[j]
    return mat1


def calc_hor_slope(mat, ratio=0.3):
    """
    Calculate the slope of horizontal lines against the horizontal axis.

    Parameters
    ----------
    mat : array_like
        2D binary array.
    ratio : float
        Used to select the ROI around the middle of an image.

    Returns
    -------
    float
        Horizontal slope of the grid.
    """
    coarse_range = 30.0  # Degree
    radi = np.pi / 180.0
    mat = np.int16(clear_border(_select_roi(mat, ratio)))
    (height, width) = mat.shape
    list_angle = 90.0 + np.arange(-coarse_range, coarse_range + 1.0)
    projections = radon(np.float32(mat), theta=list_angle, circle=False)
    list_max = np.amax(projections, axis=0)
    best_angle = -(list_angle[np.argmax(list_max)] - 90.0)
    dist_error = 0.5 * width * (np.tan(radi) / np.cos(best_angle * radi))
    mat_label, num_dots = ndi.label(mat)
    list_index = np.arange(1, num_dots + 1)
    list_cent = np.asarray(
        center_of_mass(mat, labels=mat_label, index=list_index))
    list_cent = - list_cent  # For coordinate consistency
    mean_x = np.mean(list_cent[:, 1])
    mean_y = np.mean(list_cent[:, 0])
    index_mid_dot = np.argsort(np.sqrt((mean_x - list_cent[:, 1]) ** 2
                                       + (mean_y - list_cent[:, 0]) ** 2))[0]
    used_dot = list_cent[index_mid_dot]
    line_slope = np.tan(best_angle * radi)
    list_tmp = np.sqrt(
        np.ones(num_dots, dtype=np.float32) * line_slope ** 2 + 1.0)
    list_tmp2 = used_dot[0] * np.ones(
        num_dots, dtype=np.float32) - line_slope * used_dot[1]
    list_dist = np.abs(
        line_slope * list_cent[:, 1] - list_cent[:, 0] + list_tmp2) / list_tmp
    dots_selected = np.asarray(
        [dot for i, dot in enumerate(list_cent) if list_dist[i] < dist_error])
    if len(dots_selected) > 1:
        slope = np.polyfit(dots_selected[:, 1], dots_selected[:, 0], 1)[0]
    else:
        slope = line_slope
    return slope


def calc_ver_slope(mat, ratio=0.3):
    """
    Calculate the slope of vertical lines against the vertical axis.

    Parameters
    ----------
    mat : array_like
        2D binary array.
    ratio : float
        Used to select the ROI around the middle of a image.

    Returns
    -------
    float
        Vertical slope of the grid.
    """
    coarse_range = 30.0
    radi = np.pi / 180.0
    mat = np.int16(clear_border(_select_roi(mat, ratio)))
    (height, width) = mat.shape
    list_angle = np.arange(-coarse_range, coarse_range + 1.0)
    projections = radon(np.float32(mat), theta=list_angle, circle=False)
    list_max = np.amax(projections, axis=0)
    best_angle = (list_angle[np.argmax(list_max)])
    dist_error = 0.5 * width * np.tan(radi) / np.cos(best_angle * radi)
    mat_label, num_dots = ndi.label(mat)
    list_index = np.arange(1, num_dots + 1)
    list_cent = np.fliplr(
        np.asarray(center_of_mass(mat, labels=mat_label, index=list_index)))
    mean_x = np.mean(list_cent[:, 1])
    mean_y = np.mean(list_cent[:, 0])
    index_mid_dot = np.argsort(np.sqrt((mean_x - list_cent[:, 1]) ** 2
                                       + (mean_y - list_cent[:, 0]) ** 2))[0]
    used_dot = list_cent[index_mid_dot]
    line_slope = np.tan(best_angle * radi)
    list_tmp = np.sqrt(
        np.ones(num_dots, dtype=np.float32) * line_slope ** 2 + 1.0)
    list_tmp2 = used_dot[0] * np.ones(
        num_dots, dtype=np.float32) - line_slope * used_dot[1]
    list_dist = np.abs(
        line_slope * list_cent[:, 1] - list_cent[:, 0] + list_tmp2) / list_tmp
    dots_selected = np.asarray(
        [dot for i, dot in enumerate(list_cent) if list_dist[i] < dist_error])
    if len(dots_selected) > 1:
        slope = np.polyfit(dots_selected[:, 1], dots_selected[:, 0], 1)[0]
    else:
        slope = line_slope
    return slope


def _check_dot_on_line(dot1, dot2, slope, dot_dist, ratio, num_dot_miss):
    """
    Check if dot1 and dot2 belong to the same group.

    Parameters
    ----------
    dot1 : list of float
        Coordinate of dot1.
    dot2 : list of float
        Coordinate of dot2.
    slope : float
        Slope of the line of dots.
    dot_dist : float
        Median distance of two nearest dots.
    ratio : float
        Acceptable variation.
    num_dot_miss : int
        Acceptable missing dots between dot1 and dot2.

    Returns
    -------
    bool
    """
    check = False
    dist_error = ratio * dot_dist
    search_dist = num_dot_miss * dot_dist
    if len(dot1) == 2 and len(dot2) == 2:
        xmin = dot1[1] - search_dist
        xmax = dot1[1] + search_dist
        if xmin < dot2[1] < xmax:
            ntemp1 = np.sqrt(slope * slope + 1.0)
            ntemp2 = dot1[0] - slope * dot1[1]
            dist_d12 = np.abs(slope * dot2[1] - dot2[0] + ntemp2) / ntemp1
            if dist_d12 < dist_error:
                check = True
    else:
        raise ValueError("Invalid input!!!")
    return check


def group_dots_hor_lines(mat, slope, dot_dist, ratio=0.3, num_dot_miss=6,
                         accepted_ratio=0.65):
    """
    Group dots into horizontal lines.

    Parameters
    ----------
    mat : array_like
        A binary image or a list of (y,x)-coordinates of points.
    slope : float
        Horizontal slope of the grid.
    dot_dist : float
        Median distance of two nearest dots.
    ratio : float
        Acceptable variation.
    num_dot_miss : int
        Acceptable missing dots between dot1 and dot2.
    accepted_ratio : float
        Use to select lines having the number of dots equal to or larger than
        the multiplication of the `accepted_ratio` and the maximum number of
        dots per line.

    Returns
    -------
    list of array_like
        List of 2D arrays. Each list is the coordinates (y, x) of dot-centroids
        belong to the same group. Length of each list may be different.
    """
    mat = np.asarray(mat)
    if mat.shape[1] != 2:
        mat_label, num_dots = ndi.label(np.int16(mat))
        list_dots = np.copy(center_of_mass(
            mat, labels=mat_label, index=np.arange(1, num_dots + 1)))
    else:
        list_dots = mat
        num_dots = len(list_dots)
    num_dots_left = num_dots
    list_dots_left = np.copy(list_dots)
    list_dots_left = list_dots_left[list_dots_left[:, 1].argsort()]
    list_lines = []
    while num_dots_left > 1:
        dot1 = list_dots_left[0]
        dots_selected = np.asarray([dot1])
        pos_get = [0]
        for i in range(1, len(list_dots_left)):
            dot2 = list_dots_left[i]
            check = _check_dot_on_line(
                dot1, dot2, slope, dot_dist, ratio, num_dot_miss)
            if check:
                dot1 = dot2
                dots_selected = np.vstack((dots_selected, dot2))
                pos_get.append(i)
        list_pos = np.arange(0, len(list_dots_left), dtype=np.int32)
        pos_get = np.asarray(pos_get, dtype=np.int32)
        list_pos_left = np.asarray(
            [pos for pos in list_pos if pos not in pos_get], dtype=np.int32)
        list_dots_left = list_dots_left[list_pos_left]
        num_dots_left = len(list_dots_left)
        if len(dots_selected) > 1:
            list_lines.append(dots_selected)
    list_len = [len(i) for i in list_lines]
    len_accepted = np.int16(accepted_ratio * np.max(list_len))
    lines_selected = [line for line in list_lines if len(line) > len_accepted]
    lines_selected = sorted(
        lines_selected, key=lambda list_: np.mean(list_[:, 0]))
    return lines_selected


def group_dots_ver_lines(mat, slope, dot_dist, ratio=0.3, num_dot_miss=6,
                         accepted_ratio=0.65):
    """
    Group dots into vertical lines.

    Parameters
    ----------
    mat : array_like
        A binary image or a list of (y,x)-coordinates of points.
    slope : float
        Vertical slope of the grid.
    dot_dist : float
        Median distance of two nearest dots.
    ratio : float
        Acceptable variation.
    num_dot_miss : int
        Acceptable missing dots between dot1 and dot2.
    accepted_ratio : float
        Use to select lines having the number of dots equal to or larger than
        the multiplication of the `accepted_ratio` and the maximum number of
        dots per line.

    Returns
    -------
    list of array_like
        List of 2D arrays. Each list is the coordinates (y, x) of dot-centroids
        belong to the same group. Length of each list may be different.
    """
    mat = np.asarray(mat)
    if mat.shape[1] != 2:
        mat_label, num_dots = ndi.label(np.int16(mat))
        list_dots = np.copy(center_of_mass(
            mat, labels=mat_label, index=np.arange(1, num_dots + 1)))
    else:
        list_dots = mat
        num_dots = len(list_dots)
    list_dots = np.fliplr(list_dots)  # Swap the coordinates
    num_dots_left = num_dots
    list_dots_left = np.copy(list_dots)
    list_dots_left = list_dots_left[list_dots_left[:, 1].argsort()]
    list_lines = []
    while num_dots_left > 1:
        dot1 = list_dots_left[0]
        dots_selected = np.asarray([dot1])
        pos_get = [0]
        for i in range(1, len(list_dots_left)):
            dot2 = list_dots_left[i]
            check = _check_dot_on_line(
                dot1, dot2, slope, dot_dist, ratio, num_dot_miss)
            if check:
                dot1 = dot2
                dots_selected = np.vstack((dots_selected, dot2))
                pos_get.append(i)

        list_pos = np.arange(0, len(list_dots_left), dtype=np.int32)
        pos_get = np.asarray(pos_get, dtype=np.int32)
        list_pos_left = np.asarray(
            [pos for pos in list_pos if pos not in pos_get], dtype=np.int32)
        list_dots_left = list_dots_left[list_pos_left]
        num_dots_left = len(list_dots_left)
        if len(dots_selected) > 1:
            dots_selected = np.fliplr(dots_selected)  # Swap back
            list_lines.append(dots_selected)
    list_length = [len(i) for i in list_lines]
    len_accepted = np.int16(accepted_ratio * np.max(list_length))
    lines_selected = [line for line in list_lines if len(line) > len_accepted]
    lines_selected = sorted(
        lines_selected, key=lambda list_: np.mean(list_[:, 1]))
    return lines_selected


def remove_residual_dots_hor(list_lines, slope, residual=2.5):
    """
    Remove dots having distances larger than a certain value from fitted
    horizontal parabolas.

    Parameters
    ----------
    list_lines : list of array_like
        List of the coordinates of dot-centroids on horizontal lines.
    slope : float
        Horizontal slope of the grid.
    residual : float
        Acceptable distance in pixel unit between a dot and a fitted parabola.

    Returns
    -------
    list of array_like
        List of 2D arrays. Each list is the coordinates (y, x) of dot-centroids
        belong to the same group. Length of each list may be different.
    """
    list_lines2 = []
    for i, list_ in enumerate(list_lines):
        (a2, a1, a0) = np.polyfit(list_[:, 1], list_[:, 0], 2)
        error = np.abs(
            ((a2 * list_[:, 1] ** 2 + a1 * list_[:, 1] + a0)
             - list_[:, 0]) * np.cos(np.arctan(slope)))
        dots_left = np.asarray(
            [dot for i, dot in enumerate(list_) if error[i] < residual])
        list_lines2.append(dots_left)
    return list_lines2


def remove_residual_dots_ver(list_lines, slope, residual=2.5):
    """
    Remove dots having distances larger than a certain value from fitted
    vertical parabolas.

    Parameters
    ---------
    list_lines : list of float
        List of the coordinates of the dot-centroids on the vertical lines.
    slope : float
        Slope of the vertical line.
    residual : float
        Acceptable distance in pixel unit between the dot and the fitted
        parabola.

    Returns
    -------
    list of float
        List of 2D array. Each list is the coordinates (y, x) of dot-centroids
        belong to the same group. Length of each list may be different.
    """
    list_lines2 = []
    for i, list_ in enumerate(list_lines):
        list_ = np.fliplr(list_)  # Swap the coordinates
        (a2, a1, a0) = np.polyfit(list_[:, 1], list_[:, 0], 2)
        error = np.abs(
            ((a2 * list_[:, 1] ** 2 + a1 * list_[:, 1] + a0)
             - list_[:, 0]) * np.cos(np.arctan(slope)))
        dots_left = np.asarray(
            [dot for i, dot in enumerate(list_) if error[i] < residual])
        dots_left = np.fliplr(dots_left)  # Swap back
        list_lines2.append(dots_left)
    return list_lines2


def _make_circle_mask(width, ratio):
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
