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
# E-mail: 
# Publication date: 10th July 2018
# ============================================================================
# Contributors:
# ============================================================================

"""
Module of pre-processing methods:

-   Normalize, binarize an image.
-   Determine the median dot-size, median distance between two nearest dots,
    and the slopes of grid-lines of a dot-pattern image.
-   Remove non-dot objects or misplaced dots.
-   Group points into horizontal lines and vertical lines.
-   Calculate a threshold value for binarizing.
-   Create a 2D mask with parabolic boundary.
-   Remove points outside a parabolic mask in a 2D space.
-   Extract dot-centroids of a dot-pattern from a binary or grayscale image.
-   Group points on strongly curved lines into horizontal lines and vertical
    lines based on polynomial fit.
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


def _apply_fft_filter(mat, sigma, pad, mode='reflect'):
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
    mat = np.pad(mat, ((pad, pad), (pad, pad)), mode=mode)
    (height, width) = mat.shape
    window = _make_window(height, width, sigma)
    xlist = np.arange(0, width)
    ylist = np.arange(0, height)
    x, y = np.meshgrid(xlist, ylist)
    matsign = np.power(-1.0, x + y)
    mat = np.real(ifft2(fft2(mat * matsign) * window) * matsign)
    return mat[pad:height - pad, pad:width - pad]


def normalization_fft(mat, sigma=10, pad=100, mode='reflect'):
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
    mode : str
        Padding mode.

    Returns
    -------
    array_like
        2D array. Corrected background image.
    """
    mat_bck = _apply_fft_filter(mat, sigma, pad, mode=mode)
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
        mat_roi = mat[depad_hei:height - depad_hei,
                      depad_wid:width - depad_wid]
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
    list_sum = ndi.sum(mat, labels=mat_label, index=list_index)
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
    list_dist = np.asarray([np.sort(
        np.sqrt((dot[0] - list_cent[:, 0]) ** 2 + (
                    dot[1] - list_cent[:, 1]) ** 2))[1:4]
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
        List of 2D arrays. Each array is (y, x)-coordinates of points belong
        to the same group. Length of each list may be different.
    """
    mat = np.asarray(mat)
    if mat.shape[-1] > 2:
        mat_label, num_dots = ndi.label(np.int16(mat))
        list_dots = np.copy(center_of_mass(
            mat, labels=mat_label, index=np.arange(1, num_dots + 1)))
    else:
        list_dots = mat
        num_dots = len(list_dots)
        if num_dots == 0:
            raise ValueError("Input is empty!!!")
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
    len_accepted = int(accepted_ratio * np.max(list_len))
    lines_selected = [line for line in list_lines if len(line) > len_accepted]
    lines_selected = sorted(
        lines_selected, key=lambda list_: np.mean(list_[:, 0]))
    return lines_selected


def group_dots_ver_lines(mat, slope, dot_dist, ratio=0.3, num_dot_miss=6,
                         accepted_ratio=0.75):
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
        List of 2D arrays. Each array is (y,x)-coordinates of points belong
        to the same group. Length of each list may be different.
    """
    mat = np.asarray(mat)
    if mat.shape[-1] > 2:
        mat_label, num_dots = ndi.label(np.int16(mat))
        list_dots = np.copy(center_of_mass(
            mat, labels=mat_label, index=np.arange(1, num_dots + 1)))
    else:
        list_dots = mat
        num_dots = len(list_dots)
        if num_dots == 0:
            raise ValueError("Input is empty!!!")
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
    len_accepted = int(accepted_ratio * np.max(list_length))
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
        if len(dots_left) > 0:
            list_lines2.append(dots_left)
    if len(list_lines2) == 0:
        raise ValueError("No dots left. Check the input or residual !!!")
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
        if len(dots_left) > 0:
            dots_left = np.fliplr(dots_left)  # Swap back
            list_lines2.append(dots_left)
    if len(list_lines2) == 0:
        raise ValueError("No dots left. Check the input or residual !!!")
    return list_lines2


def calculate_threshold(mat, bgr="bright", snr=2.0):
    """
    Calculate a threshold value based on Algorithm 4 in Ref. [1].

    Parameters
    ----------
    mat : array_like
        2D array.
    bgr : {"bright", "dark"}
        To indicate the brightness of the background against image features.
    snr : float
        Ratio (>1.0) used to separate image features against noise. Greater is
        less sensitive.

    Returns
    -------
    float
        Threshold value.

    References
    ----------
        https://doi.org/10.1364/OE.26.028396
    """
    size = max(mat.shape)
    list_sort = np.sort(np.ndarray.flatten(mat))
    list_dsp = ndi.zoom(list_sort, 1.0 * size/len(list_sort), mode='nearest')
    npoint = len(list_dsp)
    xlist = np.arange(0, npoint, 1.0)
    ndrop = int(0.25 * npoint)
    (slope, intercept) = np.polyfit(xlist[ndrop:-ndrop - 1],
                                    list_dsp[ndrop:-ndrop - 1], 1)[:2]
    y_end = intercept + slope * xlist[-1]
    noise_level = np.abs(y_end - intercept)
    if bgr == "bright":
        threshold = intercept - noise_level * snr * 0.5
    else:
        threshold = y_end + noise_level * snr * 0.5
    return threshold


def make_parabola_mask(height, width, hor_curviness=0.3, ver_curviness=0.3,
                       hor_margin=100, ver_margin=100, rotate=0.0):
    """
    Create a 2D mask with parabolic boundary.
    The function generates a mask where the boundary is defined by parabolic
    curves on four sides of a 2D array. The mask can be rotated by a
    specified angle.

    Parameters
    ----------
    height : int
        Height of the mask.
    width : int
        Width of the mask.
    hor_curviness : float, optional
        Horizontal curvature factor. Controls the curvature of the
        parabolas along the left and right boundaries. Larger is stronger.
    ver_curviness : float, optional
        Vertical curvature factor. Controls the curvature of the parabolas
        along the top and bottom boundaries. Larger is stronger.
    hor_margin : int or tuple of int, optional
        Horizontal margins (left and right of an image).
    ver_margin : int or tuple of int, optional
        Vertical margins (top and bottom of an image).
    rotate : float, optional
        Angle (degrees) to rotate the mask.

    Returns
    -------
    mask : array_like
        2D array which values are 1.0 inside the mask and 0.0 outside.
    """
    if isinstance(ver_margin, tuple) or isinstance(ver_margin, list):
        top_margin = ver_margin[0]
        bot_margin = ver_margin[-1]
    else:
        top_margin = bot_margin = ver_margin
    if isinstance(hor_margin, tuple) or isinstance(hor_margin, list):
        left_margin = hor_margin[0]
        right_margin = hor_margin[-1]
    else:
        left_margin = right_margin = hor_margin
    if (left_margin + right_margin) > width:
        raise ValueError("Invalid horizontal margin!!!")
    if (top_margin + bot_margin) > height:
        raise ValueError("Invalid vertical margin!!!")
    y, x = np.ogrid[:height, :width]
    parabola_y = (ver_curviness / width) * (
            x - width / 2) ** 2 + top_margin
    mask_top = np.float32(y > parabola_y)
    parabola_y = -(ver_curviness / width) * (
            x - width / 2) ** 2 + height - bot_margin
    mask_bot = np.float32(y < parabola_y)
    parabola_x = (hor_curviness / height) * (
            y - height / 2) ** 2 + left_margin
    mask_left = np.float32(x > parabola_x)
    parabola_x = -(hor_curviness / height) * (
            y - height / 2) ** 2 + width - right_margin
    mask_right = np.float32(x < parabola_x)
    mask = mask_bot * mask_top * mask_left * mask_right
    if rotate != 0.0:
        mask = np.round(ndi.rotate(mask, rotate, reshape=False))
    return np.float32(mask)


def remove_points_using_parabola_mask(points, height, width, hor_curviness=0.3,
                                      ver_curviness=0.3, hor_margin=100,
                                      ver_margin=100, rotate=0.0):
    """
    Remove points outside a parabolic mask in a 2D space.
    The mask is defined by its dimensions, curvature, margins, and an optional
    rotation angle.

    Parameters
    ----------
    points : array_like
       Array of shape (N, 2) of (y, x)-coordinates of points.
    height : int
        Height of the mask.
    width : int
        Width of the mask.
    hor_curviness : float, optional
        Horizontal curvature factor. Controls the curvature of the
        parabolas along the left and right boundaries. Larger is stronger.
    ver_curviness : float, optional
        Vertical curvature factor. Controls the curvature of the parabolas
        along the top and bottom boundaries. Larger is stronger.
    hor_margin : int or tuple of int, optional
        Horizontal margins (left and right of an image).
    ver_margin : int or tuple of int, optional
        Vertical margins (top and bottom of an image).
    rotate : float, optional
        Angle (degrees) to rotate the mask.

    Returns
    -------
    array_like
        Filtered points. Array of shape (M, 2) of (y, x)-coordinates of points.
    """
    mask = make_parabola_mask(height, width, hor_curviness=hor_curviness,
                              ver_curviness=ver_curviness,
                              hor_margin=hor_margin, ver_margin=ver_margin,
                              rotate=rotate)
    y_cors, x_cors = np.int32(points[:, 0]), np.int32(points[:, 1])
    valid_indices = ((y_cors >= 0) & (y_cors < height) &
                     (x_cors >= 0) & (x_cors < width) &
                     (mask[y_cors, x_cors] == 1.0))
    return points[valid_indices]


def get_points_dot_pattern(mat, binarize=True, ratio=0.3, thres=None):
    """
    Extract dot-centroids of a dot-pattern from a binary or grayscale image.

    Parameters
    ----------
    mat : array_like
        2D array of the dot pattern. If `binary` is False, the input must be
        a binary image (values 0 and 1 only).
    binarize : bool, optional
        To select if the input is a grayscale image that needs to be binarized.
    ratio : float
        Used to select the ROI around the middle of the image for calculating
        threshold.
    thres : float, optional
        Threshold for binarizing. Automatically calculated if None.

    Returns
    -------
    array_like
       Array of shape (N, 2) of the (y, x)-coordinates of dots' center.
    """
    if binarize:
        mat = binarization(mat, ratio=ratio, thres=thres)
    else:
        if np.max(mat) != 1.0 or np.min(mat) != 0.0:
            raise ValueError("Input not a binary image, "
                             "e.i. maximum_value=1 and minimum value=0!!!")
    mat_label, num_dots = ndi.label(np.int16(mat))
    points = ndi.center_of_mass(mat, labels=mat_label,
                                index=np.arange(1, num_dots + 1))
    return np.asarray(points)


def rotate_points(points, angle, degree_unit=True):
    """
    Rotate a set of 2D points by a specified angle.

    Parameters
    ----------
    points : array_like
        Array of shape (N, 2) of (y, x)-coordinates of points.
    angle : float
        Rotation angle in degree or radian. Positive values rotate
        counterclockwise.
    degree_unit : bool, optional
        To select the unit of the input angle.

    Returns
    -------
    array_like
        Array of shape (N, 2) of rotated (y, x)-coordinates of the points.
    """
    points = np.asarray(points)
    if degree_unit:
        angle = np.deg2rad(angle)
    x, y = points[:, 1], points[:, 0]
    xr = x * np.cos(angle) - y * np.sin(angle)
    yr = x * np.sin(angle) + y * np.cos(angle)
    return np.column_stack((yr, xr))


def remove_subset_points(selected_points, points):
    """
    Remove a subset points from a set of points.

    Parameters
    ----------
    selected_points : array_like
        Array of shape (N, 2) of (y, x)-coordinates of points to exclude.
    points : array_like
        Array of shape (M, 2) of (y, x)-coordinates of points.

    Returns
    -------
    array_like
        Array of shape (K, 2) of points after the removal.
    """
    elements_set = set(map(tuple, selected_points))
    filtered_list = [pair for pair in points if
                     tuple(pair) not in elements_set]
    return np.asarray(filtered_list)


def _get_nearby_hor_points(current_points, points, residual, order=2):
    """
    Find nearby horizontal points based on polynomial fit.
    This function gets points that are within a specified residual
    distance from a fitted curve of current points.

    Parameters
    ----------
    current_points : array_like
        Array of shape (N, 2) of (y, x)-coordinates of current points.
    points : array_like
        Array of shape (M, 2) of (y, x)-coordinates of points to evaluate.
    residual : float
        Maximum allowable distance between a point and the fitted curve.
    order : int, optional
        Polynomial order

    Returns
    -------
    nearby_points : array_like
        Selected points.
    """
    cx, cy = current_points[:, 1], current_points[:, 0]
    poly_fit = np.poly1d(np.polyfit(cx, cy, int(order)))
    x, y = points[:, 1], points[:, 0]
    distances = np.abs(y - poly_fit(x))
    nearby_points = points[distances <= residual]
    return nearby_points


def _get_nearby_hor_points_iter(initial_points, points, x_left, x_right,
                                search_dist, residual, overlap_ratio=0.5,
                                order=2):
    """
    Iteratively find nearby horizontal points based on polynomial fit.
    This function iteratively gets points within a horizontal search range
    and appends those to the selected points.

    Parameters
    ----------
    initial_points : array_like
        Array of shape (N, 2) of (y, x)-coordinates of initial points.
    points : array_like
        Array of shape (M, 2) of (y, x)-coordinates of points to evaluate.
    x_left : float
        Initial x-coordinate to start the leftward search.
    x_right : float
        Initial x-coordinate to start the rightward search.
    search_dist : float
        Distance to expand the horizontal search range in each iteration.
    residual : float
        Maximum allowable distance between a point and the fitted curve.
    overlap_ratio : float, optional
        To specify the overlap of the searching range between searching steps.
    order : int, optional
        Polynomial order

    Returns
    -------
    selected_points : array_like
        (y, x)-coordinates of selected points.
    """
    overlap_ratio = np.clip(overlap_ratio, 0.0, 1.0)
    overlap = search_dist * overlap_ratio
    xr_curr = x_right
    xr_next = xr_curr + search_dist
    xl_curr = x_left
    xl_next = xl_curr - search_dist
    x_list = points[:, 1]
    selected_points = initial_points
    check = True
    while check:
        xr_next1, xr_curr1 = xr_next + overlap, xr_curr - overlap
        xl_next1, xl_curr1 = xl_next - overlap, xl_curr + overlap
        indices = np.where(((xr_next1 >= x_list) & (x_list > xr_curr1)) |
                           ((xl_next1 <= x_list) & (x_list < xl_curr1)))[0]
        if len(indices) > 0:
            nearby_points = _get_nearby_hor_points(selected_points,
                                                   points[indices], residual,
                                                   order=order)
            if len(nearby_points) > 0:
                selected_points = np.vstack([selected_points, nearby_points])
                selected_points = np.unique(selected_points, axis=0)
            else:
                check = False
        else:
            check = False
        xr_curr = xr_next
        xr_next = xr_curr + search_dist
        xl_curr = xl_next
        xl_next = xl_curr - search_dist
    return selected_points


def group_dots_hor_lines_based_polyfit(points, slope, line_dist, ratio=0.1,
                                       num_dot_miss=3, accepted_ratio=0.65,
                                       overlap_ratio=0.5, order=2):
    """
    Group dots into horizontal lines.

    Parameters
    ----------
    points : array_like
        Array of shape (N, 2) of (y, x)-coordinates of points.
    slope : float
        Horizontal slope of the grid.
    line_dist : float
        Nominal distance between two lines.
    ratio : float
        Acceptable variation.
    num_dot_miss : int
        Acceptable missing dots between dot1 and dot2.
    accepted_ratio : float
        Use to select lines having the number of dots equal to or larger than
        the multiplication of the `accepted_ratio` and the maximum number of
        dots per line.
    overlap_ratio : float, optional
        To specify the overlap of the searching range between searching steps.
    order : int, optional
        Polynomial order.

    Returns
    -------
    list of array_like
        List of 2D arrays. Each array is (y, x)-coordinates of points belong
        to the same group. Length of each list may be different.
    """
    angle = -np.arctan(slope)
    num_points = len(points)
    if num_points == 0:
        raise ValueError("Input is empty!!!")
    points = rotate_points(points, angle, degree_unit=False)
    points = points[points[:, 1].argsort()]
    x_list = points[:, 1]
    x_min, x_max = x_list[0], x_list[-1]
    x_mid = (x_min + x_max) * 0.5
    num_dot_miss = np.clip(num_dot_miss, 1, num_points)
    search_dist = num_dot_miss * line_dist + 0.5 * line_dist
    x_start = np.clip(x_mid - search_dist, x_min, x_max)
    x_stop = np.clip(x_mid + search_dist, x_min, x_max)
    idx_list = np.where((x_list >= x_start) & (x_list <= x_stop))[0]
    list_lines = []
    if len(idx_list) > 0:
        selected_points = points[idx_list]
        grouped_points = group_dots_hor_lines(
            selected_points, 0.0, line_dist, ratio=ratio,
            num_dot_miss=num_dot_miss, accepted_ratio=accepted_ratio)
        if len(grouped_points) > 0:
            residual = ratio * line_dist
            for current_points in grouped_points:
                if len(current_points) > 2:
                    selected_points = current_points
                    x_left = selected_points[0, 1]
                    x_right = selected_points[-1, 1]
                    selected_points = _get_nearby_hor_points_iter(
                        selected_points, points, x_left, x_right, search_dist,
                        residual=residual, overlap_ratio=overlap_ratio,
                        order=order)
                if len(selected_points) > 2:
                    selected_points = rotate_points(selected_points, -angle,
                                                    degree_unit=False)
                    selected_points = selected_points[
                        selected_points[:, 1].argsort()]
                    list_lines.append(selected_points)
    list_len = [len(i) for i in list_lines]
    len_accepted = np.int16(accepted_ratio * np.max(list_len))
    lines_selected = [line for line in list_lines if len(line) > len_accepted]
    y_vals = [np.median(line[:, 0]) for line in lines_selected]
    ids = np.where(np.abs(np.diff(y_vals)) > 0.1 * line_dist)[0]
    if len(ids) > 0:
        ids = np.insert(ids + 1, 0, 0)
        lines_tmp = [line for idx, line in enumerate(lines_selected) if
                     idx in ids]
        lines_selected = lines_tmp
    lines_selected = sorted(
        lines_selected, key=lambda list_: np.mean(list_[:, 0]))
    return lines_selected


def _get_nearby_ver_points(current_points, points, residual, order=2):
    """
    Find nearby vertical points based on polynomial fit.
    This function gets points that are within a specified residual distance
    from a fitted curve of current points.

    Parameters
    ----------
    current_points : array_like
        Array of shape (N, 2) of (y, x)-coordinates of current points.
    points : array_like
        Array of shape (M, 2) of (y, x)-coordinates of points to evaluate.
    residual : float
        Maximum allowable distance between a point and the fitted curve.
    order : int, optional
        polynomial order

    Returns
    -------
    nearby_points : array_like
        Selected points.
    """
    cx, cy = current_points[:, 1], current_points[:, 0]
    poly_fit = np.poly1d(np.polyfit(cy, cx, int(order)))
    x, y = points[:, 1], points[:, 0]
    distances = np.abs(x - poly_fit(y))
    nearby_points = points[distances <= residual]
    return nearby_points


def _get_nearby_ver_points_iter(initial_points, points, y_left, y_right,
                                search_dist, residual, overlap_ratio=0.5,
                                order=2):
    """
    Iteratively find nearby vertical points based on polynomial fit.
    This function iteratively gets points within a horizontal search range
    and appends those to the selected points.

    Parameters
    ----------
    initial_points : array_like
        Array of shape (N, 2) of (y, x)-coordinates of initial points.
    points : array_like
        Array of shape (M, 2) of (y, x)-coordinates of points to evaluate.
    x_left : float
        Initial x-coordinate to start the leftward search.
    x_right : float
        Initial x-coordinate to start the rightward search.
    search_dist : float
        Distance to expand the horizontal search range in each iteration.
    residual : float
        Maximum allowable distance between a point and the fitted curve.
    overlap_ratio : float, optional
        To specify the overlap of the searching range between searching steps.
    order : int, optional
        polynomial order

    Returns
    -------
    selected_points : array_like
        (y, x)-coordinates of selected points.
    """
    overlap_ratio = np.clip(overlap_ratio, 0.0, 1.0)
    overlap = search_dist * overlap_ratio
    yr_curr = y_right
    yr_next = yr_curr + search_dist
    yl_curr = y_left
    yl_next = yl_curr - search_dist
    y_list = points[:, 0]
    selected_points = initial_points
    check = True
    while check:
        yr_next1, yr_curr1 = yr_next + overlap, yr_curr - overlap
        yl_next1, yl_curr1 = yl_next - overlap, yl_curr + overlap
        indices = np.where(((yr_next1 >= y_list) & (y_list > yr_curr1)) |
                           ((yl_next1 <= y_list) & (y_list < yl_curr1)))[0]
        if len(indices) > 0:
            nearby_points = _get_nearby_ver_points(selected_points,
                                                   points[indices], residual,
                                                   order=order)
            if len(nearby_points) > 0:
                selected_points = np.vstack([selected_points, nearby_points])
                selected_points = np.unique(selected_points, axis=0)
            else:
                check = False
        else:
            check = False
        yr_curr = yr_next
        yr_next = yr_curr + search_dist
        yl_curr = yl_next
        yl_next = yl_curr - search_dist
    return selected_points


def group_dots_ver_lines_based_polyfit(points, slope, line_dist, ratio=0.1,
                                       num_dot_miss=3, accepted_ratio=0.65,
                                       overlap_ratio=0.5, order=2):
    """
    Group dots into vertical lines.

    Parameters
    ----------
    points : array_like
        Array of shape (N, 2) of (y, x)-coordinates of points.
    slope : float
        Vertical slope of the grid.
    line_dist : float
        Nominal distance between two lines.
    ratio : float
        Acceptable variation.
    num_dot_miss : int
        Acceptable missing dots between dot1 and dot2.
    accepted_ratio : float
        Use to select lines having the number of dots equal to or larger than
        the multiplication of the `accepted_ratio` and the maximum number of
        dots per line.
    overlap_ratio : float, optional
        To specify the overlap of the searching range between searching steps.
    order : int, optional
        Polynomial order.

    Returns
    -------
    list of array_like
        List of 2D arrays. Each array is (y, x)-coordinates of points belong
        to the same group. Length of each list may be different.
    """
    angle = np.arctan(slope)
    num_points = len(points)
    if num_points == 0:
        raise ValueError("Input is empty!!!")
    points = rotate_points(points, angle, degree_unit=False)
    points = points[points[:, 0].argsort()]
    y_list = points[:, 0]
    y_min, y_max = y_list[0], y_list[-1]
    y_mid = (y_min + y_max) * 0.5
    num_dot_miss = np.clip(num_dot_miss, 1, num_points)
    search_dist = num_dot_miss * line_dist + 0.5 * line_dist
    y_start = np.clip(y_mid - search_dist, y_min, y_max)
    y_stop = np.clip(y_mid + search_dist, y_min, y_max)
    idx_list = np.where((y_list >= y_start) & (y_list <= y_stop))[0]
    list_lines = []
    if len(idx_list) > 0:
        selected_points = points[idx_list]
        grouped_points = group_dots_ver_lines(
            selected_points, 0.0, line_dist, ratio=ratio,
            num_dot_miss=num_dot_miss, accepted_ratio=accepted_ratio)
        if len(grouped_points) > 0:
            residual = ratio * line_dist
            for current_points in grouped_points:
                if len(current_points) > 2:
                    selected_points = current_points
                    y_left = selected_points[0, 0]
                    y_right = selected_points[-1, 0]
                    selected_points = _get_nearby_ver_points_iter(
                        selected_points, points, y_left, y_right, search_dist,
                        residual, overlap_ratio=overlap_ratio, order=order)
                if len(selected_points) > 2:
                    selected_points = rotate_points(selected_points, -angle,
                                                    degree_unit=False)
                    selected_points = selected_points[
                        selected_points[:, 0].argsort()]
                    list_lines.append(selected_points)
    list_len = [len(i) for i in list_lines]
    len_accepted = np.int16(accepted_ratio * np.max(list_len))
    lines_selected = [line for line in list_lines if len(line) > len_accepted]
    x_vals = [np.median(line[:, 1]) for line in lines_selected]
    ids = np.where( np.abs(np.diff(x_vals)) > 0.1 * line_dist)[0]
    if len(ids) > 0:
        ids = np.insert(ids + 1, 0, 0)
        lines_tmp = [line for idx, line in enumerate(lines_selected) if
                     idx in ids]
        lines_selected = lines_tmp
    lines_selected = sorted(
        lines_selected, key=lambda list_: np.mean(list_[:, 1]))
    return lines_selected
