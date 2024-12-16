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
# Publication date: 21 November 2021
# ============================================================================
# Contributors:
# ============================================================================

"""
Module of pre-processing methods for handling line-pattern images and
chessboard images:

-   Determine the slopes and distances between lines in horizontal and
    vertical directions.
-   Convert a chessboard image to a line-pattern image.
-   Extract line profiles across an image.
-   Get points belonging to a line by locating local extrema.
-   Get points belonging to the edges of a chessboard image by converting
    a line profile into a list of slopes, calculated using a linear fit of
    neighboring points.
-   Select good points using Gaussian peak fitting.
"""

import warnings
import numpy as np
import scipy.ndimage as ndi
from skimage.transform import radon
from scipy.optimize import curve_fit, OptimizeWarning
import discorpy.prep.preprocessing as prep


def locate_subpixel_point(list_point, option="min"):
    """
    Locate the extremum point of a 1D array with subpixel accuracy.

    Parameters
    ----------
    list_point : array_like
        1D array.
    option : {"min", "max"}
        To locate the minimum point or the maximum point.

    Returns
    -------
    float
        Subpixel position of the extremum point.
    """
    num_point = len(list_point)
    a, b, c = np.polyfit(np.arange(num_point), list_point, 2)
    if option == "min":
        pos = np.argmin(list_point)
    else:
        pos = np.argmax(list_point)
    if a != 0.0:
        num = - b / (2 * a)
        if (num >= 0) and (num < num_point):
            pos = num
    return pos


def __gauss_function(x, a, b, c, d):
    y = a * np.exp(-np.power((x - c) / (2 * b ** 2), 2)) + d
    return y


def _get_gauss_peak_fit(list_data):
    """
    Fit a Gaussian to 1d-array and return fitted result, peak position, and a
    successful flag.
    """
    npoint = len(list_data)
    list_x = np.arange(npoint) - npoint // 2
    check = False
    try:
        init_guess = [1.0, 1.0, 0.0, 0.0]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=OptimizeWarning)
            parameters = curve_fit(__gauss_function, list_x, list_data,
                                   p0=init_guess)[0]
            a = parameters[0]
            b, c = parameters[1], parameters[2]
            d = parameters[3]
            fit_data = __gauss_function(list_x, a, b, c, d)
            check = True
    except:
        fit_data = list_data
        c, d = list_x[0], 1.0
    return fit_data, c, d, check


def select_good_peaks(list_data, peaks, tol=0.2, radius=11, sigma=0,
                      use_offset=True):
    """
    Select good peaks from the input data based on Gaussian fitting and
    tolerance criteria.

    Parameters
    ----------
    list_data : array_like
        1d-array.
    peaks : list of int
        Indices of candidate peaks in the data.
    tol : float, optional
        Tolerance for peak fitting accuracy.
    radius : int, optional
        Radius around each peak to consider for fitting.
    sigma : float, optional
        Standard deviation for Gaussian smoothing.
    use_offset : bool, optional
        Use fitted offset value for judging if True.

    Returns
    -------
    list of int
        Indices of the selected good peaks.
    """
    if sigma > 0:
        list_data = ndi.gaussian_filter1d(list_data, sigma)
    good_peaks = []
    npoint = len(list_data)
    half_radius = radius // 2
    for p in peaks:
        start = max(0, p - radius)
        stop = min(npoint, p + radius + 1)
        if (stop - start) > 3:
            list_sub = list_data[start:stop]
            std = np.std(list_sub)
            if std != 0.0:
                list_norm = (list_sub - np.min(list_sub)) / std
                fit_data, del_x, offset, check = _get_gauss_peak_fit(list_norm)
                num = np.percentile(np.abs(fit_data - list_norm), 80)
                if check and abs(del_x) < half_radius and num < tol:
                    if use_offset:
                        if abs(offset) < tol:
                            good_peaks.append(p)
                    else:
                        good_peaks.append(p)
    return np.asarray(good_peaks)


def sliding_window_slope(list_data, size=3, norm=True):
    """
    Compute the absolute slopes of a linear fit within a sliding window
    across the input data, normalized by the mean slope if norm is True.

    Parameters
    ----------
    list_data : array-like
        1d-array.
    size : int, optional
        Size of the sliding window for the linear fit.
    norm : bool, optional
        Normalize the result by the mean slope if True.

    Returns
    -------
    array_like
        Normalized absolute slopes for each data point.
    """
    npoint = len(list_data)
    if npoint < 3:
        raise ValueError("Data size must be larger than 2")
    size = int(np.clip(size, 3, npoint))
    if size % 2 == 0:
        size += 1
    radius = size // 2
    padded_data = np.pad(list_data, (radius, radius), 'edge')
    x_list = np.arange(size)
    slopes = np.zeros_like(list_data)
    for i in range(len(list_data)):
        local_data = padded_data[i:i + size]
        slopes[i] = np.polyfit(x_list, local_data, 1)[0]
    slopes = np.abs(slopes)
    if norm is True:
        nmean = np.mean(slopes)
        if nmean != 0.0:
            slopes = slopes / nmean
    return slopes


def get_local_extrema_points(list_data, option="min", radius=7, sensitive=0.1,
                             denoise=True, norm=True, subpixel=True,
                             select_peaks=False, **kwargs):
    """
    Get a list of local extremum points from a 1D array.

    Parameters
    ----------
    list_data : array_like
        1D array.
    option : {"min", "max"}
        To get minimum points or maximum points
    radius : int
        Search radius. Used to locate extremum points.
    sensitive : float
        To detect extremum points against random noise. Smaller is more
        sensitive.
    denoise : bool, optional
        Applying a smoothing filter if True.
    norm : bool, optional
        Apply background normalization to the array.
    subpixel : bool, optional
        Locate points with subpixel accuracy.
    select_peaks : bool, optional
        To select good points based on Gaussian fitting.
    **kwargs : optional
        Additional parameters for the method 'select_good_peaks'

    Returns
    -------
    array_like
        1D array. Positions of local extremum points.

    """
    list_data = np.copy(list_data)
    if denoise is True:
        list_data = ndi.gaussian_filter(list_data, 3)
    if option == "max":
        list_data = np.max(list_data) - list_data
    num_point = len(list_data)
    radius = int(np.clip(radius, 1, num_point // 4))
    if norm is True:
        xlist = np.arange(num_point)
        mat_comb = np.asarray(np.vstack((xlist, list_data)))
        mat_sort = mat_comb[:, mat_comb[1, :].argsort()]
        list_sort = mat_sort[1]
        ndrop = int(0.25 * num_point)
        (a1, a0) = np.polyfit(xlist[ndrop:-ndrop - 1],
                              list_sort[ndrop:-ndrop - 1], 1)[:2]
        list_fit = a1 * xlist + a0
        l_thres, u_thres = a0, a1 * xlist[-1] + a0
        list_sort[(list_fit >= l_thres) & (list_fit <= u_thres)] = list_fit[
            (list_fit >= l_thres) & (list_fit <= u_thres)]
        mat_sort[1] = list_sort
        nmean = np.mean(np.abs(list_fit))
        backgr = mat_sort[:, mat_sort[0, :].argsort()][1]
        list_data = np.divide(list_data, backgr,
                              out=nmean * np.ones_like(list_data),
                              where=backgr != 0)
    points = []
    for i in range(radius, num_point - radius - 1, 1):
        val, pos = list_data[i], i
        list_sort = np.sort(list_data[i - radius:i + radius + 1])
        num1 = list_sort[0] - val
        nmean = np.mean(list_sort[-radius:])
        num2 = np.abs((val - nmean) / nmean) if nmean != 0 else 0.0
        if num1 == 0.0 and num2 > sensitive:
            points.append(pos)
    if select_peaks is True:
        points = select_good_peaks(np.max(list_data) - list_data, points,
                                   radius=radius, **kwargs)
    if subpixel is True:
        points_sub = []
        if len(points) > 0:
            for i in points:
                pos_sub = i - 1 + locate_subpixel_point(list_data[i - 1:i + 2],
                                                        option="min")
                points_sub.append(pos_sub)
            return np.asarray(points_sub)
    return np.asarray(points)


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


def calc_slope_distance_hor_lines(mat, ratio=0.3, search_range=30.0, radius=9,
                                  sensitive=0.1, bgr="bright", denoise=True,
                                  norm=True, subpixel=True, chessboard=False,
                                  select_peaks=False, **kwargs):
    """
    Calculate the representative distance between horizontal lines and the
    representative slope of these lines using the ROI around the middle of a
    line-pattern image or a chessboard image.

    Parameters
    ----------
    mat : array_like
        2D array.
    ratio : float
        Used to select the ROI around the middle of an image.
    search_range : float
        Search range in Degree to determine the slope of lines.
    radius : int
        Search radius. Used to locate lines.
    sensitive : float
        To detect lines against random noise. Smaller is more sensitive.
    bgr : {"bright", "dark"}
        Specify the brightness of the background against the lines.
    denoise : bool, optional
        Apply a smoothing filter if True.
    norm : bool, optional
        Apply background normalization to the array.
    subpixel : bool, optional
        Locate points with subpixel accuracy.
    chessboard : bool, optional
        If True, converts the input chessboard image to a line-pattern image.
    select_peaks : bool, optional
        To select good points based on Gaussian fitting.
    **kwargs : optional
        Additional parameters for the method 'select_good_peaks'

    Returns
    -------
    slope : float
        Slope of horizontal lines in Radian.
    distance : float
        Distance between horizontal lines.
    """
    if chessboard is True:
        mat = convert_chessboard_to_linepattern(mat)
    if denoise is True:
        mat = ndi.gaussian_filter(mat, 3)
    mat_roi = prep._select_roi(mat, ratio, square=True)
    if bgr == "bright":
        mat_roi = np.max(mat_roi) - mat_roi
    angle_coarse = 90.0 + np.arange(-search_range, search_range + 1.0)
    mask = _make_circle_mask(mat_roi.shape[0], 0.92)
    sinogram1 = radon(mat_roi * mask, theta=angle_coarse, circle=True)
    list_max1 = np.amax(sinogram1, axis=0)
    pos_max1 = np.argmax(list_max1)
    best_angle1 = angle_coarse[pos_max1]
    angle_fine = np.arange(best_angle1 - 1.0, best_angle1 + 1.05, 0.05)
    sinogram2 = radon(mat_roi * mask, theta=angle_fine, circle=True)
    list_max2 = np.amax(sinogram2, axis=0)
    pos_max2 = np.argmax(list_max2)
    best_angle2 = -(angle_fine[pos_max2] - 90)
    slope = np.tan(best_angle2 * np.pi / 180.0)
    list_ext_point = get_local_extrema_points(sinogram2[:, pos_max2],
                                              option="max", radius=radius,
                                              denoise=denoise, norm=norm,
                                              subpixel=subpixel,
                                              sensitive=sensitive,
                                              select_peaks=select_peaks,
                                              **kwargs)
    if len(list_ext_point) > 3:
        distance = np.median(np.abs(np.diff(list_ext_point)))
    else:
        distance = np.mean(np.abs(np.diff(list_ext_point)))
    return slope, distance


def calc_slope_distance_ver_lines(mat, ratio=0.3, search_range=30.0, radius=9,
                                  sensitive=0.1, bgr="bright", denoise=True,
                                  norm=True, subpixel=True, chessboard=False,
                                  select_peaks=False, **kwargs):
    """
    Calculate the representative distance between vertical lines and the
    representative slope of these lines using the ROI around the middle of a
    line-pattern image or a chessboard image.

    Parameters
    ----------
    mat : array_like
        2D array.
    ratio : float
        Used to select the ROI around the middle of an image.
    search_range : float
        Search range in Degree to determine the slope of lines.
    radius : int
        Search radius. Used to locate lines.
    sensitive : float
        To detect lines against random noise. Smaller is more sensitive.
    bgr : {"bright", "dark"}
        Specify the brightness of the background against the lines.
    denoise : bool, optional
        Applying a smoothing filter if True.
    subpixel : bool, optional
        Locate points with subpixel accuracy.
    chessboard : bool, optional
        If True, converts the input chessboard image to a line-pattern image.
    select_peaks : bool, optional
        To select good points based on Gaussian fitting.
    **kwargs : optional
        Additional parameters for the method 'select_good_peaks'

    Returns
    -------
    slope : float
        Slope of vertical lines in Radian.
    distance : float
        Distance between vertical lines.
    """
    if chessboard is True:
        mat = convert_chessboard_to_linepattern(mat)
    if denoise is True:
        mat = ndi.gaussian_filter(mat, 3)
    mat_roi = prep._select_roi(mat, ratio, square=True)
    if bgr == "bright":
        mat_roi = np.max(mat_roi) - mat_roi
    angle_coarse = np.arange(-search_range, search_range + 1.0)
    mask = _make_circle_mask(mat_roi.shape[0], 0.92)
    sinogram1 = radon(mat_roi * mask, theta=angle_coarse, circle=True)
    list_max1 = np.amax(sinogram1, axis=0)
    pos_max1 = np.argmax(list_max1)
    best_angle1 = angle_coarse[pos_max1]
    angle_fine = np.arange(best_angle1 - 1.0, best_angle1 + 1.05, 0.05)
    sinogram2 = radon(mat_roi * mask, theta=angle_fine, circle=True)
    list_max2 = np.amax(sinogram2, axis=0)
    pos_max2 = np.argmax(list_max2)
    best_angle2 = angle_fine[pos_max2]
    slope = np.tan(best_angle2 * np.pi / 180.0)
    list_ext_point = get_local_extrema_points(sinogram2[:, pos_max2],
                                              option="max", radius=radius,
                                              denoise=denoise, norm=norm,
                                              subpixel=subpixel,
                                              sensitive=sensitive,
                                              select_peaks=select_peaks,
                                              **kwargs)
    if len(list_ext_point) > 3:
        distance = np.median(np.abs(np.diff(list_ext_point)))
    else:
        distance = np.mean(np.abs(np.diff(list_ext_point)))
    return slope, distance


def _calc_index_range(height, width, angle_deg, direction):
    """
    Calculate extractable range of tilted line-profile. Positive angle is
    counterclockwise.

    Parameters
    ----------
    height : int
        Height of the image.
    width : int
        Width of the image.
    angle_deg : float
        Tilted angle in Degree.
    direction : {"horizontal", "vertical"}
        Direction of line-profile.

    Returns
    -------
    min_idx : int
        Minimum index of lines.
    max_idx : int
        Maximum index of lines.
    """
    angle = angle_deg * np.pi / 180.0
    if direction == "horizontal":
        if np.abs(angle_deg) == 90.0:
            raise ValueError("If the input angle is around 90-degree, use "
                             "the 'vertical' option and update the angle to "
                             "around 0-degree instead!!!")
        else:
            if angle_deg > 0:
                min_idx = int(np.ceil(width * np.tan(angle)))
                max_idx = height - 1
            else:
                min_idx = 0
                max_idx = height - 1 - int(
                    np.floor(width * np.tan(np.abs(angle))))
            if (min_idx < 0) or (min_idx >= height) or (max_idx < 0) or (
                    max_idx >= height):
                raise ValueError("Row index is out of range, please select "
                                 "the direction correctly !!!")
    else:
        if np.abs(angle_deg) == 90.0:
            raise ValueError("If the input angle is around 90-degree, use "
                             "the 'horizontal' option and update the angle to "
                             "around 0-degree instead!!!")
        else:
            if angle_deg > 0:
                min_idx = 0
                max_idx = width - 1 - int(np.ceil(height * np.tan(angle)))
            else:
                min_idx = int(np.floor(height * np.tan(np.abs(angle))))
                max_idx = width - 1
            if (min_idx < 0) or (min_idx >= width) or (max_idx < 0) or (
                    max_idx >= width):
                raise ValueError("Column index is out of range, please select "
                                 "the direction correctly !!!")
    return min_idx, max_idx


def get_tilted_profile(mat, index, angle_deg, direction):
    """
    Get the intensity-profile along a tilted line across an image. Positive
    angle is counterclockwise.

    Parameters
    ----------
    mat : array_like
        2D array.
    index : int
        Index of the line.
    angle_deg : float
        Tilted angle in Degree.
    direction : {"horizontal", "vertical"}
        Direction of line-profile.

    Returns
    -------
    xlist : array_like
        1D array. x-positions of points on the line.
    ylist : array_like
        1D array. y-positions of points on the line.
    profile : array_like
        1D array. Intensities of points on the line.
    """
    if mat.ndim != 2:
        raise ValueError("Input must be a 2D array !!!")
    (height, width) = mat.shape
    (min_idx, max_idx) = _calc_index_range(height, width, angle_deg, direction)
    angle = angle_deg * np.pi / 180.0
    if (index < min_idx) or (index > max_idx):
        raise ValueError("Input index is out of possible range: "
                         "[{0}, {1}]".format(min_idx, max_idx))
    if direction == "horizontal":
        rlist = np.linspace(0, np.floor(width / np.cos(angle)), width)
        xlist = rlist * np.cos(angle)
        ylist = rlist * np.sin(-angle)
        xlist = np.clip(xlist, 0, width - 1)
        ylist = np.clip(index + ylist, 0, height - 1)
        ymin = np.int16(np.floor(np.amin(ylist)))
        ymax = np.int16(np.ceil(np.amax(ylist))) + 1
        indices = ylist - ymin, xlist
        profile = ndi.map_coordinates(mat[ymin:ymax, :], indices, order=3,
                                      mode='nearest')
    else:
        rlist = np.linspace(0, np.floor(height / np.cos(angle)), height)
        ylist = rlist * np.cos(angle)
        xlist = rlist * np.sin(angle)
        xlist = np.clip(index + xlist, 0, width - 1)
        ylist = np.clip(ylist, 0, height - 1)
        xmin = np.int16(np.floor(np.amin(xlist)))
        xmax = np.int16(np.ceil(np.amax(xlist))) + 1
        indices = ylist, xlist - xmin
        profile = ndi.map_coordinates(mat[:, xmin:xmax], indices, order=3,
                                      mode='nearest')
    return xlist, ylist, profile


def convert_chessboard_to_linepattern(mat, smooth=True, bgr="bright",
                                      sigma=3):
    """
    Convert a chessboard image to a line-pattern image.

    Parameters
    ----------
    mat : array_like
        2D array.
    smooth : bool, optional
        Apply a gaussian smoothing filter if True.
    bgr : {'bright', 'dark'}
        Select the background of the output image.
    sigma : int
        Sigma of the Gaussian window, if smooth is True.

    Returns
    -------
    array_like
        Line-pattern image.
    """
    if smooth is True:
        mat = ndi.gaussian_filter(mat, sigma, mode="nearest")
    mat_line = np.mean(np.abs(np.gradient(mat)), axis=0)
    if smooth is True:
        mat_line = np.pad(mat_line[4:-4, 4:-4], 4, mode="edge")
    else:
        mat_line = np.pad(mat_line[2:-2, 2:-2], 2, mode="edge")
    if bgr == "bright":
        mat_line = np.max(mat_line) - mat_line
    mat_line = mat_line / np.mean(np.abs(mat_line))
    return mat_line


def get_cross_points_hor_lines(mat, slope_ver, dist_ver, ratio=0.3, norm=True,
                               offset=0, bgr="bright", radius=11,
                               sensitive=0.1, denoise=True, subpixel=True,
                               chessboard=False, select_peaks=False, **kwargs):
    """
    Get points on horizontal lines of a line-pattern image, or a chessboard
    image, by intersecting with a list of generated vertical-lines.

    Parameters
    ----------
    mat : array_like
        2D array.
    slope_ver : float
        Representative slope of vertical lines, can be calculated by using the
        'calc_slope_distance_ver_lines' method.
    dist_ver : float
        Representative distance between vertical lines, can be calculated by
        using the 'calc_slope_distance_ver_lines' method.
    ratio : float
        To adjust the distance (=ratio * dist_ver) between generated lines
        to create more or fewer lines.
    norm : bool, optional
        Apply background normalization to the array.
    offset : int
        Starting index of generated lines.
    bgr : {"bright", "dark"}
        Specify the brightness of the background relative to the lines.
    radius : int
        Search radius. Used to locate extremum points.
    sensitive : float
        To detect extremum points against random noise. Smaller is more
        sensitive.
    denoise : bool, optional
        Applying a smoothing filter if True.
    subpixel : bool, optional
        Locate points with subpixel accuracy.
    chessboard : bool, optional
        If True, cross points are located by finding local maxima of slopes
        of local linear fits of the generated lines.
    select_peaks : bool, optional
        To select good points based on Gaussian fitting.
    **kwargs : optional
        Additional parameters for the method 'select_good_peaks'

    Returns
    -------
    array_like
        List of (y,x)-coordinates of points.
    """
    (height, width) = mat.shape
    if bgr == "bright":
        mat = np.max(mat) - mat
    if norm is True:
        mat = prep.normalization_fft(mat, 5)
    if denoise is True:
        mat = ndi.gaussian_filter(mat, 3)
    angle = np.arctan(slope_ver)
    min_row, max_row = _calc_index_range(height, width, np.rad2deg(angle),
                                         direction="vertical")
    offset = int(np.clip(offset, 0, min(height, width) // 3))
    list_points = []
    for i in np.arange(min_row + offset, max_row - offset, ratio * dist_ver):
        xlist, ylist, profile = get_tilted_profile(mat, i, np.rad2deg(angle),
                                                   direction="vertical")
        scale = np.sqrt((xlist[-1] - xlist[0]) ** 2
                        + (ylist[-1] - ylist[0]) ** 2) / (height - 1)
        if chessboard is True:
            profile = sliding_window_slope(profile, size=3)
        rlist = get_local_extrema_points(profile, option="max", radius=radius,
                                         sensitive=sensitive,
                                         denoise=not denoise, norm=not norm,
                                         subpixel=subpixel,
                                         select_peaks=select_peaks,
                                         **kwargs) * scale
        xlist1 = rlist * np.sin(angle) + xlist[0]
        ylist1 = rlist * np.cos(angle) + ylist[0]
        list_points.extend(np.asarray(list(zip(ylist1, xlist1))))
    return np.asarray(list_points)


def get_cross_points_ver_lines(mat, slope_hor, dist_hor, ratio=0.3, norm=True,
                               offset=0, bgr="bright", radius=11,
                               sensitive=0.1, denoise=True, subpixel=True,
                               chessboard=False, select_peaks=False, **kwargs):
    """
    Get points on vertical lines of a line-pattern image, or a chessboard
    image, by intersecting with a list of generated horizontal-lines.

    Parameters
    ----------
    mat : array_like
        2D array.
    slope_hor : float
        Representative slope of horizontal lines, can be calculated by using
        the 'calc_slope_distance_hor_lines' method.
    dist_hor : float
        Representative distance between horizontal lines, can be calculated by
        using the 'calc_slope_distance_hor_lines' method.
    ratio : float
        To adjust the distance (=ratio * dist_hor) between generated lines
        to create more or fewer lines.
    norm : bool, optional
        Apply background normalization to the array.
    offset : int
        Starting index of generated lines.
    bgr : {"bright", "dark"}
        Specify the brightness of the background against the lines.
    radius : int
        Search radius. Used to locate extremum points.
    sensitive : float
        To detect extremum points against random noise. Smaller is more
        sensitive.
    denoise : bool, optional
        Applying a smoothing filter if True.
    subpixel : bool, optional
        Locate points with subpixel accuracy.
    chessboard : bool, optional
        If True, cross points are located by finding local maxima of slopes
        of local linear fits of the generated lines.
    select_peaks : bool, optional
        To select good points based on Gaussian fitting.
    **kwargs : optional
        Additional parameters for the method 'select_good_peaks'

    Returns
    -------
    array_like
        List of (y,x)-coordinates of points.
    """
    (height, width) = mat.shape
    if bgr == "bright":
        mat = np.max(mat) - mat
    if norm is True:
        mat = prep.normalization_fft(mat, 5)
    if denoise is True:
        mat = ndi.gaussian_filter(mat, 3)
    angle = np.arctan(slope_hor)
    min_col, max_col = _calc_index_range(height, width, -np.rad2deg(angle),
                                         direction="horizontal")
    offset = int(np.clip(offset, 0, min(height, width) // 8))
    list_points = []
    for i in np.arange(min_col + offset, max_col - offset, ratio * dist_hor):
        xlist, ylist, profile = get_tilted_profile(mat, i, -np.rad2deg(angle),
                                                   direction="horizontal")
        scale = np.sqrt((xlist[-1] - xlist[0]) ** 2
                        + (ylist[-1] - ylist[0]) ** 2) / (width - 1)
        if chessboard is True:
            profile = sliding_window_slope(profile, size=3)
        rlist = get_local_extrema_points(profile, option="max", radius=radius,
                                         sensitive=sensitive,
                                         denoise=not denoise, norm=not norm,
                                         subpixel=subpixel,
                                         select_peaks=select_peaks,
                                         **kwargs) * scale
        xlist1 = rlist * np.cos(angle) + xlist[0]
        ylist1 = rlist * np.sin(angle) + ylist[0]
        list_points.extend(np.asarray(list(zip(ylist1, xlist1))))
    return np.asarray(list_points)
