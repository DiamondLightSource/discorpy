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
# Description: Python implementation of the author's methods of
# distortion correction, Nghia T. Vo et al "Radial lens distortion
# correction with sub-pixel accuracy for X-ray micro-tomography"
# Optics Express 23, 32859-32868 (2015), https://doi.org/10.1364/OE.23.032859
# Publication date: 10th July 2018
# ============================================================================
# Contributors:
# ============================================================================

"""
Module of processing methods:

- Fit lines of dots to parabolas, find the center of distortion.
- Calculate undistorted intercepts of gridlines.
- Calculate distortion coefficients of the backward model, the forward model,
  and the backward-from-forward model.
- Correct perspective distortion affecting curve lines.
- Generate non-perspective points or lines from perspective points or lines.
- Calculate perspective coefficients.

"""
import numpy as np
from scipy import optimize


def _para_fit_hor(list_lines, xcenter, ycenter):
    """
    Fit horizontal lines of dots to parabolas.

    Parameters
    ----------
    list_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.

    Returns
    -------
    list_coef : list of 1D arrays
        List of the coefficients of each parabola (y=ax**2+bx+c).
    list_slines : list of 2D arrays
        List of the shifted (y,x)-coordinates of dot-centroids on each line.
    """
    num_line = len(list_lines)
    list_coef = np.zeros((num_line, 3), dtype=np.float32)
    list_slines = []
    for i, iline in enumerate(list_lines):
        line = np.asarray(iline)
        list_coef[i] = np.asarray(np.polyfit(line[:, 1] - xcenter,
                                             line[:, 0] - ycenter, 2))
        list_temp = np.asarray(
            [(dot[0] - ycenter, dot[1] - xcenter) for dot in line])
        list_slines.append(list_temp)
    return list_coef, list_slines


def _para_fit_ver(list_lines, xcenter, ycenter):
    """
    Fit vertical lines of dots to parabolas.

    Parameters
    ----------
    list_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.

    Returns
    -------
    list_coef : list of 1D arrays
        List of the coefficients of each parabola (x=ay**2+by+c).
    list_slines : list of 2D arrays
        List of the shifted (y,x)-coordinates of dot-centroids on each line.
    """
    num_line = len(list_lines)
    list_coef = np.zeros((num_line, 3), dtype=np.float32)
    list_slines = []
    for i, iline in enumerate(list_lines):
        line = np.asarray(iline)
        list_coef[i] = np.asarray(
            np.polyfit(line[:, 0] - ycenter, line[:, 1] - xcenter, 2))
        list_temp = np.asarray(
            [(dot[0] - ycenter, dot[1] - xcenter) for dot in line])
        list_slines.append(list_temp)
    return list_coef, list_slines


def find_cod_coarse(list_hor_lines, list_ver_lines):
    """
    Coarse estimation of the center of distortion.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each vertical line.

    Returns
    -------
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    """
    (list_coef_hor, list_hor_lines) = _para_fit_hor(list_hor_lines, 0.0, 0.0)
    (list_coef_ver, list_ver_lines) = _para_fit_ver(list_ver_lines, 0.0, 0.0)
    pos_hor = np.argmax(np.abs(np.diff(np.sign(list_coef_hor[:, 0])))) + 1
    pos_ver = np.argmax(np.abs(np.diff(np.sign(list_coef_ver[:, 0])))) + 1
    ycenter0 = (list_coef_hor[pos_hor - 1, 2] + list_coef_hor[
        pos_hor, 2]) * 0.5
    xcenter0 = (list_coef_ver[pos_ver - 1, 2] + list_coef_ver[
        pos_ver, 2]) * 0.5
    slope_hor = (list_coef_hor[pos_hor - 1, 1] + list_coef_hor[
        pos_hor, 1]) * 0.5
    slope_ver = (list_coef_ver[pos_ver - 1, 1] + list_coef_ver[
        pos_ver, 1]) * 0.5
    ycenter = (ycenter0 + xcenter0 * slope_hor) / (1.0 - slope_hor * slope_ver)
    xcenter = (xcenter0 + ycenter0 * slope_ver) / (1.0 - slope_hor * slope_ver)
    return xcenter, ycenter


def _func_dist(x, a, b, c):
    """
    Function for finding the minimum distance.
    """
    return x ** 2 + (a * x ** 2 + b * x + c) ** 2


def _calc_error(list_coef_hor, list_coef_ver):
    """
    Calculate a metric of measuring how close fitted lines to the coordinate
    origin by: locating points on each parabola having the minimum distance
    to the origin, applying linear fits to these points, adding intercepts of
    the fits.

    Parameters
    ----------
    list_coef_hor : list of 1D arrays
        Coefficients of parabolic fits of horizontal lines.
    list_coef_ver : list of 1D arrays
        Coefficients of parabolic fits of vertical lines.

    Returns
    -------
    float
    """
    num_hline = len(list_coef_hor)
    num_vline = len(list_coef_ver)
    list_hpoint = np.zeros((num_hline, 2), dtype=np.float32)
    for i, coefs in enumerate(list_coef_hor):
        minimum = optimize.minimize(_func_dist, 0.0, args=tuple(coefs))
        xm = minimum.x[0]
        ym = coefs[0] * xm ** 2 + coefs[1] * xm + coefs[2]
        list_hpoint[i, 0] = xm
        list_hpoint[i, 1] = ym
    list_vpoint = np.zeros((num_vline, 2), dtype=np.float32)
    for i, coefs in enumerate(list_coef_ver):
        minimum = optimize.minimize(_func_dist, 0.0, args=tuple(coefs))
        ym = minimum.x[0]
        xm = coefs[0] * ym ** 2 + coefs[1] * ym + coefs[2]
        list_vpoint[i, 0] = ym
        list_vpoint[i, 1] = xm
    error_h = np.polyfit(list_hpoint[:, 0], list_hpoint[:, 1], 1)[-1]
    error_v = np.polyfit(list_vpoint[:, 0], list_vpoint[:, 1], 1)[-1]
    return np.abs(error_h) + np.abs(error_v)


def _calc_metric(list_hor_lines, list_ver_lines, xcenter, ycenter,
                 list_xshift, list_yshift):
    """
    Calculate a metric for determining the best center of distortion.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each vertical line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_xshift : list of float
        List of x-offsets from the x-center.
    list_yshift : list of float
        List of y-offsets from the y-center.

    Returns
    -------
    xshift : float
        Shift in x-direction from the x-center.
    yshift : float
        Shift in y-direction from the y-center.
    """
    (list_coef_hor, list_hor_lines) = _para_fit_hor(
        list_hor_lines, xcenter, ycenter)
    (list_coef_ver, list_ver_lines) = _para_fit_ver(
        list_ver_lines, xcenter, ycenter)
    pos_hor = np.argmin(np.abs(list_coef_hor[:, 2]))
    pos_ver = np.argmin(np.abs(list_coef_ver[:, 2]))
    mat_metric = np.zeros(
        (len(list_xshift), len(list_yshift)), dtype=np.float32)
    num_hline = len(list_hor_lines)
    num_vline = len(list_ver_lines)
    numuse = min(5, num_hline // 2 - 1, num_vline // 2 - 1)
    (posh1, posh2) = (
        max(0, pos_hor - numuse), min(num_hline, pos_hor + numuse + 1))
    (posv1, posv2) = (
        max(0, pos_ver - numuse), min(num_vline, pos_ver + numuse + 1))
    for j, pos_x in enumerate(list_xshift):
        for i, pos_y in enumerate(list_yshift):
            (list_coef_hor, _) = _para_fit_hor(
                list_hor_lines[posh1:posh2], pos_x, pos_y)
            (list_coef_ver, _) = _para_fit_ver(
                list_ver_lines[posv1:posv2], pos_x, pos_y)
            mat_metric[i, j] = _calc_error(list_coef_hor, list_coef_ver)
    min_pos = (np.unravel_index(mat_metric.argmin(), mat_metric.shape))
    xshift = list_xshift[min_pos[1]]
    yshift = list_yshift[min_pos[0]]
    return xshift, yshift


def find_cod_fine(list_hor_lines, list_ver_lines, xcenter, ycenter, dot_dist):
    """
    Find the best center of distortion (CoD) by searching around the coarse
    estimation of the CoD.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each vertical line.
    xcenter : float
        Coarse estimation of the CoD in x-direction.
    ycenter : float
        Coarse estimation of the CoD in y-direction.
    dot_dist : float
        Median distance of two nearest dots.

    Returns
    -------
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    """
    step0 = 2.0
    list_xshift = np.arange(-dot_dist, dot_dist + step0, step0)
    list_yshift = list_xshift
    (xshift, yshift) = _calc_metric(
        list_hor_lines, list_ver_lines, xcenter, ycenter, list_xshift,
        list_yshift)
    xcenter1 = xcenter + xshift
    ycenter1 = ycenter + yshift
    step = 0.5
    list_xshift = np.arange(-step0, step0 + step, step)
    list_yshift = list_xshift
    (xshift, yshift) = _calc_metric(
        list_hor_lines, list_ver_lines, xcenter1, ycenter1, list_xshift,
        list_yshift)
    xcenter2 = xcenter1 + xshift
    ycenter2 = ycenter1 + yshift
    return xcenter2, ycenter2


def _check_missing_lines(list_coef_hor, list_coef_ver, threshold=0.3):
    """
    Check if there are missing lines

    Parameters
    ----------
    list_coef_hor : list of 1D arrays
        Coefficients of parabolic fits of horizontal lines.
    list_coef_ver : list of 1D arrays
        Coefficients of parabolic fits of vertical lines.
    threshold : float
        To determine if there are missing lines. Larger is less sensitive.

    Returns
    -------
    bool
    """
    check = False
    list_dist_hor = np.abs(np.diff(list_coef_hor[:, 2]))
    list_dist_ver = np.abs(np.diff(list_coef_ver[:, 2]))
    list_hindex = np.arange(len(list_dist_hor))
    list_vindex = np.arange(len(list_dist_ver))
    hfact = np.polyfit(list_hindex, list_dist_hor, 2)
    vfact = np.polyfit(list_vindex, list_dist_ver, 2)
    list_fit_hor = hfact[0] * list_hindex ** 2 + \
                   hfact[1] * list_hindex + hfact[2]
    list_fit_ver = vfact[0] * list_vindex ** 2 + \
                   vfact[1] * list_vindex + vfact[2]
    herror = np.max(np.abs((list_dist_hor - list_fit_hor) / list_fit_hor))
    verror = np.max(np.abs((list_dist_ver - list_fit_ver) / list_fit_ver))
    if (herror > threshold) or (verror > threshold):
        check = True
    return check


def _func_opt(d0, c0, indexc0, *list_inter):
    """
    Function for finding the optimum undistorted distance for radial
    distortion correction.
    """
    return np.sum(
        np.asarray([(np.sign(c) * np.abs(i - indexc0) * d0 + c0 - c) ** 2
                    for i, c in enumerate(list_inter)]))


def _optimize_intercept(dist_hv, pos_hv, list_inter):
    """
    Find the optimum undistorted distance for radial-distortion correction.
    """
    list_arg = [list_inter[pos_hv], pos_hv]
    list_arg.extend(list_inter)
    minimum = optimize.minimize(_func_opt, dist_hv, args=tuple(list_arg))
    return minimum.x[0]


def _calc_undistor_intercept(list_hor_lines, list_ver_lines, xcenter, ycenter,
                             optimizing=False, threshold=0.3):
    """
    Calculate the intercepts of undistorted lines.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each vertical line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    optimizing : bool, optional
        Apply optimization if True.
    threshold : float
        To determine if there are missing lines. Larger is less sensitive.

    Returns
    -------
    list_hor_uc : list of floats
        Intercepts of undistorted horizontal lines.
    list_ver_uc : list of floats
        Intercepts of undistorted vertical lines.
    """
    (list_coef_hor, list_hor_lines) = _para_fit_hor(
        list_hor_lines, xcenter, ycenter)
    (list_coef_ver, list_ver_lines) = _para_fit_ver(
        list_ver_lines, xcenter, ycenter)
    check = _check_missing_lines(list_coef_hor, list_coef_ver,
                                 threshold=threshold)
    if check:
        msg = "!!! Parameters of the methods of grouping dots need to be " \
              "adjusted !!!\n!!! Check if there are missing lines or adjust " \
              "the threshold value !!!"
        raise ValueError(msg)
    pos_hor = np.argmin(np.abs(list_coef_hor[:, 2]))
    pos_ver = np.argmin(np.abs(list_coef_ver[:, 2]))
    num_hline = len(list_hor_lines)
    num_vline = len(list_ver_lines)
    num_use = min(3, num_hline // 2 - 1, num_vline // 2 - 1)
    (posh1, posh2) = (
        max(0, pos_hor - num_use), min(num_hline, pos_hor + num_use + 1))
    (posv1, posv2) = (
        max(0, pos_ver - num_use), min(num_vline, pos_ver + num_use + 1))
    dist_hor = np.mean(np.abs(np.diff(list_coef_hor[posh1: posh2, 2])))
    dist_ver = np.mean(np.abs(np.diff(list_coef_ver[posv1: posv2, 2])))
    if optimizing is True:
        dist_hor = _optimize_intercept(dist_hor, pos_hor, list_coef_hor[:, 2])
        dist_ver = _optimize_intercept(dist_ver, pos_ver, list_coef_ver[:, 2])
    list_hor_uc = np.zeros(num_hline, dtype=np.float32)
    list_ver_uc = np.zeros(num_vline, dtype=np.float32)
    for i in range(num_hline):
        dist = np.abs(i - pos_hor) * dist_hor
        list_hor_uc[i] = np.sign(list_coef_hor[i, 2]) * dist + list_coef_hor[
            pos_hor, 2]
    for i in range(num_vline):
        dist = np.abs(i - pos_ver) * dist_ver
        list_ver_uc[i] = np.sign(list_coef_ver[i, 2]) * dist + list_coef_ver[
            pos_ver, 2]
    return list_hor_uc, list_ver_uc


def calc_coef_backward(list_hor_lines, list_ver_lines, xcenter, ycenter,
                       num_fact, optimizing=False, threshold=0.3):
    """
    Calculate the distortion coefficients of a backward mode.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each vertical line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    num_fact : int
        Number of the factors of polynomial.
    optimizing : bool, optional
        Apply optimization if True.
    threshold : float
        To determine if there are missing lines. Larger is less sensitive.

    Returns
    -------
    list_fact : list of float
        Coefficients of the polynomial.
    """
    num_fact = np.int16(np.clip(num_fact, 1, None))
    (list_hor_uc, list_ver_uc) = _calc_undistor_intercept(
        list_hor_lines, list_ver_lines, xcenter, ycenter,
        optimizing=optimizing, threshold=threshold)
    (list_coef_hor, list_hor_lines) = _para_fit_hor(
        list_hor_lines, xcenter, ycenter)
    (list_coef_ver, list_ver_lines) = _para_fit_ver(
        list_ver_lines, xcenter, ycenter)
    Amatrix = []
    Bmatrix = []
    list_expo = np.arange(num_fact, dtype=np.int16)
    for i, line in enumerate(list_hor_lines):
        (a_coef, _, c_coef) = np.float64(list_coef_hor[i])
        uc_coef = np.float64(list_hor_uc[i])
        for _, point in enumerate(line):
            xd = np.float64(point[1])
            yd = np.float64(point[0])
            rd = np.sqrt(xd * xd + yd * yd)
            Fb = (a_coef * xd * xd + c_coef) / uc_coef
            Amatrix.append(np.power(rd / Fb, list_expo))
            Bmatrix.append(Fb)
    for i, line in enumerate(list_ver_lines):
        (a_coef, _, c_coef) = np.float64(list_coef_ver[i])
        uc_coef = np.float64(list_ver_uc[i])
        for _, point in enumerate(line):
            xd = np.float64(point[1])
            yd = np.float64(point[0])
            rd = np.sqrt(xd * xd + yd * yd)
            Fb = (a_coef * yd * yd + c_coef) / uc_coef
            Amatrix.append(np.power(rd / Fb, list_expo))
            Bmatrix.append(Fb)
    Amatrix = np.asarray(Amatrix, dtype=np.float64)
    Bmatrix = np.asarray(Bmatrix, dtype=np.float64)
    list_fact = np.linalg.lstsq(Amatrix, Bmatrix, rcond=1e-64)[0]
    return list_fact


def calc_coef_forward(list_hor_lines, list_ver_lines, xcenter, ycenter,
                      num_fact, optimizing=False, threshold=0.3):
    """
    Calculate the distortion coefficients of a forward mode.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each vertical line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    num_fact : int
        Number of the factors of polynomial.
    optimizing : bool, optional
        Apply optimization if True.
    threshold : float
        To determine if there are missing lines. Larger is less sensitive.

    Returns
    -------
    list_fact : list of float
        Coefficients of the polynomial.
    """
    num_fact = np.int16(np.clip(num_fact, 1, None))
    (list_hor_uc, list_ver_uc) = _calc_undistor_intercept(
        list_hor_lines, list_ver_lines, xcenter, ycenter,
        optimizing=optimizing, threshold=threshold)
    (list_coef_hor, list_hor_lines) = _para_fit_hor(
        list_hor_lines, xcenter, ycenter)
    (list_coef_ver, list_ver_lines) = _para_fit_ver(
        list_ver_lines, xcenter, ycenter)
    list_expo = np.arange(num_fact, dtype=np.int16)
    Amatrix = []
    Bmatrix = []
    for i, line in enumerate(list_hor_lines):
        (a_coef, _, c_coef) = np.float64(list_coef_hor[i])
        uc_coef = np.float64(list_hor_uc[i])
        if uc_coef != 0.0:
            for _, point in enumerate(line):
                xd = np.float64(point[1])
                yd = np.float64(point[0])
                rd = np.sqrt(xd * xd + yd * yd)
                Fb = uc_coef / (a_coef * xd * xd + c_coef)
                if Fb != 0.0:
                    Amatrix.append(np.power(rd, list_expo))
                    Bmatrix.append(Fb)
    for i, line in enumerate(list_ver_lines):
        (a_coef, _, c_coef) = np.float64(list_coef_ver[i])
        uc_coef = np.float64(list_ver_uc[i])
        if uc_coef != 0.0:
            for _, point in enumerate(line):
                xd = np.float64(point[1])
                yd = np.float64(point[0])
                rd = np.sqrt(xd * xd + yd * yd)
                Fb = uc_coef / (a_coef * yd * yd + c_coef)
                if Fb != 0.0:
                    Amatrix.append(np.power(rd, list_expo))
                    Bmatrix.append(Fb)
    Amatrix = np.asarray(Amatrix, dtype=np.float64)
    Bmatrix = np.asarray(Bmatrix, dtype=np.float64)
    list_fact = np.linalg.lstsq(Amatrix, Bmatrix, rcond=1e-64)[0]
    return list_fact


def calc_coef_backward_from_forward(list_hor_lines, list_ver_lines, xcenter,
                                    ycenter, num_fact, optimizing=False,
                                    threshold=0.3):
    """
    Calculate the distortion coefficients of a backward mode from a forward
    model.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the (y,x)-coordinates of dot-centroids on each vertical line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    num_fact : int
        Number of the factors of polynomial.
    optimizing : bool, optional
        Apply optimization if True.
    threshold : float
        To determine if there are missing lines. Larger is less sensitive.

    Returns
    -------
    list_ffact : list of floats
        Polynomial coefficients of the forward model.
    list_bfact : list of floats
        Polynomial coefficients of the backward model.
    """
    num_fact = np.int16(np.clip(num_fact, 1, None))
    list_ffact = np.float64(
        calc_coef_forward(list_hor_lines, list_ver_lines, xcenter, ycenter,
                          num_fact, optimizing=optimizing,
                          threshold=threshold))
    (_, list_hor_lines) = _para_fit_hor(list_hor_lines, xcenter, ycenter)
    (_, list_ver_lines) = _para_fit_ver(list_ver_lines, xcenter, ycenter)
    list_expo = np.arange(num_fact, dtype=np.int16)
    Amatrix = []
    Bmatrix = []
    for _, line in enumerate(list_hor_lines):
        for _, point in enumerate(line):
            xd = np.float64(point[1])
            yd = np.float64(point[0])
            rd = np.sqrt(xd * xd + yd * yd)
            ffactor = np.float64(np.sum(list_ffact * np.power(rd, list_expo)))
            if ffactor != 0.0:
                Fb = 1 / ffactor
                ru = ffactor * rd
                Amatrix.append(np.power(ru, list_expo))
                Bmatrix.append(Fb)
    for _, line in enumerate(list_ver_lines):
        for _, point in enumerate(line):
            xd = np.float64(point[1])
            yd = np.float64(point[0])
            rd = np.sqrt(xd * xd + yd * yd)
            ffactor = np.float64(np.sum(list_ffact * np.power(rd, list_expo)))
            if ffactor != 0.0:
                Fb = 1 / ffactor
                ru = ffactor * rd
                Amatrix.append(np.power(ru, list_expo))
                Bmatrix.append(Fb)
    Amatrix = np.asarray(Amatrix, dtype=np.float64)
    Bmatrix = np.asarray(Bmatrix, dtype=np.float64)
    list_bfact = np.linalg.lstsq(Amatrix, Bmatrix, rcond=1e-64)[0]
    return list_ffact, list_bfact


def transform_coef_backward_and_forward(list_fact, mapping="backward",
                                        ref_points=None):
    """
    Transform polynomial coefficients of a radial distortion model between
    forward mapping and backward mapping.

    Parameters
    ----------
    list_fact : list of floats
        Polynomial coefficients of the radial distortion model.
    mapping : {'backward', 'forward'}
        Transformation direction.
    ref_points : list of 1D-arrays, optional
        List of the (y,x)-coordinates of points used for the transformation.
        Generated if None given.

    Returns
    -------
    list of floats
        Polynomial coefficients of the reversed model.
    """
    if ref_points is None:
        ref_points = [[i, j] for i in np.arange(-1000, 1000, 50) for j in
                      np.arange(-1000, 1000, 50)]
    else:
        num_points = len(ref_points)
        if num_points < len(list_fact):
            raise ValueError("Number of reference-points must be equal or "
                             "larger than the number of coefficients!!!")
    Amatrix = []
    Bmatrix = []
    list_expo = np.arange(len(list_fact), dtype=np.int16)
    if mapping == "forward":
        for point in ref_points:
            xu = np.float64(point[1])
            yu = np.float64(point[0])
            ru = np.sqrt(xu * xu + yu * yu)
            factor = np.float64(
                np.sum(list_fact * np.power(ru, list_expo)))
            if factor != 0.0:
                Fb = 1 / factor
                rd = factor * ru
                Amatrix.append(np.power(rd, list_expo))
                Bmatrix.append(Fb)
    else:
        for point in ref_points:
            xd = np.float64(point[1])
            yd = np.float64(point[0])
            rd = np.sqrt(xd * xd + yd * yd)
            factor = np.float64(
                np.sum(list_fact * np.power(rd, list_expo)))
            if factor != 0.0:
                Fb = 1 / factor
                ru = factor * rd
                Amatrix.append(np.power(ru, list_expo))
                Bmatrix.append(Fb)
    Amatrix = np.asarray(Amatrix, dtype=np.float64)
    Bmatrix = np.asarray(Bmatrix, dtype=np.float64)
    trans_fact = np.linalg.lstsq(Amatrix, Bmatrix, rcond=1e-64)[0]
    return trans_fact


def find_cod_bailey(list_hor_lines, list_ver_lines, iteration=2):
    """
    Find the center of distortion (COD) using the Bailey's approach (Ref. [1]).

    Parameters
    ----------
    list_hor_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each horizontal line.
    list_ver_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each vertical line.

    Returns
    -------
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.

    References
    ----------
    [1].. https://www-ist.massey.ac.nz/dbailey/sprg/pdfs/2002_IVCNZ_59.pdf
    """
    (xcenter, ycenter) = find_cod_coarse(list_hor_lines, list_ver_lines)
    list_coef_hor = _para_fit_hor(list_hor_lines, xcenter, ycenter)[0]
    list_coef_ver = _para_fit_ver(list_ver_lines, xcenter, ycenter)[0]
    a1, b1 = np.polyfit(list_coef_hor[:, 2], list_coef_hor[:, 0], 1)[0:2]
    a2, b2 = np.polyfit(list_coef_ver[:, 2], list_coef_ver[:, 0], 1)[0:2]
    xcenter = xcenter - b2 / a2
    ycenter = ycenter - b1 / a1
    for i in range(iteration):
        list_coef_hor = _para_fit_hor(list_hor_lines, xcenter, ycenter)[0]
        list_coef_ver = _para_fit_ver(list_ver_lines, xcenter, ycenter)[0]
        a1, b1 = np.polyfit(list_coef_hor[:, 2], list_coef_hor[:, 0], 1)[0:2]
        a2, b2 = np.polyfit(list_coef_ver[:, 2], list_coef_ver[:, 0], 1)[0:2]
        xcenter = xcenter - b2 / a2
        ycenter = ycenter - b1 / a1
    return xcenter, ycenter


def _generate_non_perspective_parabola_coef(list_hor_lines, list_ver_lines):
    """
    Correct the deviation of fitted parabola coefficients of each line caused
    by perspective distortion. Note that the resulting coefficients are
    referred to a different origin-coordinate instead of (0, 0).

    Parameters
    ----------
    list_hor_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each horizontal line.
    list_ver_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each vertical line.

    Returns
    -------
    list_coef_hor : list of 1D-arrays
        List of the corrected coefficients for horizontal lines.
    list_coef_ver : list of 1D-arrays
        List of the corrected coefficients for vertical lines.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    """
    num_hline, num_vline = len(list_hor_lines), len(list_ver_lines)
    xcenter, ycenter = find_cod_bailey(list_hor_lines, list_ver_lines)
    list_coef_hor = _para_fit_hor(list_hor_lines, xcenter, ycenter)[0]
    list_coef_ver = _para_fit_ver(list_ver_lines, xcenter, ycenter)[0]
    ah, bh = np.polyfit(list_coef_hor[:, 2], list_coef_hor[:, 1], 1)[0:2]
    av, bv = np.polyfit(list_coef_ver[:, 2], -list_coef_ver[:, 1], 1)[0:2]
    if np.abs(ah - av) >= 0.001:
        b0 = (ah * bv - av * bh) / (ah - av)
    else:
        b0 = (bh + bv) * 0.5
    list_coef_hor[:, 1] = b0 * np.ones(num_hline)
    list_coef_ver[:, 1] = -b0 * np.ones(num_vline)
    pos_hor = np.argmax(np.abs(np.diff(np.sign(list_coef_hor[:, 0])))) + 1
    pos_ver = np.argmax(np.abs(np.diff(np.sign(list_coef_ver[:, 0])))) + 1
    num_use = min(3, num_hline // 2 - 1, num_vline // 2 - 1)
    (posh1, posh2) = (
        max(0, pos_hor - num_use), min(num_hline, pos_hor + num_use + 1))
    (posv1, posv2) = (
        max(0, pos_ver - num_use), min(num_vline, pos_ver + num_use + 1))
    dist_hor = np.mean(np.abs(np.diff(list_coef_hor[posh1: posh2, 2])))
    dist_ver = np.mean(np.abs(np.diff(list_coef_ver[posv1: posv2, 2])))
    if dist_hor > dist_ver:
        list_coef_ver[:, 2] = list_coef_ver[:, 2] * dist_hor / dist_ver
        list_coef_ver[:, 0] = list_coef_ver[:, 0] * dist_hor / dist_ver
    else:
        list_coef_hor[:, 2] = list_coef_hor[:, 2] * dist_ver / dist_hor
        list_coef_hor[:, 0] = list_coef_hor[:, 0] * dist_ver / dist_hor
    return list_coef_hor, list_coef_ver, xcenter, ycenter


def _find_cross_point_between_parabolas(para_coef_hor, para_coef_ver):
    """
    Find a cross point between two parabolas.

    Parameters
    ----------
    para_coef_hor : array_like
        Coefficients of a horizontal parabola (y=ax**2+bx+c).
    para_coef_ver : array_like
        Coefficients of a vertical parabola (x=ay**2+by+c).

    Returns
    -------
    x, y : floats
        Coordinate of the cross point.
    """
    a1, b1, c1 = para_coef_hor[0:3]
    a2, b2, c2 = para_coef_ver[0:3]
    # coefs = [a1 ** 2 * a2, 2 * a1 * a2 * b1,
    #          a2 * b1 ** 2 + a1 * b2 + 2 * a1 * a2 * c1,
    #          -1 + b1 * b2 + 2 * a2 * b1 * c1,
    #          b2 * c1 + a2 * c1 ** 2 + c2]
    # xvals = np.float32(np.real(np.roots(coefs)))
    # if len(xvals) == 0:
    #     raise ValueError("Can't find a cross point between two parabolas")
    # if len(xvals) > 1:
    #     x = xvals[np.argmin(np.abs(xvals - c2))]
    # else:
    #     x = xvals[0]
    # y = a1 * x ** 2 + b1 * x + c1

    def equations(p):
        x, y = p
        return a1 * x ** 2 + b1 * x + c1 - y, a2 * y ** 2 + b2 * y + c2 - x

    initial_guess = (0, 0)
    solution = optimize.fsolve(equations, initial_guess)

    x, y = solution
    return x, y


def regenerate_grid_points_parabola(list_hor_lines, list_ver_lines,
                                    perspective=True):
    """
    Regenerating grid points by finding cross points between horizontal lines
    and vertical lines using their parabola coefficients.

    Parameters
    ----------
    list_hor_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each horizontal line.
    list_ver_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each vertical line.
    perspective : bool, optional
        Apply perspective correction if True.

    Returns
    -------
    new_hor_lines : list of 2D-arrays
        List of the updated (y,x)-coordinates of points on each horizontal
        line.
    new_ver_lines : list of 2D-arrays
        List of the updated (y,x)-coordinates of points on each vertical line.
    """
    if perspective is True:
        results = _generate_non_perspective_parabola_coef(list_hor_lines,
                                                          list_ver_lines)
        list_coef_hor, list_coef_ver, xcenter, ycenter = results
    else:
        xcenter, ycenter = find_cod_bailey(list_hor_lines, list_ver_lines)
        list_coef_hor = _para_fit_hor(list_hor_lines, xcenter, ycenter)[0]
        list_coef_ver = _para_fit_ver(list_ver_lines, xcenter, ycenter)[0]
    num_hline, num_vline = len(list_coef_hor), len(list_coef_ver)
    new_hor_lines = np.zeros((num_hline, num_vline, 2), dtype=np.float32)
    new_ver_lines = np.zeros((num_vline, num_hline, 2), dtype=np.float32)
    for i in range(num_hline):
        for j in range(num_vline):
            x, y = _find_cross_point_between_parabolas(list_coef_hor[i],
                                                       list_coef_ver[j])
            new_hor_lines[i, j] = np.asarray([y + ycenter, x + xcenter])
            new_ver_lines[j, i] = np.asarray([y + ycenter, x + xcenter])
    return new_hor_lines, new_ver_lines


def _generate_linear_coef(list_hor_lines, list_ver_lines, xcenter=0.0,
                          ycenter=0.0):
    """
    Get linear coefficients of horizontal and vertical lines from linear fit.

    Parameters
    ----------
    list_hor_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each horizontal line.
    list_ver_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each vertical line.
    xcenter : float
        X-origin of the coordinate system.
    ycenter : float
        Y-origin of the coordinate system.

    Returns
    -------
    list_coef_hor : list of 1D-arrays
        List of the linear coefficients for horizontal lines.
    list_coef_ver : list of 1D-arrays
        List of the linear coefficients for vertical lines.
    """
    num_hline, num_vline = len(list_hor_lines), len(list_ver_lines)
    list_coef_hor = np.zeros((num_hline, 2), dtype=np.float32)
    list_coef_ver = np.zeros((num_vline, 2), dtype=np.float32)
    for i in range(num_hline):
        list_coef_hor[i] = np.polyfit(list_hor_lines[i][:, 1] - xcenter,
                                      list_hor_lines[i][:, 0] - ycenter, 1)
    for i in range(num_vline):
        list_coef_ver[i] = np.polyfit(list_ver_lines[i][:, 0] - ycenter,
                                      list_ver_lines[i][:, 1] - xcenter, 1)
    return list_coef_hor, list_coef_ver


def _find_cross_point_between_lines(line_coef_hor, line_coef_ver):
    """
    Find a cross point between two lines.

    Parameters
    ----------
    line_coef_hor : array_like
        Coefficients of a horizontal line (y=ax+b).
    line_coef_ver : array_like
        Coefficients of a vertical line (x=ay+b).

    Returns
    -------
    x, y : floats
        Coordinate of the cross point.
    """
    a1, b1 = line_coef_hor
    a2, b2 = line_coef_ver
    y = (a1 * b2 + b1) / (1.0 - a1 * a2)
    x = a2 * y + b2
    return x, y


def _func_opt_pers(d0, c0, index_c0, *list_inter):
    """
    Function for finding the optimum undistorted distance for
    perspective-distortion correction.
    """
    return np.sum(
        np.asarray([((i - index_c0) * d0 + c0 - c) ** 2
                    for i, c in enumerate(list_inter)]))


def _optimize_intercept_perspective(dist_hv, pos_hv, list_inter):
    """
    Find the optimum undistorted distance for perspective-distortion
    correction.
    """
    list_arg = [list_inter[pos_hv], pos_hv]
    list_arg.extend(list_inter)
    minimum = optimize.minimize(_func_opt_pers, dist_hv, args=tuple(list_arg))
    return minimum.x[0]


def _calc_undistor_intercept_perspective(list_hor_lines, list_ver_lines,
                                         equal_dist=True, scale="mean",
                                         optimizing=True):
    """
    Calculate the intercepts of undistorted lines from perspective distortion.

    Parameters
    ----------
    list_hor_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each horizontal line.
    list_ver_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each vertical line.
    equal_dist : bool
        Use the condition that lines are equidistant if True.
    scale : {'mean', 'median', 'min', 'max', float}
        Scale option for the undistorted grid.
    optimizing : bool
        Apply optimization for finding line-distance if True.

    Returns
    -------
    u_intercept_hor : array_like
        1D array. List of undistorted intercepts of the horizontal lines.
    u_intercept_ver : array_like
        1D array. List of undistorted intercepts of the vertical lines.
    """
    list_coef_hor, list_coef_ver = _generate_linear_coef(list_hor_lines,
                                                         list_ver_lines)
    num_hline, num_vline = len(list_hor_lines), len(list_ver_lines)
    pos_hor, pos_ver = num_hline // 2, num_vline // 2
    num_use = min(np.clip(num_hline // 2 - 1, 1, None),
                  np.clip(num_vline // 2 - 1, 1, None))
    (posh1, posh2) = (max(0, pos_hor - num_use),
                      min(num_hline, pos_hor + num_use + 1))
    (posv1, posv2) = (max(0, pos_ver - num_use),
                      min(num_vline, pos_ver + num_use + 1))
    if scale == "max":
        dist_hor = np.max(np.abs(np.diff(list_coef_hor[posh1: posh2, 1])))
        dist_ver = np.max(np.abs(np.diff(list_coef_ver[posv1: posv2, 1])))
    elif scale == "min":
        dist_hor = np.min(np.abs(np.diff(list_coef_hor[posh1: posh2, 1])))
        dist_ver = np.min(np.abs(np.diff(list_coef_ver[posv1: posv2, 1])))
    elif scale == "median":
        dist_hor = np.median(np.abs(np.diff(list_coef_hor[posh1: posh2, 1])))
        dist_ver = np.median(np.abs(np.diff(list_coef_ver[posv1: posv2, 1])))
    else:
        dist_hor = np.mean(np.abs(np.diff(list_coef_hor[posh1: posh2, 1])))
        dist_ver = np.mean(np.abs(np.diff(list_coef_ver[posv1: posv2, 1])))
        if isinstance(scale, float):
            dist_hor = scale * dist_hor
            dist_ver = scale * dist_ver
    if optimizing is True:
        dist_hor = _optimize_intercept_perspective(dist_hor, pos_hor,
                                                   list_coef_hor[:, 1])
        dist_ver = _optimize_intercept_perspective(dist_ver, pos_ver,
                                                   list_coef_ver[:, 1])
    if equal_dist is True:
        if scale == "max":
            dist = max(dist_hor, dist_ver)
        elif scale == "min":
            dist = min(dist_hor, dist_ver)
        else:
            dist = (dist_hor + dist_ver) * 0.5
        dist_hor = dist_ver = dist
    u_intercept_hor = np.zeros(num_hline, dtype=np.float32)
    u_intercept_ver = np.zeros(num_vline, dtype=np.float32)
    for i in range(num_hline):
        dist = (i - pos_hor) * dist_hor
        u_intercept_hor[i] = dist + list_coef_hor[pos_hor, 1]
    for i in range(num_vline):
        dist = (i - pos_ver) * dist_ver
        u_intercept_ver[i] = dist + list_coef_ver[pos_ver, 1]
    return u_intercept_hor, u_intercept_ver


def regenerate_grid_points_linear(list_hor_lines, list_ver_lines):
    """
    Regenerating grid points by finding cross points between horizontal lines
    and vertical lines using their linear coefficients.

    Parameters
    ----------
    list_hor_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each horizontal line.
    list_ver_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each vertical line.

    Returns
    -------
    new_hor_lines : list of 2D-arrays
        List of the updated (y,x)-coordinates of points on each horizontal
        line.
    new_ver_lines : list of 2D-arrays
        List of the updated (y,x)-coordinates of points on each vertical line.
    """
    num_hline, num_vline = len(list_hor_lines), len(list_ver_lines)
    list_coef_hor, list_coef_ver = _generate_linear_coef(list_hor_lines,
                                                         list_ver_lines)
    new_hor_lines = np.zeros((num_hline, num_vline, 2), dtype=np.float32)
    new_ver_lines = np.zeros((num_vline, num_hline, 2), dtype=np.float32)
    for i in range(num_hline):
        for j in range(num_vline):
            x, y = _find_cross_point_between_lines(list_coef_hor[i],
                                                   list_coef_ver[j])
            new_hor_lines[i, j] = np.asarray([y, x])
            new_ver_lines[j, i] = np.asarray([y, x])
    return new_hor_lines, new_ver_lines


def generate_undistorted_perspective_lines(list_hor_lines, list_ver_lines,
                                           equal_dist=True, scale="mean",
                                           optimizing=True):
    """
    Generate undistorted lines from perspective lines.

    Parameters
    ----------
    list_hor_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each horizontal line.
    list_ver_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each vertical line.
    equal_dist : bool
        Use the condition that lines are equidistant if True.
    scale : {'mean', 'median', 'min', 'max', float}
        Scale option for the undistorted grid.
    optimizing : bool
        Apply optimization for finding line-distance if True.

    Returns
    -------
    list_uhor_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on undistorted horizontal
        lines.
    list_uver_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on undistorted vertical lines.
    """
    num_hline, num_vline = len(list_hor_lines), len(list_ver_lines)
    list_coef_hor, list_coef_ver = _generate_linear_coef(list_hor_lines,
                                                         list_ver_lines)
    ah, bh = np.polyfit(list_coef_hor[:, 1], list_coef_hor[:, 0], 1)[0:2]
    av, bv = np.polyfit(list_coef_ver[:, 1], -list_coef_ver[:, 0], 1)[0:2]
    if np.abs(ah - av) >= 0.0001:
        a0 = (ah * bv - av * bh) / (ah - av)
    else:
        a0 = (bh + bv) * 0.5
    list_coef_uhor = np.copy(list_coef_hor)
    list_coef_uver = np.copy(list_coef_ver)
    list_coef_uhor[:, 0] = a0 * np.ones(num_hline)
    list_coef_uver[:, 0] = -a0 * np.ones(num_vline)
    results = _calc_undistor_intercept_perspective(list_hor_lines,
                                                   list_ver_lines, equal_dist,
                                                   scale, optimizing)
    list_coef_uhor[:, 1] = results[0]
    list_coef_uver[:, 1] = results[1]
    list_uhor_lines = np.zeros((num_hline, num_vline, 2), dtype=np.float32)
    list_uver_lines = np.zeros((num_vline, num_hline, 2), dtype=np.float32)
    for i in range(num_hline):
        for j in range(num_vline):
            x, y = _find_cross_point_between_lines(list_coef_uhor[i],
                                                   list_coef_uver[j])
            list_uhor_lines[i, j] = np.asarray([y, x])
            list_uver_lines[j, i] = np.asarray([y, x])
    return list_uhor_lines, list_uver_lines


def generate_source_target_perspective_points(list_hor_lines, list_ver_lines,
                                              equal_dist=True, scale="mean",
                                              optimizing=True):
    """
    Generate source points (distorted) and target points (undistorted).

    Parameters
    ----------
    list_hor_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each horizontal line.
    list_ver_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on each vertical line.
    equal_dist : bool
        Use the condition that lines are equidistant if True.
    scale : {'mean', 'median', 'min', 'max', float}
        Scale option for the undistorted grid.
    optimizing : bool
        Apply optimization for finding line-distance if True.

    Returns
    -------
    source_points : list of 1D-arrays
        List of the (y,x)-coordinates of distorted points.
    target_points : list of 1D-arrays
        List of the (y,x)-coordinates of undistorted points.
    """
    list_hor_slines, list_ver_slines = regenerate_grid_points_linear(
        list_hor_lines, list_ver_lines)
    list_hor_tlines, _ = generate_undistorted_perspective_lines(
        list_hor_slines, list_ver_slines, equal_dist, scale, optimizing)
    source_points = []
    target_points = []
    for i in range(len(list_hor_slines)):
        for j in range(len(list_ver_slines)):
            source_points.append(list_hor_slines[i, j])
            target_points.append(list_hor_tlines[i, j])
    return np.asarray(source_points), np.asarray(target_points)


def generate_4_source_target_perspective_points(points, input_order="yx",
                                                equal_dist=False,
                                                scale="mean"):
    """
    Generate 4 rectangular points corresponding to 4 perspective-distorted
    points.

    Parameters
    ----------
    points : list of 1D-arrays
        List of the coordinates of 4 perspective-distorted points.
    input_order : {'yx', 'xy'}
        Order of the coordinates of input-points.
    equal_dist : bool
        Use the condition that the rectangular making of 4-points is square if
        True.
    scale : {'mean', 'min', 'max', float}
        Scale option for the undistorted points.

    Returns
    -------
    source_points : list of 1D-arrays
        List of the (y,x)-coordinates of distorted points.
    target_points : list of 1D-arrays
        List of the (y,x)-coordinates of undistorted points.
    """
    points = np.asarray(points, dtype=np.float32)
    if input_order == "xy":
        points = np.fliplr(points)
    if len(points) != 4:
        raise ValueError("Input must be a list of 4 points!!!")
    list_sort = points[points[:, 0].argsort()]
    p12 = list_sort[0:2]
    p12 = p12[p12[:, 1].argsort()]
    ((y1, x1), (y2, x2)) = p12
    p34 = list_sort[-2:]
    p34 = p34[p34[:, 1].argsort()]
    ((y3, x3), (y4, x4)) = p34
    source_points = np.asarray([[y1, x1], [y2, x2], [y3, x3], [y4, x4]])
    a12 = (y1 - y2) / (x1 - x2)
    b12 = y1 - a12 * x1
    a34 = (y3 - y4) / (x3 - x4)
    b34 = y3 - a34 * x3
    ah, bh = (a12 + a34) * 0.5, (b12 + b34) * 0.5
    a13 = (x1 - x3) / (y1 - y3)
    b13 = x1 - a13 * y1
    a24 = (x2 - x4) / (y2 - y4)
    b24 = x2 - a24 * y2
    av, bv = (a13 + a24) * 0.5, (b13 + b24) * 0.5
    a0 = np.sign(ah) * (np.abs(ah) + np.abs(av)) * 0.5
    dist12 = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    dist13 = np.sqrt((x1 - x3) ** 2 + (y1 - y3) ** 2)
    dist24 = np.sqrt((x2 - x4) ** 2 + (y2 - y4) ** 2)
    dist34 = np.sqrt((x3 - x4) ** 2 + (y3 - y4) ** 2)
    if scale == "max":
        dist_h = max(dist12, dist34)
        dist_v = max(dist13, dist24)
        if equal_dist is True:
            dist_h = dist_v = max(dist_v, dist_h)
    elif scale == "min":
        dist_h = min(dist12, dist34)
        dist_v = min(dist13, dist24)
        if equal_dist is True:
            dist_h = dist_v = min(dist_v, dist_h)
    else:
        dist_h = (dist12 + dist34) * 0.5
        dist_v = (dist13 + dist24) * 0.5
        if isinstance(scale, float):
            dist_h = dist_h * scale
            dist_v = dist_v * scale
        if equal_dist is True:
            dist_h = dist_v = (dist_v + dist_h) * 0.5
    dist_h, dist_v = dist_h * 0.5, dist_v * 0.5
    b1 = bh - np.abs(dist_v / np.cos(np.arctan(a0)))
    b2 = bh + np.abs(dist_v / np.cos(np.arctan(a0)))
    b3 = bv - np.abs(dist_h / np.cos(np.arctan(a0)))
    b4 = bv + np.abs(dist_h / np.cos(np.arctan(a0)))
    y1 = (a0 * b3 + b1) / (1.0 + a0 ** 2)
    x1 = -a0 * y1 + b3
    y2 = (a0 * b4 + b1) / (1.0 + a0 ** 2)
    x2 = -a0 * y2 + b4
    y3 = (a0 * b3 + b2) / (1.0 + a0 ** 2)
    x3 = -a0 * y3 + b3
    y4 = (a0 * b4 + b2) / (1.0 + a0 ** 2)
    x4 = -a0 * y4 + b4
    target_points = np.asarray([[y1, x1], [y2, x2], [y3, x3], [y4, x4]])
    return source_points, target_points


def calc_perspective_coefficients(source_points, target_points,
                                  mapping="backward"):
    """
    Calculate perspective coefficients of a matrix to map from source points
    to target points (Ref. [1]). Note that the coordinate of a point are in
    (y,x)-order. This is to be consistent with other functions in the module.

    Parameters
    ----------
    source_points : array_like
        List of the (y,x)-coordinates of distorted points.
    target_points : array_like
        List of the (y,x)-coordinates of undistorted points.
    mapping : {'backward', 'forward'}
        To select mapping direction.

    Returns
    -------
    array_like
        1D array of 8 coefficients.

    References
    ----------
    [1].. https://doi.org/10.1016/S0262-8856(98)00183-8

    """
    if mapping == "forward":
        s_points = np.fliplr(np.asarray(source_points))
        t_points = np.fliplr(np.asarray(target_points))
    else:
        s_points = np.fliplr(np.asarray(target_points))
        t_points = np.fliplr(np.asarray(source_points))
    Amatrix = []
    for p1, p2 in zip(s_points, t_points):
        Amatrix.append(
            [p1[0], p1[1], 1, 0, 0, 0, -p2[0] * p1[0], -p2[0] * p1[1]])
        Amatrix.append(
            [0, 0, 0, p1[0], p1[1], 1, -p2[1] * p1[0], -p2[1] * p1[1]])
    Amatrix = np.asarray(Amatrix, dtype=np.float64)
    Bmatrix = np.transpose(
        np.ndarray.flatten(np.asarray(t_points, dtype=np.float64)))
    list_coef = np.linalg.lstsq(Amatrix, Bmatrix, rcond=1e-64)[0]
    return list_coef


def update_center(list_lines, xcenter, ycenter):
    """
    Update the coordinate-center of points on lines.

    Parameters
    ----------
    list_lines : list of 2D-arrays
        List of the (y,x)-coordinates of points on lines.
    xcenter : float
        X-origin of the coordinate system.
    ycenter : float
        Y-origin of the coordinate system.

    Returns
    -------
        list of 2D-arrays.
    """
    updated_lines = []
    for i, iline in enumerate(list_lines):
        line = np.asarray(iline)
        list_temp = np.asarray(
            [(dot[0] + ycenter, dot[1] + xcenter) for dot in line])
        updated_lines.append(list_temp)
    return updated_lines
