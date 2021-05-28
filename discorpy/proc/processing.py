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
Module of processing methods:
- Fit lines of dots to parabolas, find the center of distortion.
- Calculate undistorted intercepts of gridlines.
- Calculate distortion coefficients of the backward model, the forward model,
and the backward-from-forward model.
"""
import numpy as np
from scipy import optimize


def _para_fit_hor(list_lines, xcenter, ycenter):
    """
    Fit horizontal lines of dots to parabolas.

    Parameters
    ----------
    list_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.

    Returns
    -------
    list_coef : list of 1D arrays
        List of the coefficients of each parabola.
    list_slines : list of 2D arrays
        List of the shifted coordinates of dot-centroids on each line.
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
        List of the coordinates of dot-centroids on each line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.

    Returns
    -------
    list_coef : list of 1D arrays
        List of the coefficients of each parabola.
    list_slines : list of 2D arrays
        List of the shifted coordinates of dot-centroids on each line.
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
        List of the coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each vertical line.

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
    ycenter0 = (list_coef_hor[pos_hor - 1, 2] + list_coef_hor[pos_hor, 2]) * 0.5
    xcenter0 = (list_coef_ver[pos_ver - 1, 2] + list_coef_ver[pos_ver, 2]) * 0.5
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
        List of the coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each vertical line.
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
        List of the coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each vertical line.
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


def _check_missing_lines(list_coef_hor, list_coef_ver):
    """
    Check if there are missing lines

    Parameters
    ----------
    list_coef_hor : list of 1D arrays
        Coefficients of parabolic fits of horizontal lines.
    list_coef_ver : list of 1D arrays
        Coefficients of parabolic fits of vertical lines.

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
    if (herror > 0.3) or (verror > 0.3):
        check = True
    return check


def _func_opt(d0, c0, indexc0, *list_inter):
    """
    Function for finding the optimum undistorted distance.
    """
    return np.sum(
        np.asarray([(np.sign(c) * np.abs(i - indexc0) * d0 + c0 - c) ** 2
                    for i, c in enumerate(list_inter)]))


def _optimize_intercept(dist_hv, pos_hv, list_inter):
    """
    Find the optimum undistorted distance.
    """
    list_arg = [list_inter[pos_hv], pos_hv]
    list_arg.extend(list_inter)
    minimum = optimize.minimize(_func_opt, dist_hv, args=tuple(list_arg))
    return minimum.x[0]


def _calc_undistor_intercept(list_hor_lines, list_ver_lines, xcenter, ycenter,
                             optimizing=False):
    """
    Calculate the intercepts of undistorted lines.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each vertical line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    optimizing : bool, optional
        Apply optimization if True.

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
    check = _check_missing_lines(list_coef_hor, list_coef_ver)
    if check:
        print("!!! ERROR !!!")
        print("Parameters of the methods of grouping dots need to be adjusted")
        raise ValueError("There're missing lines, algorithm will not work!!!")
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
                       num_fact):
    """
    Calculate the distortion coefficients of a backward mode.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each vertical line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    num_fact : int
        Number of the factors of polynomial.

    Returns
    -------
    list_fact : list of float
        Coefficients of the polynomial.
    """
    num_fact = np.int16(np.clip(num_fact, 1, None))
    (list_hor_uc, list_ver_uc) = _calc_undistor_intercept(
        list_hor_lines, list_ver_lines, xcenter, ycenter)
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
                      num_fact):
    """
    Calculate the distortion coefficients of a forward mode.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each vertical line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    num_fact : int
        Number of the factors of polynomial.

    Returns
    -------
    list_fact : list of float
        Coefficients of the polynomial.
    """
    num_fact = np.int16(np.clip(num_fact, 1, None))
    (list_hor_uc, list_ver_uc) = _calc_undistor_intercept(
        list_hor_lines, list_ver_lines, xcenter, ycenter)
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
                                    ycenter, num_fact):
    """
    Calculate the distortion coefficients of a backward mode from a forward
    model.

    Parameters
    ----------
    list_hor_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each horizontal line.
    list_ver_lines : list of 2D arrays
        List of the coordinates of dot-centroids on each vertical line.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    num_fact : int
        Number of the factors of polynomial.

    Returns
    -------
    list_ffact : list of floats
        Polynomial coefficients of the forward model.
    list_bfact : list of floats
        Polynomial coefficients of the backward model.
    """
    num_fact = np.int16(np.clip(num_fact, 1, 5))
    list_ffact = np.float64(
        calc_coef_forward(list_hor_lines, list_ver_lines, xcenter, ycenter,
                          num_fact))
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
