# ============================================================================
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
# ============================================================================
# Contributors:
# ============================================================================

"""
Tests for methods in linepattern.py
"""

import unittest
import numpy as np
import scipy.ndimage as ndi
import discorpy.prep.linepattern as lipa


class LinepatternMethods(unittest.TestCase):

    def setUp(self):
        self.eps = 10 ** (-6)
        hei, wid, pad, step = 128, 128, 4, 20
        mat = np.zeros((hei, wid), dtype=np.float32)
        num_hline = 0
        for i in range(step + pad, hei - pad, step):
            mat[i - 2:i + 3, step + pad - 2:wid - pad - step + 3] = 1.0
            num_hline = num_hline + 1
        num_vline = 0
        for i in range(step + pad, wid - pad, step):
            mat[step + pad - 2:hei - step - pad + 3, i - 2:i + 3] = 1.0
            num_vline = num_vline + 1
        mat_lines = ndi.gaussian_filter(1.0 - 0.2 * mat, 1)
        np.random.seed(1)
        self.mat = mat_lines + 0.05 * np.random.rand(hei, wid)
        self.dist = step
        self.num_hline, self.num_vline = num_hline, num_vline

    def test_calc_slope_distance_hor_lines(self):
        slope, distance = lipa.calc_slope_distance_hor_lines(self.mat,
                                                             ratio=0.8,
                                                             radius=4,
                                                             denoise=False,
                                                             norm=False,
                                                             subpixel=False)
        self.assertTrue(np.abs(slope) < self.eps and
                        np.abs(distance - self.dist) <= 1.0)

    def test_calc_slope_distance_ver_lines(self):
        slope, distance = lipa.calc_slope_distance_ver_lines(self.mat,
                                                             ratio=0.8,
                                                             radius=4,
                                                             denoise=False,
                                                             norm=False,
                                                             subpixel=False)
        self.assertTrue(np.abs(slope) < self.eps and
                        np.abs(distance - self.dist) <= 1.0)

    def test_get_cross_points_hor_lines(self):
        slope_ver, dist_ver = lipa.calc_slope_distance_ver_lines(self.mat,
                                                                 ratio=0.5,
                                                                 radius=4,
                                                                 denoise=False,
                                                                 norm=False)
        list_points = lipa.get_cross_points_hor_lines(self.mat, slope_ver,
                                                      dist_ver, bgr="bright",
                                                      radius=4, ratio=0.5,
                                                      denoise=True, norm=True,
                                                      offset=0)
        list_data = np.abs(np.diff(np.sort(list_points[:, 0])))
        num_line = len(lipa.get_local_extrema_points(list_data, option="max",
                                                     radius=4, denoise=False,
                                                     norm=False,
                                                     subpixel=False))
        self.assertTrue(num_line == self.num_hline - 1)

    def test_get_cross_points_ver_lines(self):
        slope_hor, dist_hor = lipa.calc_slope_distance_hor_lines(self.mat,
                                                                 ratio=0.5,
                                                                 radius=4,
                                                                 denoise=False,
                                                                 norm=False)
        list_points = lipa.get_cross_points_ver_lines(self.mat, slope_hor,
                                                      dist_hor, bgr="bright",
                                                      radius=4, ratio=0.5,
                                                      denoise=True, norm=True,
                                                      offset=0)
        list_data = np.abs(np.diff(np.sort(list_points[:, 1])))
        num_line = len(lipa.get_local_extrema_points(list_data, option="max",
                                                     radius=4, denoise=False,
                                                     norm=False,
                                                     subpixel=False))
        self.assertTrue(num_line == self.num_vline - 1)
