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

    def test_get_local_extrema_points(self):
        f_alias = lipa.get_local_extrema_points
        size = 800
        np.random.seed(1)
        list_data = np.ones(size)
        num_point = 0
        for i in range(10, size - 10, 50):
            list_data[i - 4: i + 4] = 0.0
            num_point += 1
        list_data = list_data + 0.2 * np.random.rand(size)
        list_data = ndi.gaussian_filter1d(list_data, 2)
        points = f_alias(list_data, option="min", radius=7, sensitive=0.2,
                         denoise=False, norm=True)
        list_mp = list_data[np.int16(points)]
        thres = 0.2
        self.assertTrue(num_point == len(points) and np.min(list_mp) < thres
                        and np.max(list_mp) < thres)

        list_data = np.float32(1.0 - list_data)
        thres = 1.0 - thres
        points = f_alias(list_data, option="max", radius=7, sensitive=0.2,
                         denoise=False, norm=True)
        list_mp = list_data[np.int16(points)]
        self.assertTrue(num_point == len(points) and np.min(list_mp) > thres
                        and np.max(list_mp) > thres)

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

    @staticmethod
    def __make_chessboard(hei, wid, step):
        mat = np.ones((hei, wid), dtype=np.float32)
        num = 0
        for i in range(0, hei, step):
            if num % 2 == 0:
                num1 = 0
                for j in range(0, wid, step):
                    if num1 % 2 == 0:
                        mat[i: i + step, j: j + step] = 1.0
                    else:
                        mat[i: i + step, j: j + step] = 0.0
                    num1 = num1 + 1
            else:
                num1 = 0
                for j in range(0, wid, step):
                    if num1 % 2 == 0:
                        mat[i: i + step, j: j + step] = 0.0
                    else:
                        mat[i: i + step, j: j + step] = 1.0
                    num1 = num1 + 1
            num = num + 1
        return mat

    def test_convert_chessboard_to_linepattern(self):
        f_alias = lipa.convert_chessboard_to_linepattern
        chessboard = self.__make_chessboard(90, 120, 30)
        np.random.seed(1)
        chessboard = chessboard + 0.4 * np.random.rand(90, 120)
        line_pattern = f_alias(chessboard, smooth=True, bgr="bright")

        line1 = ndi.gaussian_filter1d(line_pattern[10], 3)
        line2 = ndi.gaussian_filter1d(line_pattern[:, 10], 3)
        points1 = lipa.get_local_extrema_points(line1, radius=7, sensitive=0.2,
                                                denoise=False, norm=False)
        points2 = lipa.get_local_extrema_points(line2, radius=7, sensitive=0.2,
                                                denoise=False, norm=False)
        num1, num2 = len(points1), len(points2)
        self.assertTrue(num1 == 3 and num2 == 2)

        line_pattern = f_alias(chessboard, smooth=False, bgr="dark")
        line1 = ndi.gaussian_filter1d(line_pattern[10], 3)
        line2 = ndi.gaussian_filter1d(line_pattern[:, 10], 3)
        points1 = lipa.get_local_extrema_points(line1, option="max", radius=7,
                                                sensitive=0.2, denoise=False)
        points2 = lipa.get_local_extrema_points(line2, option="max", radius=7,
                                                sensitive=0.2, denoise=False)
        num1, num2 = len(points1), len(points2)
        self.assertTrue(num1 == 3 and num2 == 2)

    def test_get_tilted_profile(self):
        chessboard = self.__make_chessboard(90, 120, 30)
        np.random.seed(1)
        chessboard = chessboard + 0.2 * np.random.rand(90, 120)
        line_pattern = lipa.convert_chessboard_to_linepattern(chessboard,
                                                              smooth=True,
                                                              bgr="bright")
        line1 = lipa.get_tilted_profile(line_pattern, 22, 10, "horizontal")[-1]
        line2 = lipa.get_tilted_profile(line_pattern, 22, -10, "vertical")[-1]
        points1 = lipa.get_local_extrema_points(line1, option="min", radius=7,
                                                sensitive=0.2, denoise=True)
        points2 = lipa.get_local_extrema_points(line2, option="min", radius=7,
                                                sensitive=0.2, denoise=True)
        num1, num2 = len(points1), len(points2)
        self.assertTrue(num1 == 3 and num2 == 2)

