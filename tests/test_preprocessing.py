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
# E-mail: 
# ============================================================================
# Contributors:
# ============================================================================

"""
Tests for methods in preprocessing.py
"""

import os
import unittest
import numpy as np
import scipy.ndimage as ndi
import discorpy.losa.loadersaver as losa
import discorpy.prep.preprocessing as prep


class PreprocessingMethods(unittest.TestCase):

    def setUp(self):
        self.var = 0.05
        sigma = 30
        (self.hei, self.wid) = (64, 64)
        (ycen, xcen) = (self.hei // 2, self.wid // 2)
        y, x = np.ogrid[-ycen:self.hei - ycen, -xcen:self.wid - xcen]
        num = 2.0 * sigma * sigma
        self.bck = np.exp(-(x * x / num + y * y / num))
        mat = np.zeros((self.hei, self.wid), dtype=np.float32)
        mat[7:self.hei:10, 7:self.wid:10] = 1
        self.num_dots = int(np.sum(mat))
        self.mat_dots = np.float32(ndi.binary_dilation(mat, iterations=2))
        self.rem_dot = 6

    def test_normalization(self):
        mat_nor = prep.normalization(self.bck, 3)
        std_val = np.std(mat_nor)
        self.assertTrue(std_val <= self.var)

    def test_normalization_fft(self):
        mat_nor = prep.normalization_fft(self.bck, sigma=5, pad=10)
        std_val = np.std(mat_nor)
        self.assertTrue(std_val <= self.var)

    def test_binarization(self):
        bck = 0.5 * np.random.rand(self.hei, self.wid)
        mat_bin = prep.binarization(self.mat_dots + bck, denoise=False)
        num_dots2 = ndi.label(mat_bin)[-1]
        self.assertTrue(self.num_dots == num_dots2)

    def test_check_num_dots(self):
        bck = 0.5 * np.random.rand(self.hei, self.wid)
        mat_bin = prep.binarization(self.mat_dots + bck, denoise=False)
        check = prep.check_num_dots(mat_bin)
        self.assertTrue(not check)

    def test_calc_size_distance(self):
        dot_size, dot_dist = prep.calc_size_distance(self.mat_dots)
        dot_size = int(dot_size)
        dot_dist = int(dot_dist)
        self.assertTrue(dot_size == 13 and dot_dist == 10)

    def test_select_dots_based_size(self):
        dot_size, dot_dist = prep.calc_size_distance(self.mat_dots)
        mat = np.copy(self.mat_dots)
        mat_label, _ = ndi.label(mat)
        list_dots = ndi.find_objects(mat_label)
        mat1 = np.zeros_like(mat)
        for i, j in enumerate(list_dots):
            mat1[j] = mat[j]
            if i < self.rem_dot:
                mat1[j] = ndi.binary_erosion(mat[j], iterations=2)
        mat2 = prep.select_dots_based_size(mat1, dot_size, 0.1)
        num_dots2 = ndi.label(mat2)[-1]
        self.assertTrue(num_dots2 == (self.num_dots - self.rem_dot))

    def test_select_dots_based_ratio(self):
        mat = np.zeros((64, 64), dtype=np.float32)
        mat[7:64:10, 7:64:10] = np.float32(1.0)
        mat[7, 7 + 1:64:10] = np.float32(1.0)
        mat[7, 7 + 2:64:10] = np.float32(1.0)
        mat[7, 7 - 1:64:10] = np.float32(1.0)
        mat[7, 7 - 2:64:10] = np.float32(1.0)
        mat = ndi.binary_dilation(mat, iterations=2)
        mat1 = prep.select_dots_based_ratio(mat, 0.05)
        num_dots2 = ndi.label(mat1)[-1]
        self.assertTrue(num_dots2 == (self.num_dots - self.rem_dot))

    def test_select_dots_based_distance(self):
        mat = np.zeros((64, 64), dtype=np.float32)
        mat[7:64:10, 7:64:10] = np.float32(1.0)
        mat[7 + 5, 7 + 5] = np.float32(1.0)
        mat = ndi.binary_dilation(mat, iterations=2)
        dot_dist = prep.calc_size_distance(mat)[-1]
        mat1 = prep.select_dots_based_distance(mat, dot_dist, ratio=0.05)
        num_dots2 = ndi.label(mat1)[-1]
        self.assertTrue(num_dots2 == self.num_dots)

    def test_calc_hor_slope(self):
        mat_rot = np.int16(
            np.ceil(ndi.rotate(self.mat_dots, -3.0, reshape=False, order=1)))
        hor_slope = prep.calc_hor_slope(mat_rot, ratio=1.0)
        angle = np.rad2deg(np.arctan(hor_slope))
        self.assertTrue(np.abs(angle - 3.0) <= 0.2)

    def test_calc_ver_slope(self):
        mat_rot = np.int16(
            np.ceil(ndi.rotate(self.mat_dots, -3.0, reshape=False, order=1)))
        ver_slope = prep.calc_ver_slope(mat_rot, ratio=1.0)
        angle = np.rad2deg(np.arctan(ver_slope))
        self.assertTrue(np.abs(angle + 3.0) <= 0.2)

    def test_group_dots_hor_lines(self):
        dot_dist = prep.calc_size_distance(self.mat_dots, ratio=0.9)[1]
        hor_slope = prep.calc_hor_slope(self.mat_dots, ratio=1.0)
        list_lines = prep.group_dots_hor_lines(self.mat_dots, hor_slope,
                                               dot_dist, ratio=0.1,
                                               num_dot_miss=3,
                                               accepted_ratio=0.9)
        num = np.sum(np.asarray([len(line) for line in list_lines]))
        self.assertTrue(num == self.num_dots)

    def test_group_dots_ver_lines(self):
        dot_dist = prep.calc_size_distance(self.mat_dots, ratio=0.9)[1]
        ver_slope = prep.calc_ver_slope(self.mat_dots, ratio=1.0)
        list_lines = prep.group_dots_ver_lines(self.mat_dots, ver_slope,
                                               dot_dist, ratio=0.1,
                                               num_dot_miss=3,
                                               accepted_ratio=0.9)
        num = np.sum(np.asarray([len(line) for line in list_lines]))
        self.assertTrue(num == self.num_dots)

    def test_remove_residual_dots_hor(self):
        mat1 = np.copy(self.mat_dots)
        mat1[9:11, 42:44] = np.float32(1.0)
        list_lines = prep.group_dots_hor_lines(mat1, 0.0, 10.0, ratio=0.3,
                                               num_dot_miss=3,
                                               accepted_ratio=0.8)
        num1 = np.sum(np.asarray([len(line) for line in list_lines]))
        list_lines2 = prep.remove_residual_dots_hor(list_lines, 0.0, 1.5)
        num2 = np.sum(np.asarray([len(line) for line in list_lines2]))
        self.assertTrue(num1 == num2 + 1)

    def test_remove_residual_dots_ver(self):
        mat1 = np.copy(self.mat_dots)
        mat1[42:44, 9:11] = np.float32(1.0)
        list_lines = prep.group_dots_ver_lines(mat1, 0.0, 10.0, ratio=0.3,
                                               num_dot_miss=3,
                                               accepted_ratio=0.8)
        num1 = np.sum(np.asarray([len(line) for line in list_lines]))
        list_lines2 = prep.remove_residual_dots_ver(list_lines, 0.0, 1.5)
        num2 = np.sum(np.asarray([len(line) for line in list_lines2]))
        self.assertTrue(num1 == num2 + 1)

    def test_calculate_threshold(self):
        mat = 0.2 * np.ones((64, 64))
        mat[16:30, 30: 40] = 1.0
        mat = mat + 0.2 * np.random.rand(64,64)
        threshold = prep.calculate_threshold(mat, bgr="dark")
        self.assertTrue(threshold > 0.5)
        mat = np.max(mat) - mat
        threshold = prep.calculate_threshold(mat, bgr="bright")
        self.assertTrue(threshold > 0.5)

    def test_make_parabola_mask(self):
        height, width, margin = 60, 80, 10
        mask1 = prep.make_parabola_mask(height, width, hor_curviness=0.3,
                                        ver_curviness=0.3, hor_margin=10,
                                        ver_margin=10, rotate=0.0)
        self.assertEqual(mask1.shape, (height, width))
        self.assertTrue(np.min(mask1)==0.0 and np.max(mask1)==1.0)
        num_col = len(np.where(np.mean(mask1, axis=0) > 0)[0])
        self.assertTrue(num_col != margin)
        num_row = len(np.where(np.mean(mask1, axis=1) > 0)[0])
        self.assertTrue(num_row != margin)
        self.assertRaises(ValueError, prep.make_parabola_mask,
                          height, width, 40, 40)
        rotate_angle = 45.0
        mask_rotated = prep.make_parabola_mask(height, width,
                                               hor_margin=margin,
                                               ver_margin=margin,
                                               rotate=rotate_angle)
        self.assertEqual(mask_rotated.shape, (height, width))
        self.assertTrue(np.any(mask_rotated == 1.0))

    def test_remove_points_using_parabola_mask(self):
        f_alias = prep.remove_points_using_parabola_mask
        height, width, margin = 60, 80, 5
        points = np.array([[25, 25], [30, 35], [40, 60]], dtype=np.float32)
        filtered_points = f_alias(points, height, width,
                                  hor_curviness=0.1, ver_curviness=0.1,
                                  hor_margin=margin, ver_margin=margin)
        self.assertTrue(len(filtered_points) == len(points))

        points_outside = np.array([[0, 0], [59, 79], [59, 79]],
                                  dtype=np.float32)
        filtered_points = f_alias(points_outside, height, width,
                                  hor_curviness=0.1, ver_curviness=0.1,
                                  hor_margin=margin, ver_margin=margin)
        self.assertEqual(len(filtered_points), 0)

        mixed_points = np.array([[0, 0], [30, 30], [40, 7]], dtype=np.float32)
        filtered_points = f_alias(mixed_points, height, width,
                                  hor_curviness=0.1, ver_curviness=0.1,
                                  hor_margin=margin, ver_margin=margin)
        self.assertTrue(len(filtered_points) < len(mixed_points))

    def test_get_points_dot_pattern(self):
        points = prep.get_points_dot_pattern(self.mat_dots, binarize=False)
        self.assertTrue(len(points) == self.num_dots)

        mat_noise = self.mat_dots + 0.2 * np.random.rand(self.hei, self.wid)
        points = prep.get_points_dot_pattern(mat_noise, binarize=True)
        self.assertTrue(len(points) == self.num_dots)

        with self.assertRaises(ValueError):
            prep.get_points_dot_pattern(mat_noise, binarize=False)

    def test_rotate_points(self):
        points = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])
        rotated_points = prep.rotate_points(points, 90)
        expected_points = np.array([[0, -1], [1, 0], [0, 1], [-1, 0]])
        np.testing.assert_almost_equal(rotated_points, expected_points,
                                       decimal=6)

        points = np.array([[1, 0], [0, 1]])
        rotated_points = prep.rotate_points(points, 45)
        expected_points = np.array([[0.7071, -0.7071], [0.7071, 0.7071]])
        np.testing.assert_almost_equal(rotated_points, expected_points,
                                       decimal=4)

        points = np.array([[1, 0]])
        rotated_points = prep.rotate_points(points, np.pi / 2,
                                            degree_unit=False)
        expected_points = np.array([[0, -1]])
        np.testing.assert_almost_equal(rotated_points, expected_points,
                                       decimal=6)

    def test_remove_subset_points(self):
        points = np.array([[1, 2], [3, 4], [5, 6]])
        selected_points = np.array([[3, 4]])
        filtered_points = prep.remove_subset_points(selected_points, points)
        expected_points = np.array([[1, 2], [5, 6]])
        np.testing.assert_array_equal(filtered_points, expected_points)

        points = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
        selected_points = np.array([[1, 2], [7, 8]])
        filtered_points = prep.remove_subset_points(selected_points, points)
        expected_points = np.array([[3, 4], [5, 6]])
        np.testing.assert_array_equal(filtered_points, expected_points)

        points = np.array([[1, 2], [3, 4], [5, 6]])
        selected_points = np.array([[7, 8]])
        filtered_points = prep.remove_subset_points(selected_points, points)
        np.testing.assert_array_equal(filtered_points, points)

    def test_group_dots_based_polyfit(self):
        num_hor_line = 27
        num_ver_line = 37
        current_dir = os.path.dirname(__file__)
        data_path = os.path.join(current_dir, "data_for_test",
                                 "data_for_grouping.pkl")
        data = losa.load_python_list(data_path)
        slope_hor, dist_hor = data[0]
        slope_ver, dist_ver = data[1]
        points = np.asarray(data[2])
        list_hor_lines = prep.group_dots_hor_lines_based_polyfit(points,
                                                                 slope_hor,
                                                                 dist_hor,
                                                                 order=2)
        list_ver_lines = prep.group_dots_ver_lines_based_polyfit(points,
                                                                 slope_ver,
                                                                 dist_ver,
                                                                 order=2)
        self.assertTrue(len(list_hor_lines) == num_hor_line)
        self.assertTrue(len(list_ver_lines) == num_ver_line)
        self.assertTrue(len(list_hor_lines[0]) == num_ver_line)
        self.assertTrue(len(list_ver_lines[0]) == num_hor_line)
