# ============================================================================
# ============================================================================
# Copyright (c) 2023 Diamond Light Source Ltd. All rights reserved.
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
Tests for methods in utility.py
"""

import unittest
import numpy as np
import discorpy.util.utility as util
import discorpy.proc.processing as proc

class UtilityMethods(unittest.TestCase):

    def test_make_circle_mask(self):
        width, ratio = 10, 0.5
        mask = util.make_circle_mask(width, ratio)
        self.assertTrue(mask.shape == (width, width) and np.mean(mask) > 0.0)

    def test_make_dot_pattern(self):
        height, width = 10, 20
        dot_distance, dot_size, margin = 2, 3, 1
        pattern = util.make_dot_pattern(height, width, dot_distance, dot_size,
                                        margin)
        self.assertTrue(pattern.shape == (height, width) and
                        np.mean(pattern) > 0.0)

    def test_make_line_pattern(self):
        height, width = 10, 20
        line_distance, line_size, margin = 2, 3, 1
        pattern = util.make_line_pattern(height, width, line_distance,
                                         line_size, margin)
        self.assertTrue(pattern.shape == (height, width) and
                        np.mean(pattern) > 0.0)

    def test_make_chessboard(self):
        height, width = 10, 20
        size, margin = 2, 1
        margin_grayscale = 0.5
        chessboard = util.make_chessboard(height, width, size, margin,
                                          margin_grayscale)
        self.assertEqual(chessboard.shape, (height, width))
        self.assertTrue(np.mean(chessboard) > 0.0)

    def test_find_point_to_point(self):
        points = (3, 5)
        xcenter = ycenter = 4.0
        list_fact = [1.0, 0.5]
        output_order = "xy"
        xo, yo = util.find_point_to_point(points, xcenter, ycenter,
                                          list_fact, output_order)
        self.assertTrue(np.abs(xo - 5.7) < 0.1)
        self.assertTrue(np.abs(yo - 2.3) < 0.1)

    def test_unwarp_color_image_backward(self):
        height, width, channel = 40, 60, 3
        mat = np.ones((height, width), dtype=np.float32)
        mat[:5] = 0
        mat[-5:] = 0
        mat[:, :5] = 0
        mat[:, -5:] = 0
        xcenter = 30.0
        ycenter = 20.0
        list_fact = [1.0, 0.5]
        order, mode, pad, pad_mode = 1, "constant", False, "constant"
        corrected_mat = util.unwarp_color_image_backward(mat, xcenter,
                                                         ycenter, list_fact,
                                                         order, mode, pad,
                                                         pad_mode)
        self.assertEqual(corrected_mat.shape, (height, width))
        self.assertTrue(np.mean(corrected_mat) < 1.0)

        pad = True
        list_fact = [1.0, 0.1, 0.01]
        ref = [[i - ycenter, j - xcenter] for i in range(0, height, 10)
               for j in range(0, width, 10)]
        list_tfact = proc.transform_coef_backward_and_forward(list_fact,
                                                              ref_points=ref)
        corrected_mat = util.unwarp_color_image_backward(mat, xcenter,
                                                         ycenter, list_tfact,
                                                         order, mode, pad,
                                                         pad_mode)
        self.assertTrue(corrected_mat.shape != (height, width))

        pad = 10
        corrected_mat = util.unwarp_color_image_backward(mat, xcenter,
                                                         ycenter, list_tfact,
                                                         order, mode, pad,
                                                         pad_mode)
        self.assertTrue(corrected_mat.shape == (height + 20, width + 20))

        pad = (10, 20, 5, 5)
        corrected_mat = util.unwarp_color_image_backward(mat, xcenter,
                                                         ycenter, list_tfact,
                                                         order, mode, pad,
                                                         pad_mode)
        self.assertTrue(corrected_mat.shape == (height + 30, width + 10))

        pad = False
        mat = np.ones((height, width, channel), dtype=np.float32)
        mat[:5, :, :] = 0
        mat[-5:, :, :] = 0
        mat[:, :5, :] = 0
        mat[:, -5:, :] = 0
        corrected_mat = util.unwarp_color_image_backward(mat, xcenter,
                                                         ycenter, list_fact,
                                                         order, mode, pad,
                                                         pad_mode)
        self.assertEqual(corrected_mat.shape, (height, width, channel))
        self.assertTrue(np.mean(corrected_mat) < 1.0)

        pad = 10
        mat = np.ones((height, width, channel), dtype=np.float32)
        mat[:5, :, :] = 0
        mat[-5:, :, :] = 0
        mat[:, :5, :] = 0
        mat[:, -5:, :] = 0
        corrected_mat = util.unwarp_color_image_backward(mat, xcenter,
                                                         ycenter, list_fact,
                                                         order, mode, pad,
                                                         pad_mode)
        self.assertEqual(corrected_mat.shape,
                         (height + 20, width + 20, channel))
        self.assertTrue(np.mean(corrected_mat) < 1.0)
