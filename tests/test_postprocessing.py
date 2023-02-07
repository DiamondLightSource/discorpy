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
Tests for methods in postprocessing.py
"""

import unittest
import numpy as np
import scipy.ndimage as ndi
import discorpy.post.postprocessing as post
import discorpy.proc.processing as proc


class PostprocessingMethods(unittest.TestCase):

    def setUp(self):
        self.x0, self.y0 = (33.5, 35.5)
        height, width = (64, 64)
        dot_dist = 2.0
        list_fact = [1.0, -2.0 * 10 ** (-3)]
        self.list_ffact = [1.0, 2.0 * 10 ** (-3)]
        list_lines = [[[height - y, x] for x in np.arange(1, width, dot_dist)]
                      for y in np.arange(1, height, dot_dist)]
        list_dlines = []
        for line in list_lines:
            line1 = np.asarray(line)
            xu_list = line1[:, 1] - self.x0
            yu_list = line1[:, 0] - self.y0
            ru_list = np.sqrt(xu_list ** 2 + yu_list ** 2)
            flist = np.sum(np.asarray(
                [factor * ru_list ** i for i, factor in enumerate(list_fact)]),
                axis=0)
            xd_list = self.x0 + xu_list * flist
            yd_list = self.y0 + yu_list * flist
            list_dlines.append(np.asarray(list(zip(yd_list, xd_list))))
        self.list_lines = [np.asarray(line) for line in list_lines]
        self.list_dlines = list_dlines
        self.list_fact = list_fact
        self.hei, self.wid = (height, width)

    def test_unwarp_line_forward(self):
        list_clines = post.unwarp_line_forward(self.list_dlines, self.x0,
                                               self.y0, self.list_ffact)
        error = np.max(np.asarray(
            [np.max(np.abs(line - self.list_lines[i])) for (i, line) in
             enumerate(list_clines)]))
        self.assertTrue(error <= 1.0)

    def test_unwarp_line_backward(self):
        list_clines = post.unwarp_line_backward(self.list_dlines, self.x0,
                                                self.y0, self.list_fact)
        error = np.max(np.asarray(
            [np.max(np.abs(line - self.list_lines[i])) for (i, line) in
             enumerate(list_clines)]))
        self.assertTrue(error <= 1.0)

    def test_unwarp_image_backward(self):
        x0, y0 = (self.wid // 2, self.hei // 2)
        list_fact = [1.0, 3.0 * 10 ** (-3)]
        mat = np.zeros((self.hei, self.wid), dtype=np.float32)
        mat[4:-3, 4:-3] = 1.0
        mat_warp = post.unwarp_image_backward(mat, x0, y0, list_fact)
        vals = np.mean(mat_warp, axis=0)[11:-10]
        pos = len(vals) // 2
        self.assertTrue(vals[0] < vals[pos] and vals[-1] < vals[pos])

    def test_unwarp_image_forward(self):
        x0, y0 = (self.wid // 2, self.hei // 2)
        list_fact = [1.0, -6.0 * 10 ** (-3)]
        mat = np.zeros((self.hei, self.wid), dtype=np.float32)
        mat[4:-3, 4:-3] = 1.0
        mat_warp = ndi.gaussian_filter(
            post.unwarp_image_forward(mat, x0, y0, list_fact), 2)
        vals = np.mean(mat_warp, axis=0)[11:-10]
        pos = len(vals) // 2
        self.assertTrue(vals[0] < vals[pos] and vals[-1] < vals[pos])

    def test_unwarp_slice_backward(self):
        x0, y0 = (self.wid // 2, self.hei // 2)
        list_fact = [1.0, 3.0 * 10 ** (-3)]
        mat = np.zeros((self.hei, self.wid), dtype=np.float32)
        mat[:, 6:-8:8] = 1.0
        mat = np.float32(ndi.binary_dilation(np.int16(mat), iterations=1))
        mat3d = np.zeros((10, self.hei, self.wid), dtype=np.float32)
        mat3d[:] = mat
        slice_warp = post.unwarp_slice_backward(mat3d, x0, y0, list_fact, y0)
        error = np.max(mat3d[:, y0, :] - slice_warp)
        self.assertTrue(error > 0.1)

    def test_unwarp_chunk_slice_backward(self):
        x0, y0 = (self.wid // 2, self.hei // 2)
        list_fact = [1.0, 3.0 * 10 ** (-3)]
        mat = np.zeros((self.hei, self.wid), dtype=np.float32)
        mat[:, 6:-8:8] = 1.0
        mat = np.float32(ndi.binary_dilation(np.int16(mat), iterations=1))
        mat3d = np.zeros((10, self.hei, self.wid), dtype=np.float32)
        mat3d[:] = mat
        chunk_warp = post.unwarp_chunk_slices_backward(mat3d, x0, y0,
                                                       list_fact, y0 - 5,
                                                       y0 + 5)
        error1 = np.max(mat3d[:, y0 - 5, :] - chunk_warp[:, 0, :])
        error2 = np.max(mat3d[:, y0 + 5, :] - chunk_warp[:, -1, :])
        self.assertTrue(error1 > 0.1 and error2 > 0.1)

    def test_calc_residual_hor(self):
        list_clines = post.unwarp_line_forward(self.list_dlines, self.x0,
                                               self.y0, self.list_ffact)
        list_res = post.calc_residual_hor(list_clines, self.x0, self.y0)
        check = post.check_distortion(list_res)
        error = np.max(list_res[:, 1])
        self.assertTrue(error < 0.5 and check is False)

    def test_calc_residual_ver(self):
        dot_dist = 2.0
        list_fact = [1.0, -2.0 * 10 ** (-2)]
        list_ver_lines = [
            [[self.hei - y, x] for y in np.arange(1, self.hei, dot_dist)] for x
            in np.arange(1, self.wid, dot_dist)]
        list_ver_dlines = []
        for line in list_ver_lines:
            line1 = np.asarray(line)
            xu_list = line1[:, 1] - self.x0
            yu_list = line1[:, 0] - self.y0
            ru_list = np.sqrt(xu_list ** 2 + yu_list ** 2)
            flist = np.sum(np.asarray(
                [factor * ru_list ** i for i, factor in enumerate(list_fact)]),
                axis=0)
            xd_list = self.x0 + xu_list * flist
            yd_list = self.y0 + yu_list * flist
            list_ver_dlines.append(np.asarray(list(zip(yd_list, xd_list))))

        list_ver_clines = post.unwarp_line_backward(list_ver_dlines, self.x0,
                                                    self.y0, list_fact)
        list_res = post.calc_residual_ver(list_ver_clines, self.x0, self.y0)
        check = post.check_distortion(list_res)
        error = np.max(list_res[:, 1])
        self.assertTrue(error > 1.0 and check is True)

    def test_correct_perspective_line(self):
        hor_line1 = np.asarray([[10.0 + 2 * i / 32, i] for i in range(32)])
        hor_line2 = np.asarray([[26.0 - 5 * i / 32, i] for i in range(32)])
        list_hor_lines = [hor_line1, hor_line2]
        ver_line1 = np.asarray([[i, 0.0 + 3 * i / 32] for i in range(32)])
        ver_line2 = np.asarray([[i, 20.0 - 3 * i / 32] for i in range(32)])
        list_ver_lines = [ver_line1, ver_line2]

        f_alias = proc.generate_source_target_perspective_points
        source_points, target_points = f_alias(list_hor_lines, list_ver_lines,
                                               equal_dist=False, scale="mean",
                                               optimizing=False)

        f_alias1 = proc.calc_perspective_coefficients
        pers_fcoef = f_alias1(source_points, target_points, mapping="forward")
        pers_bcoef = f_alias1(source_points, target_points, mapping="backward")

        f_alias2 = post.correct_perspective_line
        list_hor_clines = f_alias2(list_hor_lines, pers_fcoef)
        list_ver_clines = f_alias2(list_ver_lines, pers_fcoef)
        list_hor_dlines = f_alias2(list_hor_clines, pers_bcoef)
        list_ver_dlines = f_alias2(list_ver_clines, pers_bcoef)

        hor_dline1 = list_hor_dlines[0]
        hor_dline2 = list_hor_dlines[1]
        ver_dline1 = list_ver_dlines[0]
        ver_dline2 = list_ver_dlines[1]

        h_slope1a = np.polyfit(hor_line1[:, 1], hor_line1[:, 0], 1)[0]
        h_slope1b = np.polyfit(hor_dline1[:, 1], hor_dline1[:, 0], 1)[0]
        h_slope2a = np.polyfit(hor_line2[:, 1], hor_line2[:, 0], 1)[0]
        h_slope2b = np.polyfit(hor_dline2[:, 1], hor_dline2[:, 0], 1)[0]

        v_slope1a = np.polyfit(ver_line1[:, 0], ver_line1[:, 1], 1)[0]
        v_slope1b = np.polyfit(ver_dline1[:, 0], ver_dline1[:, 1], 1)[0]
        v_slope2a = np.polyfit(ver_line2[:, 0], ver_line2[:, 1], 1)[0]
        v_slope2b = np.polyfit(ver_dline2[:, 0], ver_dline2[:, 1], 1)[0]

        error1 = np.abs(h_slope1a - h_slope1b)
        error2 = np.abs(h_slope2a - h_slope2b)
        error3 = np.abs(v_slope1a - v_slope1b)
        error4 = np.abs(v_slope2a - v_slope2b)
        val = 1.0e-6
        self.assertTrue(error1 < val and error2 < val and error3 < val
                        and error4 < val)

    def test_correct_perspective_image(self):
        hor_line1 = np.asarray([[10.0 + 2 * i / 32, i] for i in range(64)])
        hor_line2 = np.asarray([[60.0 - 5 * i / 32, i] for i in range(64)])
        list_hor_lines = [hor_line1, hor_line2]
        ver_line1 = np.asarray([[i, 5.0 + 3 * i / 32] for i in range(64)])
        ver_line2 = np.asarray([[i, 60.0 - 3 * i / 32] for i in range(64)])
        list_ver_lines = [ver_line1, ver_line2]

        f_alias = proc.generate_source_target_perspective_points
        source_points, target_points = f_alias(list_hor_lines, list_ver_lines,
                                               equal_dist=False, scale="mean",
                                               optimizing=False)

        f_alias1 = proc.calc_perspective_coefficients
        pers_fcoef = f_alias1(source_points, target_points, mapping="forward")
        pers_bcoef = f_alias1(source_points, target_points, mapping="backward")

        mat = np.zeros((64, 64), dtype=np.float32)
        pos1, pos2 = 10, 45
        mat[pos1:pos1 + 1] = np.float32(1.0)
        mat[pos2:pos2 + 1] = np.float32(0.5)
        list0 = np.mean(mat, axis=1)
        num0 = len(list0[list0 > 0])
        max_pos0 = np.argmax(list0)

        mat_cor1 = post.correct_perspective_image(mat, pers_bcoef)
        list1 = np.mean(mat_cor1, axis=1)
        num1 = len(list1[list1 > 0])
        self.assertTrue(num1 > num0)

        mat_cor2 = post.correct_perspective_image(mat_cor1, pers_fcoef)
        list2 = np.mean(mat_cor2, axis=1)
        max_pos2 = np.argmax(list2)
        self.assertTrue(max_pos0 == max_pos2)
