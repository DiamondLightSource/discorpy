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

    def test_unwarp_slice_backward(self):
        x0, y0 = (self.wid // 2, self.hei // 2)
        list_fact = [1.0, 3.0 * 10 ** (-3)]
        mat = np.zeros((self.hei, self.wid), dtype=np.float32)
        mat[:, 6:-8:8] = 1.0
        mat = np.float32(ndi.binary_dilation(np.int16(mat), iterations=1))
        mat3d = np.zeros((10, self.hei, self.wid), dtype=np.float32)
        mat3d[:] = mat
        chunk_warp = post.unwarp_chunk_slices_backward(mat3d, x0, y0, list_fact,
                                                       y0 - 5, y0 + 5)
        error1 = np.max(mat3d[:, y0 - 5, :] - chunk_warp[:, 0, :])
        error2 = np.max(mat3d[:, y0 + 5, :] - chunk_warp[:, -1, :])
        self.assertTrue(error1 > 0.1 and error2 > 0.1)

    def test_calc_residual_hor(self):
        list_clines = post.unwarp_line_forward(self.list_dlines, self.x0,
                                               self.y0, self.list_ffact)
        list_res = post.calc_residual_hor(list_clines, self.x0, self.y0)
        error = np.max(list_res[:, 1])
        self.assertTrue(error < 0.5)

    def test_calc_residual_ver(self):
        dot_dist = 2.0
        list_fact = [1.0, -2.0 * 10 ** (-3)]
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
        error = np.max(list_res[:, 1])
        self.assertTrue(error < 0.5)
