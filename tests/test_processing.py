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
Tests for methods in processing.py
"""

import unittest
import numpy as np
import discorpy.proc.processing as proc


class ProcessingMethods(unittest.TestCase):

    def setUp(self):
        self.x0 = 33.0
        self.y0 = 35.0
        height = 64
        width = 64
        dot_dist = 2.0
        list_fact = [1.0, -2.0 * 10 ** (-3)]
        list_hor_lines = [
            [[height - y, x] for x in np.arange(1, width, dot_dist)] for y
            in np.arange(1, height, dot_dist)]
        self.list_hor_dlines = []
        for line in list_hor_lines:
            line1 = np.asarray(line)
            xu_list = line1[:, 1] - self.x0
            yu_list = line1[:, 0] - self.y0
            ru_list = np.sqrt(xu_list ** 2 + yu_list ** 2)
            flist = np.sum(np.asarray(
                [factor * ru_list ** i for i, factor in enumerate(list_fact)]),
                axis=0)
            xd_list = self.x0 + xu_list * flist
            yd_list = self.y0 + yu_list * flist
            self.list_hor_dlines.append(
                np.asarray(list(zip(yd_list, xd_list))))
        list_ver_lines = [
            [[height - y, x] for y in np.arange(1, height, dot_dist)] for x
            in np.arange(1, width, dot_dist)]
        self.list_ver_dlines = []
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
            self.list_ver_dlines.append(
                np.asarray(list(zip(yd_list, xd_list))))
        self.dot_dist = dot_dist
        self.list_fact = list_fact

    def test_find_cod_coarse(self):
        x_cod, y_cod = proc.find_cod_coarse(self.list_hor_dlines,
                                            self.list_ver_dlines)
        self.assertTrue(
            (np.abs(x_cod - self.x0) < self.dot_dist) and (
                    np.abs(y_cod - self.y0) < self.dot_dist))

    def test_find_cod_fine(self):
        x_cod, y_cod = proc.find_cod_coarse(self.list_hor_dlines,
                                            self.list_ver_dlines)
        x_cod, y_cod = proc.find_cod_fine(self.list_hor_dlines,
                                          self.list_ver_dlines, x_cod, y_cod,
                                          self.dot_dist)
        self.assertTrue(isinstance(x_cod, float) and isinstance(y_cod, float))

    def test_calc_coef_backward(self):
        x_cod, y_cod = proc.find_cod_coarse(self.list_hor_dlines,
                                            self.list_ver_dlines)
        list_fact = proc.calc_coef_backward(self.list_hor_dlines,
                                            self.list_ver_dlines,
                                            x_cod, y_cod, 2)
        error1 = np.abs((list_fact[0] - self.list_fact[0]) / self.list_fact[0])
        error2 = np.abs((list_fact[1] - self.list_fact[1]) / self.list_fact[1])
        self.assertTrue(error1 < 0.1 and error2 < 0.1)

        list_fact = proc.calc_coef_backward(self.list_hor_dlines,
                                            self.list_ver_dlines,
                                            x_cod, y_cod, 2, optimizing=True)
        error1 = np.abs((list_fact[0] - self.list_fact[0]) / self.list_fact[0])
        error2 = np.abs((list_fact[1] - self.list_fact[1]) / self.list_fact[1])
        self.assertTrue(error1 < 0.1 and error2 < 0.15)

    def test_calc_coef_fordward(self):
        x_cod, y_cod = proc.find_cod_coarse(self.list_hor_dlines,
                                            self.list_ver_dlines)
        list_fact = proc.calc_coef_forward(self.list_hor_dlines,
                                           self.list_ver_dlines,
                                           x_cod, y_cod, 2)
        error1 = np.abs((list_fact[0] - self.list_fact[0]) / self.list_fact[0])
        error2 = np.abs((list_fact[1] + self.list_fact[1]) / self.list_fact[1])
        self.assertTrue(error1 < 0.1 and error2 < 0.2)

    def test_calc_coef_backward_from_forward(self):
        x_cod, y_cod = proc.find_cod_coarse(self.list_hor_dlines,
                                            self.list_ver_dlines)
        list_ffact, list_bfact = proc.calc_coef_backward_from_forward(
            self.list_hor_dlines, self.list_ver_dlines, x_cod, y_cod, 2)
        error1 = np.abs(
            (list_ffact[0] - self.list_fact[0]) / self.list_fact[0])
        error2 = np.abs(
            (list_ffact[1] + self.list_fact[1]) / self.list_fact[1])
        error3 = np.abs(
            (list_bfact[0] - self.list_fact[0]) / self.list_fact[0])
        error4 = np.abs(
            (list_bfact[1] - self.list_fact[1]) / self.list_fact[1])
        self.assertTrue(
            error1 < 0.1 and error2 < 0.2 and error3 < 0.1 and error4 < 0.2)

    def test_find_cod_bailey(self):
        x_cod, y_cod = proc.find_cod_bailey(self.list_hor_dlines,
                                            self.list_ver_dlines)
        self.assertTrue((np.abs(x_cod - self.x0) < 1.0) and
                        (np.abs(y_cod - self.y0) < 1.0))

    def test_regenerate_grid_points_parabola(self):
        list_hline1, list_vline1 = proc.regenerate_grid_points_parabola(
            self.list_hor_dlines, self.list_ver_dlines, perspective=True)
        list_hline2, list_vline2 = proc.regenerate_grid_points_parabola(
            self.list_hor_dlines, self.list_ver_dlines, perspective=False)
        num_hpoint1 = np.sum(np.asarray([len(line) for line in list_hline1]))
        num_vpoint1 = np.sum(np.asarray([len(line) for line in list_vline1]))
        num_hpoint2 = np.sum(np.asarray([len(line) for line in list_hline2]))
        num_vpoint2 = np.sum(np.asarray([len(line) for line in list_vline2]))
        self.assertTrue(num_vpoint1 == num_vpoint2 and
                        num_hpoint1 == num_hpoint2 and
                        num_hpoint2 == num_vpoint1)

    def test_regenerate_grid_points_linear(self):
        list_hline, list_vline = proc.regenerate_grid_points_linear(
            self.list_hor_dlines, self.list_ver_dlines)
        num_hpoint = np.sum(np.asarray([len(line) for line in list_hline]))
        num_vpoint = np.sum(np.asarray([len(line) for line in list_vline]))
        self.assertTrue(num_vpoint == num_hpoint)

    def test_generate_undistorted_perspective_lines(self):
        f_alias = proc.generate_undistorted_perspective_lines
        uhor_lines = f_alias(self.list_hor_dlines, self.list_ver_dlines,
            equal_dist=True, optimizing=False)[0]
        num_hpoint1 = np.sum(np.asarray([len(line) for line in uhor_lines]))
        uhor_lines = f_alias(self.list_hor_dlines, self.list_ver_dlines,
            equal_dist=False, optimizing=True)[0]
        num_hpoint2 = np.sum(np.asarray([len(line) for line in uhor_lines]))
        self.assertTrue(num_hpoint1 == num_hpoint2)

        uhor_lines = f_alias(self.list_hor_dlines, self.list_ver_dlines,
                             scale="max")[0]
        num_hpoint1 = np.sum(np.asarray([len(line) for line in uhor_lines]))
        uhor_lines = f_alias(self.list_hor_dlines, self.list_ver_dlines,
                             scale="min")[0]
        num_hpoint2 = np.sum(np.asarray([len(line) for line in uhor_lines]))
        uhor_lines = f_alias(self.list_hor_dlines, self.list_ver_dlines,
                             scale="median")[0]
        num_hpoint3 = np.sum(np.asarray([len(line) for line in uhor_lines]))
        uhor_lines = f_alias(self.list_hor_dlines, self.list_ver_dlines,
                             scale=1.0)[0]
        num_hpoint4 = np.sum(np.asarray([len(line) for line in uhor_lines]))

        self.assertTrue(num_hpoint1 == num_hpoint2
                        and num_hpoint1 == num_hpoint3
                        and num_hpoint1 == num_hpoint4)

    def test_generate_source_target_perspective_points(self):
        num_points = np.sum(
            np.asarray([len(line) for line in self.list_hor_dlines]))
        s_points, t_points = proc.generate_source_target_perspective_points(
            self.list_hor_dlines, self.list_ver_dlines)
        self.assertTrue(num_points == len(s_points) and
                        num_points == len(t_points))

    def test_generate_4_source_target_perspective_points(self):
        f_alias = proc.generate_4_source_target_perspective_points
        points = [[5, 5], [6, 50], [40, 7], [45, 57]]
        t_points0 = [[3.58143506, 2.58661269], [7.83739762, 50.02633148],
                     [40.77223206, -0.74988769], [45.02819462, 46.6898311]]

        s_points, t_points = f_alias(points, scale="mean", equal_dist=False)
        num2 = np.mean(np.abs(
            np.float32(t_points) - np.asarray(t_points0, dtype=np.float32)))
        self.assertTrue(num2 <= 1.0e-6)

        s_points, t_points = f_alias(points, scale="max", equal_dist=True)
        self.assertTrue(len(s_points) == 4 and len(t_points) == 4)

        s_points, t_points = f_alias(points, scale="min", equal_dist=True)
        self.assertTrue(len(s_points) == 4 and len(t_points) == 4)

        s_points, t_points = f_alias(points, scale="median", equal_dist=True)
        self.assertTrue(len(s_points) == 4 and len(t_points) == 4)

        s_points, t_points = f_alias(points, scale=1.0, equal_dist=True)
        self.assertTrue(len(s_points) == 4 and len(t_points) == 4)

        s_points2, _ = f_alias(points, scale=1.0, equal_dist=True,
                               input_order="xy")
        self.assertTrue(
            np.mean(np.abs(s_points[:, 0] - s_points2[:, 0])) > 1e-6)

    def test_calc_perspective_coefficients(self):
        s_points = [[5, 5], [6, 50], [40, 7], [45, 57]]
        t_points = [[3.58143506, 2.58661269], [7.83739762, 50.02633148],
                    [40.77223206, -0.74988769], [45.02819462, 46.6898311]]
        backward_coef = proc.calc_perspective_coefficients(s_points, t_points,
                                                           mapping="backward")
        b_coef0 = [8.31034232e-01, 1.11425384e-01, 2.38551326e+00,
                   -6.50926172e-02, 8.30299316e-01, 2.12884603e+00,
                   -1.67982946e-03, -2.46465092e-03]

        forward_coef = proc.calc_perspective_coefficients(s_points, t_points,
                                                          mapping="forward")
        f_coef0 = [1.19832778e+00, -1.68236843e-01, -2.50047647e+00,
                   8.82260677e-02, 1.19760396e+00, -2.75997890e+00,
                   2.23043277e-03, 2.66906651e-03]
        num1 = np.mean(np.abs(backward_coef - np.asarray(b_coef0)))
        num2 = np.mean(np.abs(forward_coef - np.asarray(f_coef0)))
        self.assertTrue(num1 <= 1.0e-6 and num2 <= 1.0e-6)

    def test_update_center(self):
        lines1 = np.asarray(
            [[[1, 2], [1, 6], [1, 10]], [[3, 2], [3, 6], [3, 10]]])
        xc, yc = 5, 6
        lines2 = proc.update_center(lines1, xc, yc)
        results = np.concatenate((lines2[0], lines2[1]),
                                 axis=0) - np.concatenate(
            (lines1[0], lines1[1]), axis=0)
        num1 = np.abs(np.mean(results[:, 0]) - yc)
        num2 = np.abs(np.mean(results[:, 1]) - xc)
        self.assertTrue(num1 == 0.0 and num2 == 0.0)

    def test_transform_coef_backward_and_forward(self):
        ffacts1 = np.asarray([1.0, -2.0e-3, 5.0e-6])
        bfacts1 = np.asarray([1.0, 2.0e-3, 2.0e-6])
        points = [[i, j] for i in range(30) for j in range(30)]
        bfacts2 = proc.transform_coef_backward_and_forward(ffacts1,
                                                           mapping="backward",
                                                           ref_points=points)
        ffacts2 = proc.transform_coef_backward_and_forward(bfacts2,
                                                           mapping="forward",
                                                           ref_points=points)
        num1 = np.mean(np.abs(ffacts2 - ffacts1))
        num2 = np.mean(np.abs(bfacts2 - bfacts1))
        self.assertTrue(num1 <= 1.0e-3 and num2 <= 1.0e-3)
