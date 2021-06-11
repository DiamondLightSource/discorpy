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
            self.list_hor_dlines.append(np.asarray(list(zip(yd_list, xd_list))))
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
            self.list_ver_dlines.append(np.asarray(list(zip(yd_list, xd_list))))
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
        error1 = np.abs((list_ffact[0] - self.list_fact[0]) / self.list_fact[0])
        error2 = np.abs((list_ffact[1] + self.list_fact[1]) / self.list_fact[1])
        error3 = np.abs((list_bfact[0] - self.list_fact[0]) / self.list_fact[0])
        error4 = np.abs((list_bfact[1] - self.list_fact[1]) / self.list_fact[1])
        self.assertTrue(
            error1 < 0.1 and error2 < 0.2 and error3 < 0.1 and error4 < 0.2)
