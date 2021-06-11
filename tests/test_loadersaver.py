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
Tests for methods in loadersaver.py
"""

import os
import sys
import unittest
import shutil
import h5py
import numpy as np
import discorpy.losa.loadersaver as losa


class LoaderSaverMethods(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not os.path.isdir("data"):
            os.makedirs("data")

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir("data"):
            shutil.rmtree("data")

    def test_load_image(self):
        file_path = "data/img.tif"
        losa.save_image(file_path, np.random.rand(64, 64))
        mat = losa.load_image(file_path)
        self.assertTrue(mat.shape == (64, 64))

    def test_load_hdf_file(self):
        file_path = "data/data.hdf"
        ifile = h5py.File(file_path, "w")
        ifile.create_dataset("entry/data", data=np.random.rand(64, 64))
        ifile.close()
        data = losa.load_hdf_file(file_path, "entry/data")[:]
        self.assertTrue(isinstance(data, np.ndarray) and data.shape == (64, 64))

    def test_load_hdf_object(self):
        file_path = "data/data1.hdf"
        ifile = h5py.File(file_path, "w")
        ifile.create_dataset("entry/data", data=np.random.rand(64, 64))
        ifile.close()
        data = losa.load_hdf_object(file_path, "entry/data")
        self.assertTrue(isinstance(data, object) and data.shape == (64, 64))

    def test_save_image(self):
        file_path = "data/img.tif"
        losa.save_image(file_path, np.random.rand(64, 64))
        self.assertTrue(os.path.isfile(file_path))

    def test_save_plot_image(self):
        ver = sys.version[:3]
        check = True
        if ver != "2.7":
            file_path = "data/plot.png"
            list_lines = [np.asarray([(i, j) for j in range(5, 64, 5)]) for i in
                         range(5, 64, 5)]
            losa.save_plot_image(file_path, list_lines, 64, 64, dpi=100)
            check = os.path.isfile(file_path)
        self.assertTrue(check)

    def test_save_residual_plot(self):
        file_path = "data/plot1.png"
        list_data = np.ones((64, 2), dtype=np.float32)
        list_data[:, 0] = 0.5 * np.random.rand(64)
        list_data[:, 1] = np.arange(64)
        losa.save_residual_plot(file_path, list_data, 64, 64, dpi=100)
        self.assertTrue(os.path.isfile(file_path))

    def test_save_hdf_file(self):
        file_path = "data/data2.hdf"
        losa.save_hdf_file(file_path, np.random.rand(64, 64))
        self.assertTrue(os.path.isfile(file_path))

    def test_open_hdf_stream(self):
        data_out = losa.open_hdf_stream("data/data.hdf", (64, 64))
        self.assertTrue(isinstance(data_out, object))

    def test_save_metadata_txt(self):
        file_path = "data/coef.txt"
        losa.save_metadata_txt(file_path, 31, 32, (1.0, 0.0))
        self.assertTrue(os.path.isfile(file_path))

    def test_load_metadata_txt(self):
        file_path = "data/coef1.txt"
        losa.save_metadata_txt(file_path, 31.0, 32.0, [1.0, 0.0])
        (x, y, facts) = losa.load_metadata_txt(file_path)
        self.assertTrue(((x == 31.0) and (y == 32.0)) and facts == [1.0, 0.0])
