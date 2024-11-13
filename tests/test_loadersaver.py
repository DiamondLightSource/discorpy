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
        if not os.path.isdir("./tmp/data"):
            os.makedirs("./tmp/data")

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir("./tmp"):
            shutil.rmtree("./tmp")

    def test_load_image(self):
        file_path = "./tmp/data/img.tif"
        losa.save_image(file_path, np.float32(np.random.rand(64, 64)))
        mat = losa.load_image(file_path)
        self.assertTrue(mat.shape == (64, 64))

        file_path = "./tmp/data/img.png"
        losa.save_image(file_path, np.ones((64, 64, 3), dtype=np.uint8))
        mat = losa.load_image(file_path, average=True)
        self.assertTrue(mat.shape == (64, 64))

        # self.assertRaises(ValueError, losa.load_image, "./tmp/data\\img.png")
        self.assertRaises(Exception, losa.load_image, "./tmp/data/img1.png")

    def test_get_hdf_information(self):
        file_path = "./tmp/data/data.hdf"
        ifile = h5py.File(file_path, "w")
        ifile.create_dataset("entry/data",
                             data=np.float32(np.random.rand(64, 64)))
        ifile.create_dataset("entry/energy", data=25.0)
        ifile.create_group("entry/position/stage")
        ifile.close()
        results = losa.get_hdf_information(file_path)
        self.assertTrue(len(results) == 3 and isinstance(results[0][0], str))

    def test_find_hdf_key(self):
        file_path = "./tmp/data/data.hdf"
        ifile = h5py.File(file_path, "w")
        key = "entry/energy"
        ifile.create_dataset(key, data=25.0)
        key2 = "entry/data"
        ifile.create_dataset(key2, data=np.random.rand(10))
        ifile.create_group("entry/stages/x_stage")
        ifile.create_group("entry/stages/y_stage")
        ifile.close()

        results = losa.find_hdf_key(file_path, key)
        self.assertTrue(str(results[0][0]) == key)

        results = losa.find_hdf_key(file_path, "data")
        self.assertTrue(results[0][0] == key2)

        results = losa.find_hdf_key(file_path, "x_stage")
        self.assertTrue(results[0][0] == "entry/stages/x_stage")

    def test_load_hdf_file(self):
        f_alias = losa.load_hdf_file
        file_path = "./tmp/data/data.hdf"
        ifile = h5py.File(file_path, "w")
        ifile.create_dataset("entry/data", data=np.random.rand(64, 64))
        ifile.close()
        data = f_alias(file_path, "entry/data")[:]
        self.assertTrue(isinstance(data, np.ndarray) and
                        data.shape == (64, 64))

        file_path = "./tmp/data/data2.hdf"
        ifile = h5py.File(file_path, "w")
        ifile.create_dataset("entry/data", data=np.random.rand(32, 64, 64))
        ifile.close()
        data = f_alias(file_path, None, index=(0, 32, 2), axis=0)
        self.assertTrue(data.shape[0] == 16)

        data = f_alias(file_path, None, index=(0, 32, 2), axis=1)
        self.assertTrue(data.shape[1] == 16)

        data = f_alias(file_path, None, index=(0, 32, 2), axis=2)
        self.assertTrue(data.shape[2] == 16)

        data = f_alias(file_path, None, index=1, axis=0)
        self.assertTrue(data.shape == (64, 64))

        data = f_alias(file_path, None, index=1, axis=1)
        self.assertTrue(data.shape == (32, 64))

        data = f_alias(file_path, None, index=1, axis=2)
        self.assertTrue(data.shape == (32, 64))

        data = f_alias(file_path, None, index=None, axis=0)
        self.assertTrue(data.shape == (32, 64, 64))

        data = f_alias(file_path, None, index=(0, 10), axis=0)
        self.assertTrue(data.shape == (10, 64, 64))

        data = f_alias(file_path, None, index=[1, 5, 7, 10], axis=0)
        self.assertTrue(data.shape == (4, 64, 64))

        data = f_alias(file_path, None, index=1, axis=1)
        self.assertTrue(data.shape == (32, 64))

        file_path = "./tmp/data/data3.hdf"
        ifile = h5py.File(file_path, "w")
        ifile.create_dataset("entry/test", data=np.random.rand(32, 64, 64))
        ifile.close()
        self.assertRaises(ValueError, f_alias, "./tmp/data/data3.hdf", None,
                          None, axis=0)
        self.assertRaises(ValueError, f_alias, "./tmp/data/data3.hdf",
                          "entry/test1", None, axis=0)

        file_path = "./tmp/data/data4.hdf"
        ifile = h5py.File(file_path, "w")
        ifile.create_dataset("entry/data", data=np.random.rand(2, 32, 64, 64))
        ifile.close()
        self.assertRaises(ValueError, f_alias, "./tmp/data/data4.hdf", None,
                          None, axis=0)

    def test_load_hdf_object(self):
        f_alias = losa.load_hdf_object
        file_path = "./tmp/data/data1.hdf"
        ifile = h5py.File(file_path, "w")
        ifile.create_dataset("entry/data", data=np.random.rand(64, 64))
        ifile.close()
        data = f_alias(file_path, "entry/data")
        self.assertTrue(isinstance(data, object) and data.shape == (64, 64))
        self.assertRaises(ValueError, f_alias, file_path, "entry/data1")

    def test_save_image(self):
        mat1 = np.float32(np.random.rand(64, 64))
        mat2 = np.float32(np.random.rand(64, 64, 3))
        file_path = "./tmp/data/tmp/img.tif"
        losa.save_image(file_path, mat1)
        self.assertTrue(os.path.isfile(file_path))

        file_path = "./tmp/data/tmp/img.tif"
        losa.save_image(file_path, mat2)
        self.assertTrue(os.path.isfile(file_path))

        file_path = "./tmp/data/tmp/img.png"
        losa.save_image(file_path, mat2)
        self.assertTrue(os.path.isfile(file_path))

        path = losa.save_image(file_path, mat1, overwrite=False)
        self.assertTrue(os.path.isfile(path))

        path1 = losa.save_image(file_path, mat1, overwrite=False)
        self.assertTrue(os.path.isfile(path1))

    def test_save_plot_image(self):
        f_alias = losa.save_plot_image
        list_lines = [np.asarray([(i, j) for j in range(5, 64, 5)]) for i
                      in range(5, 64, 5)]
        file_path = "./tmp/data/plot.png"
        f_alias(file_path, list_lines, 64, 64, dpi=100)
        self.assertTrue(os.path.isfile(file_path))
        path = f_alias(file_path, list_lines, 64, 64, dpi=100,
                       overwrite=False)
        self.assertTrue(os.path.isfile(path))

    def test_save_residual_plot(self):
        f_alias = losa.save_residual_plot
        ver = sys.version[:3]
        list_data = np.ones((64, 2), dtype=np.float32)
        file_path = "./tmp/data/plot1.png"
        list_data[:, 0] = 0.5 * np.random.rand(64)
        list_data[:, 1] = np.arange(64)
        f_alias(file_path, list_data, 64, 64, dpi=100)
        self.assertTrue(os.path.isfile(file_path))

        path = f_alias(file_path, list_data, 64, 64,
                       dpi=100, overwrite=False)
        self.assertTrue(os.path.isfile(path))

    def test_save_hdf_file(self):
        f_alias = losa.save_hdf_file
        file_path = "./tmp/data/data2.hdf"
        f_alias(file_path, np.random.rand(64, 64))
        self.assertTrue(os.path.isfile(file_path))

        path = f_alias(file_path, np.random.rand(64, 64), overwrite=False)
        self.assertTrue(path != file_path)

        file_path = "./tmp/data/data3"
        f_alias(file_path, np.random.rand(64, 64))
        self.assertTrue(os.path.isfile(file_path + ".hdf"))

    def test_open_hdf_stream(self):
        f_alias = losa.open_hdf_stream
        data_out = f_alias("./tmp/data/data.hdf", (64, 64))
        self.assertTrue(isinstance(data_out, object))

        data_out1 = f_alias("./tmp/data/data.hdf", (64, 64), overwrite=False)
        self.assertTrue(isinstance(data_out1, object))

        file_path = "./tmp/data/data2"
        data_out2 = f_alias(file_path, (64, 64))
        self.assertTrue(isinstance(data_out2, object))
        self.assertTrue(os.path.isfile(file_path + ".hdf"))

        data_out3 = f_alias("./tmp/data/data3.hdf", (64, 64),
                            options={"entry/energy": 25.0})
        self.assertTrue(isinstance(data_out3, object))

        self.assertRaises(ValueError, f_alias, "./tmp/data/data4.hdf",
                          (64, 64), options={"energy/entry/data": 25.0})

    def test_save_plot_points(self):
        f_alias = losa.save_plot_points
        list_data = np.ones((64, 2), dtype=np.float32)
        file_path = "./tmp/data/plot1.png"
        list_data[:, 0] = 0.5 * np.random.rand(64)
        list_data[:, 1] = np.arange(64)
        f_alias(file_path, list_data, 64, 64, dpi=100)
        self.assertTrue(os.path.isfile(file_path))

        path = f_alias(file_path, list_data, 64, 64,
                       dpi=100, overwrite=False)
        self.assertTrue(os.path.isfile(path))

    def test_save_metadata_txt(self):
        f_alias = losa.save_metadata_txt
        file_path = "./tmp/data/coef.txt"
        f_alias(file_path, 31, 32, (1.0, 0.0))
        self.assertTrue(os.path.isfile(file_path))

        path = f_alias(file_path, 31, 32, (1.0, 0.0), overwrite=False)
        self.assertTrue(path != file_path)

        file_path = "./tmp/data/coef1"
        f_alias(file_path, 31, 32, (1.0, 0.0))
        self.assertTrue(os.path.isfile(file_path + ".txt"))

    def test_load_metadata_txt(self):
        file_path = "./tmp/data/coef1.txt"
        losa.save_metadata_txt(file_path, 31.0, 32.0, [1.0, 0.0])
        (x, y, facts) = losa.load_metadata_txt(file_path)
        self.assertTrue(((x == 31.0) and (y == 32.0)) and facts == [1.0, 0.0])

    def test_save_metadata_json(self):
        f_alias = losa.save_metadata_json
        file_path = "./tmp/data/coef.json"
        f_alias(file_path, 31, 32, [1.0, 0.0])
        self.assertTrue(os.path.isfile(file_path))

        path = f_alias(file_path, 31, 32, [1.0, 0.0], overwrite=False)
        self.assertTrue(path != file_path)

        file_path_no_ext = "./tmp/data/coef1"
        f_alias(file_path_no_ext, 31, 32, [1.0, 0.0])
        self.assertTrue(os.path.isfile(file_path_no_ext + ".json"))

    def test_load_metadata_json(self):
        file_path = "./tmp/data/coef1.json"
        losa.save_metadata_json(file_path, 31.0, 32.0, [1.0, 0.0])
        x, y, facts = losa.load_metadata_json(file_path)
        self.assertTrue((x == 31.0) and (y == 32.0) and facts == [1.0, 0.0])
