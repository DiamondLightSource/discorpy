#============================================================================
# Copyright (c) 2018 Nghia T. Vo. All rights reserved.
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
#============================================================================
# Author: Nghia T. Vo
# E-mail: nghia.vo@diamond.ac.uk
# Description: Python implementation (2.7) of the author's methods of
# distortion correction, Nghia T. Vo et al "Radial lens distortion
# correction with sub-pixel accuracy for X-ray micro-tomography"
# Optics Express 23, 32859-32868 (2015), https://doi.org/10.1364/OE.23.032859
# Publication date: 10th July 2018
#============================================================================

import timeit
import numpy as np
import vounwarp.losa.loadersaver as io
import vounwarp.prep.preprocessing as prep
import vounwarp.proc.processing as proc
import vounwarp.post.postprocessing as post


"""
Example to show how to apply distortion correction along axis-1 instead of
axis-0 of a 3D dataset. This is useful for tomographic data where one can
generate directly corrected sinograms instead of doing two separate steps
of correcting distortion in the projection space and generating sinograms.

"""


def calc_distor_coef(mat, poly_order):
    # Pre-processing
    mat1 = prep.binarization(mat)
    (dot_size, dot_dist) = prep.calc_size_distance(mat1)
    mat1 = prep.select_dots_based_size(mat1, dot_size)
    mat1 = prep.select_dots_based_ratio(mat1)
    hor_slope = prep.calc_hor_slope(mat1)
    ver_slope = prep.calc_ver_slope(mat1)
    list_hor_lines = prep.group_dots_hor_lines(mat1, hor_slope, dot_dist)
    list_ver_lines = prep.group_dots_ver_lines(mat1, ver_slope, dot_dist)
    list_hor_lines = prep.remove_residual_dots_hor(list_hor_lines, hor_slope)
    list_ver_lines = prep.remove_residual_dots_ver(list_ver_lines, ver_slope)
    # Processing
    (xcenter, ycenter) = proc.find_cod_coarse(list_hor_lines, list_ver_lines)
    list_fact = proc.calc_coef_backward(
        list_hor_lines, list_ver_lines, xcenter, ycenter, poly_order)
    return xcenter, ycenter, list_fact

time_start = timeit.default_timer()
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Initial parameters
file_path = "../../vounwarp/data/dot_pattern_05.jpg"
output_base = "C:/home/unwarp/"
poly_order = 5  # Order of a polynomial
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Input
mat0 = io.load_image(file_path)
(height, width) = mat0.shape
# Distortion coefficients calculation
(xcenter, ycenter, list_fact) = calc_distor_coef(mat0, poly_order)
io.save_metadata_txt(
    output_base + "/coefficients_bw.txt", xcenter, ycenter, list_fact)
# Generate a 3D dataset for demonstration.
# Replace this step with a real 3D data in your codes.
mat3D = np.zeros((600, height, width), dtype=np.float32)
mat3D[:] = mat0
# Generate a chunk of corrected slices in the range of
# (height//2 - 10; height//2 + 10)
start_index = 14
stop_index = 20
corrected_slices = post.unwarp_chunk_slices_backward(
    mat3D, xcenter, ycenter, list_fact, start_index, stop_index)
for i in range(start_index, stop_index):
    name = "0000" + str(i)
    output_name = output_base + "/before/slice_" + name[-5:] + ".tif"
    io.save_image(output_name, mat3D[:, i, :])

    output_name = output_base + "/after/corrected_slice_" + name[-5:] + ".tif"
    io.save_image(output_name, corrected_slices[:, i - start_index, :])

time_stop = timeit.default_timer()
print("Calculation completes in {} second !").format(time_stop - time_start)
