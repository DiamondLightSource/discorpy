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
# Description: Python implementation of the author's methods of
# distortion correction, Nghia T. Vo et al "Radial lens distortion
# correction with sub-pixel accuracy for X-ray micro-tomography"
# Optics Express 23, 32859-32868 (2015), https://doi.org/10.1364/OE.23.032859
# Publication date: 10th July 2018
#============================================================================

import timeit
import vounwarp.losa.loadersaver as io
import vounwarp.prep.preprocessing as prep
import vounwarp.proc.processing as proc
import vounwarp.post.postprocessing as post


"""
Example to show the simplest use of the package.

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
file_path = "../data/dot_pattern_01.jpg"
output_base = "C:/correction/"
poly_order = 5  # Number of polynomial coefficients
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Input
mat0 = io.load_image(file_path)
# Distortion coefficients calculation
(xcenter, ycenter, list_fact) = calc_distor_coef(mat0, poly_order)
# Output
corrected_mat = post.unwarp_image_backward(
    mat0, xcenter, ycenter, list_fact)
io.save_image(output_base + "/corrected_image_bw.tif", corrected_mat)
io.save_metadata_txt(
    output_base + "/coefficients_bw.txt", xcenter, ycenter, list_fact)
time_stop = timeit.default_timer()
print(("Calculation completes in {} second !").format(time_stop - time_start))
