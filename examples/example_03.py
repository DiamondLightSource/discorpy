#============================================================================
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
#============================================================================
# Author: Nghia T. Vo
# E-mail: nghia.vo@diamond.ac.uk
# Description: Python implementation of the author's methods of
# distortion correction, Nghia T. Vo et al "Radial lens distortion
# correction with sub-pixel accuracy for X-ray micro-tomography"
# Optics Express 23, 32859-32868 (2015), https://doi.org/10.1364/OE.23.032859
# Publication date: 10th July 2018
#============================================================================

import numpy as np
import timeit
import vounwarp.losa.loadersaver as io
import vounwarp.prep.preprocessing as prep
import vounwarp.proc.processing as proc
import vounwarp.post.postprocessing as post

"""
Example to show how to adjust parameters of methods to process a challenging 
image. 'dot_pattern_04.jpg' in the 'data/' folder is a low quality image with
lots of blobs, misplaced dots, and missing dots, so some parameters of 
pre-procesing methods need to be adjusted.
"""

time_start = timeit.default_timer()
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
file_path = "../data/dot_pattern_04.jpg"
output_base = "E:/correction/"
num_coef = 5  # Number of polynomial coefficients
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Load an image
print("Load image: {}".format(file_path))
mat0 = io.load_image(file_path)
(height, width) = mat0.shape

# Correct non-uniform background. Use sigma = 5
mat1 = prep.normalization_fft(mat0, sigma=5)

# Binarize the image.
# For tiny dots, you may need to set the threshold manually and add a parameter
# denoise=False.
mat1 = prep.binarization(mat1, ratio=0.5, thres=None)
io.save_image(output_base + "/binarized_image.tif", mat1)

# Calculate the median dot-size and the median distance of two nearest dots.
(dot_size, dot_dist) = prep.calc_size_distance(mat1, ratio=0.3)
print("Median size of dots: {0}\nMedian distance between two dots: {1}".format(
    dot_size, dot_dist))

# Select dots with a certain range of size.
# Ratio is changed to 0.9 as the size of binary dots varies significantly.
mat1 = prep.select_dots_based_size(mat1, dot_size, ratio=0.9)
io.save_image(output_base + "/cleaned_1_image.tif", mat1)

# Select dots with a certain ratio between the major and the minor axis
# Strong distortion at the edges of the image makes round-shape dots become
# elliptical dots, so the acceptable variation is changed to be larger,
# ratio=1.0.
mat1 = prep.select_dots_based_ratio(mat1, ratio=1.0)
io.save_image(output_base + "/cleaned_2_image.tif", mat1)

# Calculate the horizontal slope and the vertical slope of the grid.
hor_slope = prep.calc_hor_slope(mat1, ratio=0.3)
ver_slope = prep.calc_ver_slope(mat1, ratio=0.3)
print("Horizontal slope: {0}\nVertical slope: {1}".format(hor_slope, ver_slope))

# Group dots into horizontal lines and vertical lines.
list_hor_lines = prep.group_dots_hor_lines(mat1, hor_slope, dot_dist, ratio=0.3,
                                           num_dot_miss=10, accepted_ratio=0.6)
list_ver_lines = prep.group_dots_ver_lines(mat1, ver_slope, dot_dist, ratio=0.3,
                                           num_dot_miss=10, accepted_ratio=0.6)

# Remove residual dots.
list_hor_lines = prep.remove_residual_dots_hor(list_hor_lines, hor_slope,
                                               residual=2.0)
list_ver_lines = prep.remove_residual_dots_ver(list_ver_lines, ver_slope,
                                               residual=2.0)
io.save_plot_image(output_base + "/group_horizontal_dots.png", list_hor_lines,
                   height, width)
io.save_plot_image(output_base + "/group_vertical_dots.png", list_ver_lines,
                   height, width)

# Following steps are straightforward:

# Check if the distortion is significant.
list_hor_data = post.calc_residual_hor(list_hor_lines, 0.0, 0.0)
io.save_residual_plot(output_base + "/residual_hor_before_correction.png",
                      list_hor_data, height, width)
list_ver_data = post.calc_residual_ver(list_ver_lines, 0.0, 0.0)
io.save_residual_plot(output_base + "/residual_ver_before_correction.png",
                      list_ver_data, height, width)
check1 = post.check_distortion(list_hor_data)
check2 = post.check_distortion(list_ver_data)
if (not check1) and (not check2):
    print("!!! Distortion is not significant !!!")

# Calculate the center of distortion. xcenter is the center from the left
# of the image. ycenter is the center from the top of the image.
(xcenter, ycenter) = proc.find_cod_coarse(list_hor_lines, list_ver_lines)
# Use fine-search if there's no perspective distortion.
# (xcenter, ycenter) = proc.find_cod_fine(
#     list_hor_lines, list_ver_lines, xcenter, ycenter, dot_dist)
print("\nCenter of distortion:\nx-center (from the left of the image): "
      "{0}\ny-center (from the top of the image): {1}\n".format(
        xcenter, ycenter))

# Calculate distortion coefficients of a backward-from-forward model.
print("Calculate distortion coefficients...")
list_ffact, list_bfact = proc.calc_coef_backward_from_forward(list_hor_lines,
                                                              list_ver_lines,
                                                              xcenter, ycenter,
                                                              num_coef)

# Apply distortion correction
corrected_mat = post.unwarp_image_backward(mat0, xcenter, ycenter, list_bfact)
io.save_image(output_base + "/corrected_image.tif", corrected_mat)
io.save_image(output_base + "/diff_corrected_image.tif",
              np.abs(corrected_mat - mat0))
io.save_metadata_txt(output_base + "/distortion_coefficients_bw.txt", xcenter,
                     ycenter, list_bfact)

# Check the correction results.
list_uhor_lines = post.unwarp_line_backward(list_hor_lines, xcenter, ycenter,
                                            list_bfact)
list_uver_lines = post.unwarp_line_backward(list_ver_lines, xcenter, ycenter,
                                            list_bfact)
io.save_plot_image(output_base + "/horizontal_dots_unwarped.png",
                   list_uhor_lines, height, width)
io.save_plot_image(output_base + "/vertical_dots_unwarped.png",
                   list_uver_lines, height, width)
list_hor_data = post.calc_residual_hor(list_uhor_lines, xcenter, ycenter)
list_ver_data = post.calc_residual_ver(list_uver_lines, xcenter, ycenter)
io.save_residual_plot(output_base + "/residual_hor_after_correction.png",
                      list_hor_data, height, width)
io.save_residual_plot(output_base + "/residual_ver_after_correction.png",
                      list_ver_data, height, width)
check1 = post.check_distortion(list_hor_data)
check2 = post.check_distortion(list_ver_data)
if check1 or check2:
    print("!!! Correction results are not at sub-pixel accuracy !!!")
time_stop = timeit.default_timer()
print("Done!!!\nRunning time is {} second!".format(time_stop - time_start))
