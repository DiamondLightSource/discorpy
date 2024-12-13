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
# Description: Usage demonstration
# Publication date: 10th July 2018
# ============================================================================
# Contributors:
# ============================================================================

import os
import sys
import timeit
import argparse
import numpy as np
import warnings

warnings.filterwarnings("ignore")
from scipy import ndimage as ndi
import discorpy.losa.loadersaver as losa
import discorpy.prep.preprocessing as prep
import discorpy.proc.processing as proc
import discorpy.post.postprocessing as post

"""
Standalone script for calculating distortion coefficients from a dot pattern.
Acceptable file formats: tif, jpg, png, hdf, or nxs.
Example of use:
python unwarp.py -i home/user/data/dot_pattern_01.tif -o home/user/correction -n 5
"""

parser = argparse.ArgumentParser(
    description="Load an image of a dot pattern (tif, jpg, png, hdf, or nxs)"
                " and calculate distortion coefficients")
parser.add_argument("-i", dest="input", help="Path to the file", required=True)
parser.add_argument("-k", dest="key", help="Key path to the dataset if input is a hdf/nxs/h5 file",
                    required=False, default=None)
parser.add_argument("-o", dest="output", help="Output folder", required=True)
parser.add_argument("-f", dest="flat",
                    help="Path to a flat-field, or use 'norm' for self-normalization",
                    required=False, default="none")
parser.add_argument("-n", dest="order",
                    help="Number of polynomial coefficients", type=int,
                    required=False, default=5)
parser.add_argument("-p", dest="perspective",
                    help="Enable perspective correction", required=False,
                    default=False)
print("******************************************************************\n")
print("                    Start the script!!!\n                            ")
print("******************************************************************\n")
args = parser.parse_args()
file_path = args.input
output_base = args.output
flat_path = args.flat
poly_order = args.order
key_path_hdf = args.key
perspective = args.perspective

time_start = timeit.default_timer()
# Load data
print("1 ---> Load file: {}".format(file_path))
_, file_ext = os.path.splitext(file_path)
if (file_ext == ".hdf") or (file_ext == ".nxs"):
    mat0 = losa.load_hdf_file(file_path, key_path=key_path_hdf, index=None,
                            axis=0)
    if len(mat0.shape) == 3:
        mat0 = np.mean(mat0, axis=0)
else:
    mat0 = losa.load_image(file_path)
(height, width) = mat0.shape

# Load flat-field
if (flat_path != "none") and (flat_path != "norm"):
    _, file_ext = os.path.splitext(flat_path)
    if (file_ext == ".hdf") or (file_ext == ".nxs"):
        flat = losa.load_hdf_file(flat_path, key_path=key_path_hdf, index=None,
                                axis=0)
        if len(flat.shape) == 3:
            flat = np.mean(flat, axis=0)
    else:
        flat = losa.load_image(flat_path)
    (height1, width1) = flat.shape
    if (height != height1) or (width != width):
        raise ValueError(
            "Shape of the image and the flat-field are not the same")
    nmean = np.mean(flat)
    flat[flat == 0.0] = nmean
    mat1 = mat0 / flat
elif flat_path == "norm":
    mat1 = prep.normalization_fft(mat0, sigma=5, pad=40)
else:
    mat1 = mat0

# Binarize
print("2 ---> Binarize and clean the dot-pattern image")
mat1 = prep.binarization(mat1, ratio=0.3, thres=None)
losa.save_image(output_base + "/binarized_image.tif", mat1)
# Calculate dot size, distance of two nearest dots
(dot_size, dot_dist) = prep.calc_size_distance(mat1, ratio=0.5)
# Select dots with a certain size
mat1 = prep.select_dots_based_size(mat1, dot_size, ratio=0.5)
# Select dots with a certain ratio between major and minor axis
mat1 = prep.select_dots_based_ratio(mat1, ratio=0.3)
losa.save_image(output_base + "/cleaned_image.tif", mat1)

# Calculate the horizontal slope and the vertical slope of the grid
print("3 ---> Calculate the slopes of lines of dots\n")
hor_slope = prep.calc_hor_slope(mat1, ratio=0.3)
ver_slope = prep.calc_ver_slope(mat1, ratio=0.3)
print("	Horizontal slope: {0}  Vertical slope: {1}\n".format(hor_slope,
                                                                ver_slope))

print("4 ---> Group dots into horizontal lines and vertical lines")
# Group dots into horizontal lines and vertical lines
list_hor_lines = prep.group_dots_hor_lines(mat1, hor_slope, dot_dist, ratio=0.3,
                                           num_dot_miss=6, accepted_ratio=0.65)
list_ver_lines = prep.group_dots_ver_lines(mat1, ver_slope, dot_dist, ratio=0.3,
                                           num_dot_miss=6, accepted_ratio=0.65)

print("5 ---> Remove misplaced dots")
# Remove residual dots.
list_hor_lines = prep.remove_residual_dots_hor(list_hor_lines, hor_slope,
                                               residual=2.0)
list_ver_lines = prep.remove_residual_dots_ver(list_ver_lines, ver_slope,
                                               residual=2.0)
losa.save_plot_image(output_base + "/group_horizontal_dots.png", list_hor_lines,
                   height, width)
losa.save_plot_image(output_base + "/group_vertical_dots.png", list_ver_lines,
                   height, width)

print("6 ---> Calculate residuals before distortion correction")
# Check if the distortion is significant.
list_hor_data = post.calc_residual_hor(list_hor_lines, 0.0, 0.0)
losa.save_residual_plot(output_base + "/residual_hor_before_correction.png",
                      list_hor_data, height, width)
list_ver_data = post.calc_residual_ver(list_ver_lines, 0.0, 0.0)
losa.save_residual_plot(output_base + "/residual_ver_before_correction.png",
                      list_ver_data, height, width)

# check1 = post.check_distortion(list_hor_data)
# check2 = post.check_distortion(list_ver_data)
# if (not check1) and (not check2):
#     print("!!! Distortion is not significant !!!")
#     sys.exit(0)

if perspective is True:
    try:
        list_hor_lines, list_ver_lines = proc.regenerate_grid_points_parabola(
            list_hor_lines, list_ver_lines, perspective=perspective)
    except AttributeError:
        raise ValueError("Perspective correction only available from "
                         "Discorpy 1.4!!!")

# Calculate center of distortion. xcenter is the center from the left
# of the image. ycenter is the center from the top of the image.
print("7 ---> Calculate the center of distortion\n")
(xcenter, ycenter) = proc.find_cod_coarse(list_hor_lines, list_ver_lines)
# Use fine-search if there's no perspective distortion
# (xcenter, ycenter) = proc.find_cod_fine(
#     list_hor_lines, list_ver_lines, xcenter, ycenter, dot_dist)
print("\nCenter of distortion:\nx-center (from the left of the image): "
      "{0}\ny-center (from the top of the image): {1}\n".format(
        xcenter, ycenter))

print("8 ---> Calculate distortion coefficients")
# Calculate distortion coefficients using the backward-from-forward model
list_ffact, list_bfact = proc.calc_coef_backward_from_forward(
    list_hor_lines, list_ver_lines, xcenter, ycenter, poly_order)

# # Calculate distortion coefficients using the backward model
# list_bfact = proc.calc_coef_backward(
#     list_hor_lines, list_ver_lines, xcenter, ycenter, poly_order)

# Apply distortion correction
print("9 ---> Apply distortion correction")
corrected_mat = post.unwarp_image_backward(mat0, xcenter, ycenter, list_bfact)
losa.save_image(output_base + "/corrected_image.tif", corrected_mat)
losa.save_image(output_base + "/original_image.tif", mat0)
losa.save_metadata_txt(output_base + "/distortion_coefficients_bw.txt", xcenter,
                     ycenter, list_bfact)
# Check the correction results
print("10---> Evaluate the correction results")
list_uhor_lines = post.unwarp_line_backward(list_hor_lines, xcenter, ycenter,
                                            list_bfact)
list_uver_lines = post.unwarp_line_backward(list_ver_lines, xcenter, ycenter,
                                            list_bfact)
losa.save_plot_image(output_base + "/horizontal_dots_unwarped.png",
                   list_uhor_lines, height, width)
losa.save_plot_image(output_base + "/vertical_dots_unwarped.png",
                   list_uver_lines, height, width)
list_hor_data = post.calc_residual_hor(list_uhor_lines, xcenter, ycenter)
list_ver_data = post.calc_residual_ver(list_uver_lines, xcenter, ycenter)
losa.save_residual_plot(output_base + "/residual_hor_after_correction.png",
                      list_hor_data, height, width)
losa.save_residual_plot(output_base + "/residual_ver_after_correction.png",
                      list_ver_data, height, width)
check1 = post.check_distortion(list_hor_data)
check2 = post.check_distortion(list_ver_data)
if check1 or check2:
    print("!!! Correction results are not at sub-pixel accuracy !!!")
time_stop = timeit.default_timer()
print("******************************************************************\n")
print("   Time cost  is: {} second \n".format(time_stop - time_start))
print("******************************************************************\n")
print("   Results are at: {} \n".format(output_base))
print("******************************************************************\n")
