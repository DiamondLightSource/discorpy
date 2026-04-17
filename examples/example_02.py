import timeit
import numpy as np
import discorpy.losa.loadersaver as losa
import discorpy.prep.preprocessing as prep
import discorpy.proc.processing as proc
import discorpy.post.postprocessing as post

"""
Example to show the simplest use of the package if a good quality 
dot-pattern image is used.
"""


def calc_distor_coef(mat, num_coef, perspective=False):
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
    if perspective is True:
        try:
            # Find center of distortion
            xcenter, ycenter = proc.find_center_based_vanishing_points_iteration(
                list_hor_lines, list_ver_lines, iteration=2)
            print(f"Center of distortion: X-center {xcenter}. Y-center {ycenter}")
            # Correct perspective distortion
            list_hor_lines, list_ver_lines = proc.correct_perspective_effect(
                list_hor_lines, list_ver_lines, xcenter, ycenter)
        except AttributeError:
            raise ValueError("Robust perspective correction only available from Discorpy 1.7!!!")
    if not perspective:
        (xcenter, ycenter) = proc.find_cod_coarse(list_hor_lines, list_ver_lines)
    list_fact = proc.calc_coef_backward(list_hor_lines, list_ver_lines, xcenter,
                                        ycenter, num_coef)
    return xcenter, ycenter, list_fact


time_start = timeit.default_timer()
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Initial parameters
file_path = "../data/dot_pattern_01.jpg"
output_base = "E:/correction/"
num_coef = 5  # Number of polynomial coefficients
perspective = False # Correct perspective distortion if True
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Input
print("Load image")
mat0 = losa.load_image(file_path)
# Calculation of distortion coefficients
print("Calculate distortion coefficients of a backward model...")
(xcenter, ycenter, list_fact) = calc_distor_coef(mat0, num_coef,
                                                 perspective=perspective)
# Output
corrected_mat = post.unwarp_image_backward(mat0, xcenter, ycenter, list_fact)
losa.save_image(output_base + "/corrected_image.tif", corrected_mat)
losa.save_image(output_base + "/difference.tif", corrected_mat - mat0)
losa.save_metadata_txt(output_base + "/coefficients_bw.txt", xcenter, ycenter,
                     list_fact)
time_stop = timeit.default_timer()
print("Done!!!\nRunning time is {} second !".format(time_stop - time_start))
