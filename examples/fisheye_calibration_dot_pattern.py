#==============================================================================
# Demonstration on how to calibrate a fisheye camera using a dot-pattern image
# =============================================================================

import numpy as np
import discorpy.losa.loadersaver as losa
import discorpy.prep.preprocessing as prep
import discorpy.proc.processing as proc
import discorpy.util.utility as util
import matplotlib.pyplot as plt


file_path = r"..\data\fisheye\GoPro8_dot_pattern.jpg"
test_file_path = r"..\data\fisheye\GoPro8_testing_image.jpg"
output_base = r"F:\tmp\fisheye_correction_dot_pattern"

print("  1-> Load image ...")
img0 = losa.load_image(file_path)
(height, width) = img0.shape
img_norm = prep.normalization_fft(img0, 10)

num_factor = 5
mat1 = prep.binarization(img_norm, ratio=0.3)
(dot_size, dist_dot) = prep.calc_size_distance(mat1, ratio=0.3)

print("  2-> Calculate slope and distance between lines...")
slope_hor = prep.calc_hor_slope(mat1, ratio=0.3)
slope_ver = prep.calc_ver_slope(mat1, ratio=0.3)
dist_hor = dist_ver = dist_dot
print(f"       Horizontal slope: {slope_hor} Distance: {dist_hor}")
print(f"       Vertical slope  : {slope_ver} Distance: {dist_ver}")
print("  3-> Extract reference-points !!!!")

list_points = prep.get_points_dot_pattern(mat1, binarize=False)
list_points_hor_lines = list_points
list_points_ver_lines = np.copy(list_points)

hor_margin = (450, 100)
ver_margin = (100, 100)
list_points_hor_lines = prep.remove_points_using_parabola_mask(list_points_hor_lines,
                                                               height, width,
                                                               hor_curviness=0.4,
                                                               ver_curviness=0.3,
                                                               hor_margin=hor_margin,
                                                               ver_margin=ver_margin)

list_points_ver_lines = prep.remove_points_using_parabola_mask(list_points_ver_lines,
                                                               height, width,
                                                               hor_curviness=0.4,
                                                               ver_curviness=0.3,
                                                               hor_margin=hor_margin,
                                                               ver_margin=ver_margin)

# mask = prep.make_parabola_mask(height, width,hor_curviness=0.4, ver_curviness=0.3,
#                                hor_margin=hor_margin, ver_margin=ver_margin)
# plt.imshow(img_norm * mask, origin="lower")
# plt.plot(list_points_hor_lines[:, 1], list_points_hor_lines[:, 0], ".", color="red")
# plt.plot(list_points_ver_lines[:, 1], list_points_ver_lines[:, 0], ".", color="blue")
# plt.show()

print("  4-> Group points into lines !!!!")
list_hor_lines = prep.group_dots_hor_lines_based_polyfit(list_points_hor_lines,
                                                         slope_hor, dist_hor,
                                                         ratio=0.1, num_dot_miss=3,
                                                         accepted_ratio=0.65, order=2)
list_ver_lines = prep.group_dots_ver_lines_based_polyfit(list_points_ver_lines,
                                                         slope_ver, dist_ver,
                                                         ratio=0.1, num_dot_miss=3,
                                                         accepted_ratio=0.65, order=2)

list_hor_lines = prep.remove_residual_dots_hor(list_hor_lines, slope_hor, 2.0)
list_ver_lines = prep.remove_residual_dots_ver(list_ver_lines, slope_ver, 2.0)

plt.imshow(img_norm, origin="lower")
for line in list_hor_lines:
    plt.plot(line[:, 1], line[:, 0], "-o", color="red")
plt.show()

plt.imshow(img_norm, origin="lower")
for line in list_ver_lines:
    plt.plot(line[:, 1], line[:, 0], "-o", color="blue")
plt.show()

# Find center of distortion
xcenter, ycenter = proc.find_center_based_vanishing_points_iteration(
    list_hor_lines, list_ver_lines, iteration=2)
print(f"Center of distortion: X-center {xcenter}. Y-center {ycenter}")

# Correct perspective distortion
corr_hor_lines, corr_ver_lines = proc.correct_perspective_effect(
    list_hor_lines, list_ver_lines, xcenter, ycenter)

# Calculate polynomial coefficients of the radial distortion
list_bfact = proc.calc_coef_backward(corr_hor_lines, corr_ver_lines, xcenter,
                                     ycenter, num_factor)
print(f"Polynomial coefficients of radial distortion: {list_bfact}")
losa.save_metadata_json(output_base + "/distortion_parameters.json", xcenter,
                        ycenter, list_bfact)

# Load calibration image as color image
img0 = losa.load_image(file_path, average=False)
img_corr = util.unwarp_color_image_backward(img0, xcenter, ycenter, list_bfact,
                                            pad=400)
###  Using OpenCV-remap backend for fast computing.
# img_corr = util.unwarp_image_backward_cv2(img0, xcenter, ycenter, list_bfact,
#                                           pad=400)
losa.save_image(output_base + "/corrected_img.jpg", img_corr)

# Load another image for correction.
img0 = losa.load_image(test_file_path, average=False)
img_corr = util.unwarp_color_image_backward(img0, xcenter, ycenter, list_bfact,
                                            pad=400)
losa.save_image(output_base + "/corrected_test_img.jpg", img_corr)
