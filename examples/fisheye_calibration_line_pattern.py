#==============================================================================
# Demonstration on how to calibrate a fisheye camera using a line-pattern image
# =============================================================================

import discorpy.losa.loadersaver as losa
import discorpy.prep.preprocessing as prep
import discorpy.prep.linepattern as lprep
import discorpy.proc.processing as proc
import discorpy.post.postprocessing as post
import discorpy.util.utility as util

file_path = "../data/fisheye/GoPro8_line_pattern.jpg"
output_base = "E:/fisheye_correction"
num_factor = 5

print("1-> Load image ...")
img0 = losa.load_image(file_path)
(height, width) = img0.shape
img_norm = prep.normalization_fft(img0, 10)

print("2-> Calculate slope and distance between lines...")
slope_hor, dist_hor = lprep.calc_slope_distance_hor_lines(img_norm, chessboard=False)
slope_ver, dist_ver = lprep.calc_slope_distance_ver_lines(img_norm, chessboard=False)
print(f"       Horizontal slope: {slope_hor} Distance: {dist_hor}")
print(f"       Vertical slope  : {slope_ver} Distance: {dist_ver}")
print("3-> Extract reference-points !!!!")

# Detect points on lines, lines are dark, background is bright.
list_points_hor_lines = lprep.get_cross_points_hor_lines(img_norm, slope_ver,
                                                         dist_ver, bgr='bright',
                                                         chessboard=False, radius=9,
                                                         sensitive=0.1)
list_points_ver_lines = lprep.get_cross_points_ver_lines(img_norm, slope_hor,
                                                         dist_hor, bgr='bright',
                                                         chessboard=False, radius=9,
                                                         sensitive=0.1)
# Optional: Remove unwanted points at image border
list_points_hor_lines = prep.remove_points_using_parabola_mask(
    list_points_hor_lines, height, width, hor_curviness=0.4, ver_curviness=0.3,
    hor_margin=(400, 300), ver_margin=(150, 200))

list_points_ver_lines = prep.remove_points_using_parabola_mask(
    list_points_ver_lines, height, width, hor_curviness=0.4, ver_curviness=0.3,
    hor_margin=(400, 300), ver_margin=(150, 200))

print("4-> Group points into lines !!!!")
list_hor_lines = prep.group_dots_hor_lines_based_polyfit(list_points_hor_lines,
                                                         slope_hor, dist_hor,
                                                         ratio=0.1, num_dot_miss=3,
                                                         accepted_ratio=0.65, order=2)
list_ver_lines = prep.group_dots_ver_lines_based_polyfit(list_points_ver_lines,
                                                         slope_ver, dist_ver,
                                                         ratio=0.1, num_dot_miss=3,
                                                         accepted_ratio=0.65, order=2)
list_hor_lines = prep.remove_residual_dots_hor(list_hor_lines, slope_hor, 3.0)
list_ver_lines = prep.remove_residual_dots_ver(list_ver_lines, slope_ver, 3.0)

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
# img_corr = util.unwarp_image_backward_cv2(img0, xcenter, ycenter,
#                                                list_bfact, pad=400)
losa.save_image(output_base + "/corrected_img.jpg", img_corr)
