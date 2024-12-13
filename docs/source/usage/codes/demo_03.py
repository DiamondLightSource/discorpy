import numpy as np
import discorpy.losa.loadersaver as losa
import discorpy.prep.preprocessing as prep
import discorpy.proc.processing as proc
import discorpy.post.postprocessing as post

# Initial parameters
file_path = "../../data/dot_pattern_04.jpg"
output_base = "E:/output_demo_03/"
num_coef = 5  # Number of polynomial coefficients
mat0 = losa.load_image(file_path)  # Load image
(height, width) = mat0.shape

# Correct non-uniform background.
mat1 = prep.normalization_fft(mat0, sigma=20)
losa.save_image(output_base + "/image_normed.tif", mat1)

# Segment dots.
threshold = prep.calculate_threshold(mat1, bgr="bright", snr=3.0)
mat1 = prep.binarization(mat1, ratio=0.5, thres=threshold)
losa.save_image(output_base + "/image_binarized.tif", mat1)

# Calculate the median dot size and distance between them.
(dot_size, dot_dist) = prep.calc_size_distance(mat1)
# Remove non-dot objects
mat1 = prep.select_dots_based_size(mat1, dot_size, ratio=0.8)
# Remove non-elliptical objects
mat1 = prep.select_dots_based_ratio(mat1, ratio=0.8)
losa.save_image(output_base + "/image_cleaned.tif", mat1)

# Calculate the slopes of horizontal lines and vertical lines.
hor_slope = prep.calc_hor_slope(mat1)
ver_slope = prep.calc_ver_slope(mat1)
print("Horizontal slope: {0}. Vertical slope: {1}".format(hor_slope, ver_slope))

# Group points into lines
list_hor_lines = prep.group_dots_hor_lines(mat1, hor_slope, dot_dist, ratio=0.3,
                                           num_dot_miss=10, accepted_ratio=0.65)
list_ver_lines = prep.group_dots_ver_lines(mat1, ver_slope, dot_dist, ratio=0.3,
                                           num_dot_miss=10, accepted_ratio=0.65)
# Remove outliers
list_hor_lines = prep.remove_residual_dots_hor(list_hor_lines, hor_slope,
                                               residual=2.0)
list_ver_lines = prep.remove_residual_dots_ver(list_ver_lines, ver_slope,
                                               residual=2.0)

# Save output for checking
losa.save_plot_image(output_base + "/horizontal_lines.png", list_hor_lines,
                   height, width)
losa.save_plot_image(output_base + "/vertical_lines.png", list_ver_lines,
                   height, width)
list_hor_data = post.calc_residual_hor(list_hor_lines, 0.0, 0.0)
list_ver_data = post.calc_residual_ver(list_ver_lines, 0.0, 0.0)
losa.save_residual_plot(output_base + "/hor_residual_before_correction.png",
                      list_hor_data, height, width)
losa.save_residual_plot(output_base + "/ver_residual_before_correction.png",
                      list_ver_data, height, width)

# Regenerate grid points after correcting the perspective effect.
list_hor_lines, list_ver_lines = proc.regenerate_grid_points_parabola(
    list_hor_lines, list_ver_lines, perspective=True)

# Calculate parameters of the radial correction model
(xcenter, ycenter) = proc.find_cod_coarse(list_hor_lines, list_ver_lines)
list_fact = proc.calc_coef_backward(list_hor_lines, list_ver_lines,
                                    xcenter, ycenter, num_coef)
losa.save_metadata_txt(output_base + "/coefficients_radial_distortion.txt",
                     xcenter, ycenter, list_fact)
print("X-center: {0}. Y-center: {1}".format(xcenter, ycenter))
print("Coefficients: {0}".format(list_fact))

# Check the correction results:
# Apply correction to the lines of points
list_uhor_lines = post.unwarp_line_backward(list_hor_lines, xcenter, ycenter,
                                            list_fact)
list_uver_lines = post.unwarp_line_backward(list_ver_lines, xcenter, ycenter,
                                            list_fact)
# Calculate the residual of the unwarpped points.
list_hor_data = post.calc_residual_hor(list_uhor_lines, xcenter, ycenter)
list_ver_data = post.calc_residual_ver(list_uver_lines, xcenter, ycenter)
# Save the results for checking
losa.save_plot_image(output_base + "/unwarpped_horizontal_lines.png",
                   list_uhor_lines, height, width)
losa.save_plot_image(output_base + "/unwarpped_vertical_lines.png",
                   list_uver_lines, height, width)
losa.save_residual_plot(output_base + "/hor_residual_after_correction.png",
                      list_hor_data, height, width)
losa.save_residual_plot(output_base + "/ver_residual_after_correction.png",
                      list_ver_data, height, width)

# Correct the image
corrected_mat = post.unwarp_image_backward(mat0, xcenter, ycenter, list_fact)
# Save results. Note that the output is 32-bit-tif.
losa.save_image(output_base + "/corrected_image.tif", corrected_mat)
losa.save_image(output_base + "/difference.tif", corrected_mat - mat0)
