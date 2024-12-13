import numpy as np
import matplotlib.pyplot as plt
import discorpy.losa.loadersaver as losa
import discorpy.prep.preprocessing as prep
import discorpy.proc.processing as proc
import discorpy.post.postprocessing as post


# Initial parameters
file_path = "../../data/dot_pattern_06.jpg"
output_base = "E:/output_demo_05/"
num_coef = 4  # Number of polynomial coefficients
mat0 = losa.load_image(file_path) # Load image
(height, width) = mat0.shape

# Normalize background
mat1 = prep.normalization_fft(mat0, sigma=20)
# Segment dots
threshold = prep.calculate_threshold(mat1, bgr="bright", snr=1.5)
mat1 = prep.binarization(mat1, thres=threshold)
losa.save_image(output_base + "/segmented_dots.jpg", mat1)

# Calculate the median dot size and distance between them.
(dot_size, dot_dist) = prep.calc_size_distance(mat1)
# Calculate the slopes of horizontal lines and vertical lines.
hor_slope = prep.calc_hor_slope(mat1)
ver_slope = prep.calc_ver_slope(mat1)
print("Horizontal slope: {0}. Vertical slope: {1}".format(hor_slope, ver_slope))

# Group dots into lines.
list_hor_lines0 = prep.group_dots_hor_lines(mat1, hor_slope, dot_dist,
                                            ratio=0.3, num_dot_miss=2,
                                            accepted_ratio=0.6)
list_ver_lines0 = prep.group_dots_ver_lines(mat1, ver_slope, dot_dist,
                                            ratio=0.3, num_dot_miss=2,
                                            accepted_ratio=0.6)

# Save output for checking
losa.save_plot_image(output_base + "/horizontal_lines.png", list_hor_lines0, height, width)
losa.save_plot_image(output_base + "/vertical_lines.png", list_ver_lines0, height, width)
list_hor_data = post.calc_residual_hor(list_hor_lines0, 0.0, 0.0)
list_ver_data = post.calc_residual_ver(list_ver_lines0, 0.0, 0.0)
losa.save_residual_plot(output_base + "/hor_residual_before_correction.png", list_hor_data,
                      height, width)
losa.save_residual_plot(output_base + "/ver_residual_before_correction.png", list_ver_data,
                      height, width)

# Optional: for checking perspective distortion
(xcen_tmp, ycen_tmp) = proc.find_cod_bailey(list_hor_lines0, list_ver_lines0)
list_hor_coef = proc._para_fit_hor(list_hor_lines0, xcen_tmp, ycen_tmp)[0]
list_ver_coef = proc._para_fit_ver(list_ver_lines0, xcen_tmp, ycen_tmp)[0]
# Optional: plot the results
plt.figure(0)
plt.plot(list_hor_coef[:, 2], list_hor_coef[:, 0], "-o")
plt.plot(list_ver_coef[:, 2], list_ver_coef[:, 0], "-o")
plt.xlabel("c-coefficient")
plt.ylabel("a-coefficient")

plt.figure(1)
plt.plot(list_hor_coef[:, 2], -list_hor_coef[:, 1], "-o")
plt.plot(list_ver_coef[:, 2], list_ver_coef[:, 1], "-o")
plt.xlabel("c-coefficient")
plt.ylabel("b-coefficient")
plt.show()
# Optional: correct parabola coefficients
hor_coef_corr, ver_coef_corr = proc._generate_non_perspective_parabola_coef(
    list_hor_lines0, list_ver_lines0)[0:2]
# Optional: plot to check the results
plt.figure(0)
plt.plot(hor_coef_corr[:, 2], hor_coef_corr[:, 0], "-o")
plt.plot(ver_coef_corr[:, 2], ver_coef_corr[:, 0], "-o")
plt.xlabel("c-coefficient")
plt.ylabel("a-coefficient")
plt.figure(1)
plt.plot(hor_coef_corr[:, 2], -hor_coef_corr[:, 1], "-o")
plt.plot(ver_coef_corr[:, 2], ver_coef_corr[:, 1], "-o")
plt.xlabel("c-coefficient")
plt.ylabel("b-coefficient")
plt.show()

# Regenerate grid points with the correction of perspective effect.
list_hor_lines1, list_ver_lines1 = proc.regenerate_grid_points_parabola(
    list_hor_lines0, list_ver_lines0, perspective=True)

# Save output for checking
losa.save_plot_image(output_base + "/horizontal_lines_regenerated.png",
                   list_hor_lines1, height, width)
losa.save_plot_image(output_base + "/vertical_lines_regenerated.png",
                   list_ver_lines1, height, width)

# Calculate parameters of the radial correction model
(xcenter, ycenter) = proc.find_cod_coarse(list_hor_lines1, list_ver_lines1)
list_fact = proc.calc_coef_backward(list_hor_lines1, list_ver_lines1,
                                    xcenter, ycenter, num_coef)
losa.save_metadata_txt(output_base + "/coefficients_radial_distortion.txt",
                     xcenter, ycenter, list_fact)
print("X-center: {0}. Y-center: {1}".format(xcenter, ycenter))
print("Coefficients: {0}".format(list_fact))

# Regenerate the lines without perspective correction for later use.
list_hor_lines2, list_ver_lines2 = proc.regenerate_grid_points_parabola(
    list_hor_lines0, list_ver_lines0, perspective=False)

# Unwarp lines using the backward model:
list_uhor_lines = post.unwarp_line_backward(list_hor_lines2, xcenter, ycenter, list_fact)
list_uver_lines = post.unwarp_line_backward(list_ver_lines2, xcenter, ycenter, list_fact)

# Optional: unwarp lines using the forward model. Note that unwarping lines
# using the backward model is based on optimization (https://zenodo.org/record/1322720)
# which may result in strong fluctuation if lines are strongly curved. In such case,
# using the forward model is more reliable.

# list_ffact = proc.calc_coef_forward(list_hor_lines1, list_ver_lines1,
#                                     xcenter, ycenter, num_coef)
# list_uhor_lines = post.unwarp_line_forward(list_hor_lines2, xcenter, ycenter,
#                                            list_ffact)
# list_uver_lines = post.unwarp_line_forward(list_ver_lines2, xcenter, ycenter,
#                                            list_ffact)

# Check the residual of unwarped lines:
list_hor_data = post.calc_residual_hor(list_uhor_lines, xcenter, ycenter)
list_ver_data = post.calc_residual_ver(list_uver_lines, xcenter, ycenter)
losa.save_residual_plot(output_base + "/hor_residual_after_correction.png", list_hor_data, height, width)
losa.save_residual_plot(output_base + "/ver_residual_after_correction.png", list_ver_data, height, width)

# Unwarp the image
mat_rad_corr = post.unwarp_image_backward(mat0, xcenter, ycenter, list_fact)
# Save results
losa.save_image(output_base + "/image_radial_corrected.jpg", mat_rad_corr)
losa.save_image(output_base + "/radial_difference.jpg", mat_rad_corr - mat0)

# -------------------------------------
# For perspective distortion correction
# -------------------------------------

# Generate source points and target points to calculate coefficients of a perspective model
source_points, target_points = proc.generate_source_target_perspective_points(list_uhor_lines, list_uver_lines,
                                                                              equal_dist=True, scale="mean",
                                                                              optimizing=False)
# # Optional: generate non-perspective grid-lines for checking
# hor_lines_pers_corr, ver_lines_pers_corr = proc.generate_undistorted_perspective_lines(list_uhor_lines, list_uver_lines,
#                                                                               equal_dist=True, scale="mean",
#                                                                               optimizing=False)

# Calculate perspective coefficients:
pers_coef = proc.calc_perspective_coefficients(source_points, target_points, mapping="backward")
image_pers_corr = post.correct_perspective_image(mat_rad_corr, pers_coef)
# Save results
np.savetxt(output_base + "/perspective_coefficients.txt", np.transpose([pers_coef]))
losa.save_image(output_base + "/image_radial_perspective_corrected.jpg", image_pers_corr)
losa.save_image(output_base + "/perspective_difference.jpg", image_pers_corr - mat_rad_corr)
