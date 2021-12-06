import numpy as np
import discorpy.losa.loadersaver as io
import discorpy.prep.preprocessing as prep
import discorpy.prep.linepattern as lprep
import discorpy.proc.processing as proc
import discorpy.post.postprocessing as post
import matplotlib.pyplot as plt
import scipy.ndimage as ndi

# Initial parameters
file_path = "../../../data/laptop_camera/chessboard.jpg"
output_base = "./for_demo_06/"
num_coef = 5  # Number of polynomial coefficients
mat0 = io.load_image(file_path) # Load image
(height, width) = mat0.shape

# Convert the chessboard image to a line-pattern image
mat1 = lprep.convert_chessboard_to_linepattern(mat0)
io.save_image(output_base + "/line_pattern_converted.jpg", mat1)

# Calculate slope and distance between lines
slope_hor, dist_hor = lprep.calc_slope_distance_hor_lines(mat1, radius=15, sensitive=0.5)
slope_ver, dist_ver = lprep.calc_slope_distance_ver_lines(mat1, radius=15, sensitive=0.5)
print("Horizontal slope: ", slope_hor, " Distance: ", dist_hor)
print("Vertical slope: ", slope_ver, " Distance: ", dist_ver)

# Extract reference-points
list_points_hor_lines = lprep.get_cross_points_hor_lines(mat1, slope_ver, dist_ver,
                                                         ratio=0.3, norm=True, offset=450,
                                                         bgr="bright", radius=15,
                                                         sensitive=1.0, denoise=True,
                                                         subpixel=True)
list_points_ver_lines = lprep.get_cross_points_ver_lines(mat1, slope_hor, dist_hor,
                                                         ratio=0.3, norm=True, offset=150,
                                                         bgr="bright", radius=15,
                                                         sensitive=1.0, denoise=True,
                                                         subpixel=True)
io.save_plot_points(output_base + "/ref_points_horizontal.png", list_points_hor_lines,
                    height, width, color="red")
io.save_plot_points(output_base + "/ref_points_vertical.png", list_points_ver_lines,
                    height, width, color="blue")

# Group points into lines
list_hor_lines = prep.group_dots_hor_lines(list_points_hor_lines, slope_hor, dist_hor,
                                           ratio=0.1, num_dot_miss=2, accepted_ratio=0.8)
list_ver_lines = prep.group_dots_ver_lines(list_points_ver_lines, slope_ver, dist_ver,
                                           ratio=0.1, num_dot_miss=2, accepted_ratio=0.8)
# Remove residual dots
list_hor_lines = prep.remove_residual_dots_hor(list_hor_lines, slope_hor, 2.0)
list_ver_lines = prep.remove_residual_dots_ver(list_ver_lines, slope_ver, 2.0)

# Save output for checking
io.save_plot_image(output_base + "/horizontal_lines.png", list_hor_lines, height, width)
io.save_plot_image(output_base + "/vertical_lines.png", list_ver_lines, height, width)
list_hor_data = post.calc_residual_hor(list_hor_lines, 0.0, 0.0)
list_ver_data = post.calc_residual_ver(list_ver_lines, 0.0, 0.0)
io.save_residual_plot(output_base + "/hor_residual_before_correction.png",
                      list_hor_data, height, width)
io.save_residual_plot(output_base + "/ver_residual_before_correction.png",
                      list_ver_data, height, width)

# Regenerate grid points after correcting the perspective effect.
list_hor_lines, list_ver_lines = proc.regenerate_grid_points_parabola(
    list_hor_lines, list_ver_lines, perspective=True)

# Calculate parameters of the radial correction model
(xcenter, ycenter) = proc.find_cod_coarse(list_hor_lines, list_ver_lines)
list_fact = proc.calc_coef_backward(list_hor_lines, list_ver_lines,
                                    xcenter, ycenter, num_coef)
io.save_metadata_txt(output_base + "/coefficients_radial_distortion.txt",
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
io.save_plot_image(output_base + "/unwarpped_horizontal_lines.png",
                   list_uhor_lines, height, width)
io.save_plot_image(output_base + "/unwarpped_vertical_lines.png",
                   list_uver_lines, height, width)
io.save_residual_plot(output_base + "/hor_residual_after_correction.png",
                      list_hor_data, height, width)
io.save_residual_plot(output_base + "/ver_residual_after_correction.png",
                      list_ver_data, height, width)

# Correct the image
corrected_mat = post.unwarp_image_backward(mat0, xcenter, ycenter, list_fact)
# Save results.
io.save_image(output_base + "/corrected_image.jpg", corrected_mat)
io.save_image(output_base + "/difference.jpg", corrected_mat - mat0)

# ----------------------------------------
# Load original color image for correction
# ------------------------------------------

# Load coefficients from previous calculation if need to
# (xcenter, ycenter, list_fact) = io.load_metadata_txt(
#     output_base + "/coefficients_radial_distortion.txt")

img = io.load_image(file_path, average=False)
img_corrected = np.copy(img)
# Unwarped each color channel of the image
for i in range(img.shape[-1]):
    img_corrected[:, :, i] = post.unwarp_image_backward(img[:, :, i], xcenter,
                                                        ycenter, list_fact)
# Save the unwarped image.
io.save_image(output_base + "/corrected_image_color.jpg", img_corrected)

# Load a test image and correct it.
img = io.load_image("../../../data/laptop_camera/test_image.jpg", average=False)
img_corrected = np.copy(img)
for i in range(img.shape[-1]):
    img_corrected[:, :, i] = post.unwarp_image_backward(img[:, :, i], xcenter,
                                                        ycenter, list_fact)
io.save_image(output_base + "/test_image_corrected.jpg", img_corrected)
