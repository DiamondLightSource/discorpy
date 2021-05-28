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
#============================================================================

"""
Example to show how to calculate parameters of a backward model from an
estimated forward model of the Hazard Cameras (Hazcams) on the underside of
NASA’s Perseverance Mars rover.
"""
import numpy as np
import discorpy.losa.loadersaver as io
import discorpy.post.postprocessing as post

# Load image
mat0 = io.load_image("Sol0_1st_color.png")
output_base = "figs/"
(height, width) = mat0.shape
mat0 = mat0 / np.max(mat0)

# Estimated forward model
xcenter = width / 2.0 + 110.0
ycenter = height / 2.0 - 20.0
list_pow = np.asarray([1.0, 10**(-4), 10**(-7), 10**(-10), 10**(-13)])
list_coef = np.asarray([1.0, 4.0, 5.0, 17.0, 3.0])
list_ffact = list_pow * list_coef

# Calculate parameters of backward model from the estimated forward model
list_hor_lines = []
for i in range(20, height-20, 50):
    list_tmp = []
    for j in range(20, width-20, 50):
        list_tmp.append([i - ycenter, j - xcenter])
    list_hor_lines.append(list_tmp)
Amatrix = []
Bmatrix = []
list_expo = np.arange(len(list_ffact), dtype=np.int16)
for _, line in enumerate(list_hor_lines):
    for _, point in enumerate(line):
        xd = np.float64(point[1])
        yd = np.float64(point[0])
        rd = np.sqrt(xd * xd + yd * yd)
        ffactor = np.float64(np.sum(list_ffact * np.power(rd, list_expo)))
        if ffactor != 0.0:
            Fb = 1 / ffactor
            ru = ffactor * rd
            Amatrix.append(np.power(ru, list_expo))
            Bmatrix.append(Fb)
Amatrix = np.asarray(Amatrix, dtype=np.float64)
Bmatrix = np.asarray(Bmatrix, dtype=np.float64)
list_bfact = np.linalg.lstsq(Amatrix, Bmatrix, rcond=1e-64)[0]

# Apply distortion correction
corrected_mat = post.unwarp_image_backward(mat0, xcenter, ycenter, list_bfact)
io.save_image(output_base + "/after.png", corrected_mat)
io.save_image(output_base + "/before.png", mat0)
io.save_metadata_txt(
    output_base + "/coefficients.txt", xcenter, ycenter, list_bfact)