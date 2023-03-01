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
# E-mail: 
#============================================================================

"""
Example to show how to guess parameters of a forward model from
an unknown camera. In this case it's from  the Hazard Cameras (Hazcams) on the
underside of NASA’s Perseverance Mars rover.
https://mars.nasa.gov/system/downloadable_items/45689_PIA24430-Perseverance's_first_full-color_look_at_Mars.png
"""


import numpy as np
import discorpy.losa.loadersaver as io
import discorpy.post.postprocessing as post

# Load image
mat0 = io.load_image("Sol0_1st_color.png")
output_base = "figs/"
(height, width) = mat0.shape
mat0 = mat0 / np.max(mat0)

# Create line pattern
line_pattern = np.zeros((height, width), dtype=np.float32)
for i in range(50, height - 50, 40):
    line_pattern[i - 1:i + 2] = 1.0

# Estimate parameters by visual inspection.
# Coarse estimation
xcenter = width / 2.0 + 110.0
ycenter = height / 2.0 - 20.0
list_pow = np.asarray([1.0, 10**(-4), 10**(-7), 10**(-10), 10**(-13)])
# Fine estimation
list_coef = np.asarray([1.0, 4.0, 5.0, 17.0, 3.0])
list_ffact = list_pow * list_coef

pad = width
mat_pad = np.pad(line_pattern, pad, mode='edge')
mat_cor = post.unwarp_image_backward(mat_pad, xcenter + pad,
                                     ycenter + pad, list_ffact)
mat_cor = mat_cor[pad:pad + height, pad:pad + width]
io.save_image(output_base + "/overlay.jpg", (mat0 + 0.5*mat_cor))
