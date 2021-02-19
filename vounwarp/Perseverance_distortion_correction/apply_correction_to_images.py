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
Example to show how to apply distortion correction to images of the Hazard Cameras
(Hazcams) on the underside of NASA’s Perseverance Mars rover.
"""
import numpy as np
import vounwarp.losa.loadersaver as io
import vounwarp.post.postprocessing as post
from PIL import Image

# Load color image
file_path = "Sol0_1st_color.png"
output_base = "figs/"
mat = np.asarray(Image.open(file_path), dtype=np.float32)
# Import distortion coefficients
(xcenter, ycenter, list_fact) = io.load_metadata_txt("figs/coefficients.txt")

for i in range(mat.shape[-1]):
    mat[:, :, i] = post.unwarp_image_backward(mat[:, :, i], xcenter, ycenter,
                                              list_fact)
io.save_image(output_base + "/Sol0_1st_color_correction.png", mat)