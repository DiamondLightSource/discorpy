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
# Description: Python implementation of the author's methods of
# distortion correction, Nghia T. Vo et al "Radial lens distortion
# correction with sub-pixel accuracy for X-ray micro-tomography"
# Optics Express 23, 32859-32868 (2015), https://doi.org/10.1364/OE.23.032859
# Publication date: 10th July 2018
#============================================================================

import timeit
import numpy as np
import discorpy.losa.loadersaver as io
import discorpy.post.postprocessing as post


"""
In practice, we may already have distortion coefficients of a lens, but the
center of distortion (CoD) may be changed due to lens exchange. For some
reasons, we can't collect a dot pattern again. Following is the example
which is applicable to tomographic data to search for the best CoD.
The idea is that we generate unwarped sinograms using estimated CoDs,
then we can do reconstruction and check the results to find the best CoD.
"""

time_start = timeit.default_timer()
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Initial parameters
file_path = "../data/coef_dot_05.txt"
output_base = "E:/correction/"
search_range = 100
step = 20
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Import distortion coefficients
(xcenter, ycenter, list_fact) = io.load_metadata_txt(file_path)

# Generate a 3D dataset for demonstration.
# Replace this step with a real 3D data in your codes.
mat0 = io.load_image("../data/dot_pattern_05.jpg")
(height, width) = mat0.shape
mat3D = np.zeros((600, height, width), dtype=np.float32)
mat3D[:] = mat0

# Generate unwarped slice with the index of 14
# at different xcenter and ycenter.
index = 14
for x_search in range(-search_range, search_range + step, step):
    for y_search in range(-search_range, search_range + step, step):
        corrected_slice = post.unwarp_slice_backward(
            mat3D, xcenter + x_search, ycenter + y_search, list_fact, index)
        # ----------------------------------------------------------
        # Do reconstruction here using other packages: tomopy, astra
        # ----------------------------------------------------------
        output_name = output_base + "/xcenter_"\
            + "{:5.2f}".format(xcenter + x_search) + "_ycenter_"\
            + "{:5.2f}".format(ycenter + y_search) + ".tif"
        io.save_image(output_name, corrected_slice)

time_stop = timeit.default_timer()
print("Running time is {} second!".format(time_stop - time_start))
