# ============================================================================
# ============================================================================
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
# ============================================================================
# Author: Nghia T. Vo
# E-mail: nghia.vo@diamond.ac.uk
# Description: Python implementation of the author's methods of
# distortion correction, Nghia T. Vo et al "Radial lens distortion
# correction with sub-pixel accuracy for X-ray micro-tomography"
# Optics Express 23, 32859-32868 (2015), https://doi.org/10.1364/OE.23.032859
# Publication date: 10th July 2018
# ============================================================================
# Contributors:
# ============================================================================

"""
Module for I/O tasks:
- Load data from an image file (tif, png, jpg) or a hdf file.
- Save a 2D array as a tif/png/jpg image or 2D, 3D array to a hdf file.
- Save a plot of data points as an image.
- Save and load metadata to and from a text file.
"""

import os
import h5py
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from collections import OrderedDict
import errno


def load_image(file_path, average=True):
    """
    Load data from an image.

    Parameters
    ----------
    file_path : str
        Path to a file.
    average : bool, optional
        Average a multi-channel image if True.

    Returns
    -------
    array_like
    """
    if "\\" in file_path:
        raise ValueError(
            "Please use a file path following the Unix convention")
    try:
        mat = np.asarray(Image.open(file_path), dtype=np.float32)
    except IOError:
        print("No such file or directory: {}".format(file_path))
        raise
    if len(mat.shape) > 2 and average is True:
        axis_m = np.argmin(mat.shape)
        mat = np.mean(mat, axis=axis_m)
    return mat


def _get_key(name, obj):
    """
    Find a key path having 'data' in a dataset. Use with Group.visititems()
    method to walk through a hdf5 tree.
    """
    wanted_key = 'data'
    if isinstance(obj, h5py.Group):
        for key, val in list(obj.items()):
            if key == wanted_key:
                if isinstance(obj[key], h5py.Dataset):
                    key_path = obj.name + "/" + key
                    return key_path


def get_hdf_information(file_path):
    """
    Get information of datasets in a hdf/nxs file.

    Parameters
    ----------
    file_path : str
        Path to the file.

    Returns
    -------
    list_key : str
        Keys to the datasets.
    list_shape : tuple of int
        Shapes of the datasets.
    list_type : str
        Types of the datasets.
    """
    ifile = h5py.File(file_path, 'r')
    keys = []
    ifile.visit(keys.append)
    list_key = []
    list_shape = []
    list_type = []
    for key in keys:
        data = ifile[key]
        if isinstance(data, h5py.Group):
            for key2, _ in list(data.items()):
                list_key.append(key + "/" + key2)
        else:
            list_key.append(data.name)
    for i, key in enumerate(list_key):
        data = ifile[list_key[i]]
        try:
            shape = data.shape
        except AttributeError:
            shape = "None"
        try:
            dtype = data.dtype
        except AttributeError:
            dtype = "None"
        if isinstance(data, list):
            if len(data) == 1:
                if not isinstance(data, np.ndarray):
                    dtype = str(list(data)[0])
                    dtype = dtype.replace("b'", "'")
        list_shape.append(shape)
        list_type.append(dtype)
    ifile.close()
    return list_key, list_shape, list_type


def find_hdf_key(file_path, pattern):
    """
    Find datasets matching the pattern in a hdf/nxs file.

    Parameters
    ----------
    file_path : str
        Path to the file.
    pattern : str
        Pattern to find the full names of the datasets.

    Returns
    -------
    list_key : str
        Keys to the datasets.
    list_shape : tuple of int
        Shapes of the datasets.
    list_type : str
        Types of the datasets.
    """
    ifile = h5py.File(file_path, 'r')
    list_key = []
    keys = []
    ifile.visit(keys.append)
    for key in keys:
        data = ifile[key]
        if isinstance(data, h5py.Group):
            for key2, _ in list(data.items()):
                list_key.append(data.name + "/" + key2)
        else:
            list_key.append(data.name)
    list_dkey = []
    list_dshape = []
    list_dtype = []
    for _, key in enumerate(list_key):
        if pattern in key:
            list_dkey.append(key)
            data = ifile[key]
            try:
                shape = data.shape
            except AttributeError:
                shape = "None"
            list_dshape.append(shape)
            try:
                dtype = data.dtype
            except AttributeError:
                dtype = "None"
            list_dtype.append(dtype)
            if isinstance(data, list):
                if len(data) == 1:
                    dtype = str(list(data)[0])
                    dtype = dtype.replace("b'", "'")
    ifile.close()
    return list_dkey, list_dshape, list_dtype


def load_hdf_file(file_path, key_path=None, index=None, axis=0):
    """
    Load data from a hdf5/nxs file.

    Parameters
    ----------
    file_path : str
        Path to a hdf/nxs file.
    key_path : str
        Key path to a dataset.
    index : int or tuple of int
        Values for slicing data. Can be integer, tuple or list,
        e.g index=(start,stop,step) or index=(slice1, slice2, slice3,slice4).
    axis : int
        Slice direction

    Returns
    -------
    array_like
        2D array or 3D array.
    """
    mat = None
    if "\\" in file_path:
        raise ValueError(
            "Please use a file path following the Unix convention")
    try:
        ifile = h5py.File(file_path, 'r')
    except IOError:
        print("Couldn't open file: {}".format(file_path))
        raise
    if key_path is None:
        key_path = ifile.visititems(_get_key)  # Find the key automatically
        if key_path is None:
            raise ValueError("Please provide the key path to the dataset!")
    check = key_path in ifile
    if not check:
        print("Couldn't open object with the key path: {}".format(key_path))
        raise ValueError("!!! Wrong key !!!")
    idata = ifile[key_path]
    shape = idata.shape
    if len(shape) < 2 or len(shape) > 3:
        raise ValueError("Require a 2D or 3D dataset!")
    if len(shape) == 2:
        mat = np.asarray(idata)
    if len(shape) == 3:
        axis = np.clip(axis, 0, 3)
        if index is None:
            mat = np.float32(idata[:, :, :])
        else:
            if isinstance(index, int):
                try:
                    if axis == 0:
                        mat = np.float32(idata[index, :, :])
                    elif axis == 1:
                        mat = np.float32(idata[:, index, :])
                    else:
                        mat = np.float32(idata[:, :, index])
                except ValueError:
                    print("Index out of range!")
            if isinstance(index, tuple) or isinstance(index, list):
                if len(index) == 3:
                    starti = index[0]
                    stopi = index[1]
                    stepi = index[2]
                    list_index = list(range(starti, stopi, stepi))
                elif len(index) == 2:
                    starti = index[0]
                    stopi = index[1]
                    list_index = list(range(starti, stopi))
                else:
                    list_index = list(index)
                try:
                    if axis == 0:
                        mat = np.float32(idata[list_index, :, :])
                    elif axis == 1:
                        mat = np.float32(idata[:, list_index, :])
                    else:
                        mat = np.float32(idata[:, :, list_index])
                except ValueError:
                    print("Index out of range!")
            if mat.shape[axis] == 1:
                mat = np.swapaxes(mat, axis, 0)[0]
            if mat.shape[axis] == 0:
                raise ValueError("Empty indices!")
    return mat


def load_hdf_object(file_path, key_path):
    """
    Load a hdf/nexus dataset as an object.

    Parameters
    ----------
    file_path : str
        Path to a hdf/nxs file.
    key_path : str
        Key path to a dataset.

    Returns
    -------
    object
        hdf/nxs object.
    """
    try:
        ifile = h5py.File(file_path, 'r')
    except IOError:
        print("Couldn't open file: {}".format(file_path))
        raise
    check = key_path in ifile
    if not check:
        print("Couldn't open object with the key path: {}".format(key_path))
        raise ValueError("!!! Wrong key !!!")
    return ifile[key_path]


def _create_folder(file_path):
    """
    Create a folder to save a file if not exists.

    Parameters
    ----------
    file_path : str
        Path to a file
    """
    file_base = os.path.dirname(file_path)
    if not os.path.exists(file_base):
        try:
            os.makedirs(file_base)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


def _create_file_name(file_path):
    """
    Create a file name to avoid overwriting.

    Parameters
    ----------
    file_path : str
        Path to a file

    Returns
    -------
    str
        Updated file path.
    """
    file_base, file_ext = os.path.splitext(file_path)
    if os.path.isfile(file_path):
        nfile = 1
        check = True
        while check:
            name_add = '0000' + str(nfile)
            file_path = file_base + "_" + name_add[-4:] + file_ext
            if os.path.isfile(file_path):
                nfile = nfile + 1
            else:
                check = False
    return file_path


def save_image(file_path, mat, overwrite=True):
    """
    Save 2D data to an image.

    Parameters
    ----------
    file_path : str
        Output file path.
    mat : array_like
        2D array.
    overwrite : bool, optional
        Overwrite an existing file if True.

    Returns
    -------
    str
        Updated file path.
    """
    if "\\" in file_path:
        raise ValueError(
            "Please use a file path following the Unix convention")
    file_base, file_ext = os.path.splitext(file_path)
    if not ((file_ext == ".tif") or (file_ext == ".tiff")):
        mat = np.uint8(255 * (mat - np.min(mat)) / (np.max(mat) - np.min(mat)))
    _create_folder(file_path)
    if not overwrite:
        file_path = _create_file_name(file_path)
    image = Image.fromarray(mat)
    try:
        image.save(file_path)
    except IOError:
        print("Couldn't write to file {}".format(file_path))
        raise
    return file_path


def save_plot_image(file_path, list_lines, height, width, overwrite=True,
                    dpi=100):
    """
    Save the plot of dot-centroids to an image. Useful to check if the dots
    are arranged properly where dots on the same line having the same color.

    Parameters
    ----------
    file_path : str
        Output file path.
    list_lines : list of array_like
        List of 2D arrays. Each list is the coordinates of dots on a line.
    height : int
        Height of the image.
    width : int
        Width of the image.
    overwrite : bool, optional
        Overwrite the existing file if True.
    dpi : int, optional
        The resolution in dots per inch.

    Returns
    -------
    str
        Updated file path.
    """
    if "\\" in file_path:
        raise ValueError(
            "Please use a file path following the Unix convention")
    _create_folder(file_path)
    if not overwrite:
        file_path = _create_file_name(file_path)
    fig = plt.figure()
    fig.set_size_inches(width / dpi, height / dpi)
    ax = plt.Axes(fig, [0., 0., 1.0, 1.0])
    ax.set_axis_off()
    fig.add_axes(ax)
    plt.axis((0, width, 0, height))
    m_size = 0.5 * min(height / dpi, width / dpi)
    for line in list_lines:
        plt.plot(line[:, 1], height - line[:, 0], '-o', markersize=m_size)
    try:
        plt.savefig(file_path, dpi=dpi)
    except IOError:
        print("Couldn't write to file {}".format(file_path))
        raise
    plt.close()
    return file_path


def save_residual_plot(file_path, list_data, height, width, overwrite=True,
                       dpi=100):
    """
    Save the plot of residual against radius to an image. Useful to check the
    accuracy of unwarping results.

    Parameters
    ----------
    file_path : str
        Output file path.
    list_data : array_like
        2D array. List of [residual, radius] of each dot.
    height : int
        Height of the output image.
    width : int
        Width of the output image.
    overwrite : bool, optional
        Overwrite the existing file if True.
    dpi : int, optional
        The resolution in dots per inch.

    Returns
    -------
    str
        Updated file path.
    """
    if "\\" in file_path:
        raise ValueError(
            "Please use a file path following the Unix convention")
    _create_folder(file_path)
    if not overwrite:
        file_path = _create_file_name(file_path)
    fig = plt.figure()
    fig.set_size_inches(width / dpi, height / dpi)
    m_size = 0.5 * min(height / dpi, width / dpi)
    plt.rc('font', size=np.int16(m_size * 3))
    plt.xlabel('Radius')
    plt.ylabel('Residual')
    plt.plot(list_data[:, 0], list_data[:, 1], '.', markersize=m_size)
    try:
        plt.savefig(file_path, dpi=dpi, bbox_inches='tight')
    except IOError:
        print("Couldn't write to file {}".format(file_path))
        raise
    plt.close()
    return file_path


def save_hdf_file(file_path, idata, key_path='entry', overwrite=True):
    """
    Write data to a hdf5 file.

    Parameters
    ----------
    file_path : str
        Output file path.
    idata : array_like
        Data to be saved.
    key_path : str
        Key path to the dataset.
    overwrite : bool, optional
        Overwrite an existing file if True.

    Returns
    -------
    str
        Updated file path.
    """
    if "\\" in file_path:
        raise ValueError(
            "Please use a file path following the Unix convention")
    file_base, file_ext = os.path.splitext(file_path)
    if not ((file_ext == '.hdf') or (file_ext == '.h5')):
        file_ext = '.hdf'
    file_path = file_base + file_ext
    _create_folder(file_path)
    if not overwrite:
        file_path = _create_file_name(file_path)
    try:
        ofile = h5py.File(file_path, 'w')
    except IOError:
        print("Couldn't write file: {}".format(file_path))
        raise
    grp = ofile.create_group(key_path)
    grp.create_dataset("data", data=idata)
    ofile.close()
    return file_path


def open_hdf_stream(file_path, data_shape, key_path='entry/data',
                    data_type='float32', overwrite=True, **options):
    """
    Open stream to write data to a hdf/nxs file with options to add metadata.

    Parameters
    ----------
    file_path : str
        Path to the file.
    data_shape : tuple of int
        Shape of the data.
    key_path : str
        Key path to the dataset.
    data_type: str
        Type of data.
    overwrite : bool
        Overwrite the existing file if True.
    options : dict, optional
        Add metadata. E.g. options={"entry/angles": angles, "entry/energy": 53}.

    Returns
    -------
    object
        hdf object.
    """
    file_base, file_ext = os.path.splitext(file_path)
    if not (file_ext == '.hdf' or file_ext == '.h5' or file_ext == ".nxs"):
        file_ext = '.hdf'
    file_path = file_base + file_ext
    _create_folder(file_path)
    if not overwrite:
        file_path = _create_file_name(file_path)
    try:
        ofile = h5py.File(file_path, 'w')
    except IOError:
        print("Couldn't write to file: {}".format(file_path))
        raise
    if len(options) != 0:
        for opt_name in options:
            opts = options[opt_name]
            for key in opts:
                if key_path in key:
                    msg = "!!! Selected key path, '{0}', can not be a child" \
                          " key-path of '{1}' !!!\n!!! Change to make sure " \
                          "they are at the same level !!!".format(key, key_path)
                    raise ValueError(msg)
                ofile.create_dataset(key, data=opts[key])
    data_out = ofile.create_dataset(key_path, data_shape, dtype=data_type)
    return data_out


def save_metadata_txt(file_path, xcenter, ycenter, list_fact, overwrite=True):
    """
    Write metadata to a text file.

    Parameters
    ----------
    file_path : str
        Output file path.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_fact : float
        1D array. Coefficients of a polynomial.
    overwrite : bool, optional
        Overwrite an existing file if True.

    Returns
    -------
    str
        Updated file path.
    """
    file_base, file_ext = os.path.splitext(file_path)
    if not ((file_ext == '.txt') or (file_ext == '.dat')):
        file_ext = '.txt'
    file_path = file_base + file_ext
    _create_folder(file_path)
    if not overwrite:
        file_path = _create_file_name(file_path)
    metadata = OrderedDict()
    metadata['xcenter'] = xcenter
    metadata['ycenter'] = ycenter
    for i, fact in enumerate(list_fact):
        kname = 'factor' + str(i)
        metadata[kname] = fact
    with open(file_path, "w") as f:
        for line in metadata:
            f.write(str(line) + " = " + str(metadata[line]))
            f.write('\n')
    return file_path


def load_metadata_txt(file_path):
    """
    Load distortion coefficients from a text file.

    Parameters
    ----------
    file_path : str
        Path to a file.

    Returns
    -------
    tuple of floats and list
        Tuple of (xcenter, ycenter, list_fact).
    """
    if ("\\" in file_path):
        raise ValueError(
            "Please use a file path following the Unix convention")
    with open(file_path, 'r') as f:
        x = f.read().splitlines()
        list_data = []
        for i in x:
            list_data.append(float(i.split()[-1]))
    xcenter = list_data[0]
    ycenter = list_data[1]
    list_fact = list_data[2:]
    return xcenter, ycenter, list_fact
