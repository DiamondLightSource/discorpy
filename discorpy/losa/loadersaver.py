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
# E-mail: 
# Publication date: 10th July 2018
# ============================================================================
# Contributors:
# ============================================================================

"""
Module for I/O tasks:

- Load data from an image file (tif, png, jpg) or a hdf file.
- Save a 2D array as a tif/png/jpg image or a 2D, 3D array to a hdf file.
- Save a plot of data points to an image.
- Save/load metadata to/from a text/json file.
- Save/load python list.

"""

import json
import pickle
import platform
from pathlib import Path
import h5py
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib import font_manager
from collections import OrderedDict


def __correct_path(file_path):
    """
    Correct escaped sequences in WinOS file path.
    """
    if isinstance(file_path, Path):
        file_path = str(file_path)
    escape_sequences = {
        '\a': r'\a',
        '\b': r'\b',
        '\f': r'\f',
        '\n': r'\n',
        '\r': r'\r',
        '\t': r'\t',
        '\v': r'\v',
        '\0': r'\0',
    }
    for char, escaped in escape_sequences.items():
        if char in file_path:
            file_path = file_path.replace(char, escaped)
    file_path = file_path.replace('\\', '/')
    return Path(file_path)


def __get_path(file_path, check_exist=True):
    """
    Get/check a file path
    """
    if platform.system() == 'Windows':
        file_path = __correct_path(file_path)
    else:
        file_path = Path(file_path)
    if check_exist:
        if not file_path.exists():
            raise ValueError(f"No such file: {file_path}")
    return file_path


def load_image(file_path, average=True):
    """
    Load data from an image.

    Parameters
    ----------
    file_path : str
        Path to a file.
    average : bool, optional
        Average a multichannel image if True.

    Returns
    -------
    array_like
    """
    try:
        mat = np.array(Image.open(__get_path(file_path)), dtype=np.float32)
    except Exception as error:
        raise ValueError(error)
    if len(mat.shape) > 2 and average is True:
        axis_m = np.argmin(mat.shape)
        mat = np.mean(mat, axis=axis_m)
    return mat


def get_hdf_information(file_path, display=False):
    """
    Get information of datasets in a hdf/nxs file.

    Parameters
    ----------
    file_path : str
        Path to the file.
    display : bool
        Print the results onto the screen if True.

    Returns
    -------
    list_key : str
        Keys to the datasets.
    list_shape : tuple of int
        Shapes of the datasets.
    list_type : str
        Types of the datasets.
    """
    hdf_object = h5py.File(__get_path(file_path), 'r')
    keys = []
    hdf_object.visit(keys.append)
    list_key, list_shape, list_type = [], [], []
    for key in keys:
        try:
            data = hdf_object[key]
            if isinstance(data, h5py.Group):
                list_tmp = list(data.items())
                if list_tmp:
                    for key2, _ in list_tmp:
                        list_key.append(key + "/" + key2)
                else:
                    list_key.append(key)
            else:
                list_key.append(data.name)
        except KeyError:
            list_key.append(key)
            pass
    for i, key in enumerate(list_key):
        shape, dtype = None, None
        try:
            data = hdf_object[list_key[i]]
            if isinstance(data, h5py.Dataset):
                shape, dtype = data.shape, data.dtype
            list_shape.append(shape)
            list_type.append(dtype)
        except KeyError:
            list_shape.append(shape)
            list_type.append(dtype)
            pass
    hdf_object.close()
    if display:
        if list_key:
            for i, key in enumerate(list_key):
                print(key + " : " + str(list_shape[i]) + " : " + str(
                    list_type[i]))
        else:
            print("Empty file !!!")
    return list_key, list_shape, list_type


def find_hdf_key(file_path, pattern, display=False):
    """
    Find datasets matching the name-pattern in a hdf/nxs file.

    Parameters
    ----------
    file_path : str
        Path to the file.
    pattern : str
        Pattern to find the full names of the datasets.
    display : bool
        Print the results onto the screen if True.

    Returns
    -------
    list_key : str
        Keys to the datasets.
    list_shape : tuple of int
        Shapes of the datasets.
    list_type : str
        Types of the datasets.
    """
    hdf_object = h5py.File(__get_path(file_path), 'r')
    list_key, keys = [], []
    hdf_object.visit(keys.append)
    for key in keys:
        try:
            data = hdf_object[key]
            if isinstance(data, h5py.Group):
                list_tmp = list(data.items())
                if list_tmp:
                    for key2, _ in list_tmp:
                        list_key.append(key + "/" + key2)
                else:
                    list_key.append(key)
            else:
                list_key.append(data.name)
        except KeyError:
            pass
    list_dkey, list_dshape, list_dtype = [], [], []
    for _, key in enumerate(list_key):
        if pattern in key:
            list_dkey.append(key)
            shape, dtype = None, None
            try:
                data = hdf_object[key]
                if isinstance(data, h5py.Dataset):
                    shape, dtype = data.shape, data.dtype
                list_dtype.append(dtype)
                list_dshape.append(shape)
            except KeyError:
                list_dtype.append(dtype)
                list_dshape.append(shape)
                pass
    hdf_object.close()
    if display:
        if list_dkey:
            for i, key in enumerate(list_dkey):
                print(key + " : " + str(list_dshape[i]) + " : " + str(
                    list_dtype[i]))
        else:
            print("Can't find datasets with keys matching the "
                  "pattern: {}".format(pattern))
    return list_dkey, list_dshape, list_dtype


def _get_key(name, obj):
    """
    Find a key path containing 'data' in a dataset. Use with Group.visititems()
    method to walk through an HDF5 tree.
    """
    if isinstance(obj, h5py.Group):
        for key, val in obj.items():
            if key == "data" and isinstance(val, h5py.Dataset):
                return f"{obj.name}/{key}"


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
    try:
        ifile = h5py.File(__get_path(file_path), 'r')
    except Exception as error:
        raise ValueError(f"Error: {error}")
    if key_path is None:
        key_path = ifile.visititems(_get_key)  # Find the key automatically
        if key_path is None:
            raise ValueError("Please provide the key path to the dataset!")
    check = key_path in ifile
    if not check:
        raise ValueError("Couldn't open object with the key path: "
                         "{}".format(key_path))
    idata = ifile[key_path]
    shape = idata.shape
    if len(shape) < 2 or len(shape) > 3:
        raise ValueError("Require a 2D or 3D dataset!")
    if len(shape) == 2:
        mat = np.asarray(idata)
    if len(shape) == 3:
        axis = int(np.clip(axis, 0, 2))
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
                    axis = np.clip(axis, 0, 1)
                except IndexError:
                    raise
            if isinstance(index, tuple) or isinstance(index, list):
                if len(index) == 3:
                    start = index[0]
                    stop = index[1]
                    step = index[2]
                    list_index = list(range(start, stop, step))
                elif len(index) == 2:
                    start = index[0]
                    stop = index[1]
                    list_index = list(range(start, stop))
                else:
                    list_index = list(index)
                try:
                    if axis == 0:
                        mat = np.float32(idata[list_index, :, :])
                    elif axis == 1:
                        mat = np.float32(idata[:, list_index, :])
                    else:
                        mat = np.float32(idata[:, :, list_index])
                except IndexError:
                    raise
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
        ifile = h5py.File(__get_path(file_path), 'r')
    except Exception as error:
        raise ValueError(f"Error: {error}")
    check = key_path in ifile
    if not check:
        raise ValueError(f"Couldn't open object with the key: {key_path}")
    return ifile[key_path]


def _create_folder(file_path):
    """
    Create a folder to save a file if not exists.

    Parameters
    ----------
    file_path : str
        Path to a file
    """
    path = Path(file_path).resolve()
    if path.suffix:
        folder_path = path.parent
    else:
        folder_path = path
    if not folder_path.exists():
        try:
            folder_path.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise ValueError(f"Can't create : {folder_path}. Error: {e}")


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
    file_path = Path(file_path)
    file_base = file_path.stem
    file_ext = file_path.suffix
    parent_dir = file_path.parent
    if file_path.exists():
        nfile = 0
        while True:
            name_add = f"_{nfile:04d}"
            new_file_name = f"{file_base}{name_add}{file_ext}"
            new_file_path = parent_dir / new_file_name
            if new_file_path.exists():
                nfile += 1
            else:
                file_path = new_file_path
                break
    return str(file_path)


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
    file_path = __get_path(file_path, check_exist=False).resolve()
    file_ext = file_path.suffix
    if not ((file_ext == ".tif") or (file_ext == ".tiff")):
        if mat.dtype != np.uint8:
            nmin, nmax = np.min(mat), np.max(mat)
            if nmax != nmin:
                mat = np.uint8(255.0 * (mat - nmin) / (nmax - nmin))
            else:
                mat = np.uint8(mat)
    else:
        if len(mat.shape) > 2:
            axis_m = np.argmin(mat.shape)
            mat = np.mean(mat, axis=axis_m)
    _create_folder(str(file_path))
    if not overwrite:
        file_path = _create_file_name(str(file_path))
    image = Image.fromarray(mat)
    try:
        image.save(file_path)
    except Exception as error:
        raise ValueError(f"Couldn't write to file: {file_path}. Error {error}")
    return file_path


def save_plot_image(file_path, list_lines, height, width, overwrite=True,
                    dpi=100):
    """
    Save the plot of points to an image. Useful to check if the points are
    arranged properly where points on the same line having the same color.

    Parameters
    ----------
    file_path : str
        Output file path.
    list_lines : list of array_like
        List of 2D arrays. Each list is the coordinates of points on a line.
    height : int
        Height of the image.
    width : int
        Width of the image.
    overwrite : bool, optional
        Overwrite the existing file if True.
    dpi : int, optional
        The resolution in points per inch.

    Returns
    -------
    str
        Updated file path.
    """
    file_path = str(__get_path(file_path, check_exist=False))
    _create_folder(file_path)
    if not overwrite:
        file_path = _create_file_name(file_path)
    fig = plt.figure(frameon=False)
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
    except Exception as error:
        raise ValueError(f"Couldn't write to file: {file_path}. Error {error}")
    plt.close()
    return file_path


def __check_font(font_family):
    """
    Check if a specific font is available in Matplotlib.

    Parameters
    ----------
    font_family : str
        Name of the font to check.

    Returns
    -------
    bool
        True if font is available, False otherwise.
    """
    try:
        font_manager.findfont(font_family, fallback_to_default=False)
        return True
    except:
        return False


def save_residual_plot(file_path, list_data, height, width, overwrite=True,
                       dpi=100, font_family='Times New Roman'):
    """
    Save the plot of residual against radius to an image. Useful to check the
    accuracy of unwarping results.

    Parameters
    ----------
    file_path : str
        Output file path.
    list_data : array_like
        2D array. List of [residual, radius] of each point.
    height : int
        Height of the output image.
    width : int
        Width of the output image.
    overwrite : bool, optional
        Overwrite the existing file if True.
    dpi : int, optional
        The resolution in points per inch.
    font_family : str, optional
        To set the font family

    Returns
    -------
    str
        Updated file path.
    """
    file_path = str(__get_path(file_path, check_exist=False))
    _create_folder(file_path)
    if not overwrite:
        file_path = _create_file_name(file_path)
    fig = plt.figure(frameon=False)
    fig.set_size_inches(width / dpi, height / dpi)
    m_size = 0.5 * min(height / dpi, width / dpi)
    plt.rc('font', size=np.int16(m_size * 4))
    if __check_font(font_family):
        plt.rcParams['font.family'] = font_family
    plt.rcParams['font.weight'] = 'bold'
    plt.xlabel('Radius', fontweight='bold')
    plt.ylabel('Residual', fontweight='bold')
    plt.plot(list_data[:, 0], list_data[:, 1], '.', markersize=m_size)
    try:
        plt.savefig(file_path, dpi=dpi, bbox_inches='tight')
    except Exception as error:
        raise ValueError(f"Couldn't write to file: {file_path}. Error {error}")
    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)
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
    file_path = __get_path(file_path, check_exist=False).resolve()
    if file_path.suffix.lower() not in {'.hdf', '.h5', '.nxs', '.hdf5'}:
        file_path = file_path.with_suffix('.hdf')
    _create_folder(str(file_path))
    if not overwrite:
        file_path = _create_file_name(str(file_path))
    try:
        ofile = h5py.File(file_path, 'w')
    except Exception as error:
        raise ValueError(f"Couldn't write to file: {file_path}. Error {error}")
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
        Add metadata. Example:
        options={"entry/angles": angles, "entry/energy": 53}.

    Returns
    -------
    object
        hdf object.
    """
    file_path = __get_path(file_path, check_exist=False).resolve()
    if file_path.suffix.lower() not in {'.hdf', '.h5', '.nxs', '.hdf5'}:
        file_path = file_path.with_suffix('.hdf')
    _create_folder(str(file_path))
    if not overwrite:
        file_path = _create_file_name(str(file_path))
    try:
        ofile = h5py.File(file_path, 'w')
    except Exception as error:
        raise ValueError(f"Couldn't write to file: {file_path}. Error {error}")
    if len(options) != 0:
        for opt_name in options:
            opts = options[opt_name]
            for key in opts:
                if key_path in key:
                    msg = "!!! Selected key path, '{0}', can not be a " \
                          "child key-path of '{1}' !!!\n!!! Change to make " \
                          "sure they are at the same level " \
                          "!!!".format(key, key_path)
                    raise ValueError(msg)
                ofile.create_dataset(key, data=opts[key])
    data_out = ofile.create_dataset(key_path, data_shape, dtype=data_type)
    return data_out


def save_plot_points(file_path, list_points, height, width, overwrite=True,
                     dpi=100, marker="o", color="blue"):
    """
    Save the plot of points to an image. Useful to check if the points are
    arranged properly where points on the same line having the same color.

    Parameters
    ----------
    file_path : str
        Output file path.
    list_points : list of 1D-array
        List of the (y-x)-coordinates of points.
    height : int
        Height of the image.
    width : int
        Width of the image.
    overwrite : bool, optional
        Overwrite the existing file if True.
    dpi : int, optional
        The resolution in points per inch.
    marker : str
        Plot marker. Full list is at:
        https://matplotlib.org/stable/api/markers_api.html
    color : str
        Marker color. Full list is at:
        https://matplotlib.org/stable/tutorials/colors/colors.html

    Returns
    -------
    str
        Updated file path.
    """
    file_path = __get_path(file_path, check_exist=False).resolve()
    _create_folder(str(file_path))
    if not overwrite:
        file_path = _create_file_name(str(file_path))
    fig = plt.figure(frameon=False)
    fig.set_size_inches(width / dpi, height / dpi)
    ax = plt.Axes(fig, [0., 0., 1.0, 1.0])
    ax.set_axis_off()
    fig.add_axes(ax)
    plt.axis((0, width, 0, height))
    m_size = 0.5 * min(height / dpi, width / dpi)
    for point in list_points:
        plt.plot(point[1], height - point[0], marker, color=color,
                 markersize=m_size)
    try:
        plt.savefig(file_path, dpi=dpi)
    except IOError:
        raise ValueError("Couldn't write to file {}".format(file_path))
    plt.close()
    return file_path


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
    file_path = __get_path(file_path, check_exist=False).resolve()
    if file_path.suffix.lower() not in {'.txt', '.dat'}:
        file_path = file_path.with_suffix('.txt')
    _create_folder(str(file_path))
    if not overwrite:
        file_path = _create_file_name(str(file_path))
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
    with open(__get_path(file_path), 'r') as f:
        x = f.read().splitlines()
        list_data = []
        for i in x:
            list_data.append(float(i.split()[-1]))
    xcenter = list_data[0]
    ycenter = list_data[1]
    list_fact = list_data[2:]
    return xcenter, ycenter, list_fact


def __numpy_encoder(obj):
    if isinstance(obj, (np.int8, np.int16, np.int32, np.int64,
                        np.uint8, np.uint16, np.uint32, np.uint64)):
        return int(obj)
    elif isinstance(obj, (np.float16, np.float32, np.float64)):
        return float(obj)
    elif isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    raise TypeError(f"Object of type '{type(obj).__name__}' "
                    f"is not JSON serializable")


def save_metadata_json(file_path, xcenter, ycenter, list_fact, overwrite=True):
    """
    Write metadata to a JSON file.

    Parameters
    ----------
    file_path : str
        Output file path.
    xcenter : float
        Center of distortion in x-direction.
    ycenter : float
        Center of distortion in y-direction.
    list_fact : list of float
        Coefficients of a polynomial.
    overwrite : bool, optional
        Overwrite an existing file if True.

    Returns
    -------
    str
        Updated file path.
    """
    file_path = __get_path(file_path, check_exist=False).resolve()
    if file_path.suffix.lower() != '.json':
        file_path = file_path.with_suffix('.json')
    _create_folder(str(file_path))
    if not overwrite:
        file_path = _create_file_name(str(file_path))
    metadata = {
        'xcenter': float(xcenter),
        'ycenter': float(ycenter),
        'list_fact': list_fact
    }
    with open(file_path, "w") as f:
        json.dump(metadata, f, indent=4, default=__numpy_encoder)
    return file_path


def load_metadata_json(file_path):
    """
    Load distortion coefficients from a JSON file.

    Parameters
    ----------
    file_path : str
        Path to a JSON file.

    Returns
    -------
    tuple of floats and list
        Tuple of (xcenter, ycenter, list_fact).
    """
    with open(__get_path(file_path), 'r') as f:
        metadata = json.load(f)
    xcenter = metadata['xcenter']
    ycenter = metadata['ycenter']
    list_fact = metadata['list_fact']
    return xcenter, ycenter, list_fact


def load_python_list(file_path):
    """
    Load a Python list from a pickle file (.pkl).

    Parameters
    ----------
    file_path : str
        Path to the pickle file.

    Returns
    -------
    list
        The Python list.
    """
    with open(__get_path(file_path), 'rb') as f:
        loaded_data = pickle.load(f)
    return loaded_data


def save_python_list(file_path, python_list, overwrite=True):
    """
    Write python list to a pickle file (.pkl).

    Parameters
    ----------
    file_path : str
        Output file path.
    python_list : list
        Python list.
    overwrite : bool, optional
        Overwrite an existing file if True.

    Returns
    -------
    str
        Updated file path.
    """
    file_path = __get_path(file_path, check_exist=False).resolve()
    if file_path.suffix.lower() != '.pkl':
        file_path = file_path.with_suffix('.pkl')
    _create_folder(str(file_path))
    if not overwrite:
        file_path = _create_file_name(str(file_path))
    with open(file_path, 'wb') as f:
        pickle.dump(python_list, f)
    return file_path


def find_file(path):
    """
    Search file

    Parameters
    ----------
    path : str
        Path and pattern to find files.

    Returns
    -------
    str or list of str
        List of files.
    """
    path = __correct_path(path)
    file_paths = list(path.parent.glob(path.name))
    if not file_paths:
        raise FileNotFoundError(f"No files found matching: {path}")
    return sorted([file.as_posix() for file in file_paths])
