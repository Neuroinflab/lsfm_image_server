#!/usr/bin/env python
# -*- coding: utf-8 -*

import os
import sys
import json
import tqdm
import h5py
import copy
import errno
import imageio
import inspect
import h5py_cache

import tifffile as tf
import numpy as np
import nibabel as nib
import SimpleITK as sitk
import datetime as dt

import dump_metadata as dm

from utils import PathUtil
from utils import InputImageType
from lxml import etree
from collections import namedtuple
from scipy.misc import imsave

import logging

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)

M_ORIGIN = 'origin'
M_SHAPE = 'shape'
M_VOXEL_SIZE = 'spacing'

RAS_INVERSE_DIRECTION = np.array([-1, -1, 1])
RAS_DIRECTION = np.array([1, 1, -1])

INVERSE_AFFINE_SIGN = np.array([[1, 1, -1, -1],
                                [1, 1, -1, -1],
                                [-1, -1, 1, 1],
                                [1, 1, 1, 1]])

SIGN_LPI_RAS_T = np.array([1, 1, -1])
SIGN_LPI_RAS_A = np.array([[-1, 1, -1], [1, -1, -1], [1, 1, -1]])

SIGN_RAS_T = np.array([-1, -1, -1])
SIGN_RAS_A = np.array([[-1, 1, -1],
                       [1, -1, -1],
                       [1, 1, -1]])
DIRECTION_RAS = (1.0, 0.0, 0.0,
                 0.0, 1.0, 0.0,
                 0.0, 0.0, -1.0)

SIGN_LPI_T = np.array([-1, -1, 1])
SIGN_LPI_A = np.array([[1, 1, -1], [1, 1, -1], [-1, -1, 1]])
DIRECTION_LPI = (-1.0, 0.0, 0.0,
                 0.0, -1.0, 0.0,
                 0.0, 0.0, 1.0)

TransformTuple = namedtuple('TransformTuple', 'order name type invert')


class ImageProxy(object):
    def __init__(self,
                 channel_name,
                 specimen_id,
                 multifile_prefix,
                 metadata,
                 xml_file_name,
                 RAM_limit,
                 orientation='RAS',
                 is_multichannel=False,
                 is_segmentation=False):

        self.channel_name = channel_name
        self.xml_file_name = xml_file_name
        self.specimen_id = specimen_id
        self.multifile_prefix = multifile_prefix
        self.metadata = metadata
        self.RAM_limit = RAM_limit
        self.is_multichannel = is_multichannel
        self.is_segmentation = is_segmentation
        self.orientation = orientation

        self.pixel_type = self._get_pixel_type()
        self.direction = (-1.0, 0.0, 0.0,
                          0.0, -1.0, 0.0,
                          0.0, 0.0, 1.0)

        self.is_stream = False
        self.is_nifti = False

        self.data_shape = None
        self.voxel_size = None

    @property
    def nifti_affine(self):
        return None

    @property
    def nifti_header(self):
        return None

    @property
    def stack_size(self):
        return None

    @property
    def stack_shape(self):
        return None

    @property
    def overhead_stack_shape(self):
        return None

    @property
    def image_shape(self):
        return self.data_shape[:3]

    @property
    def data_type(self):
        if self.metadata.bits_per_sample == 16:
            return np.int16
        elif self.metadata.bits_per_sample == 8:
            return np.int8
        else:
            raise ValueError("Unknown data type")

    @property
    def plane_size(self):
        return self.data_shape[1] * self.data_shape[2] * self.metadata.bits_per_sample / 8. / (1000 ** 3)

    @property
    def size(self):
        return self.plane_size * self.data_shape[0]

    @property
    def affine(self):
        return np.eye(4, 4) * np.hstack([self.voxel_size, [1.]])

    @property
    def num_of_stacks(self):
        if self.is_stream and self.stack_size < self.data_shape[0]:
            return self.data_shape[0] // self.stack_size
        else:
            return 1

    @property
    def interpolation_type(self):
        if self.is_segmentation:
            return sitk.sitkNearestNeighbor
        else:
            return sitk.sitkLinear

    @property
    def cast_interpolation_type(self):
        if self.is_segmentation:
            return sitk.sitkUInt16
        elif self.is_multichannel:
            return sitk.sitkVectorFloat32
        else:
            return sitk.sitkFloat32

    @property
    def cast_output_type(self):
        if self.is_segmentation:
            return sitk.sitkUInt16
        elif self.is_multichannel:
            return sitk.sitkVectorUInt16
        else:
            return sitk.sitkUInt16

    def _read_image(self, img_path):
        raise NotImplementedError

    def stream_data(self):
        raise NotImplementedError

    def _get_pixel_type(self):
        '''image_io = itk.ImageIOFactory.CreateImageIO(str(self.metadata.file_path),
                                                    itk.ImageIOFactory.ReadMode)
        image_io.SetFileName(str(self.metadata.file_path))
        return image_io.GetPixelType()'''
        if self.is_segmentation:
            return sitk.sitkUInt16

        if self.is_multichannel:
            return sitk.__getattribute__('sitkVectorFloat' + str(self.metadata.bits_per_sample))

        return sitk.__getattribute__('sitkUInt' + str(self.metadata.bits_per_sample))

    def _get_sitk_pixel_type(self):
        image_io = sitk.ReadImage(str(self.metadata.file_path))
        return image_io.GetPixelID()

    def _get_transpose(self):
        target_space = np.array(list('RAS'))
        reverse_axes_dict = dict(zip('RASLPI', 'LPIRAS'))
        source_space = np.array(list(self.orientation))

        transpose = []
        flip = []
        for i, letter in enumerate(target_space):
            if letter == source_space[i]:
                transpose.append(i)
            elif letter == reverse_axes_dict[source_space[i]]:
                transpose.append(i)
            else:
                if letter in source_space:
                    to = np.argwhere(letter == source_space)[0][0]
                    transpose.append(to)
                else:
                    reverse_source_space = np.array([reverse_axes_dict[l] for l in source_space])
                    to = np.argwhere(letter == reverse_source_space)[0][0]
                    transpose.append(to)

        source_space = source_space[transpose]

        for i, letter in enumerate(source_space):
            if letter == target_space[i]:
                flip.append(1)
            elif letter == reverse_axes_dict[target_space[i]]:
                flip.append(-1)
            else:
                logger.debug("Something went wrong with reorienting RAI codes")

        return np.array(transpose), np.array(flip)

    def __repr__(self):
        return "The plane size is: {:.2f} GB, \nThe whole size is: {:.2f} GB" \
               "\nStack size is: {} planes \nIs stream: {} \nAffine: \n{}\n" \
               "Pixel type: {}\n".format(self.plane_size,
                                         self.size,
                                         self.stack_size,
                                         self.is_stream,
                                         self.affine,
                                         self.pixel_type)

    @staticmethod
    def reverse_axes(img_data, reverse=[False, False, False]):
        if reverse[0]:
            img_data = img_data[::-1, :, :]
        if reverse[1]:
            img_data = img_data[:, ::-1, :]
        if reverse[2]:
            img_data = img_data[:, :, ::-1]

        return img_data

    @staticmethod
    def get_available_ram():
        assert os.name == 'posix'
        available_ram = float(os.popen("free -m").readlines()[1].split()[-1]) / 1000.
        logger.debug("Available RAM: {:.2f}".format(available_ram))
        return available_ram

    @staticmethod
    def get_image_proxy_class(meta):
        """

        :type meta: dm.ImageMetaData
        """
        logger.debug(meta.file_type)
        if meta.file_type == InputImageType.MultipleTiff.name:
            return StreamableTiffProxy
        elif meta.file_type == InputImageType.OmeTiff.name:
            return StreamableOMEProxy
        elif meta.file_type == InputImageType.Nifti.name:
            return NiftiProxy
        elif meta.file_type == InputImageType.ImageJTiff.name or meta.file_type == InputImageType.Tiff.name:
            return NonStreamableTiffProxy

        else:
            raise ValueError("Unsupported file format specified in image metadata")


class StreamableTiffProxy(ImageProxy):
    """
    This is currently implemented for multifile tiff data, _read_image and
    stream_data methods need to be reimplemented for other file types, or for single
    file OME Tiff.
    """

    def __init__(self, *args, **kwargs):

        super(StreamableTiffProxy, self).__init__(*args, **kwargs)

        try:
            self.voxel_size = np.array([self.metadata.voxel_size_z,
                                        self.metadata.voxel_size_y,
                                        self.metadata.voxel_size_x], dtype=np.float64)

            self.data_shape = np.array([self.metadata.image_size_z,
                                        self.metadata.image_size_y,
                                        self.metadata.image_size_x], dtype=np.int32)

            self.file_dir, self.file_name = os.path.split(self.metadata.file_path)

        except AttributeError:
            logger.error("Attribute not found in metadata object", exc_info=True)
            raise

        self.is_stream = True

    @property
    def pattern(self):
        return os.path.join(self.file_dir, self.multifile_prefix)

    @property
    def stack_size(self):
        stack_size = int(self.RAM_limit // self.plane_size)
        if stack_size > self.data_shape[0]:
            return self.data_shape[0]
        else:
            return stack_size

    @property
    def origin(self):
        return np.array([0., 0., 0.])

    @property
    def stack_shape(self):

        return np.array([self.stack_size, self.data_shape[1], self.data_shape[2]],
                        dtype=np.uint16)

    @property
    def overhead_stack_shape(self):
        return np.array([self.data_shape[0] - (self.stack_size * self.num_of_stacks),
                         self.data_shape[1],
                         self.data_shape[2]],
                        dtype=np.uint16)

    def _craft_img_path(self, idx):
        return self.pattern % idx

    def _read_image(self, img_path):
        return tf.imread(img_path, key=0)

    def stream_data(self):
        assert ImageProxy.get_available_ram() > self.stack_size * self.plane_size

        current_idx = 0
        while current_idx < self.data_shape[0]:
            logger.debug("Current streaming index: {}".format(current_idx))
            if current_idx + self.stack_size < self.data_shape[0]:
                next_indices = np.arange(current_idx, current_idx + self.stack_size)
            else:
                next_indices = np.arange(current_idx, self.data_shape[0])
            try:
                logger.debug("Next stream slice: {}".format(next_indices))
                data = self._read_image(self._craft_img_path(next_indices[0]))
                for idx in tqdm.tqdm(next_indices[1:]):
                    next_data = self._read_image(self._craft_img_path(idx))
                    data = np.dstack(([data, next_data]))
            except IOError:
                logger.error("Could not open file", exc_info=True)
                raise

            current_idx = current_idx + self.stack_size

            yield data.transpose(2, 0, 1)

    def __repr__(self):
        return "The plane size is: {:.2f} GB, \nThe whole size is: {:.2f} GB" \
               "\nStack size is: {} planes \nIs stream: {} \nAffine: \n{}" \
               "Stack shape: {}\n Overhead stack shape: {}\n".format(self.plane_size,
                                                                     self.size,
                                                                     self.stack_size,
                                                                     self.is_stream,
                                                                     self.affine,
                                                                     self.stack_shape,
                                                                     self.overhead_stack_shape)


class StreamableOMEProxy(StreamableTiffProxy):
    def __init__(self, *args, **kwargs):

        super(StreamableOMEProxy, self).__init__(*args, **kwargs)

        try:
            self.voxel_size = np.array([self.metadata.voxel_size_z,
                                        self.metadata.voxel_size_y,
                                        self.metadata.voxel_size_x], dtype=np.float64)

            self.data_shape = np.array([self.metadata.image_size_z,
                                        self.metadata.image_size_y,
                                        self.metadata.image_size_x], dtype=np.int32)

            self.file_dir, self.file_name = os.path.split(self.metadata.file_path)

        except AttributeError:
            logger.error("Attribute not found in metadata object", exc_info=True)
            raise

    @property
    def pattern(self):
        return None

    def _read_image(self, img_path):
        return tf.TiffFile(img_path)

    def stream_data(self):
        assert ImageProxy.get_available_ram() > self.stack_size * self.plane_size

        try:
            tiff_file = self._read_image(self.metadata.file_path)
        except IOError:
            logger.error("Could not open OME file.", exc_info=True)
            raise

        current_idx = 0
        while current_idx < self.data_shape[0]:
            next_indices = slice(current_idx, current_idx + self.stack_size)
            if current_idx + self.stack_size > self.data_shape[0]:
                next_indices = slice(current_idx, current_idx + self.overhead_stack_shape[0])
            data = tiff_file.asarray(next_indices, 0)[:next_indices.stop - next_indices.start]
            current_idx += self.stack_size
            yield data


class NonStreamableTiffProxy(ImageProxy):
    def __init__(self, *args, **kwargs):

        super(NonStreamableTiffProxy, self).__init__(*args, **kwargs)

        try:
            self.voxel_size = np.array([self.metadata.voxel_size_z,
                                        self.metadata.voxel_size_y,
                                        self.metadata.voxel_size_x], dtype=np.float64)

            self.data_shape = np.array([self.metadata.image_size_z,
                                        self.metadata.image_size_y,
                                        self.metadata.image_size_x], dtype=np.int32)

            self.file_dir, self.file_name = os.path.split(self.metadata.file_path)

        except AttributeError:
            logger.error("Attribute not found in metadata object", exc_info=True)
            raise

        self.is_stream = False

        try:
            self.image = self._read_image(self.metadata.file_path)
        except IOError:
            logger.error("Could not open file", exc_info=True)
            raise

    @property
    def stack_size(self):
        return self.data_shape[0]

    def _read_image(self, img_path):
        return tf.TiffFile(img_path)

    def stream_data(self):
        assert ImageProxy.get_available_ram() > self.size * 2
        data = [self.image.asarray()]
        return data


class NiftiProxy(ImageProxy):
    def __init__(self, *args, **kwargs):

        super(NiftiProxy, self).__init__(*args, **kwargs)

        try:
            self.voxel_size = np.array([self.metadata.voxel_size_x,
                                        self.metadata.voxel_size_y,
                                        self.metadata.voxel_size_z], dtype=np.float64)

            self.data_shape = np.array([self.metadata.image_size_x,
                                        self.metadata.image_size_y,
                                        self.metadata.image_size_z], dtype=np.int32)

            self.file_dir, self.file_name = os.path.split(self.metadata.file_path)
        except AttributeError:
            logger.error("Attribute not found in metadata object", exc_info=True)
            raise

        self.is_stream = False
        self.is_nifti = True

        try:
            self.image = self._read_image(self.metadata.file_path)
        except IOError:
            logger.error("Could not open file", exc_info=True)
            raise

        assert ImageProxy.get_available_ram() > self.size
        self.data = self.image.get_data()
        if len(self.data.shape) == 5:
            self.data = self.data[:, :, :, 0, :]

        self.data_shape = np.array(self.data.shape)

        logger.debug("Pixel type for channel {}: {}".format(self.channel_name, self.pixel_type))

    @property
    def stack_size(self):
        return self.data_shape[0]

    @property
    def nifti_affine(self):
        return self.image.affine

    @property
    def nifti_header(self):
        return self.image.header

    @property
    def data_type(self):
        return self.data.dtype

    @property
    def affine(self):
        affine = self.image.affine
        affine[:3, 3] *= np.array([-1, -1, 1])
        return affine

    def _read_image(self, img_path):
        return nib.load(img_path)

    def stream_data(self):
        return [self.data]


def logging_decor(func):
    def func_wrapper(self, *args, **kwargs):
        options = dict()
        options.update(kwargs)
        args_name = inspect.getargspec(func)[0][1:]

        for i, arg in enumerate(args):
            if hasattr(arg, '__dict__'):
                options.update(arg.__dict__)
            else:
                options.update({args_name[i]: arg})

        str_options = '\n'.join("%s=%r" % (key, val) for (key, val) in options.iteritems())

        with open('/proc/{}/cmdline'.format(os.getpid())) as c:
            str_options = '\n'.join([str_options, "cmdline={}".format(' '.join(c.read().split(b'\x00'))
)])

        str_options = np.string_(str_options)
        self.log_operation(func.__name__, str_options)

        return func(self, *args, **kwargs)
    return func_wrapper


class LightMicroscopyHDF(object):
    def __init__(self, file_path, access_mode='a'):

        self.file_path = file_path
        self.channels = dict()

        try:
            # self.h5_file = h5py.File(self.file_path, 'a', libver='latest')
            self.h5_file = h5py_cache.File(self.file_path, access_mode, chunk_cache_mem_size=1024 ** 3, libver='latest')
        except IOError as e:
            if e.errno == errno.EACCES:
                self.h5_file = h5py_cache.File(self.file_path, 'r', chunk_cache_mem_size=1024 ** 3,
                                               libver='latest')
            else:
                logger.error("File not found or couldn't open file", exc_info=True)
                raise

        self.root = self.h5_file['/']
        self.bdv = BigDataViewer(h5_file=self.h5_file)

        if self.number_of_channels is None:
            self.root.attrs['LSFM_setups_count'] = 0
        else:
            self._populate_channels()

    def __repr__(self):
        metadata = '\n'.join("%s=%r" % (key, val) for (key, val) in self.extended_metadata.iteritems())
        channel_info = '\n'.join(str(channel) for channel in self.channels.values())
        return '\n---------\n'.join([metadata, channel_info])

    @property
    def extended_metadata(self):
        if PathUtil.get_extended_metadata_path() not in self.h5_file:
            self.h5_file.create_group(PathUtil.get_extended_metadata_path())

        return self.h5_file[PathUtil.get_extended_metadata_path()].attrs

    @logging_decor
    def add_metadata(self, m_key, m_value):
        self.extended_metadata[m_key] = m_value
        self.h5_file.flush()

    def log_operation(self, cmd_name, cmd_line):
        if self.h5_file.mode == 'r':
            logger.info('File in read mode only, cannot log operations!')
            return
        log_path = PathUtil.get_log_path(dt.datetime.now(), cmd_name)
        logger.debug('log path: {:} \n command-name: {:}\n command: {:}'.format(log_path,
                                                                                cmd_name,
                                                                                cmd_line))
        self.h5_file.create_dataset(name=log_path, data=cmd_line)

    def close(self):
        self.h5_file.close()

    def open(self):
        self.h5_file = h5py_cache.File(self.file_path, 'a', chunk_cache_mem_size=1024 ** 3, libver='latest')

    @property
    def channel_names(self):
        if '/LSFM' in self.h5_file:
            return self.h5_file['/LSFM'].keys()
        else:
            return "object has no channels"

    @property
    def number_of_channels(self):
        if self.root.attrs.__contains__('LSFM_setups_count'):
            return self.root.attrs['LSFM_setups_count']
        else:
            return None

    @number_of_channels.setter
    def number_of_channels(self, value):
        self.root.attrs['LSFM_setups_count'] += value

    def _populate_channels(self):
        for channel_name in self.channel_names:
            channel = HDFChannel(self.h5_file, channel_name)
            self.channels[channel_name] = channel

    def get_channel(self, _channel_name):
        """

        :type _channel_name: str
        :rtype :HDFChannel
        """
        return self.channels[_channel_name]

    @logging_decor
    def export_image(self, export_cmd):
        """

        :type export_cmd: ExportCmd
        """

        channel = self.get_channel(export_cmd.channel_name)
        return channel.export_image(export_cmd)

    @logging_decor
    def export_slices(self, export_slices_cmd):

        channel = self.get_channel(export_slices_cmd.channel_name)
        # if export_slices_cmd.resample:
            # export_slices_cmd.ref_channel = self.get_channel(export_slices_cmd.ref_channel)

        return channel.export_slices(export_slices_cmd)

    @logging_decor
    def write_channel(self, image_proxy, create_bdv=True):

        data_sets = []

        def create_channel_datasets(_channel):
            """

            :type _channel: ProxyChannel
            """
            for pyramid_level in _channel.levels:
                try:
                    data_set = self.h5_file.create_dataset(pyramid_level.path,
                                                           pyramid_level.shape,
                                                           chunks=pyramid_level.chunks,
                                                           dtype=channel.image.data_type,
                                                           compression='gzip',
                                                           compression_opts=6)
                    data_sets.append(data_set)
                except RuntimeError:
                    logger.error("Dataset probably already present in file", exc_info=True)

        def write_metadata(_channel):
            channel_metadata = {'origin': [pl.origin for pl in _channel.levels],
                                'spacing': [pl.voxel_size for pl in _channel.levels],
                                'shape': [pl.shape for pl in _channel.levels],
                                'is_multichannel': _channel.image.is_multichannel,
                                'is_segmentation': _channel.image.is_segmentation,
                                'is_nifti': _channel.image.is_nifti,
                                'pixel_type': _channel.image.pixel_type}

            self.h5_file[_channel.image_path].attrs.update(channel_metadata)
            self.number_of_channels = 1

        def write_data(_channel):
            """

            :type _channel: ProxyChannel
            """

            data_idx = 0
            for data in _channel.image.stream_data():
                meta_image = MetaImage.meta_image_from_channel(data, _channel)
                # xkcd
                # meta_image.set_direction_RAS()
                if _channel.image.is_stream:
                    logger.debug("Processing data : {}/{}".format(data_idx + 1,
                                                                  _channel.image.num_of_stacks))
                for level_num, pyramid_level in enumerate(_channel.levels):

                    if level_num == 0:
                        writing_slice = slice(0, None)

                        if _channel.image.is_stream:
                            writing_slice = slice(data_idx * pyramid_level.stack_shape[0],
                                                  data_idx * pyramid_level.stack_shape[0] + data.shape[0])
                        try:
                            data_sets[level_num][writing_slice, :, :] = data
                        except TypeError:
                            logger.error("Probably wrong data shape", exc_info=True)
                            raise
                    else:

                        # scale_factors = _channel.resolutions[level_num - 1] / _channel.resolutions[level_num]
                        scale_factors = _channel.relative_resolution_scales[level_num][:len(_channel.image.image_shape)]
                        logger.debug("Scale factors for resampling: {}".format(scale_factors))
                        meta_image = ImageProcessor.resample_image(meta_image, scale_factors)
                        # xkcd
                        # meta_image.set_direction_RAS()
                        writing_slice = slice(0, None)

                        if _channel.image.is_stream:
                            writing_slice = slice(data_idx * pyramid_level.stack_shape[0],
                                                  data_idx * pyramid_level.stack_shape[0] +
                                                  meta_image.data_shape[0])
                        try:
                            data_sets[level_num][writing_slice, :, :] = meta_image.data
                        except TypeError:
                            logger.error("Probably wrong data shape", exc_info=True)
                            raise

                    if data_idx == 0:
                        pyramid_level.origin = meta_image.origin
                        pyramid_level.voxel_size = meta_image.voxel_size
                        pyramid_level.affine = meta_image.affine
                data_idx += 1

                write_metadata(_channel)

        def create_big_data_viewer_links(_channel):
            """

            :type _channel: ProxyChannel
            """

            if np.all(_channel.image.data_shape[0] <= _channel.image.data_shape):
                resolutions = _channel.resolutions[:, ::-1]
                subdivisions = _channel.subdivisions[:, ::-1]
                voxel_size = _channel.image.voxel_size[::-1]
                image_shape = _channel.image.data_shape[::-1]
            else:
                resolutions = _channel.resolutions
                subdivisions = _channel.subdivisions
                voxel_size = _channel.image.voxel_size
                image_shape = _channel.image.data_shape

            # let's check if datasets already present
            if _channel.resolutions_path in self.h5_file or _channel.subdivisions_path in self.h5_file:
                raise ValueError("Big Data Viewer setups already present in file!")
                # should we delete and overwrite them?

            self.h5_file.create_dataset(_channel.resolutions_path,
                                        _channel.resolutions.shape,
                                        dtype=np.float64)

            self.h5_file.create_dataset(_channel.subdivisions_path,
                                        _channel.subdivisions.shape,
                                        dtype=np.int16)

            self.h5_file[_channel.resolutions_path][...] = resolutions
            self.h5_file[_channel.subdivisions_path][...] = subdivisions

            for pyramid_level in _channel.levels:
                self.h5_file[pyramid_level.bdv_path] = h5py.SoftLink('/' + pyramid_level.path)

            affine = nib.affines.from_matvec(np.eye(3, 3) * np.array(voxel_size),
                                             np.array([0, 0, 0]))

            self.bdv.add_setup(image_shape, voxel_size, affine)
            self.bdv.create_xml_from_setups(_channel.image.xml_file_name, self.file_path)

        channel = ProxyChannel(image_proxy, self.bdv.number_of_bdv_channels)
        create_channel_datasets(channel)
        write_data(channel)
        if create_bdv and not channel.image.is_multichannel and not channel.image.is_segmentation:
            create_big_data_viewer_links(channel)

        self.channels[image_proxy.channel_name] = HDFChannel(self.h5_file, image_proxy.channel_name)


class Affine(object):
    def __init__(self, affine_name, affine_path=None):
        if affine_path is not None:
            try:
                self.itk_affine = sitk.ReadTransform(affine_path)
                self.itk_inv_affine = self.itk_affine.GetInverse()
            except RuntimeError:
                raise RuntimeError("Affine file format is most likely not compliant"
                                   "with itk")
            self.nifti_affine = Affine.get_affine_matrix(affine_path)
            self.inv_nifti_affine = np.linalg.inv(self.nifti_affine)
            self.affine_name = affine_name

    def set_affine(self, nifti_affine, itk_params, itk_fixed_params):
        self.nifti_affine = nifti_affine
        self.inv_nifti_affine = np.linalg.inv(self.nifti_affine)
        self.itk_affine = sitk.AffineTransform(3)
        self.itk_affine.SetParameters(itk_params)
        self.itk_affine.SetFixedParameters(itk_fixed_params)
        self.itk_inv_affine = self.itk_affine.GetInverse()

    def write_to_itk_file(self, file_path):
        sitk.WriteTransform(self.itk_affine, file_path)

    def get_direction_matrix(self):
        voxel_sizes = nib.affines.voxel_sizes(self.nifti_affine)
        dir_matrix = self.nifti_affine[:3, :3].copy()
        dir_matrix /= voxel_sizes
        dir_matrix *= RAS_INVERSE_DIRECTION[:, np.newaxis]
        return dir_matrix

    @staticmethod
    def compute_offset(matrix, center, translation):
        '''
        Computes translation vector relative to center of rotation

        Parameters
        ----------
        param matrix: 3D affine matrix
        type matrix: np.ndarray of float, (3, 3)

        param center: 3D vector defining center of rotation
        type center: ndarray of float, (3,)

        param translation: original 3D translation vector
        type translation: ndarray of float (3,)

        Returns
        -------
        rtype: ndarray of float (3,)

        based on:
           https://github.com/hinerm/ITK/blob/master/Modules/Core/Transform/
           include/itkMatrixOffsetTransformBase.hxx
        '''
        offset = translation + center
        offset -= matrix.dot(center)
        return offset

    @staticmethod
    def get_affine_matrix(fname):
        '''
        Reads from a file an itk matrix defining affine transformation and converts
        it to a nifti-compatible format

        Parameters
        ----------
        param fname: file name (txt)
        type fname: str

        Returns
        -------
        rtype: ndarray of float (4, 4)
        '''
        transform = sitk.ReadTransform(fname)

        params = np.array(transform.GetParameters())
        affine = params[:9].reshape(3, 3)
        translation = params[9:]

        center = np.array(transform.GetFixedParameters())
        translation = Affine.compute_offset(affine, center, translation)
        translation = translation.reshape(-1, 1) * [[-1], [-1], [1]]
        sign_flip = np.array([[1, 1, -1], [1, 1, -1], [-1, -1, 1]])
        affine *= sign_flip

        return np.vstack([np.hstack([affine, translation]),
                          [0., 0., 0., 1.]])

    @staticmethod
    def get_direction_matrix_from_itk_affine(itk_affine):

        params = np.array(itk_affine.GetParameters())
        affine = params[:9].reshape(3, 3)
        translation = params[9:]
        center = np.array(itk_affine.GetFixedParameters())
        translation = Affine.compute_offset(affine, center, translation)

        affine = nib.affines.from_matvec(affine, translation)

        voxel_sizes = nib.affines.voxel_sizes(affine)
        dir_matrix = affine[:3, :3].copy()
        dir_matrix /= voxel_sizes
        # dir_matrix *= INVERSE_AFFINE_SIGN[:3, :3]
        # dir_matrix *= RAS_DIRECTION[:, np.newaxis]
        return dir_matrix

    def __repr__(self):
        return "affine: {}".format(self.affine_name)


TransposeTuple = namedtuple('TransposeTuple', 'transpose, flip, inv_transpose, inv_flip,'
                            'ltranspose, lflip, linv_transpose, linv_flip')


class ImageExporter(object):
    def __init__(self, export_cmd, pyramid_level):
        """
        :type pyramid_level:
        :type export_cmd: ExportCmd
        :param export_cmd:
        :param pyramid_level:
        """
        self.export_cmd = export_cmd
        self.pyramid_level = pyramid_level

        self.transpose_tuple = TransposeTuple(*self.get_axes_transpose(ndim=len(self.pyramid_level.shape)))
        self.shape = np.array(pyramid_level.shape)[self.transpose_tuple.transpose]
        self.origin = np.abs(np.array(pyramid_level.origin[self.transpose_tuple.transpose])) * np.array([1, 1, -1.])
        self.voxel_size = np.array(pyramid_level.voxel_size[self.transpose_tuple.transpose])
        self.affine_vp, self.affine_pv = ImageProcessor.get_affines(self.voxel_size,
                                                                    self.origin)
        self.ct = CompositeTransform(self.export_cmd.list_of_transforms)

    def get_axes_transpose(self, ndim=3):
        transpose, flip = np.arange(ndim), np.ones(ndim)
        inv_transpose, inv_flip = np.arange(ndim), np.ones(ndim)
        transpose[:3], flip[:3] = ImageProcessor.get_transpose(self.export_cmd.input_orientation)
        inv_transpose[:3] = np.argsort(transpose[:3])
        inv_flip[:3] = flip[inv_transpose[:3]] < 0

        return transpose[:3], flip[:3], inv_transpose[:3], inv_flip[:3], transpose, flip, inv_transpose, inv_flip

    def _compute_indices(self):

        self.start_os_mm = self.export_cmd.phys_origin
        self.chunk_size_in_voxels = (self.export_cmd.phys_size // self.voxel_size).astype(np.int)
        self.end_os_mm = np.array(self.export_cmd.phys_origin) + \
                         np.array(self.export_cmd.phys_size) * np.array([1, 1, -1.])

        #self.start_is_mm = self.ct.composite.TransformPoint(self.start_os_mm)
        self.start_is_mm = np.abs(self.ct.composite.TransformPoint(self.start_os_mm)) * np.array([1, 1, -1.])
        self.start_index = np.abs(np.around(self.affine_pv.TransformPoint(self.start_is_mm)).astype(np.int))

        self.end_index = self.start_index + self.chunk_size_in_voxels
        self.start_index = np.array(self.start_index)[self.transpose_tuple.inv_transpose]
        self.end_index = np.array(self.end_index)[self.transpose_tuple.inv_transpose]

    def _prepare_chunk(self):
        self.chunk = self.pyramid_level.get_level_chunk_data(self.start_index, self.end_index,
                                                             self.transpose_tuple.inv_flip)

    def _transpose_chunk(self):

        self.chunk = self.chunk.transpose(self.transpose_tuple.ltranspose)
        self.chunk = ImageProcessor.reverse_axes(self.chunk, self.transpose_tuple.flip < 0)

    def _resample_chunk(self):

        logger.debug("_resample_chunk START_IS_MM {}".format(self.start_is_mm))

        self.chunk_affine = nib.affines.from_matvec(np.diag(self.voxel_size) * SIGN_RAS_A, self.start_is_mm)
        self.meta_image = MetaImage(self.chunk, self.chunk_affine,
                                    self.pyramid_level.pixel_type,
                                    DIRECTION_RAS,
                                    self.pyramid_level.is_segmentation,
                                    self.pyramid_level.is_nifti)

        if self.export_cmd.resample_image:
            self.scale_factors = self.voxel_size / self.export_cmd.output_resolution
        else:
            self.scale_factors = np.array([1., 1., 1.])

        self.meta_image = ImageProcessor.resample_image(self.meta_image, self.scale_factors)

    def _transform_chunk(self):

        self.chunk = ImageProcessor.apply_transformations(self.meta_image, self.ct)

    def process(self):
        self._compute_indices()
        self._prepare_chunk()
        self._transpose_chunk()
        self._resample_chunk()
        self._transform_chunk()

        return self.chunk


class ImageRegionExporter(ImageExporter):

    def _compute_indices(self):
        self.start, self.end = ImageProcessor.get_bounding_box(self.export_cmd.segmentation,
                                                               self.export_cmd.region_id)

        self.start_is_mm = self.export_cmd.segmentation.TransformIndexToPhysicalPoint(self.start)
        self.end_is_mm = self.export_cmd.segmentation.TransformIndexToPhysicalPoint(self.end)
        self.start_index = np.abs(np.around(self.affine_pv.TransformPoint(self.start_is_mm)).astype(np.int))
        self.end_index = np.abs(np.around(self.affine_pv.TransformPoint(self.end_is_mm)).astype(np.int))
        self.start_index = np.array(self.start_index)[self.transpose_tuple.inv_transpose]
        self.end_index = np.array(self.end_index)[self.transpose_tuple.inv_transpose]


class ImageWholeExporter(ImageExporter):

    def _compute_indices(self):
        if self.export_cmd.phys_origin is not None:
            self.start_os_mm = self.export_cmd.phys_origin
            self.start_is_mm = self.ct.composite.TransformPoint(self.start_os_mm)
        else:
            self.start_is_mm = self.origin
            self.start_os_mm = self.ct.affine_composite.GetInverse().TransformPoint(self.start_is_mm)

        self.start_index = np.abs(np.around(self.affine_pv.TransformPoint(self.start_is_mm)).astype(np.int))
        self.chunk_size_in_voxels = self.shape
        self.end_index = self.start_index + self.chunk_size_in_voxels
        self.start_index = np.array(self.start_index)[self.transpose_tuple.inv_transpose]
        self.end_index = np.array(self.end_index)[self.transpose_tuple.inv_transpose]


class HDFChannel(object):
    def __init__(self, h5_file, channel_name, orientation='RAI'):

        self.channel_name = channel_name
        logger.debug("Channel name: {}".format(channel_name))
        self.orientation = orientation
        self.h5_file = h5_file
        self.image_path = PathUtil.get_lsfm_image_path(self.channel_name)
        if self.image_path not in self.h5_file:
            raise ValueError("Channel {} not found!".format(self.channel_name))

        self.pyramid_levels = []
        self.image_node = self.h5_file[self.image_path]
        self.p_origins = self._get_image_meta_data(M_ORIGIN)
        self.p_shapes = self._get_image_meta_data(M_SHAPE)
        self.p_voxel_sizes = self._get_image_meta_data(M_VOXEL_SIZE)
        self.n_levels = len(self.p_shapes)
        try:
            self.is_multichannel = bool(int(self._get_image_meta_data('is_multichannel')))
            self.is_segmentation = bool(int(self._get_image_meta_data('is_segmentation')))
            self.pixel_type = int(self._get_image_meta_data('pixel_type'))
        except KeyError:
            self.is_multichannel = False
            self.is_segmentation = False
            self.pixel_type = 3

        if self.is_segmentation:
            logger.debug("{} ".format(self.is_segmentation))
            logger.debug("{} is segmentation".format(self.channel_name))
        self._populate_pyramid_levels()

    def _get_image_data_path(self, _level=0):
        return PathUtil.get_lsfm_image_cells_path(self.channel_name, _level)

    def _get_image_meta_data(self, meta_name):
        return np.array(self.image_node.attrs[meta_name])

    def _populate_pyramid_levels(self):
        for i_level in range(self.n_levels):
            p_level = HDFChannel.PyramidLevel(self.h5_file)
            p_level.id = i_level
            p_level.pixel_type = self.pixel_type
            p_level.is_segmentation = self.is_segmentation
            p_level.is_multichannel = self.is_multichannel
            p_level.is_nifti = self.is_multichannel or self.is_segmentation #TODO
            p_level.path = self._get_image_data_path(i_level)
            p_level.origin = self.p_origins[i_level]
            p_level.shape = self.p_shapes[i_level]
            p_level.voxel_size = self.p_voxel_sizes[i_level]
            p_level.physical_size = p_level.shape[:len(p_level.voxel_size)] * p_level.voxel_size
            p_level.affine = nib.affines.from_matvec(np.diag(p_level.voxel_size),
                                                     p_level.origin)

            self.pyramid_levels.append(p_level)

    def log_operation(self, cmd_name, cmd_line):
        if self.h5_file.mode == 'r':
            logger.info('File in read mode only, cannot log operations!')
            return
        log_path = PathUtil.get_log_path(dt.datetime.now(), cmd_name)
        logger.debug('log path: {:} \n command-name: {:}\n command: {:}'.format(log_path,
                                                                                cmd_name,
                                                                                cmd_line))
        self.h5_file.create_dataset(name=log_path, data=cmd_line)

    @logging_decor
    def write_affine(self, affine_name, affine_file_path):
        # type: (str, str) -> None
        affine = Affine(affine_path=affine_file_path, affine_name=affine_name)
        affine_path = PathUtil.get_lsfm_affine_path(self.channel_name, affine_name)
        affine_itk_path = PathUtil.get_lsfm_affine_itk_path(self.channel_name, affine_name)

        try:
            self.h5_file[affine_path] = affine.nifti_affine
            self.h5_file[affine_itk_path[0]] = affine.itk_affine.GetParameters()
            self.h5_file[affine_itk_path[1]] = affine.itk_affine.GetFixedParameters()
            logger.debug("Affine {} successfully written".format(affine_name))
        except RuntimeError:
            logger.error("Affine with that name ({}) is already present!".format(affine_name))

    def read_affine(self, affine_name):
        """

        :type affine_name: str
        :rtype :Affine
        """
        affine = Affine(affine_name=affine_name)
        affine_path = PathUtil.get_lsfm_affine_path(self.channel_name, affine_name)
        affine_itk_path = PathUtil.get_lsfm_affine_itk_path(self.channel_name, affine_name)

        try:
            affine_nifti = np.array(self.h5_file[affine_path], dtype=np.float64)
            affine_itk_params = self.h5_file[affine_itk_path[0]]
            affine_itk_fixed_params = self.h5_file[affine_itk_path[1]]
        except RuntimeError:
            logger.error("Affine {} not found!".format(affine_name))
            raise
        except KeyError:
            logger.error("Affine {} not found!".format(affine_name))
            sys.exit()

        affine.set_affine(affine_nifti, affine_itk_params, affine_itk_fixed_params)
        return affine

    def read_displacement_field(self, df_name):
        try:
            df_channel = HDFChannel(self.h5_file, df_name)
        except ValueError:
            logger.error("Displacement Field {} not found!".format(df_name))
            sys.exit()

        df_level = df_channel.pyramid_levels[0]
        cmd = ExportCmd(None, None, None, 'RAS', 0, [], None, None, None, None, None, None)
        #df_meta_image = df_level.get_meta_image(cmd)
        #df = df_meta_image.get_sitk_image()
        df = df_level.get_meta_image(cmd)
        df = sitk.Cast(df, sitk.sitkVectorFloat64)
        # sitk.WriteImage(df, '../results/exp_4/inner_df.nii.gz')
        dft = sitk.DisplacementFieldTransform(df)
        return dft

    def read_segmentation(self, seg_name):
        try:
            seg_channel = HDFChannel(self.h5_file, seg_name)
        except ValueError:
            logger.error("Segmentation {} not found!".format(seg_name))
            sys.exit()

        seg_level = seg_channel.pyramid_levels[0]
        cmd = ExportCmd(None, None, None, 'RAS', 0, [], None, None, None, None, None, None)
        #seg_meta_image = seg_level.get_meta_image(cmd)
        #seg = seg_meta_image.get_sitk_image()
        seg = seg_level.get_meta_image(cmd)
        seg = sitk.Cast(seg, sitk.sitkInt16)
        return seg

    def read_ref_channel(self, ref_name):
        try:
            ref_channel = HDFChannel(self.h5_file, ref_name)
        except ValueError:
            logger.error("Reference channel {} not found!".format(ref_name))
            sys.exit()

        ref_channel_level = ref_channel.pyramid_levels[0]
        ref_channel_level.is_nifti = True
        cmd = ExportCmd(None, None, None, 'RAS', 0, [], None, None, None, None, None, None)
        ref_meta_image = ref_channel_level.get_meta_image(cmd)
        return ref_meta_image

    def export_slices(self, cmd):
        """

        :type cmd: ExportSlicesCmd
        """
        if not cmd.resample:
            return self.pyramid_levels[cmd.input_resolution_level].get_slices(cmd)

        else:
            raise NotImplementedError

    def process_transforms(self, list_of_transforms):
        transforms = []
        for transform_tuple in list_of_transforms:
            if transform_tuple.type == 'affine':
                affine = self.read_affine(transform_tuple.name)
                if transform_tuple.invert:
                    affine = affine.itk_inv_affine
                else:
                    affine = affine.itk_affine
                transform = Transform(transform_tuple.type, transform_tuple.order,
                                      transform_tuple.name, affine)
                transforms.append(transform)
            elif transform_tuple.type == 'df':
                df = self.read_displacement_field(transform_tuple.name)
                transform = Transform(transform_tuple.type, transform_tuple.order,
                                      transform_tuple.name, df)
                transforms.append(transform)
            else:
                raise NotImplementedError("Wrong transformation type")

        return transforms

    def export_grid_of_chunks(self, export_cmd, level):
        """

        :param export_cmd:
        :param level:
        :return:
        """
        def get_whole_image_grid():
            logger.debug("Exporting whole image as chunks")
            transpose, flip = ImageProcessor.get_transpose(export_cmd.input_orientation)
            whole_phys_size = self.pyramid_levels[level].physical_size[transpose]
            grid_size = np.ceil(whole_phys_size / (export_cmd.phys_size - export_cmd.overlap_mm))
            grid_size = grid_size.astype(np.int32)
            phys_origin = np.abs(np.array(self.pyramid_levels[level].origin)[transpose]) * [1., 1., -1.]
            logger.debug("transpose for grid size determination: {}".format(transpose))
            logger.debug("grid size: {}".format(grid_size))
            logger.debug('phys origin: {}'.format(phys_origin))
            return phys_origin, grid_size

        def get_physical_coordinates():
            rows, cols, sls = np.indices([n_x, n_y, n_z])
            grid_indices = np.vstack([rows.flatten(), cols.flatten(), sls.flatten()]).T
            phys_affine = np.diag((np.array(export_cmd.phys_size) - export_cmd.overlap_mm) * [1., 1., -1.])
            phys_coords = np.apply_along_axis(lambda x: np.dot(phys_affine, x), 1, grid_indices)
            phys_coords += phys_origin
            return grid_indices, phys_coords

        def get_chunk_cmd():
            export_chunk_cmd = copy.copy(export_cmd)
            export_chunk_cmd.segmentation_name = None
            export_chunk_cmd.region_id = None
            export_chunk_cmd.grid_size = None
            export_chunk_cmd.overlap_mm = None
            return export_chunk_cmd

        def export_chunk(i, origin):
            logger.debug("chunk with index {} at origin {}".format(grid_indices[i], origin))
            if export_cmd.output_path:
                export_chunk_cmd.output_path = export_cmd.output_path + "_{}_{}_{}.nii.gz".format(
                    grid_indices[i][0],
                    grid_indices[i][1],
                    grid_indices[i][2])
            else:
                export_chunk_cmd.output_path = None

            export_chunk_cmd.phys_origin = origin
            return self.pyramid_levels[level].get_meta_image(export_chunk_cmd)

        logger.debug("Exporting grid of chunks")
        if np.all(export_cmd.grid_size == [0, 0, 0]):
            phys_origin, grid_size = get_whole_image_grid()
            n_x, n_y, n_z = grid_size
            logger.debug("will export whole image as chunks")
        else:
            n_x, n_y, n_z = export_cmd.grid_size
            phys_origin = export_cmd.phys_origin

        grid_indices, phys_coords = get_physical_coordinates()
        export_chunk_cmd = get_chunk_cmd()

        if export_cmd.output_path:
            for i, origin in enumerate(phys_coords):
                export_chunk(i, origin)
            sys.exit()
        else:
            logger.debug("Returning bunch of objects")
            return (export_chunk(i, origin) for (i, origin) in enumerate(phys_coords))

    def export_image(self, export_cmd):
        """
        Export a subset of image volume based on ExportCommand object.

        :type export_cmd: ExportCmd
        """

        if export_cmd.list_of_transforms:
            export_cmd.list_of_transforms = self.process_transforms(export_cmd.list_of_transforms)

        if export_cmd.input_resolution_level is None:
            level = self.find_suitable_spacing_level(export_cmd.output_resolution)
        else:
            if export_cmd.input_resolution_level not in np.arange(self.n_levels):
                raise ValueError("Error: wrong resolution level value")
            level = export_cmd.input_resolution_level

        if export_cmd.export_region:
            segmentation = self.read_segmentation(export_cmd.segmentation_name)
            export_cmd.segmentation = segmentation

        if export_cmd.grid_of_chunks:
            return self.export_grid_of_chunks(export_cmd, level)

        return self.pyramid_levels[level].get_meta_image(export_cmd)  # .save_nifti(export_cmd.output_path)

    def find_suitable_spacing_level(self, resolution):
        if resolution is None:
            raise ValueError("You must provide resolution or pyramid level for export!")
        spacings = np.array(self.p_voxel_sizes, dtype=np.float)
        resolution = np.array(resolution, dtype=np.float)
        msk = map(np.all, spacings <= resolution)
        distances = map(np.linalg.norm, (spacings[msk] - resolution))
        if len(distances) == 0:
            return 0
        else:
            return np.argmin(distances)

    def __repr__(self):
        """

        :return: string representation of channel: affines and displacement fields written,
                 information concerning levels of the pyramid of resolutions
        """

        try:
            affines = self.h5_file[PathUtil.get_lsfm_affines_group_path(self.channel_name)].keys()
            affines = "\n".join(affines)
        except KeyError:
            affines = "No affines were found for this channel."

        try:
            displacement_fields = self.h5_file[PathUtil.get_lsfm_dfs_path(self.channel_name)].keys()
            displacement_fields = "\n".join(displacement_fields)
        except KeyError:
            displacement_fields = "No displacement fields were found for this channel"

        pyramid_levels = "\n".join(str(pl) for pl in self.pyramid_levels)

        r = "Channel:\n {} \n" \
            "\n\nAffines:\n {} \n" \
            "\n\nDisplacement Fields:\n {} \n" \
            "\n\nPyramid levels:\n {}".format(self.channel_name,
                                               affines,
                                               displacement_fields,
                                               pyramid_levels)

        return r

    class PyramidLevel(object):

        def __init__(self, h5_file):
            self.h5_file = h5_file
            self.data = None
            self.pixel_type = None
            self.is_segmentation = None
            self.is_nifti = None
            self.path = None
            self.origin = None
            self.shape = None
            self.voxel_size = None
            self.physical_size = None
            self.affine = None

        def get_level_data(self):
            if self.data is None:
                self.data = self.h5_file[self.path][...]
            return self.data

        @staticmethod
        def get_slice(start, end, flip=False, step=1):
            if flip:
                start = -start
                if start == 0:
                    start = None
                return slice(-end, start, step)
            else:
                return slice(start, end, step)

        @staticmethod
        def batch_slices(slc, num=15):
            num_slabs = (slc.stop - slc.start) / slc.step // num
            for i in range(num_slabs):
                yield slice(slc.start + i * num, slc.start + (i + 1) * num * slc.step, slc.step)
            yield slice(slc.start + num_slabs * num * slc.step, slc.stop, slc.step)

        def get_resampled_slices(self, cmd):
            """
            :type cmd: ExportSlicesCmd
            :param cmd:
            :return:
            """
            logger.debug("input orientation {}".format(cmd.ref_orientation))
            transpose, flip = ImageProcessor.get_transpose(cmd.ref_orientation)
            inv_transpose = np.argsort(transpose)
            logger.debug("transpose: {}\n flip: {}".format(transpose, flip))
            target_spacing = cmd.ref_channel.pyramid_levels[cmd.ref_level].voxel_size[transpose]
            logger.debug("target_spacing: {}".format(target_spacing))
            logger.debug("source spacing: {}".format(self.voxel_size))

            scale_factors = self.voxel_size / target_spacing
            logger.debug("scaling_factors: {}".format(scale_factors))

            target_size = cmd.ref_channel.pyramid_levels[cmd.ref_level].shape[transpose]
            logger.debug("target size: {}".format(target_size))
           #source_size = self.shape[transpose]
            logger.debug("source image size: {}".format(self.shape))

            cmd_tmp = ExportCmd(None, None, None, 'RAS', 0, [], None, None, None, None)
            seg_meta_image = self.get_meta_image(cmd_tmp)
            seg = seg_meta_image.get_sitk_image()
            seg = sitk.Cast(seg, sitk.sitkInt16)

            logger.debug("seg size: {}".format(seg.GetSize()))

            roi_filter = sitk.RegionOfInterestImageFilter()

            ## try to make this axis-agnostic
            if cmd.start is None:
                cmd.start = 0
            if cmd.stop is None:
                cmd.stop = self.shape[cmd.axis]

            start_index = [0, 0, 0]
            start_index[cmd.axis] = cmd.start

            slab_size = np.copy(self.shape)
            planes_to_resample = 10
            slab_size[cmd.axis] = planes_to_resample  # cmd.step
            print(float(self.shape[cmd.axis]))
            output_planes_per_input_planes = target_size[cmd.axis] / float(self.shape[cmd.axis])

            logger.debug("output planes per input planes: {}".format(output_planes_per_input_planes))
            roi_filter.SetSize(slab_size)
            plane_no = 0
            for slab_start in xrange(cmd.start, cmd.stop - cmd.step, cmd.step):
                start_index[cmd.axis] = slab_start
                logger.debug("start index: {}".format(start_index))
                logger.debug("slab size: {}".format(slab_size))
                roi_filter.SetIndex(start_index)
                roi_slab = roi_filter.Execute(seg)
                resampled_slab = ImageProcessor.resample_sitk(roi_slab, scale_factors, sitk.sitkNearestNeighbor)
                slab_data = sitk.GetArrayFromImage(resampled_slab)
                logger.debug("slab shape: {}".format(slab_data.shape))
                # some_plane = slab_data[:, 5, :]
                # imsave(cmd.output_path % slab_start + cmd.output_ext, some_plane)
                current_plane = start_index[cmd.axis] * output_planes_per_input_planes
                logger.debug("current plane: {}".format(current_plane))

                # chosen_planes = slab_data[:, ::cmd.step, :]
                chosen_planes = slab_data
                for i in xrange(chosen_planes.shape[cmd.axis]):
                    # print("---> "+ str(current_plane + i))
                    if round(current_plane + i) % cmd.step == 0:
                        # print(current_plane + i)

                        imsave(cmd.output_path % round(current_plane + i) + cmd.output_ext,
                               chosen_planes[:, i, :])

                        # print(cmd.output_path % plane_no + cmd.output_ext)
                        # print(current_plane + i * cmd.step)
                        # imsave(cmd.output_path % round(current_plane + i * cmd.step) + cmd.output_ext,
                        #       chosen_planes[:, i, :])
                        # imsave(cmd.output_path % plane_no + cmd.output_ext, chosen_planes[:, i, :])
                        # plane_no += 1

            ##



            '''index = [0, 0, 0]
            size = [180, 20, 239]

            roi_filter.SetIndex(index)
            roi_filter.SetSize(size)

            roi_slab = roi_filter.Execute(seg)

            resampled_slab = ImageProcessor.resample_sitk(roi_slab, scale_factors, sitk.sitkLinear)
            #sitk.WriteImage(resampled_slab, '../results/slab_smally.nii.gz')
            slab_data = sitk.GetArrayFromImage(resampled_slab)
            logger.debug("slab shape: {}".format(slab_data.shape))
            some_plane = slab_data[:, 15, :]
            imsave('../results/seg_plane_15.tif', some_plane)'''


            # compute scaling factors
            # get meta image for a slab
            # resample slab, get slices, output images

        def get_slices(self, cmd):
            """

            :type cmd: ExportSlicesCmd
            """
            if cmd.output_path:
                if not os.path.exists(os.path.dirname(cmd.output_path)):
                    os.makedirs(os.path.dirname(cmd.output_path))

            transpose, flip = ImageProcessor.get_transpose(cmd.input_orientation)
            logger.debug("transpose: {}\n flip: {}".format(transpose, flip))
            inv_transpose = np.argsort(transpose)
            logger.debug("inverse transpose: {}".format(inv_transpose))
            inv_flip = flip[inv_transpose] < 0
            axis = transpose[cmd.axis]
            axis_inverse = flip[cmd.axis] < 0
            logger.debug("inverse main axis: {}".format(axis_inverse))
            axis_h, axis_w = np.delete(transpose, cmd.axis, 0)

            max_w, max_h, max_z = self.shape[axis_w], self.shape[axis_h], self.shape[axis]
            logger.debug("max width: {} max height: {} max z: {}".format(max_w, max_h, max_z))

            if cmd.start is None:
                cmd.start = 0

            if cmd.stop is None:
                cmd.stop = max_z

            if cmd.roi_ox is None:
                cmd.roi_ox = 0

            if cmd.roi_oy is None:
                cmd.roi_oy = 0

            if cmd.roi_sx is None:
                cmd.roi_sx = max_w

            if cmd.roi_sy is None:
                cmd.roi_sy = max_h

            # sl_main = self.get_slice(cmd.start, cmd.stop, axis_inverse, step=cmd.step)
            sl_main = self.get_slice(cmd.start, cmd.stop, False, step=cmd.step)
            sl_w = self.get_slice(cmd.roi_ox, cmd.roi_ox + cmd.roi_sx, inv_flip[axis_w])
            sl_h = self.get_slice(cmd.roi_oy, cmd.roi_oy + cmd.roi_sy, inv_flip[axis_h])

            logger.debug("sl_main: {}".format(sl_main))

            logger.debug("voxel size h: {0:.4f}".format(self.voxel_size[axis_h]))
            logger.debug("voxel size w: {0:.4f}".format(self.voxel_size[axis_w]))

            plane_count = 0
            sl_dict = dict()

            sl_dict[axis_w] = sl_w
            sl_dict[axis_h] = sl_h

            list_of_batch_slices = list(self.batch_slices(sl_main))
            #if axis_inverse:
            #list_of_batch_slices = list_of_batch_slices[::-1]
            print(list_of_batch_slices)

            # for main_slice in self.batch_slices(sl_main):
            for main_slice in list_of_batch_slices:
                # main_slice = slice(main_slice.start * -1, main_slice.stop * -1, main_slice.step)
                logger.debug("main slice: {}".format(main_slice))
                sl_dict[axis] = main_slice
                logger.debug("sl_dict: {}".format(sl_dict))
                data = self.h5_file[self.path][sl_dict[0], sl_dict[1], sl_dict[2]]
                data = data.transpose([axis, axis_h, axis_w])
                if axis_inverse:
                    data = data[::-1]
                logger.debug("batch data shape: {}".format(data.shape))

                for plane in data:
                    if cmd.output_path:
                        imageio.imwrite(cmd.output_path % plane_count + cmd.output_ext, plane)
                    else:
                        yield plane
                    plane_count += 1

            return

        def get_level_chunk_data(self, start_idx, end_idx, flip):

            slx = self.get_slice(start_idx[0], end_idx[0], flip[0])
            sly = self.get_slice(start_idx[1], end_idx[1], flip[1])
            slz = self.get_slice(start_idx[2], end_idx[2], flip[2])

            logger.debug("slx: {}".format(slx))
            logger.debug("sly: {}".format(sly))
            logger.debug("slz: {}".format(slz))

            expected_shape = end_idx - start_idx

            '''expected_shape = np.array([slx.stop - slx.start,
                                       sly.stop - sly.start,
                                       slz.stop - slz.start])'''

            '''expected_shape = np.array([len(xrange(*slx.indices(max(self.shape[0], slx.stop)))),
                                       len(xrange(*sly.indices(max(self.shape[1], sly.stop)))),
                                       len(xrange(*slz.indices(max(self.shape[2], sly.stop))))])'''

            logger.debug("Pyramid level shape: {}".format(self.shape))
            logger.debug("Expected chunk shape from source: {}".format(expected_shape))

            data = self.h5_file[self.path][slx, sly, slz]
            logger.debug("actual shape of data: {}".format(data.shape))

            shape_difference = expected_shape - data.shape[:3]
            logger.debug("Difference in shape: {}".format(shape_difference))

            data_slx = self.get_slice(0, data.shape[0])
            data_sly = self.get_slice(0, data.shape[1])
            data_slz = self.get_slice(0, data.shape[2])

            if np.abs(slx.start) > self.shape[0]:
                data_slx = self.get_slice(0, data.shape[0], flip=True)
            if np.abs(sly.start) > self.shape[1]:
                data_sly = self.get_slice(0, data.shape[1], flip=True)
            if np.abs(slz.start) > self.shape[2]:
                data_slz = self.get_slice(0, data.shape[2], flip=True)

            n_dims = data.ndim - 3
            proper_shape = tuple(expected_shape) + (data.shape[-n_dims:]) if n_dims else tuple(expected_shape)
            logger.debug("proper chunk shape: {}".format(proper_shape))
            proper_chunk = np.zeros(proper_shape, dtype=data.dtype)

            # proper_chunk[:data.shape[0], :data.shape[1], :data.shape[2]] = data
            logger.debug("data_slx: {}".format(data_slx))
            logger.debug("data_sly: {}".format(data_sly))
            logger.debug("data_slz: {}".format(data_slz))
            proper_chunk[data_slx, data_sly, data_slz, ...] = data

            return proper_chunk

        def get_meta_image(self, export_cmd):

            if export_cmd.whole_image:
                img_exp = ImageWholeExporter(export_cmd, self)
            else:
                if export_cmd.export_region:
                    img_exp = ImageRegionExporter(export_cmd, self)
                else:
                    img_exp = ImageExporter(export_cmd, self)

            output_img = img_exp.process()

            if export_cmd.output_path:
                sitk.WriteImage(output_img, export_cmd.output_path)
                return True
            else:
                return output_img

        def old_get_meta_image(self, export_cmd):
            """

            :type export_cmd: ExportCmd
            """

            # we want just a chunk
            if not export_cmd.whole_image:
                # if no need for resampling:
                if not export_cmd.resample_image:
                    # get the transpose
                    transpose, flip = ImageProcessor.get_transpose(export_cmd.input_orientation)
                    inv_transpose = np.argsort(transpose)
                    inv_flip = flip[inv_transpose] < 0
                    # create a dummy image
                    d_shape = tuple(np.array(self.shape)[transpose])
                    # direction fix
                    d_origin = tuple(np.array(self.origin[transpose]) * np.array([1, -1, -1]))
                    d_voxel_size = tuple(np.array(self.voxel_size[transpose]))
                    logger.debug('size in voxel space: {}'.format(export_cmd.phys_size // d_voxel_size))
                    logger.debug("origin_transposed: {}".format(d_origin))
                    logger.debug("voxel_size: {}".format(d_voxel_size))

                    affine_vp, affine_pv = ImageProcessor.get_affines(d_voxel_size, d_origin)
                    start_os_mm = tuple(export_cmd.phys_origin)

                    chunk_size_in_voxels = (export_cmd.phys_size // d_voxel_size).astype(np.int)

                    # fix directionality (signs) for this
                    end_os_mm = tuple(np.array(export_cmd.phys_origin) +
                                      np.array(export_cmd.phys_size) * np.array([1, 1, -1]))

                    logger.debug("start: {} end: {}".format(start_os_mm, end_os_mm))

                    # if there are no transforms
                    if not export_cmd.list_of_transforms:
                        # start_index = dummy_img.TransformPhysicalPointToIndex(start_os_mm)
                        # logger.debug('start_index from dummy: {}'.format(start_index))
                        start_index = np.around(affine_pv.TransformPoint(start_os_mm)).astype(np.int)
                        logger.debug('start_index from affine: {}'.format(affine_pv.TransformPoint(start_os_mm)))

                        # start_os_mm = dummy_img.TransformIndexToPhysicalPoint(start_index)
                        start_os_mm = affine_vp.TransformPoint(start_index)
                        start_index = np.abs(start_index)
                        logger.debug("start_os_mm transformed back from index {}".format(start_os_mm))

                        # end index should be recomputed by voxel space
                        # end_index = dummy_img.TransformPhysicalPointToIndex(end_os_mm)
                        end_index = start_index + chunk_size_in_voxels
                        start_index = tuple(np.array(start_index)[inv_transpose])
                        end_index = tuple(np.array(end_index)[inv_transpose])
                        logger.debug("start: {} end: {}".format(start_index, end_index))

                        chunk = self.get_level_chunk_data(start_index,
                                                          end_index,
                                                          inv_flip)
                        logger.debug("chunk dimensions: {}".format(chunk.shape))
                        chunk = chunk.transpose(transpose)
                        chunk = ImageProcessor.reverse_axes(chunk, flip < 0)

                        img_data = sitk.GetImageFromArray(chunk.transpose())
                        img_data.SetSpacing(d_voxel_size)
                        img_data.SetOrigin(start_os_mm)
                        img_data.SetDirection(DIRECTION_RAS)
                        sitk.WriteImage(img_data, export_cmd.output_path)

                    # if there are transforms
                    else:

                        # dummy_img = ImageProcessor.create_dummy_image(d_shape, DIRECTION_RAS,
                        #                                              d_voxel_size, d_origin,
                        #                                              self.pixel_type)

                        logger.debug("level shape: {}".format(self.shape))
                        composite_transform = CompositeTransform(export_cmd.list_of_transforms)
                        ct = composite_transform.composite

                        # lets go back to input space
                        start_is_mm = ct.TransformPoint(start_os_mm)
                        start_index = np.abs(np.around(affine_pv.TransformPoint(start_is_mm)).astype(np.int))
                        logger.debug('start_index from affine: {}'.format(start_index))
                        end_index = start_index + chunk_size_in_voxels
                        logger.debug('end_index from chunk size: {}'.format(end_index))
                        # end_is_mm = ct.TransformPoint(end_os_mm) this was taking too small chunks
                        end_is_mm = tuple(start_is_mm + np.array(export_cmd.phys_size) * np.array([1, 1, -1]))
                        end_index_from_phys = np.abs(np.around(affine_pv.TransformPoint(end_is_mm)).astype(np.int))
                        logger.debug('END INDEX FROM PHYS: {}'.format(end_index_from_phys))
                        logger.debug("start_is: {} end_is: {}".format(start_is_mm, end_is_mm))

                        # grab the chunk
                        # start_index = dummy_img.TransformPhysicalPointToIndex(start_is_mm)
                        # end_index = dummy_img.TransformPhysicalPointToIndex(end_is_mm)
                        start_index = tuple(np.array(start_index)[inv_transpose])
                        end_index = tuple(np.array(end_index)[inv_transpose])
                        logger.debug("start: {} end: {}".format(start_index, end_index))

                        chunk = self.get_level_chunk_data(start_index,
                                                          end_index,
                                                          inv_flip)

                        # here add resampling to output space

                        # if we want trasformed data:
                        chunk = chunk.transpose(transpose)
                        chunk = ImageProcessor.reverse_axes(chunk, flip < 0)
                        chunk_affine = nib.affines.from_matvec(np.diag(d_voxel_size) * SIGN_RAS_A,
                                                               start_is_mm)
                        # chunk_affine[:3, 3] = start_is_mm
                        logger.debug("chunk shape: {}".format(chunk.shape))
                        logger.debug("chunk affine: {}".format(chunk_affine))
                        meta_image = MetaImage(chunk, chunk_affine,
                                               self.pixel_type, DIRECTION_RAS, self.is_segmentation)
                        # meta_image.reorient(export_cmd.input_orientation)
                        # meta_image.set_direction_RAS()

                        out = ImageProcessor.apply_transformations(meta_image,
                                                                   meta_image,
                                                                   composite_transform)
                        sitk.WriteImage(out, export_cmd.output_path)
                        # end trasformateion


                        # for untransformed data:
                        '''
                        chunk = chunk.transpose(transpose)
                        chunk = ImageProcessor.reverse_axes(chunk, flip < 0)
                        img_data = sitk.GetImageFromArray(chunk.transpose())
                        img_data.SetSpacing(d_voxel_size)
                        img_data.SetOrigin(start_is_mm)
                        img_data.SetDirection(DIRECTION_RAS)
                        sitk.WriteImage(img_data, export_cmd.output_path)'''

                else:  # if resampling needed

                    transpose, flip = ImageProcessor.get_transpose(export_cmd.input_orientation)
                    inv_transpose = np.argsort(transpose)
                    inv_flip = flip[inv_transpose] < 0

                    d_shape = tuple(np.array(self.shape)[transpose])
                    # direction fix
                    d_origin = tuple(np.array(self.origin[transpose]) * np.array([1, -1, -1]))
                    d_voxel_size = tuple(np.array(self.voxel_size[transpose]))
                    logger.debug('size in voxel space: {}'.format(export_cmd.phys_size // d_voxel_size))
                    logger.debug("origin_transposed: {}".format(d_origin))
                    logger.debug("voxel_size: {}".format(d_voxel_size))

                    affine_vp, affine_pv = ImageProcessor.get_affines(d_voxel_size, d_origin)
                    start_os_mm = tuple(export_cmd.phys_origin)

                    chunk_size_in_voxels = (export_cmd.phys_size // d_voxel_size).astype(np.int)

                    # fix directionality (signs) for this
                    end_os_mm = tuple(np.array(export_cmd.phys_origin) +
                                      np.array(export_cmd.phys_size) * np.array([1, 1, -1]))

                    logger.debug("start: {} end: {}".format(start_os_mm, end_os_mm))

                    scale_factors = d_voxel_size / export_cmd.output_resolution
                    logger.debug("scale_factors: {}".format(scale_factors))

                    logger.debug("level shape: {}".format(self.shape))
                    composite_transform = CompositeTransform(export_cmd.list_of_transforms)
                    ct = composite_transform.composite

                    # lets go back to input space
                    start_is_mm = ct.TransformPoint(start_os_mm)
                    start_index = np.abs(np.around(affine_pv.TransformPoint(start_is_mm)).astype(np.int))
                    logger.debug('start_index from affine: {}'.format(start_index))
                    end_index = start_index + chunk_size_in_voxels
                    logger.debug('end_index from chunk size: {}'.format(start_index))
                    # end_is_mm = ct.TransformPoint(end_os_mm) this was taking too small chunks
                    end_is_mm = tuple(start_is_mm + np.array(export_cmd.phys_size) * np.array([1, 1, -1]))
                    logger.debug("start_is: {} end_is: {}".format(start_is_mm, end_is_mm))

                    # grab the chunk
                    # start_index = dummy_img.TransformPhysicalPointToIndex(start_is_mm)
                    # end_index = dummy_img.TransformPhysicalPointToIndex(end_is_mm)
                    start_index = tuple(np.array(start_index)[inv_transpose])
                    end_index = tuple(np.array(end_index)[inv_transpose])
                    logger.debug("start: {} end: {}".format(start_index, end_index))

                    chunk = self.get_level_chunk_data(start_index,
                                                      end_index,
                                                      inv_flip)

                    # here add resampling to output space

                    # if we want trasformed data:
                    chunk = chunk.transpose(transpose)
                    chunk = ImageProcessor.reverse_axes(chunk, flip < 0)
                    chunk_affine = nib.affines.from_matvec(np.diag(d_voxel_size) * SIGN_RAS_A,
                                                           start_is_mm)
                    # chunk_affine[:3, 3] = start_is_mm
                    logger.debug("chunk shape: {}".format(chunk.shape))
                    logger.debug("chunk affine: {}".format(chunk_affine))
                    meta_image = MetaImage(chunk, chunk_affine,
                                           self.pixel_type, DIRECTION_RAS, self.is_segmentation, self.is_nifti)

                    meta_image = ImageProcessor.resample_image(meta_image, scale_factors)

                    # meta_image.reorient(export_cmd.input_orientation)
                    # meta_image.set_direction_RAS()

                    out = ImageProcessor.apply_transformations(meta_image,
                                                               meta_image,
                                                               composite_transform)
                    sitk.WriteImage(out, export_cmd.output_path)

            if export_cmd.export_region:

                if not export_cmd.resample_image:
                    # get boundig box for region in segmentation
                    start, end = ImageProcessor.get_bounding_box(export_cmd.segmentation,
                                                                 export_cmd.region_id)

                    start_mm = export_cmd.segmentation.TransformIndexToPhysicalPoint(start)
                    end_mm = export_cmd.segmentation.TransformIndexToPhysicalPoint(end)
                    # get transpose and flip:
                    transpose, flip = ImageProcessor.get_transpose(export_cmd.input_orientation)
                    inv_transpose = np.argsort(transpose)

                    # get dummy image of current level
                    d_shape = tuple(np.array(self.shape)[transpose])
                    d_origin = tuple(np.array(self.origin[transpose]))
                    d_voxel_size = tuple(np.array(self.voxel_size[transpose]))
                    dummy_img = ImageProcessor.create_dummy_image(d_shape, DIRECTION_RAS,
                                                                  d_voxel_size, d_origin,
                                                                  self.pixel_type)

                    start_index = dummy_img.TransformPhysicalPointToIndex(start_mm)
                    end_index = dummy_img.TransformPhysicalPointToIndex(end_mm)
                    logger.debug("start: {} end: {}".format(start_index, end_index))
                    inv_flip = flip[inv_transpose] < 0

                    start_index = tuple(np.array(start_index)[inv_transpose])
                    end_index = tuple(np.array(end_index)[inv_transpose])

                    chunk = self.get_level_chunk_data(start_index,
                                                      end_index,
                                                      inv_flip)
                    logger.debug("chunk dimensions: {}".format(chunk.shape))
                    chunk = chunk.transpose(transpose)
                    chunk = ImageProcessor.reverse_axes(chunk, flip < 0)

                    img_data = sitk.GetImageFromArray(chunk.transpose())
                    img_data.SetSpacing(d_voxel_size)
                    img_data.SetOrigin(start_mm)
                    img_data.SetDirection(DIRECTION_RAS)

                    if export_cmd.list_of_transforms:
                        composite_transform = CompositeTransform(export_cmd.list_of_transforms)
                        ct = composite_transform.composite

                        chunk_affine = nib.affines.from_matvec(np.diag(d_voxel_size) * SIGN_RAS_A,
                                                               start_mm)
                        # chunk_affine[:3, 3] = start_is_mm
                        logger.debug("chunk shape: {}".format(chunk.shape))
                        logger.debug("chunk affine: {}".format(chunk_affine))
                        meta_image = MetaImage(chunk, chunk_affine,
                                               self.pixel_type, DIRECTION_RAS, self.is_segmentation, self.is_nifti)

                        img_data = ImageProcessor.apply_transformations(meta_image,
                                                                        meta_image,
                                                                        composite_transform)

                    # now extract the region # TODO parametrize to specify if chunk or actual region preferred

                    temp_seg = sitk.ReadImage('/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4/segmentation_left_hemisphere_mm.nii.gz')

                    region = ImageProcessor.extract_labeled_region(export_cmd.region_id,
                                                                   img_data,
                                                                   temp_seg)



                    #sitk.WriteImage(img_data, export_cmd.output_path)
                    sitk.WriteImage(region, export_cmd.output_path)
                    return

            # if export_cmd indicates whole image
            if export_cmd.whole_image:
                # if export_cmd indicates no resampling needed
                if not export_cmd.resample_image:
                    # get the level data, reorient, and save

                    meta_image = MetaImage(self.get_level_data(), self.affine,
                                           self.pixel_type, DIRECTION_LPI, self.is_segmentation, self.is_nifti)
                    meta_image.reorient(export_cmd.input_orientation)
                    meta_image.set_direction_RAS()
                    if export_cmd.output_path is not None:
                        meta_image.save_nifti(export_cmd.output_path)
                    else:
                        return meta_image
                else:
                    # if resampling needed
                    # get the level data, resample, reorient and save
                    meta_image = MetaImage(self.get_level_data(), self.affine,
                                           self.pixel_type, DIRECTION_LPI, self.is_segmentation, self.is_nifti)
                    meta_image.reorient(export_cmd.input_orientation)
                    scale_factors = meta_image.voxel_size / export_cmd.output_resolution
                    meta_image = ImageProcessor.resample_image(meta_image, scale_factors)
                    meta_image.set_direction_RAS()
                    if export_cmd.output_path is not None:
                        if not export_cmd.list_of_transforms:
                            if export_cmd.export_region:
                                logger.debug("Exporting region")
                                return meta_image.get_sitk_image()
                            else:
                                meta_image.save_nifti(export_cmd.output_path)
                    else:
                        return meta_image

                if export_cmd.list_of_transforms:
                    composite_transform = CompositeTransform(export_cmd.list_of_transforms)

                    # check the expected shape of input data for transformations
                    transpose, flip = ImageProcessor.get_transpose(export_cmd.input_orientation)
                    inv_transpose = np.argsort(transpose)
                    inv_flip = flip[inv_transpose] < 0

                    d_shape = tuple(np.array(self.shape)[transpose])
                    # direction fix
                    d_origin = tuple(np.array(self.origin[transpose]) * np.array([1, -1, -1]))
                    d_voxel_size = tuple(np.array(self.voxel_size[transpose]))
                    d_phys_size = np.array(self.physical_size[transpose])
                    logger.debug('size in voxel space: {}'.format(self.shape))
                    logger.debug('physical size: {}'.format(d_phys_size))
                    logger.debug("origin_transposed: {}".format(d_origin))
                    logger.debug("voxel_size: {}".format(d_voxel_size))

                    affine_vp, affine_pv = ImageProcessor.get_affines(d_voxel_size, d_origin)

                    ct = composite_transform.affine_composite
                    ct_inv = ct.GetInverse()
                    start_os_mm = ct_inv.TransformPoint(d_origin)
                    end_is_mm = d_origin + d_phys_size
                    end_os_mm = ct_inv.TransformPoint(end_is_mm)

                    logger.debug('Transposed origin (mm): {} \n Origin in output space: {} \n End in input space: {}\n'
                                 'End in output space: {}'.format(d_origin, start_os_mm, end_is_mm, end_os_mm))

                    #start_index = np.abs(np.around(affine_pv.TransformPoint(start_is_mm)).astype(np.int))
                    #logger.debug('start_index from affine: {}'.format(start_index))
                    # end check

                    ref_meta_image = export_cmd.ref_channel
                    # ref_meta_image.reorient('RAS')
                    # ref_meta_image.set_direction_RAS()

                    out = ImageProcessor.apply_transformations_with_reference(meta_image,
                                                                              ref_meta_image,
                                                                              composite_transform)
                    sitk.WriteImage(out, export_cmd.output_path)



                    # return MetaImage(self.get_level_data(), self.affine,
                    #                 self.pixel_type, DIRECTION_LPI, self.is_segmentation)

        def __repr__(self):
            return "\nLevel : {}" \
                   "\nLevel shape:      {: <7} {: <7} {: <7}" \
                   "\nLevel spacing:    {:.5f} {:.5f} {:.5f}" \
                   "\nLevel origin:     {:+.4f} {:+.4f} {:+.4f}" \
                   "\nPhysical size [mm]: {:.5f} {:.5f} {:.5f}\n".format(self.id,
                                                                         self.shape[0],
                                                                         self.shape[1],
                                                                         self.shape[2],
                                                                         self.voxel_size[0],
                                                                         self.voxel_size[1],
                                                                         self.voxel_size[2],
                                                                         self.origin[0],
                                                                         self.origin[1],
                                                                         self.origin[2],
                                                                         self.physical_size[0],
                                                                         self.physical_size[1],
                                                                         self.physical_size[2])


class ExportSlicesCmd(object):
    def __init__(self, channel_name, input_orientation, input_resolution_level, axis,
                 output_path=None, extract_roi=None, slicing_range=None,
                 ref_channel=None, ref_level=None, ref_orientation=None):

        if slicing_range is None:
            slicing_range = [None, None, 1]
        if extract_roi is None:
            extract_roi = [None, None, None, None]

        self.axis = axis
        self.channel_name = channel_name
        self.start, self.stop, self.step = np.array(slicing_range)
        self.output_path = output_path
        if self.output_path:
            self.output_path, self.output_ext = os.path.splitext(output_path)
        self.roi_ox, self.roi_oy, self.roi_sx, self.roi_sy = np.array(extract_roi)
        self.ref_channel = ref_channel
        self.ref_level = ref_level
        self.ref_orientation = ref_orientation
        self.input_orientation = input_orientation
        self.input_resolution_level = input_resolution_level
        if self.ref_channel is not None:
            self.input_orientation = 'RAS'

    @property
    def resample(self):
        if self.ref_channel is not None:
            return True


class ExportCmd(object):
    def __init__(self, channel_name, output_path, output_resolution, input_orientation,
                 input_resolution_level, list_of_transforms, phys_origin, phys_size,
                 segmentation_name, region_id, grid_size, overlap_mm):
        # type: (str, str, list, str, int, list, list, list) -> None

        self.channel_name = channel_name
        self.output_path = output_path
        if self.output_path is not None:
            if not os.path.exists(os.path.dirname(self.output_path)):
                os.makedirs(os.path.dirname(self.output_path))

        self.output_resolution = output_resolution
        if output_resolution:
            if len(output_resolution) != 3:
                raise ValueError("Wrong shape of output resolution array")
            self.output_resolution = np.array(output_resolution, np.float64)
        self.input_orientation = input_orientation
        self.input_resolution_level = input_resolution_level
        self.list_of_transforms = list_of_transforms
        self.phys_origin = phys_origin
        if phys_origin is not None:
            self.phys_origin = np.array(phys_origin, np.float64)
        self.phys_size = phys_size
        if phys_size is not None:
            self.phys_size = np.array(phys_size, np.float64)

        self.segmentation_name = segmentation_name
        self.region_id = region_id
        self.segmentation = None

        self.grid_size = grid_size
        if grid_size is not None:
            self.grid_size = np.array(grid_size, dtype=np.int32)

        self.overlap_mm = overlap_mm

    @property
    def grid_of_chunks(self):
        if self.grid_size is not None and self.overlap_mm is not None:
            return True
        else:
            return False

    @property
    def whole_image(self):
        # if self.phys_origin is None and self.phys_size is None and not self.export_region:
        if self.phys_size is None and not self.export_region:
            return True
        else:
            return False

    @property
    def resample_image(self):
        if self.input_resolution_level is None:
            if self.output_resolution is not None:
                return True
        else:
            return False

    @property
    def export_region(self):
        if self.segmentation_name is not None and self.region_id is not None:
            return True
        else:
            return False


class Transform(object):
    def __init__(self, transform_type, order, transform_name, transform):
        self.order = order
        self.transform_type = transform_type
        self.transform = transform
        self.transform_name = transform_name


class CompositeTransform(object):
    def __init__(self, list_of_transforms):
        self.transforms = list_of_transforms

    def add_transform(self, transform):
        self.transforms.append(transform)

    def sort_transforms(self):
        self.transforms.sort(key=lambda x: x.order)

    @property
    def composite(self):
        self.sort_transforms()
        ct = sitk.Transform(self.transforms[0].transform)
        for transform in self.transforms[1:]:
            ct.AddTransform(transform.transform)

        return ct

    @property
    def affine_composite(self):
        act = None
        self.sort_transforms()
        for transform in self.transforms:
            if transform.transform_type == 'affine':
                if act is None:
                    act = sitk.Transform(transform.transform)
                else:
                    act.AddTransform(transform.transform)

        return act


class ProxyChannel(object):
    def __init__(self, image_proxy, num_bdv_setup):
        """

        :type image_proxy: ImageProxy
        """
        self.image = image_proxy
        self.num_bdv_setup = num_bdv_setup
        self.image_path = PathUtil.get_lsfm_image_path(self.image.channel_name)
        self.resolutions_path = PathUtil.get_res_path(self.num_bdv_setup)
        self.subdivisions_path = PathUtil.get_sub_path(self.num_bdv_setup)
        self.setup_path = PathUtil.get_setup_path(self.num_bdv_setup)

        self.resolutions, self.subdivisions = BigDataViewer.propose_mipmaps(self.image.voxel_size,
                                                                            self.image.image_shape)

        if self.image.is_multichannel:
            additional_dimensions = len(self.image.data_shape) - 3
            self.resolutions = np.hstack([self.resolutions, np.ones((self.resolutions.shape[0],
                                                                     additional_dimensions))])
            self.subdivisions = np.hstack([self.subdivisions, np.ones((self.subdivisions.shape[0],
                                                                       additional_dimensions))])

        logger.debug("resolutions:\n {}".format(self.resolutions))
        logger.debug("subdivisions:\n {}".format(self.subdivisions))
        logger.debug("resolution scales: \n {}".format(self.resolution_scales))

        self.levels = []
        self._compute_pyramid_levels()

        for pl in self.levels:
            logger.debug(pl)

    @property
    def num_pyramid_levels(self):
        return self.resolutions.shape[0]

    @property
    def resolution_scales(self):
        if self.num_pyramid_levels > 1:
            return self.resolutions[0] / self.resolutions
        else:
            return None

    @property
    def relative_resolution_scales(self):
        if self.num_pyramid_levels > 1:
            rrs = np.vstack([np.ones(self.resolutions.shape[1]),
                             np.divide(self.resolutions[:-1],
                                       self.resolutions[1:])])
            '''if self.image.is_multichannel:
                additional_dimensions = len(self.image.data_shape) - 3
                return np.hstack([rrs, np.ones((rrs.shape[0], additional_dimensions))])'''
            return rrs
        else:
            return None

    def _compute_pyramid_levels(self):

        def compute_level_shape(num_level, atr='stack_shape'):
            pyramid_level = self.levels[num_level]
            shape = np.array(self.levels[num_level - 1].__getattribute__(atr) *
                             self.relative_resolution_scales[num_level], dtype=np.uint16)

            shape[0] = int(round(self.levels[num_level - 1].__getattribute__(atr)[0] *
                                 self.relative_resolution_scales[num_level][0]))

            overhead_shape = None
            stack_shape = shape.copy()

            if self.image.is_stream:
                overhead_shape = np.uint16(self.levels[0].overhead_stack_shape *
                                           self.relative_resolution_scales[level])

                overhead_shape[0] = int(round(self.levels[0].overhead_stack_shape[0] *
                                              self.relative_resolution_scales[level][0]))

                shape[0] = shape[0] * self.image.num_of_stacks
                shape[0] += overhead_shape[0]

            logger.debug("Level: {}\n"
                         "SHAPE: {}\n"
                         "STACK SHAPE: {}\n"
                         "OVERHEAD STACKSHAPE: {}\n"
                         "N of stacks: {}".format(num_level, shape, stack_shape, overhead_shape,
                                                  self.image.num_of_stacks))

            pyramid_level.shape = shape
            pyramid_level.stack_shape = stack_shape
            pyramid_level.overhead_stack_shape = overhead_shape

        pyramid_level = ProxyChannel.PyramidLevel()
        pyramid_level.id = 0
        pyramid_level.shape = self.image.data_shape
        pyramid_level.spacing = self.image.voxel_size
        pyramid_level.stack_shape = self.image.stack_shape
        pyramid_level.overhead_stack_shape = self.image.overhead_stack_shape
        pyramid_level.path = PathUtil.get_lsfm_image_cells_path(self.image.channel_name, 0)
        pyramid_level.bdv_path = PathUtil.get_cells_path(0, self.num_bdv_setup, 0)
        pyramid_level.chunks = tuple(self.subdivisions[0])
        self.levels.append(pyramid_level)

        if self.num_pyramid_levels == 1:
            logger.debug("Image won't be resampled (resolution too low)")
            return

        shape_attr = 'shape'
        if self.image.is_stream:
            shape_attr = 'stack_shape'

        for level in xrange(1, self.num_pyramid_levels):
            pyramid_level = ProxyChannel.PyramidLevel()
            pyramid_level.path = PathUtil.get_lsfm_image_cells_path(self.image.channel_name,
                                                                    level)
            pyramid_level.bdv_path = PathUtil.get_cells_path(0, self.num_bdv_setup, level)
            pyramid_level.id = level
            pyramid_level.spacing = self.image.voxel_size * 1. / \
                                    self.resolution_scales[level][:len(self.image.voxel_size)]

            pyramid_level.chunks = tuple(self.subdivisions[level])
            self.levels.append(pyramid_level)
            compute_level_shape(level, shape_attr)

    class PyramidLevel(object):
        def __repr__(self):
            return "stack shape: {} \n" \
                   "overhead stack shape: {}\n" \
                   "level shape: {}\n" \
                   "level spacing: {}\n".format(self.stack_shape,
                                                self.overhead_stack_shape,
                                                self.shape,
                                                self.spacing)


class SitkDataType(object):
    def __init__(self, pixel_type, is_segmentation=False):
        self.is_segmentation = is_segmentation
        self.pixel_type = pixel_type

        if self.pixel_type in [sitk.sitkVectorUInt8, sitk.sitkVectorInt8, sitk.sitkVectorUInt16,
                               sitk.sitkVectorInt16, sitk.sitkVectorUInt32, sitk.sitkVectorInt32,
                               sitk.sitkVectorUInt64, sitk.sitkVectorInt64, sitk.sitkVectorFloat32,
                               sitk.sitkVectorFloat64]:
            self.component_type = 'vector'
            self.interpolation_type = sitk.sitkLinear
            self.interpolation_pixel = sitk.sitkVectorFloat32

        elif self.pixel_type in [sitk.sitkUInt8, sitk.sitkInt8, sitk.sitkUInt16,
                                 sitk.sitkInt16, sitk.sitkUInt32, sitk.sitkInt32,
                                 sitk.sitkUInt64, sitk.sitkInt64, sitk.sitkFloat32,
                                 sitk.sitkFloat64]:

            if not self.is_segmentation:
                self.component_type = 'scalar'
                self.interpolation_type = sitk.sitkLinear
                self.interpolation_pixel = sitk.sitkFloat32
            else:
                self.component_type = 'scalar'
                self.interpolation_type = sitk.sitkNearestNeighbor
                self.interpolation_pixel = self.pixel_type

        else:
            raise ValueError("Unknown pixel type: {}".format(self.pixel_type))

    def __repr__(self):
        return "Pixel type: {}\n" \
               "Component type: {}\n" \
               "Interpolation type: {}\n" \
               "Interpolation pixel type: {}\n" \
               "Is segmentation: {}".format(self.pixel_type,
                                            self.component_type,
                                            self.interpolation_type,
                                            self.interpolation_pixel,
                                            self.is_segmentation)


class MetaImage(object):
    def __init__(self,
                 data,
                 affine,
                 pixel_type,
                 direction,
                 is_segmentation=False,
                 is_nifti=False):

        if len(data.shape) > 4:
            data = data[:, :, :, 0, :]

        self.data = data
        self.data_shape = data.shape
        self.image_shape = self.data_shape[:3]
        self.is_segmentation = is_segmentation
        self.is_nifti = is_nifti

        self.direction = direction
        self.voxel_size = np.array(nib.affines.voxel_sizes(affine), dtype=np.float64)
        self.affine = affine
        self.origin = self.affine[:3, 3]
        self.pixel_type = pixel_type
        self.sitk_data = SitkDataType(self.pixel_type, self.is_segmentation)
        self.is_segmentation = is_segmentation
        self.transpose = np.arange(len(self.data_shape))
        self.transpose[:3] = [2, 1, 0]

        logger.debug("Origin when creating meta image: {}".format(self.origin))

    def get_sitk_image(self):
        img_data = self.data.copy()
        transpose_data = [2, 1, 0]
        if self.sitk_data.component_type == 'vector':
            if len(img_data.shape) == 4:
                # img_data = img_data[:, :, :, 0, :]
                transpose_data = np.arange(len(img_data.shape))
                transpose_data[:3] = [2, 1, 0]

        img_data = img_data.transpose(self.transpose)
        img_data = sitk.GetImageFromArray(img_data)
        img_data.SetSpacing(self.voxel_size)
        img_data.SetOrigin(self.origin)
        img_data.SetDirection(self.direction)

        logger.debug("Set spacing {} \n origin {}\n and direction {}".format(self.voxel_size,
                                                                             self.origin,
                                                                             self.direction))

        return img_data

    def save_nifti(self, output_path):
        logger.debug("Origin when saving: {}".format(self.origin))
        logger.debug("affine when saving: {}".format(self.affine))
        nifti_img = nib.Nifti1Image(self.data, self.affine)
        nifti_img.header['xyzt_units'] = 2
        nifti_img.header['sform_code'] = 2
        nifti_img.header['qform_code'] = 0
        nifti_img.header['pixdim'][4:] = 0
        nib.save(nifti_img, output_path)

    def reorient(self, source_orientation):
        if source_orientation == 'RAS' or self.is_nifti:
            return
        transpose, flip = ImageProcessor.get_transpose(source_orientation)
        logger.debug("transpose: {}\n flip: {}".format(transpose, flip))
        _transpose = np.arange(len(self.data_shape))
        _transpose[:3] = transpose
        self.data = self.data.transpose(_transpose)
        self.reverse_axes(flip < 0)

        self.voxel_size = self.voxel_size[transpose]
        self.origin = np.abs(self.origin[transpose]) * SIGN_LPI_T
        affine = np.abs(np.diag(self.voxel_size)) * SIGN_LPI_A
        self.affine = nib.affines.from_matvec(affine, self.origin)
        logger.debug("Origin after reorienting: {}".format(self.origin))

    def reverse_axes(self, reverse=[False, False, False]):
        if reverse[0]:
            self.data = self.data[::-1, :, :]
        if reverse[1]:
            self.data = self.data[:, ::-1, :]
        if reverse[2]:
            self.data = self.data[:, :, ::-1]

    def set_direction_LPI(self):
        if self.is_nifti:
            return
        self.direction = DIRECTION_LPI
        self.affine[:3, :3] = np.abs(self.affine[:3, :3]) * SIGN_LPI_A
        self.affine[:3, 3] = np.abs(self.affine[:3, 3]) * SIGN_LPI_T

    def set_direction_RAS(self):

        self.direction = DIRECTION_RAS
        self.affine[:3, :3] = np.abs(self.affine[:3, :3]) * SIGN_RAS_A
        if not self.is_nifti:
            self.affine[:3, 3] = np.abs(self.affine[:3, 3]) * SIGN_RAS_T
            self.origin = self.affine[:3, 3]
        logger.debug("Origin after setting RAS direction: {}".format(self.origin))

    @staticmethod
    def meta_image_from_channel(data, _channel):
        """

        :type data: numpy.ndarray
        :type _channel: ProxyChannel
        """
        return MetaImage(data, _channel.image.affine,
                         _channel.image.pixel_type, _channel.image.direction,
                         _channel.image.is_segmentation)


class MetaMultiChannelImage(MetaImage):
    def __init__(self, *args, **kwargs):
        super(MetaMultiChannelImage, self).__init__(*args, **kwargs)


class ImageProcessor(object):

    @staticmethod
    def reverse_axes(data, reverse=[False, False, False]):
        if reverse[0]:
            data = data[::-1, :, :]
        if reverse[1]:
            data = data[:, ::-1, :]
        if reverse[2]:
            data = data[:, :, ::-1]
        return data

    @staticmethod
    def get_bounding_box(segmentation, reg_id):
        label_stats_filter = sitk.LabelStatisticsImageFilter()
        label_stats_filter.Execute(segmentation, segmentation)
        x1, x2, y1, y2, z1, z2 = label_stats_filter.GetBoundingBox(reg_id)
        start_point = (x1, y1, z1)
        end_point = (x2, y2, z2)

        return start_point, end_point

    @staticmethod
    def create_dummy_image(size, direction, spacing, origin, pixel_id):
        logger.debug("size: {}".format(size))
        img = sitk.Image(int(size[0]), int(size[1]), int(size[2]), int(pixel_id))
        img.SetDirection(direction)
        img.SetSpacing(spacing)
        img.SetOrigin(origin)

        return img

    @staticmethod
    def trim_image_to_label_region(img, labels, reg_id, margin=0):
        """
        Trims intensity image to bounding box determined by label (region_id)
        with given margin value.
        :param img:
        :param labels:
        :param reg_id:
        :param margin: how many voxels of margin should be added
        :return:
        """

        label_stats_filter = sitk.LabelStatisticsImageFilter()
        roi_filter = sitk.RegionOfInterestImageFilter()

        label_stats_filter.Execute(img, labels)
        x1, x2, y1, y2, z1, z2 = label_stats_filter.GetBoundingBox(reg_id)
        index = [x1, y1, z1]
        size = [x2 - x1, y2 - y1, z2 - z1]

        # check if margin within bounds of img
        img_size = np.array(img.GetSize())
        m_index = np.array(index) - margin
        m_size = np.array(size) + (2 * margin)

        if np.any(m_index < 0):
            print "margin for segmentation reduced at index, YOU ARE ON THE EDGE"
            m_index[np.argwhere(m_index < 0)] = 0
            reduce_size = np.array(index) - m_index
            m_size -= reduce_size

        if np.any((m_index + m_size) > img_size):
            print "margin for segmentation reduced at upper bound, YOU ARE ON THE EDGE"
            out_of_bounds = np.argwhere((m_index + m_size) > img_size)
            m_size[out_of_bounds] = img_size[out_of_bounds]

        index = list(m_index)
        size = list(m_size)

        roi_filter.SetIndex(index)
        roi_filter.SetSize(size)
        out_img = roi_filter.Execute(img)

        return out_img

    @staticmethod
    def extract_labeled_region(reg_id, img, labels, margin=0):
        """
        Extracts region of interest from img identified by region_id, which corresponds to
        label intensity value in labels'
        :param reg_id: unit16 indicating designated region
        :param img: sitk.Image
        :param labels: sitk.Image of segmentation
        :return: sitk.Image
        """

        labels = sitk.Cast(labels, sitk.sitkInt16)
        threshold_filter = sitk.ThresholdImageFilter()
        rescale_filter = sitk.RescaleIntensityImageFilter()
        resample_filter = sitk.ResampleImageFilter()
        multiply_filter = sitk.MultiplyImageFilter()

        resample_filter.SetReferenceImage(img)
        resample_filter.SetInterpolator(sitk.sitkNearestNeighbor)
        labels = resample_filter.Execute(labels)

        # trim image to bounding box around selected region in segmentation
        # TODO: furher optimize this to send a chunk of image to this function
        img = ImageProcessor.trim_image_to_label_region(img, labels, reg_id, margin=margin)
        # trim segmentation to some reasonable bounding box
        labels = ImageProcessor.trim_image_to_label_region(labels, labels, reg_id, margin=margin)

        threshold_filter.SetLower(reg_id)
        threshold_filter.SetUpper(reg_id)
        labels = threshold_filter.Execute(labels)

        rescale_filter.SetOutputMaximum(1)
        labels = rescale_filter.Execute(labels)
        labels = sitk.Cast(labels, img.GetPixelID())

        region = multiply_filter.Execute(img, labels)

        return region

    @staticmethod
    def compute_offset(sitk_image, spacing):
        offset = np.array(sitk_image.GetDirection()).reshape(3, 3) * spacing

        for i in range(offset.shape[0]):
            offset[i] *= 0.5

        return offset

    @staticmethod
    def apply_transformations(meta_image, meta_ref_image, composite_transform):

        #input_image = sitk.ReadImage('../results/exp_4/001_10_mm.nii.gz')
        input_image = meta_image.get_sitk_image()
        #sitk.WriteImage(input_image, '../results/exp_4/chunk_to_transform.nii.gz')
        #output_origin = composite_transform.composite.TransformPoint(input_image.GetOrigin())
        output_origin = composite_transform.affine_composite.GetInverse().TransformPoint(input_image.GetOrigin())

        #output_origin = composite_transform.composite.GetInverse().TransformPoint(input_image.GetOrigin())
        #output_origin = np.array(output_origin) * np.array([-1, -1, 1])

        output_size = np.array(input_image.GetSize()) #+ 2 * np.abs(np.round(output_origin / np.array(input_image.GetSpacing())))

        logger.debug("output_size: {}".format(output_size))
        logger.debug("input spacing: {}".format(input_image.GetSpacing()))

        #
        #ref_image = sitk.ReadImage('../results/exp_4/template_both_hemispheres_mm.nii.gz')

        resampler = sitk.ResampleImageFilter()
        resampler.SetInterpolator(sitk.sitkLinear)
        resampler.SetDefaultPixelValue(0)

        resampler.SetOutputOrigin(output_origin)
        resampler.SetOutputSpacing(input_image.GetSpacing())
        resampler.SetSize(map(int, output_size))
        resampler.SetOutputDirection(input_image.GetDirection())


        #resampler.SetReferenceImage(ref_image)
        resampler.SetTransform(composite_transform.composite)

        out = resampler.Execute(input_image)
        #sitk.WriteImage(out, '../results/exp_4/test_from_inner_data.nii.gz')
        return out

    @staticmethod
    def old_apply_transformations(meta_image, meta_ref_image, composite_transform):

        """
        :type meta_image: MetaImage
        :type meta_ref_image: MetaImage
        :type composite_transform: CompositeTransform
        """

        input_data = meta_image.get_sitk_image()
        logger.debug("Input origin: {}".format(input_data.GetOrigin()))

        affine_composite = composite_transform.affine_composite
        if affine_composite is not None:
            out_direction = tuple(Affine.get_direction_matrix_from_itk_affine(affine_composite.GetInverse()).flatten())
            logger.debug("Output direction: {}".format(out_direction))
        else:
            out_origin = input_data.GetOrigin()
            out_direction = input_data.GetDirection()

        #ref_data = meta_ref_image.get_sitk_image()
        ref_data = sitk.ReadImage('/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/04_new_avrgt_25_right_hemisphere.nii')
        # ref_data = sitk.ReadImage('/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/_cfos_25um.nii.gz')
        #ref_data = sitk.ReadImage('/home/sbednarek/DEV/lsfm/results/cfos/cfos_autofluo25um.nii.gz')

        scale_factors = np.array(ref_data.GetSpacing()) / np.array(input_data.GetSpacing())
        logger.debug("Scale factors from ref image: {}".format(scale_factors))

        output_size = np.int16(np.floor(np.array(ref_data.GetSize()) * scale_factors))

        logger.debug("Size from ref image: {}".format(output_size))

        resampler = sitk.ResampleImageFilter()
        resampler.SetInterpolator(meta_image.sitk_data.interpolation_type)
        resampler.SetDefaultPixelValue(0)
        resampler.SetTransform(composite_transform.composite)
        resampler.SetSize((map(int, output_size)))
        resampler.SetOutputDirection(input_data.GetDirection())
        resampler.SetOutputSpacing(input_data.GetSpacing())
        resampler.SetOutputOrigin(ref_data.GetOrigin())

        #resampler.SetReferenceImage(ref_data)

        output_data = resampler.Execute(input_data)
        logger.debug("Output origin: {}".format(output_data.GetOrigin()))

        return output_data

    @staticmethod
    def resample_sitk(sitk_image, scale_factors, interpolation_type):

        scale_factors = np.array(scale_factors, dtype=np.float64)

        input_size = np.array(sitk_image.GetSize(), dtype=np.float64)
        output_size = np.int16(np.floor(input_size * scale_factors))
        output_size[0] = int(round(input_size[0] * scale_factors[0]))

        voxel_size = np.array(sitk_image.GetSpacing())

        output_spacing = voxel_size * input_size
        output_spacing = output_spacing / output_size

        input_offset = ImageProcessor.compute_offset(sitk_image, voxel_size)
        output_offset = ImageProcessor.compute_offset(sitk_image, output_spacing)

        output_origin = np.array(sitk_image.GetOrigin()) - input_offset + output_offset
        output_origin = np.diag(output_origin)

        identity_transform = sitk.AffineTransform(3)

        resampler = sitk.ResampleImageFilter()
        resampler.SetInterpolator(interpolation_type)
        resampler.SetDefaultPixelValue(0)
        resampler.SetTransform(identity_transform)
        resampler.SetSize((map(int, output_size)))
        resampler.SetOutputOrigin(tuple(output_origin))
        resampler.SetOutputSpacing(tuple(output_spacing))
        resampler.SetOutputDirection(sitk_image.GetDirection())

        output_image = resampler.Execute(sitk_image)

        return output_image

    @staticmethod
    def resample_image(meta_image, scale_factors, debug=False):

        """
        :type debug: bool
        :type scale_factors: numpy.ndarray
        :type meta_image: MetaImage
        """

        scale_factors = np.array(scale_factors, dtype=np.float64)
        sitk_image = meta_image.get_sitk_image()

        logger.debug("Interpolator type: {}".format(meta_image.sitk_data.interpolation_type))

        # input_size = np.array(meta_image.data_shape, dtype=np.float64)
        input_size = np.array(sitk_image.GetSize(), dtype=np.float64)
        output_size = np.int16(np.floor(input_size * scale_factors))
        output_size[0] = int(round(input_size[0] * scale_factors[0]))

        output_spacing = meta_image.voxel_size * input_size
        output_spacing = output_spacing / output_size

        input_offset = ImageProcessor.compute_offset(sitk_image, meta_image.voxel_size)
        output_offset = ImageProcessor.compute_offset(sitk_image, output_spacing)

        output_origin = meta_image.origin - input_offset + output_offset
        output_origin = np.diag(output_origin)

        output_affine = nib.affines.from_matvec(np.diag(output_spacing) * SIGN_RAS_A,
                                                [0., 0., 0.])
        output_affine[:3, 3] = output_origin

        logger.debug("Resampling output affine: {}".format(output_affine))

        identity_transform = sitk.AffineTransform(3)
        sitk_image = sitk.Cast(sitk_image, meta_image.sitk_data.interpolation_pixel)

        logger.debug("Output size: {}".format(output_size))

        resampler = sitk.ResampleImageFilter()
        resampler.SetInterpolator(meta_image.sitk_data.interpolation_type)
        resampler.SetDefaultPixelValue(0)
        resampler.SetTransform(identity_transform)
        resampler.SetSize((map(int, output_size)))
        resampler.SetOutputOrigin(tuple(output_origin))
        resampler.SetOutputSpacing(tuple(output_spacing))
        resampler.SetOutputDirection(meta_image.direction)

        if debug:
            resampler.DebugOn()

        output_image = resampler.Execute(sitk_image)
        output_image = sitk.Cast(output_image, meta_image.pixel_type)

        # sitk.WriteImage(output_image, '../results/sitk_img.nii.gz')

        output_image = sitk.GetArrayFromImage(output_image)
        output_image = output_image.transpose(meta_image.transpose)

        out_meta_image = MetaImage(output_image, output_affine, meta_image.pixel_type,
                                   meta_image.direction, meta_image.is_segmentation)

        # return MetaImage(output_image, output_affine, meta_image.pixel_type,
        #                 meta_image.direction, meta_image.is_segmentation)

        return out_meta_image

    @staticmethod
    def prev_resample_image(meta_image, scale_factors, debug=False):

        """
        :type debug: bool
        :type scale_factors: numpy.ndarray
        :type meta_image: MetaImage
        """

        scale_factors = np.array(scale_factors, dtype=np.float64)
        sitk_image = meta_image.get_sitk_image()

        logger.debug("Interpolator type: {}".format(meta_image.sitk_data.interpolation_type))

        # input_size = np.array(meta_image.data_shape, dtype=np.float64)
        input_size = np.array(sitk_image.GetSize(), dtype=np.float64)
        output_size = np.int16(np.floor(input_size * scale_factors))
        output_size[0] = int(round(input_size[0] * scale_factors[0]))

        output_spacing = meta_image.voxel_size * input_size
        output_spacing = output_spacing / output_size

        input_offset = ImageProcessor.compute_offset(sitk_image, meta_image.voxel_size)
        output_offset = ImageProcessor.compute_offset(sitk_image, output_spacing)

        output_origin = meta_image.origin - input_offset + output_offset
        output_origin = np.diag(output_origin)

        output_affine = nib.affines.from_matvec(np.diag(output_spacing) * SIGN_RAS_A,
                                                [0., 0., 0.])
        output_affine[:3, 3] = output_origin

        logger.debug("Resampling output affine: {}".format(output_affine))

        identity_transform = sitk.AffineTransform(3)
        sitk_image = sitk.Cast(sitk_image, meta_image.sitk_data.interpolation_pixel)

        logger.debug("Output size: {}".format(output_size))

        resampler = sitk.ResampleImageFilter()
        resampler.SetInterpolator(meta_image.sitk_data.interpolation_type)
        resampler.SetDefaultPixelValue(0)
        resampler.SetTransform(identity_transform)
        resampler.SetSize((map(int, output_size)))
        resampler.SetOutputOrigin(tuple(output_origin))
        resampler.SetOutputSpacing(tuple(output_spacing))
        resampler.SetOutputDirection(meta_image.direction)

        if debug:
            resampler.DebugOn()

        output_image = resampler.Execute(sitk_image)
        output_image = sitk.Cast(output_image, meta_image.pixel_type)

        # sitk.WriteImage(output_image, '../results/sitk_img.nii.gz')

        output_image = sitk.GetArrayFromImage(output_image)
        output_image = output_image.transpose(meta_image.transpose)

        out_meta_image = MetaImage(output_image, output_affine, meta_image.pixel_type,
                                   meta_image.direction, meta_image.is_segmentation)

        # return MetaImage(output_image, output_affine, meta_image.pixel_type,
        #                 meta_image.direction, meta_image.is_segmentation)

        return out_meta_image

    @staticmethod
    def get_transpose(source_space):
        target_space = np.array(list('RAS'))
        reverse_axes_dict = dict(zip('RASLPI', 'LPIRAS'))
        source_space = np.array(list(source_space))

        transpose = []
        flip = []
        for i, letter in enumerate(target_space):
            if letter == source_space[i]:
                transpose.append(i)
            elif letter == reverse_axes_dict[source_space[i]]:
                transpose.append(i)
            else:
                if letter in source_space:
                    to = np.argwhere(letter == source_space)[0][0]
                    transpose.append(to)
                else:
                    reverse_source_space = np.array([reverse_axes_dict[l] for l in source_space])
                    to = np.argwhere(letter == reverse_source_space)[0][0]
                    transpose.append(to)

        source_space = source_space[transpose]

        for i, letter in enumerate(source_space):
            if letter == target_space[i]:
                flip.append(1)
            elif letter == reverse_axes_dict[target_space[i]]:
                flip.append(-1)
            else:
                logger.debug("Something went wrong with reorienting RAI codes")

        return np.array(transpose), np.array(flip)


class BigDataViewer(object):
    def __init__(self, h5_file):

        self.h5_file = h5_file
        self.root = self.h5_file['/']
        self.setups = []

        if self.number_of_bdv_channels is None:
            self.root.attrs['BDV_setups_count'] = 0
        else:
            for i in range(self.number_of_bdv_channels):
                self.setups.append(self.h5_file[PathUtil.get_setup_path(i)])

    def add_setup(self, shape, spacing, affine):
        setup = self.h5_file[PathUtil.get_setup_path(self.number_of_bdv_channels)]
        setup.attrs['shape'] = shape
        setup.attrs['spacing'] = spacing
        setup.attrs['affine'] = affine
        self.setups.append(setup)
        self.number_of_bdv_channels = 1

    def create_xml_from_setups(self, xml_fname, fname):
        BigDataViewer.create_xml(xml_fname, fname, self.setups_shapes, self.setups_affines,
                                 self.setups_spacings)

    @property
    def number_of_bdv_channels(self):
        if self.root.attrs.__contains__('BDV_setups_count'):
            return self.root.attrs['BDV_setups_count']
        else:
            return None

    @number_of_bdv_channels.setter
    def number_of_bdv_channels(self, value):
        self.root.attrs['BDV_setups_count'] += value

    @property
    def setups_shapes(self):
        return [setup.attrs['shape'] for setup in self.setups]

    @property
    def setups_spacings(self):
        return [setup.attrs['spacing'] for setup in self.setups]

    @property
    def setups_affines(self):
        return [setup.attrs['affine'] for setup in self.setups]

    @staticmethod
    def propose_mipmaps(voxel_size, shape):
        # type: (numpy.ndarray, numpy.ndarray) -> (numpy.ndarray, numpy.ndarray)
        """Automatically determines appropriate number of mipmap levels, subsampling
        factors and chunk sizes for each level.

        Based on original BigDataViewer implementation. Chunksize is set as either
        16x16x16 or 32x32x4 depending on which one is closer to isotropic.

        Parameters
        ----------

        :param voxel_size: physical size of voxel (3D)
        :type voxel_size: np.array( float )

        :param shape: resolution of input image (3D)
        :type shape: np.array( int )

        Returns
        -------

        :rtype: tuple of two numpy arrays with subsampling factors and chunk sizes"""

        LOWEST_RES = 128
        subdiv_16_16_16 = np.array([16, 16, 16])
        subdiv_32_32_4 = np.array([[4, 32, 32],
                                   [32, 4, 32],
                                   [32, 32, 4]])

        subdiv_128_128_128 = np.array([128, 128, 128])
        subdiv_256_256_32 = np.array([[32, 256, 256],
                                      [256, 32, 256],
                                      [256, 256, 32]])
        resolutions = []
        subdivisions = []

        # keep original resolution
        res = np.array([1, 1, 1])
        voxel_scale = voxel_size / voxel_size.min()

        # compute number of levels for the mipmap pyramid
        # levels = np.int(np.ceil(np.log2(shape.max())) - np.log2(LOWEST_RES))
        while True:
            logger.debug("Mipmap resolution: {}".format(shape // res))
            resolutions.append(res)
            vscale_max = voxel_scale.max()
            shape_min = (shape // res).min()
            logger.debug("Voxel scale: {}".format(voxel_scale))

            if (4 * vscale_max / 32.) > (1. / vscale_max):
                if shape_min > 128:
                    subdivisions.append(subdiv_256_256_32[voxel_scale.argmax()])
                else:
                    subdivisions.append(subdiv_32_32_4[voxel_scale.argmax()])
            else:
                if shape_min > 128:
                    subdivisions.append(subdiv_128_128_128)
                else:
                    subdivisions.append(subdiv_16_16_16)

            res = res.copy()
            voxel_scale = voxel_scale.copy()

            res[voxel_scale < 2.0] = res[voxel_scale < 2.0] * 2
            voxel_scale[voxel_scale < 2.0] = voxel_scale[voxel_scale < 2.0] * 2

            voxel_scale = voxel_scale / voxel_scale.min()
            if (shape // res).max() < LOWEST_RES or (shape // res).min() < 32:
                break

        return np.array(resolutions, dtype=np.float64), np.array(subdivisions)

    @staticmethod
    def create_xml(xml_fname, fname, img_shapes,
                   img_affines, voxel_sizes, unit='mm', partition_fnames=None):
        '''
        Creates generic BigDataViewer xml file for single setup-single timepoint
        image

        Parameters
        ----------
        param xml_fname: path to output xml file
        type xml_fname: str

        param fname: path to the hdf5 file
        type fname: str

        param img_shape: highest resolution of the 3D image
        type img_shape: array[int, int, int]

        param img_affine: affine array defining position of the image data array in
                          a reference space
        type img_affine: 4x4 matrix of floats
        '''
        output_dir, fname = os.path.split(fname)
        sizes = [' '.join([str(x).strip() for x in img_shape])
                 for img_shape in img_shapes]
        voxel_sizes = [' '.join([str(x).strip() for x in voxel_size])
                       for voxel_size in voxel_sizes]

        np.set_printoptions(suppress=True, precision=4)

        if partition_fnames is None:
            affinev = [''.join(np.str(
                list(np.hstack(img_affine[:-1])))[1:-1].split(','))
                       for img_affine in img_affines]
        else:
            affinev = [' '.join([str(x) for x in img_affine])
                       for img_affine in img_affines]

        root = etree.Element("SpimData")
        root.attrib['version'] = '0.2'
        BasePath = etree.SubElement(root, "BasePath")
        BasePath.attrib["type"] = 'relative'
        BasePath.text = '.'
        SequenceDescription = etree.SubElement(root, "SequenceDescription")
        ViewRegistrations = etree.SubElement(root, "ViewRegistrations")

        #   SequenceDescription
        ImageLoader = etree.SubElement(SequenceDescription, "ImageLoader")
        ImageLoader.attrib['format'] = "bdv.hdf5"
        hdf5 = etree.SubElement(ImageLoader, "hdf5")
        hdf5.attrib['type'] = "relative"
        hdf5.text = fname

        if partition_fnames is not None:
            for fname in partition_fnames:
                partition = etree.SubElement(ImageLoader, "partition")
                path = etree.SubElement(partition, "path")
                path.attrib['type'] = "relative"
                path.text = fname
                timepoints = etree.SubElement(partition, "timepoints")
                timepoints.text = "0"
                setups = etree.SubElement(partition, "setups")
                setups.text = "0"

        ViewSetups = etree.SubElement(SequenceDescription, "ViewSetups")
        for i, vs in enumerate(sizes):
            ViewSetup = etree.SubElement(ViewSetups, "ViewSetup")
            vs_id = etree.SubElement(ViewSetup, "id")

            vs_id.text = str(i)  # mandatory
            vs_name = etree.SubElement(ViewSetup, "name")
            vs_name.text = "channel " + str(i)  # optional
            vs_size = etree.SubElement(ViewSetup, "size")
            vs_size.text = sizes[i]  # optional
            vs_voxelSize = etree.SubElement(ViewSetup, "voxelSize")
            vs_voxelSize_unit = etree.SubElement(vs_voxelSize, "unit")
            vs_voxelSize_unit.text = unit
            vs_voxelSize_size = etree.SubElement(vs_voxelSize, "size")
            vs_voxelSize_size.text = voxel_sizes[i]

        Timepoints = etree.SubElement(SequenceDescription, "Timepoints")
        Timepoints.attrib['type'] = "range"
        first = etree.SubElement(Timepoints, "first")
        first.text = "0"
        last = etree.SubElement(Timepoints, "last")
        last.text = "0"

        # ViewRegistrations
        for i, vr in enumerate(img_affines):
            ViewRegistration = etree.SubElement(ViewRegistrations,
                                                "ViewRegistration")
            ViewRegistration.attrib['timepoint'] = "0"
            ViewRegistration.attrib['setup'] = str(i)
            ViewTransform = etree.SubElement(ViewRegistration, "ViewTransform")
            ViewTransform.attrib['type'] = "affine"
            affine = etree.SubElement(ViewTransform, "affine")
            affine.text = affinev[i]

        doc = etree.ElementTree(root)
        doc.write(os.path.join(output_dir, xml_fname),
                  encoding="UTF-8",
                  pretty_print=True)


def test_proxy(input_path, json_path):
    meta = dm.ImageMetaData(input_path)
    with open(json_path) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)

    # ip = StreamableTiffProxy('autofluo', ' mysz_2', "auto_8_Z%03d.ome.tif", meta, 0.3)
    ip = StreamableOMEProxy('cfos', 'fos_4', None, meta, 1)
    # ip = NonStreamableTiffProxy('cfos', 'mysz', None, meta, None)
    # ip = NiftiProxy('autofluo', 'exp2', None, meta, None)
    print(ip)

    k = 0
    for data_chunk in ip.stream_data():
        if k > 0:
            return
        print data_chunk.shape
        k += 1


def test_lmdhf(image_path, json_path, hdf_path, multichannel=False):
    meta = dm.ImageMetaData(image_path)
    with open(json_path) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)

    # ip = NiftiProxy('autofluo', 'exp2', None, meta, None)
    #ip = ImageProxy.get_image_proxy_class(meta)('cfos2', 'mysz', None, meta, 'single_DIR.xml', None)
    ip = ImageProxy.get_image_proxy_class(meta)('fos', 'exp_4', 'Z%06d.tif', meta, 'exp4.xml', 4.,
                                               is_multichannel=multichannel)
    # ip = StreamableOMEProxy('cfos', 'fos_4', None, meta, 1)

    lm_test = LightMicroscopyHDF(hdf_path)
    lm_test.write_channel(ip)
    print lm_test
    # lm_test.channels[0].export_image(level=4)

    # print lm_test.number_of_channels


def test_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)
    channel = lmf.get_channel('cfos2')
    channel.write_affine('to_atlas', '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/mysz_affine.txt')
    list_of_transforms = [TransformTuple(0, 'to_atlas', 'affine', True)]
    e_cmd = ExportCmd(channel_name='cfos2', output_path='../results/_RAS_cfos2_35um.nii.gz',
                      output_resolution=[0.035, 0.035, 0.035],
                      input_orientation='SPR', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)

    lmf.export_image(e_cmd)


def write_cfos_affine(hdf_path):

    lmf = LightMicroscopyHDF(hdf_path)
    channel = lmf.get_channel('right_template')
    channel.write_affine('auto_to_template',
                         '/home/sbednarek/DEV/lsfm/resources/deformation_field_test/30_transforms/template_right_physical_to_structural_001_Affine.txt')


def test_cfos_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)


    list_of_transforms = [TransformTuple(0, 'new_auto_to_template', 'affine', True)]
    e_cmd = ExportCmd(channel_name='autofluo', output_path='../results/_cfos_15um.nii.gz',
                      output_resolution=[0.015, 0.015, 0.015],
                      input_orientation='RSP', input_resolution_level=None,
                      list_of_transforms=[],
                      phys_origin=None, phys_size=None)

    lmf.export_image(e_cmd)


def test_cfos_add_displacement_field(image_path, json_path, hdf_path):

    lmf = LightMicroscopyHDF(hdf_path)

    meta = dm.ImageMetaData(image_path)
    with open(json_path) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)
    ip = ImageProxy.get_image_proxy_class(meta)('inverse_warp_template', 'mysz', None, meta, None, 2.,
                                                is_multichannel=True)

    lmf.write_channel(ip)


def test_cfos_apply_displacement_field(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)

    list_of_transforms = [TransformTuple(0, 'inverse_warp_template', 'df', True),
                          TransformTuple(1, 'new_auto_to_template', 'affine', True)]

    # list_of_transforms = [TransformTuple(0, 'new_auto_to_template', 'affine', True)]
    e_cmd = ExportCmd(channel_name='autofluo', output_path='../results/_cfos_df_af_origin_25sum.nii.gz',
                      output_resolution=[0.025, 0.025, 0.025],
                      input_orientation='RSP', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)

    lmf.export_image(e_cmd)


def write_template_to_hdf(image_path, json_path, hdf_path):
    meta = dm.ImageMetaData(image_path)
    with open(json_path) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)

    ip = ImageProxy.get_image_proxy_class(meta)('right_template', 'mysz', None, meta, 'template.xml', 2.,
                                                is_multichannel=False, is_segmentation=False)

    lm_test = LightMicroscopyHDF(hdf_path)
    lm_test.write_channel(ip)


def write_seg_to_hdf(image_path, json_path, hdf_path):
    meta = dm.ImageMetaData(image_path)
    with open(json_path) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)

    ip = ImageProxy.get_image_proxy_class(meta)('autofluo_segmentation', 'mysz', None, meta, 'w.xml', 2.,
                                                is_multichannel=False, is_segmentation=True)

    lm_test = LightMicroscopyHDF(hdf_path)
    lm_test.write_channel(ip)


def test_template_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)

    list_of_transforms = [TransformTuple(0, 'forward_warp', 'df', False),
                          TransformTuple(1, 'auto_to_template', 'affine', False)]

    list_of_transforms = [TransformTuple(0, 'auto_to_template', 'affine', False),
                          TransformTuple(1, 'forward_warp', 'df', False)]

    list_of_transforms = []

    e_cmd = ExportCmd(channel_name='right_template', output_path='../results/template_2_auto_30um.nii.gz',
                      output_resolution=None,
                      input_orientation='RAS', input_resolution_level=1,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)
    #lmf.export_image(e_cmd)

    e_cmd = ExportCmd(channel_name='right_template', output_path='../results/template_40um.nii.gz',
                      output_resolution=[0.04, 0.04, 0.04],
                      input_orientation='RAS', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)
    lmf.export_image(e_cmd)

    e_cmd = ExportCmd(channel_name='segmentation', output_path='../results/segmentation_to_auto_30um.nii.gz',
                      output_resolution=[0.030, 0.030, 0.030],
                      input_orientation='RAS', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)

    #lmf.export_image(e_cmd)


def write_exp4_affine(hdf_path):

    lmf = LightMicroscopyHDF(hdf_path)
    channel = lmf.get_channel('cfos')
    channel.write_affine('cfos_to_auto',
                         '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4/structure_001_to_signal_001_physical_affine.txt')

    channel.write_affine('cfos_to_template',
                         '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4/template_right_physical_to_signal_001_Affine.txt')


def write_transforms_signal(hdf_path, channel_name, dir_path,
                            inverse_warp='inverse_warp_structural_001.nii.gz',
                            inverse_warp_json='inverse_warp.json',
                            signal_to_structural='structure_001_to_signal_001_physical_affine.txt',
                            signal_to_template='template_right_physical_to_signal_001_Affine.txt'):

    lmf = LightMicroscopyHDF(hdf_path)
    channel = lmf.get_channel(channel_name)

    channel.write_affine('cfos_to_auto',
                         os.path.join(dir_path, signal_to_structural))

    channel.write_affine('cfos_to_template',
                         os.path.join(dir_path, signal_to_template))

    meta = dm.ImageMetaData(os.path.join(dir_path, inverse_warp))
    with open(os.path.join(dir_path, inverse_warp_json)) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)

    ip = ImageProxy.get_image_proxy_class(meta)('inverse_warp_template', 'some_id', None, meta, None, 2.,
                                                is_multichannel=True)

    lmf.write_channel(ip)


def add_displacement_field(image_path, json_path, hdf_path):

    lmf = LightMicroscopyHDF(hdf_path)

    meta = dm.ImageMetaData(image_path)
    with open(json_path) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)
    ip = ImageProxy.get_image_proxy_class(meta)('inverse_warp_template', 'autofluo', None, meta, None, 2.,
                                                is_multichannel=True)

    lmf.write_channel(ip)


def test_exp4_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)

    #list_of_transforms = [TransformTuple(0, 'auto_to_template', 'affine', True)]
    list_of_transforms = [TransformTuple(0, 'inverse_warp_template', 'df', True),
                          TransformTuple(1, 'auto_to_template', 'affine', True)]
    e_cmd = ExportCmd(channel_name='autofluo', output_path='../results/exp_4/exp4_10um_affine_df_ref.nii.gz',
                      output_resolution=[0.010, 0.010, 0.010],
                      input_orientation='PSR', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)

    lmf.export_image(e_cmd)


def test_exp4_cfos_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)

    #list_of_transforms = [TransformTuple(0, 'auto_to_template', 'affine', True)]
    list_of_transforms = [TransformTuple(0, 'cfos_to_auto', 'affine', True),
                          TransformTuple(1, 'inverse_warp_template', 'df', True),
                          TransformTuple(2, 'cfos_to_template', 'affine', True)]
    e_cmd = ExportCmd(channel_name='fos', output_path='../results/exp_4/level_2_fos.nii.gz',
                      output_resolution=[0.01, 0.01, 0.01],
                      input_orientation='PSR', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)

    lmf.export_image(e_cmd)


def test_region_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)
    e_cmd = ExportCmd(channel_name='fos', output_path='../results/exp_4/extracted_reg_382_level_0.nii.gz',
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=0,
                      list_of_transforms=[],
                      phys_origin=None, phys_size=None,
                      segmentation_name='right_signal_segmentation', region_id=382)

    lmf.export_image(e_cmd)


def test_chunk_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)
    e_cmd = ExportCmd(channel_name='fos', output_path='../results/exp_4/chunk_level_0.nii.gz',
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=0,
                      list_of_transforms=[],
                      phys_origin=[2.84, 1.81, -4.524], phys_size=[.5, .5, .5],
                      segmentation_name=None, region_id=None)

    lmf.export_image(e_cmd)


def test_chunk_transform_export(hdf_path, channel_name='cfos', output_path='../results/check_margis_2.nii.gz'):
    lmf = LightMicroscopyHDF(hdf_path)
    list_of_transforms = [TransformTuple(0, 'cfos_to_auto', 'affine', True),
                          TransformTuple(1, 'inverse_warp_template', 'df', True),
                          TransformTuple(2, 'cfos_to_template', 'affine', True)]
    e_cmd = ExportCmd(channel_name=channel_name,
                      output_path=output_path,
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=0,
                      list_of_transforms=list_of_transforms,
                      phys_origin=[7.475, 5.975, -4.475], phys_size=[2.25, 1.25, 3.25],
                      #phys_origin=[7.475, 5.975, -4.475], phys_size=[1., 1., 1.],
                      #phys_origin=[5.325, 6.6, -1.975], phys_size=[1.7, 1.5, 1.5],
                      segmentation_name=None, region_id=None)

    lmf.export_image(e_cmd)


def test_auto_chunk_transform_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)
    list_of_transforms = [TransformTuple(0, 'inverse_warp_template', 'df', True),
                          TransformTuple(1, 'auto_to_template', 'affine', True)]
    e_cmd = ExportCmd(channel_name='autofluo', output_path='../results/exp_4/auto.nii.gz',
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=0,
                      list_of_transforms=list_of_transforms,
                      phys_origin=[7.475, 5.975, -4.475], phys_size=[2.25, 1.25, 3.25],
                      segmentation_name=None, region_id=None)

    lmf.export_image(e_cmd)


def test_slice_export(hdf_path, channel, output_path):
    lmf = LightMicroscopyHDF(hdf_path)
    e_cmd = ExportSlicesCmd(channel_name=channel, output_path=output_path, input_orientation="PSR",
                            input_resolution_level=1, slicing_range=[20, 400, 2], axis=1,
                            extract_roi=[1, 1, 500, 500],
                            ref_channel=None, ref_level=None)

    lmf.export_slices(e_cmd)


if __name__ == '__main__':
    logger.info("Testing..")

    input_path = '/media/sbednarek/4BCFEE837AD4D9DD/Diana/new_autofluo_2018_01/fos_8/Z_planes/auto_8_Z000.ome.tif'
    json_path = '/home/sbednarek/DEV/lsfm/results/fos_8_metadata.json'
    ome_input_path = '/media/sbednarek/4BCFEE837AD4D9DD/Ctr1.ome.tif'
    ome_json_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/metadata/ctr1_ome.json'
    single_input_path = '/media/sbednarek/4BCFEE837AD4D9DD/TDPrat_mapa_terasticher_downsampled.tif'
    single_json_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/metadata/TDPrat.json'
    nifti_path = '/home/sbednarek/DEV/lsfm/results/Marzena_02_2018/exp2_autofluo_25um.nii.gz'
    nifti_json_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/metadata/exp2_autofluo_25.json'
    df_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4/inverse_warp_structural_001.nii.gz'
    df_json_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4/inverse_warp_meta.json'
    mysz_path = '/home/sbednarek/DEV/lsfm/resources/mysz_mapa.ome.tif'
    mysz_json_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/metadata/mysz_mapa_ome.json'
    template_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/04_new_avrgt_25_right_hemisphere.nii'
    json_template_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/metadata/template_25.json'
    seg_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4/segmentation_in_signal-001-25_deformable.nii.gz'
    auto_seg_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4/segmentation_in_structural-001-25_deformable.nii.gz'
    auto_seg_json ='/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/metadata/segementation_in_autofluo_11.json'
    json_seg_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/metadata/segmentation_in_signal_11.json'
    df_forward_path = '/home/sbednarek/DEV/lsfm/resources/deformation_field_test/30_transforms/forward_warp_structural_001.nii.gz'
    df_forward_json_path = '../results/metadata/forward_warp.json'

    exp_4_json = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/resources/647-exp4.json'
    exp_4_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/resources/exp4_cfos/Z000000.tif'

    # test_proxy(input_path, json_path)
    # test_proxy(ome_input_path, ome_json_path)
    # test_proxy(single_input_path, single_json_path)
    # test_proxy(nifti_path, nifti_json_path)

    # test_lmdhf(single_input_path, single_json_path, '../results/single_test_DIRECTION.h5')
    # test_lmdhf(ome_input_path, ome_json_path, '../results/testing_ome2.h5')
    # test_lmdhf(nifti_path, nifti_json_path, '../results/testing_nifti.h5')
    # test_lmdhf(df_path, df_json_path, '../results/df_test2.h5', True)
    # test_lmdhf(mysz_path, mysz_json_path, '../results/myszdf_RAS.h5')

    #test_export('../results/myszdf_RAS.h5')

    # test_export('../results/single_test.h5')

    # test_cfos_export('/home/sbednarek/DEV/lsfm/results/cfos.h5')

    # test_cfos_add_displacement_field(df_path, df_json_path, '/home/sbednarek/DEV/lsfm/results/cfos.h5')
    # test_cfos_apply_displacement_field('/home/sbednarek/DEV/lsfm/results/cfos.h5')

    # write_template_to_hdf(template_path, json_template_path, '../results/template.h5')

    # test_template_export('../results/template.h5')

    # test_lmdhf(df_forward_path, df_forward_json_path, '../results/template.h5', True)

    # write_cfos_affine('../results/template.h5')

    #test_template_export('../results/template.h5')

    #test_lmdhf(exp_4_path, exp_4_json, '../results/exp_4_new.h5', False)

    #write_exp4_affine('../results/exp_4_new.h5')

    #write_exp4_affine('/home/sbednarek/DEV/lsfm/results/experimental_4.h5')

    #add_displacement_field(df_path, df_json_path, '../results/exp_4_new.h5')

    #add_displacement_field(df_path, df_json_path, '/home/sbednarek/DEV/lsfm/results/experimental_4.h5')

    #test_exp4_export('../results/exp_4_new.h5')

    #test_exp4_cfos_export('../results/exp_4_new.h5')

    #write_seg_to_hdf(seg_path, json_seg_path, '../results/exp_4_new.h5')

    #test_region_export('../results/exp_4_new.h5')

    #test_chunk_export('../results/exp_4_new.h5')

    #test_chunk_transform_export('../results/exp_4_new.h5')
    #test_auto_chunk_transform_export('../results/exp_4_new.h5')

    test_chunk_transform_export('/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4_new.h5',
                                'fos', output_path='/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4/check_margis_1.nii.gz')

    #test_chunk_transform_export('/mnt/nicl/home/pmajka/results/experimental_5.h5',
    #                            output_path='../results/experimental_5/in_cfos_whole_amygdala_native_cfos.nii.gz')

    #write_transforms_signal('/mnt/nicl/home/pmajka/results/experimental_5.h5', 'cfos',
    #                        '/home/sbednarek/DEV/lsfm_pipeline/000012/lsfm_image_server/30_transforms/')

    #test_chunk_transform_export('/mnt/nicl/home/pmajka/results/experimental_3.h5',
    #                            output_path='../results/experimental_3/in_cfos_whole_amygdala_native_cfos.nii.gz')

    #test_chunk_transform_export('/mnt/nicl/home/pmajka/results/experimental_2.h5',
    #                            output_path='../results/experimental_2/in_cfos_whole_amygdala_native_cfos.nii.gz')

    #test_chunk_transform_export('/mnt/nicl/home/pmajka/results/control_2.h5',
    #                            output_path='../results/control_2/in_cfos_whole_amygdala_native_cfos.nii.gz')


    #test_slice_export('../results/exp_4_new.h5', 'autofluo', '../results/slices/img_%04d.tif')

    #write_seg_to_hdf(auto_seg_path, auto_seg_json, '../results/exp_4_new.h5')