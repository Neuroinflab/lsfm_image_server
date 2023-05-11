#!/usr/bin/env python
# -*- coding: utf-8 -*
###############################################################################
#                                                                             #
#    LSFM image server                                                        #
#                                                                             #
#    Copyright (C) 2017-2018 Sylwia Bednarek, Piotr Majka                     #
#    (Laboratory of Neuroinformatics; Nencki Institute of Experimental        #
#    Biology of Polish Academy of Sciences)                                   #
#                                                                             #
#    This software is free software: you can redistribute it and/or modify    #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This software is distributed in the hope that it will be useful,         #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this software.  If not, see http://www.gnu.org/licenses/.     #
#                                                                             #
###############################################################################

import os
import logging

import numpy as np
import nibabel as nib
import tifffile as tf
import SimpleITK as sitk

from medpy import io as medio
from .utils import InputImageType, parallelize, read_image

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


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
                 is_segmentation=False,
                 n_cpus=15):

        self.channel_name = channel_name
        self.xml_file_name = xml_file_name
        self.specimen_id = specimen_id
        self.multifile_prefix = multifile_prefix
        self.metadata = metadata
        self.RAM_limit = RAM_limit
        self.is_multichannel = is_multichannel
        self.is_segmentation = is_segmentation
        self.orientation = orientation
        self.n_cpus = n_cpus

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
        reverse_axes_dict = dict(list(zip('RASLPI', 'LPIRAS')))
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
        parallel_read_image = parallelize(read_image, self.n_cpus)

        current_idx = 0
        while current_idx < self.data_shape[0]:
            logger.debug("Current streaming index: {}".format(current_idx))
            if current_idx + self.stack_size < self.data_shape[0]:
                next_indices = np.arange(current_idx, current_idx + self.stack_size)
            else:
                next_indices = np.arange(current_idx, self.data_shape[0])
            try:
                logger.debug("Next stream slice: {}".format(next_indices))
                paths = [self._craft_img_path(n) for n in next_indices]
                data = parallel_read_image(paths)
                data = np.dstack(np.array(sorted(data, key=lambda x: x[1]))[:, 0])
                logger.debug("done reading and sorting")
            except IOError:
                logger.error("Could not open file", exc_info=True)
                raise

            current_idx = current_idx + self.stack_size
            data = data.transpose(2, 0, 1)
            data_chunk = np.zeros([self.stack_size, self.data_shape[1], self.data_shape[2]])
            data_chunk[:data.shape[0], ...] = data
            yield data_chunk

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
            data = np.zeros([self.stack_size, self.data_shape[1], self.data_shape[2]])
            data[:next_indices.stop - next_indices.start, ...] = tiff_file.asarray(next_indices,
                                                                                   0)[:next_indices.stop -
                                                                                       next_indices.start]
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
        self.pixel_type = self._get_sitk_pixel_type()

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


