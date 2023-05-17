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
import sys
import json
import h5py
import copy
import errno
import imageio
import inspect

import numpy as np
import nibabel as nib
import SimpleITK as sitk
import datetime as dt

from . import bdv
from . import image
from . import export
from . import transforms
from . import processing
from . import dump_metadata as dm
from . import constants as const

from .utils import PathUtil, highestPowerof2
#from scipy.misc import imsave

import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


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

        str_options = '\n'.join("%s=%r" % (key, val) for (key, val) in list(options.items()))

        with open('/proc/{}/cmdline'.format(os.getpid())) as c:
            str_options = '\n'.join([str_options, "cmdline={}".format(' '.join(c.read().split(str(b'\x00','utf-8','ignore'))))])

        str_options = np.string_(str_options)
        self.log_operation(func.__name__, str_options)

        return func(self, *args, **kwargs)
    return func_wrapper


class LightMicroscopyHDF(object):
    def __init__(self, file_path, access_mode='a'):

        self.file_path = file_path
        self.channels = dict()

        try:
            if not os.path.exists(os.path.dirname(self.file_path)):
                os.makedirs(os.path.dirname(self.file_path))
            self.h5_file = h5py.File(self.file_path, access_mode, rdcc_nbytes=1024 * 1024 * 64, libver='latest')
        except IOError as e:
            if e.errno == errno.EACCES:
                self.h5_file = h5py.File(self.file_path, 'r', rdcc_nbytes=1024 * 1024 * 64, libver='latest')
            else:
                logger.error("File not found or couldn't open file", exc_info=True)
                raise

        self.root = self.h5_file['/']
        self.bdv = bdv.BigDataViewer(h5_file=self.h5_file)

        if self.number_of_channels is None:
            self.root.attrs['LSFM_setups_count'] = 0
        else:
            self._populate_channels()

    def __repr__(self):
        metadata = '\n'.join("%s=%r" % (key, val) for (key, val) in list(self.extended_metadata.items()))
        channel_info = '\n'.join(str(channel) for channel in list(self.channels.values()))
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

    @logging_decor
    def get_metadata(self, m_key):
        try:
            m_value = self.extended_metadata[m_key]
        except KeyError:
            logger.error("Key not found!")
            return

        return m_value

    def log_operation(self, cmd_name, cmd_line):
        if self.h5_file.mode == 'r':
            logger.info('File in read mode only, cannot log operations!')
            return
        log_path = PathUtil.get_log_path(dt.datetime.now(),
                                         cmd_name + str(np.random.randint(10**5)))
        logger.debug('log path: {:} \n command-name: {:}\n command: {:}'.format(log_path,
                                                                                cmd_name,
                                                                                cmd_line))
        self.h5_file.create_dataset(name=log_path, data=cmd_line)

    def get_logs(self):
        log_header, log_data = [], []

        def visit_log(name, obj):
            log_header.append(name)
            log_data.append(str(obj[...]))

        self.h5_file[PathUtil.get_logs_path()].visititems(visit_log)
        return list(zip(log_header, log_data))

    def close(self):
        self.h5_file.close()

    def open(self):
        self.h5_file = h5py_cache.File(self.file_path, 'a', chunk_cache_mem_size=1024 ** 3, libver='latest')

    @property
    def channel_names(self):
        if '/LSFM' in self.h5_file:
            return list(self.h5_file['/LSFM'].keys())
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
                                                           compression_opts=6,
                                                           shuffle=True,
                                                           fletcher32=True)
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
                meta_image = processing.MetaImage.meta_image_from_channel(data, _channel)
                if _channel.image.is_stream:
                    logger.debug("Processing data : {}/{}".format(data_idx + 1,
                                                                  _channel.image.num_of_stacks + 1))
                for level_num, pyramid_level in enumerate(_channel.levels):

                    if level_num == 0:
                        writing_slice = slice(0, None)

                        if _channel.image.is_stream:
                            writing_slice = slice(data_idx * pyramid_level.stack_shape[0],
                                                  data_idx * pyramid_level.stack_shape[0] + data.shape[0])
                        try:
                            logger.debug("Writing main level")
                            data_sets[level_num][writing_slice, :, :] = data
                            logger.debug("Done writing main level")
                        except TypeError:
                            logger.error("Probably wrong data shape", exc_info=True)
                            raise
                    else:

                        scale_factors = _channel.relative_resolution_scales[level_num][:len(_channel.image.image_shape)]
                        scale_factors[0] = float(_channel.levels[level_num].shape[0]) / \
                                           _channel.levels[level_num-1].shape[0]
                        logger.debug("Scale factors for resampling: {}".format(scale_factors))
                        meta_image = processing.resample_image(meta_image, scale_factors)
                        writing_slice = slice(0, None)

                        if _channel.image.is_stream:
                            writing_slice = slice(data_idx * pyramid_level.stack_shape[0],
                                                  data_idx * pyramid_level.stack_shape[0] +
                                                  meta_image.data_shape[0])
                        try:
                            logger.debug("Writing sublevel")
                            data_sets[level_num][writing_slice, :, :] = meta_image.data
                            logger.debug("done writing sublevel")
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
        self.p_origins = self._get_image_meta_data(const.M_ORIGIN)
        self.p_shapes = self._get_image_meta_data(const.M_SHAPE)
        self.p_voxel_sizes = self._get_image_meta_data(const.M_VOXEL_SIZE)
        self.n_levels = len(self.p_shapes)
        try:
            self.is_multichannel = bool(int(self._get_image_meta_data('is_multichannel')))
            self.is_segmentation = bool(int(self._get_image_meta_data('is_segmentation')))
            self.pixel_type = int(self._get_image_meta_data('pixel_type'))
        except KeyError:
            self.is_multichannel = False
            self.is_segmentation = False
            self.pixel_type = 3

        # if self.is_segmentation:
        #    logger.debug("{} ".format(self.is_segmentation))
        #    logger.debug("{} is segmentation".format(self.channel_name))
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
            p_level.is_nifti = self.is_multichannel or self.is_segmentation
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
        affine = transforms.Affine(affine_path=affine_file_path, affine_name=affine_name)
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
        affine = transforms.Affine(affine_name=affine_name)
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
        cmd = export.ExportCmd(None, None, None, 'RAS', 0, [], None, None, None, None, None, None, None)
        df = df_level.get_meta_image(cmd)
        df = sitk.Cast(df, sitk.sitkVectorFloat64)
        dft = sitk.DisplacementFieldTransform(df)
        return dft

    def read_segmentation(self, seg_name):
        try:
            seg_channel = HDFChannel(self.h5_file, seg_name)
        except ValueError:
            logger.error("Segmentation {} not found!".format(seg_name))
            sys.exit()

        seg_level = seg_channel.pyramid_levels[0]
        cmd = export.ExportCmd(None, None, None, 'RAS', 0, [], None, None, None, None, None, None, None)
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
        cmd = export.ExportCmd(None, None, None, 'RAS', 0, [], None, None, None, None, None, None, None)
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
        _transforms = []
        for transform_tuple in list_of_transforms:
            if transform_tuple.type == 'affine':
                affine = self.read_affine(transform_tuple.name)
                if transform_tuple.invert:
                    affine = affine.itk_inv_affine
                else:
                    affine = affine.itk_affine
                transform = transforms.Transform(transform_tuple.type, transform_tuple.order,
                                                 transform_tuple.name, affine)
                _transforms.append(transform)
            elif transform_tuple.type == 'df':
                df = self.read_displacement_field(transform_tuple.name)
                transform = transforms.Transform(transform_tuple.type, transform_tuple.order,
                                                 transform_tuple.name, df)
                _transforms.append(transform)
            else:
                raise NotImplementedError("Wrong transformation type")

        return _transforms

    def export_grid_of_chunks(self, export_cmd, level):
        """

        :param export_cmd:
        :param level:
        :return:
        """
        def get_whole_image_grid():
            logger.debug("Exporting whole image as chunks")
            transpose, flip = processing.get_transpose(export_cmd.input_orientation)
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
        msk = list(map(np.all, spacings <= resolution))
        distances = list(map(np.linalg.norm, (spacings[msk] - resolution)))
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
            affines = list(self.h5_file[PathUtil.get_lsfm_affines_group_path(self.channel_name)].keys())
            affines = "\n".join(affines)
        except KeyError:
            affines = "No affines were found for this channel."

        try:
            displacement_fields = list(self.h5_file[PathUtil.get_lsfm_dfs_path(self.channel_name)].keys())
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
            transpose, flip = processing.get_transpose(cmd.ref_orientation)
            inv_transpose = np.argsort(transpose)
            logger.debug("transpose: {}\n flip: {}".format(transpose, flip))
            target_spacing = cmd.ref_channel.pyramid_levels[cmd.ref_level].voxel_size[transpose]
            logger.debug("target_spacing: {}".format(target_spacing))
            logger.debug("source spacing: {}".format(self.voxel_size))

            scale_factors = self.voxel_size / target_spacing
            logger.debug("scaling_factors: {}".format(scale_factors))

            target_size = cmd.ref_channel.pyramid_levels[cmd.ref_level].shape[transpose]
            logger.debug("target size: {}".format(target_size))
            # source_size = self.shape[transpose]
            logger.debug("source image size: {}".format(self.shape))

            cmd_tmp = export.ExportCmd(None, None, None, 'RAS', 0, [], None, None, None, None)
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
            print((float(self.shape[cmd.axis])))
            output_planes_per_input_planes = target_size[cmd.axis] / float(self.shape[cmd.axis])

            logger.debug("output planes per input planes: {}".format(output_planes_per_input_planes))
            roi_filter.SetSize(slab_size)
            plane_no = 0
            for slab_start in range(cmd.start, cmd.stop - cmd.step, cmd.step):
                start_index[cmd.axis] = slab_start
                logger.debug("start index: {}".format(start_index))
                logger.debug("slab size: {}".format(slab_size))
                roi_filter.SetIndex(start_index)
                roi_slab = roi_filter.Execute(seg)
                resampled_slab = processing.resample_sitk(roi_slab, scale_factors, sitk.sitkNearestNeighbor)
                slab_data = sitk.GetArrayFromImage(resampled_slab)
                logger.debug("slab shape: {}".format(slab_data.shape))
                # some_plane = slab_data[:, 5, :]
                # imsave(cmd.output_path % slab_start + cmd.output_ext, some_plane)
                current_plane = start_index[cmd.axis] * output_planes_per_input_planes
                logger.debug("current plane: {}".format(current_plane))

                # chosen_planes = slab_data[:, ::cmd.step, :]
                chosen_planes = slab_data
                for i in range(chosen_planes.shape[cmd.axis]):
                    # print("---> "+ str(current_plane + i))
                    if round(current_plane + i) % cmd.step == 0:
                        # print(current_plane + i)

                        imageio.imwrite(cmd.output_path % round(current_plane + i) + cmd.output_ext,
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

            transpose, flip = processing.get_transpose(cmd.input_orientation)
            logger.debug("transpose: {}\n flip: {}".format(transpose, flip))
            inv_transpose = np.argsort(transpose)
            logger.debug("inverse transpose: {}".format(inv_transpose))
            inv_flip = flip[inv_transpose] < 0
            axis = transpose[cmd.axis]
            # axis_inverse = flip[cmd.axis] < 0
            # logger.debug("inverse main axis: {}".format(axis_inverse))
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
            # if axis_inverse:
            # list_of_batch_slices = list_of_batch_slices[::-1]
            logger.debug("List of batch slices: {}".format(list_of_batch_slices))

            for main_slice in list_of_batch_slices:
                # main_slice = slice(main_slice.start * -1, main_slice.stop * -1, main_slice.step)
                logger.debug("main slice: {}".format(main_slice))
                sl_dict[axis] = main_slice
                logger.debug("sl_dict: {}".format(sl_dict))
                data = self.h5_file[self.path][sl_dict[0], sl_dict[1], sl_dict[2]]
                data = data.transpose([axis, axis_h, axis_w])
                # if axis_inverse:
                    # data = data[::-1]
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

            logger.debug("data_slx: {}".format(data_slx))
            logger.debug("data_sly: {}".format(data_sly))
            logger.debug("data_slz: {}".format(data_slz))
            proper_chunk[data_slx, data_sly, data_slz, ...] = data

            return proper_chunk

        def get_meta_image(self, export_cmd):

            if export_cmd.whole_image:
                img_exp = export.ImageWholeExporter(export_cmd, self)
            else:
                if export_cmd.export_region:
                    img_exp = export.ImageRegionExporter(export_cmd, self)
                else:
                    img_exp = export.ImageExporter(export_cmd, self)

            output_img = img_exp.process()

            if export_cmd.output_path:
                sitk.WriteImage(output_img, export_cmd.output_path)
                return True
            else:
                return output_img

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
        self.max_subdivision = 1024

        self.resolutions, self.subdivisions = bdv.BigDataViewer.propose_mipmaps(self.image.voxel_size,
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

            stack_shape = shape.copy()

            if self.image.is_stream:
                shape[0] = shape[0] * (self.image.num_of_stacks + 1)

            logger.debug("Level: {}\n"
                         "SHAPE: {}\n"
                         "STACK SHAPE: {}\n"
                         "N of stacks: {}".format(num_level, shape, stack_shape,
                                                  self.image.num_of_stacks))

            pyramid_level.shape = shape
            pyramid_level.stack_shape = stack_shape
            pyramid_level.chunks = tuple([min(self.max_subdivision,
                                         highestPowerof2(int(x))) for x in pyramid_level.stack_shape])

        pyramid_level = ProxyChannel.PyramidLevel()
        pyramid_level.id = 0
        pyramid_level.shape = self.image.data_shape.copy()
        if self.image.is_stream:
            pyramid_level.shape[0] = self.image.stack_shape[0] * (self.image.num_of_stacks + 1)
        pyramid_level.spacing = self.image.voxel_size
        pyramid_level.stack_shape = self.image.stack_shape
        pyramid_level.path = PathUtil.get_lsfm_image_cells_path(self.image.channel_name, 0)
        pyramid_level.bdv_path = PathUtil.get_cells_path(0, self.num_bdv_setup, 0)
        logger.debug("0 level shape: {}".format(self.image.stack_shape) )
        pyramid_level.chunks = tuple([min(self.max_subdivision,
                                         highestPowerof2(int(x))) for x in pyramid_level.stack_shape])
        self.subdivisions[0] = pyramid_level.chunks
        self.levels.append(pyramid_level)

        if self.num_pyramid_levels == 1:
            logger.debug("Image won't be resampled (resolution too low)")
            return

        shape_attr = 'shape'
        if self.image.is_stream:
            shape_attr = 'stack_shape'

        for level in range(1, self.num_pyramid_levels):
            pyramid_level = ProxyChannel.PyramidLevel()
            pyramid_level.path = PathUtil.get_lsfm_image_cells_path(self.image.channel_name,
                                                                    level)
            pyramid_level.bdv_path = PathUtil.get_cells_path(0, self.num_bdv_setup, level)
            pyramid_level.id = level
            pyramid_level.spacing = self.image.voxel_size * 1. / \
                                    self.resolution_scales[level][:len(self.image.voxel_size)]
            self.levels.append(pyramid_level)
            compute_level_shape(level, shape_attr)
            self.subdivisions[level] = pyramid_level.chunks
            

    class PyramidLevel(object):
        def __repr__(self):
            return "stack shape: {} \n" \
                   "level shape: {}\n" \
                   "level chunks: {}\n" \
                   "level spacing: {}\n".format(self.stack_shape,
                                                self.shape,
                                                self.chunks,
                                                self.spacing)


def test_proxy(input_path, json_path):
    meta = dm.ImageMetaData(input_path)
    with open(json_path) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)

    # ip = StreamableTiffProxy('autofluo', ' mysz_2', "auto_8_Z%03d.ome.tif", meta, 0.3)
    #ip = image.StreamableOMEProxy('cfos', 'fos_4', None, meta, 1)
    ip = image.StreamableNrrdProxy('N11', 'test_mgre', metadata=meta,
                                    xml_file_name="test_mgre.xml", RAM_limit=2.,
                                    multifile_prefix=None)
    # ip = NonStreamableTiffProxy('cfos', 'mysz', None, meta, None)
    # ip = NiftiProxy('autofluo', 'exp2', None, meta, None)
    print(ip)
    page =0

    k = 0
    for data_chunk in ip.stream_data():
        page += data_chunk.shape[0]
        if page >= 4:
            np.save("/data/sbednarek/pnas23/lsfm_sample2.npy", data_chunk, allow_pickle=False)
            exit()
        print((data_chunk.shape))
        print(page)
        k += 1


def test_lmdhf(image_path, json_path, hdf_path, multichannel=False):
    meta = dm.ImageMetaData(image_path)
    with open(json_path) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)

    # ip = NiftiProxy('autofluo', 'exp2', None, meta, None)
    # ip = ImageProxy.get_image_proxy_class(meta)('cfos2', 'mysz', None, meta, 'single_DIR.xml', None)
    # ip = ImageProxy.get_image_proxy_class(meta)('cfos', 'mysz', 'Z%06d.tif', meta, 'check_roi_bigger.xml', 4.,
    #                                           is_multichannel=multichannel)

    #ip = image.ImageProxy.get_image_proxy_class(meta)('autofluo', 'N11', None, meta, 'N11.xml', 2.,
     #                                           is_multichannel=multichannel)
    ip = image.ImageProxy.get_image_proxy_class(meta)('M4D', 'N08', None, meta, 'N08.xml', 2.,
                                                is_multichannel=multichannel)
    #ip = image.ImageProxy.get_image_proxy_class(meta)('rd', 'N05', None, meta, 'N05.xml', 2.,
    #                                            is_multichannel=multichannel)
    
    #ip = image.ImageProxy.get_image_proxy_class(meta)('cfos', 'N71', 'Z%06d.tif', meta, 'N71.xml', 2.,
    #                                            is_multichannel=multichannel)
    
    # ip = StreamableOMEProxy('cfos', 'fos_4', None, meta, 1)

    lm_test = LightMicroscopyHDF(hdf_path)
    lm_test.write_channel(ip)
    print(lm_test)
    # lm_test.channels[0].export_image(level=4)

    # print lm_test.number_of_channels


def test_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)
    channel = lmf.get_channel('cfos2')
    channel.write_affine('to_atlas', '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/mysz_affine.txt')
    list_of_transforms = [transforms.TransformTuple(0, 'to_atlas', 'affine', True)]
    e_cmd = export.ExportCmd(channel_name='cfos2', output_path='../results/_RAS_cfos2_35um.nii.gz',
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

    list_of_transforms = [transforms.TransformTuple(0, 'new_auto_to_template', 'affine', True)]
    e_cmd = export.ExportCmd(channel_name='autofluo', output_path='../results/_cfos_15um.nii.gz',
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
    ip = image.ImageProxy.get_image_proxy_class(meta)('inverse_warp_template', 'mysz', None, meta, None, 2.,
                                                is_multichannel=True)

    lmf.write_channel(ip)


def test_cfos_apply_displacement_field(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)

    list_of_transforms = [transforms.TransformTuple(0, 'inverse_warp_template', 'df', True),
                          transforms.TransformTuple(1, 'new_auto_to_template', 'affine', True)]

    # list_of_transforms = [TransformTuple(0, 'new_auto_to_template', 'affine', True)]
    e_cmd = export.ExportCmd(channel_name='autofluo', output_path='../results/_cfos_df_af_origin_25sum.nii.gz',
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

    ip = image.ImageProxy.get_image_proxy_class(meta)('right_template', 'mysz', None, meta, 'template.xml', 2.,
                                                is_multichannel=False, is_segmentation=False)

    lm_test = LightMicroscopyHDF(hdf_path)
    lm_test.write_channel(ip)


def write_seg_to_hdf(image_path, json_path, hdf_path):
    meta = dm.ImageMetaData(image_path)
    with open(json_path) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)

    ip = image.ImageProxy.get_image_proxy_class(meta)('autofluo_segmentation', 'mysz', None, meta, 'w.xml', 2.,
                                                is_multichannel=False, is_segmentation=True)

    lm_test = LightMicroscopyHDF(hdf_path)
    lm_test.write_channel(ip)


def test_template_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)

    list_of_transforms = [transforms.TransformTuple(0, 'forward_warp', 'df', False),
                          transforms.TransformTuple(1, 'auto_to_template', 'affine', False)]

    list_of_transforms = [transforms.TransformTuple(0, 'auto_to_template', 'affine', False),
                          transforms.TransformTuple(1, 'forward_warp', 'df', False)]

    list_of_transforms = []

    e_cmd = export.ExportCmd(channel_name='right_template', output_path='../results/template_2_auto_30um.nii.gz',
                      output_resolution=None,
                      input_orientation='RAS', input_resolution_level=1,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)
    # lmf.export_image(e_cmd)

    e_cmd = export.ExportCmd(channel_name='right_template', output_path='../results/template_40um.nii.gz',
                      output_resolution=[0.04, 0.04, 0.04],
                      input_orientation='RAS', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)
    lmf.export_image(e_cmd)

    e_cmd = export.ExportCmd(channel_name='segmentation', output_path='../results/segmentation_to_auto_30um.nii.gz',
                      output_resolution=[0.030, 0.030, 0.030],
                      input_orientation='RAS', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)

    # lmf.export_image(e_cmd)


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

    ip = image.ImageProxy.get_image_proxy_class(meta)('inverse_warp_template', 'some_id', None, meta, None, 2.,
                                                is_multichannel=True)

    lmf.write_channel(ip)


def add_displacement_field(image_path, json_path, hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)

    meta = dm.ImageMetaData(image_path)
    with open(json_path) as fp:
        json_meta = json.load(fp)

    meta.update(json_meta)
    ip = image.ImageProxy.get_image_proxy_class(meta)('inverse_warp_template', 'autofluo', None, meta, None, 2.,
                                                is_multichannel=True)

    lmf.write_channel(ip)


def test_exp4_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)

    # list_of_transforms = [TransformTuple(0, 'auto_to_template', 'affine', True)]
    list_of_transforms = [transforms.TransformTuple(0, 'inverse_warp_template', 'df', True),
                          transforms.TransformTuple(1, 'auto_to_template', 'affine', True)]
    e_cmd = export.ExportCmd(channel_name='autofluo', output_path='../results/exp_4/exp4_10um_affine_df_ref.nii.gz',
                      output_resolution=[0.010, 0.010, 0.010],
                      input_orientation='PSR', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)

    lmf.export_image(e_cmd)


def test_exp4_cfos_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)

    # list_of_transforms = [TransformTuple(0, 'auto_to_template', 'affine', True)]
    list_of_transforms = [transforms.TransformTuple(0, 'cfos_to_auto', 'affine', True),
                          transforms.TransformTuple(1, 'inverse_warp_template', 'df', True),
                          transforms.TransformTuple(2, 'cfos_to_template', 'affine', True)]
    e_cmd = export.ExportCmd(channel_name='fos', output_path='../results/exp_4/level_2_fos.nii.gz',
                      output_resolution=[0.01, 0.01, 0.01],
                      input_orientation='PSR', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None)

    lmf.export_image(e_cmd)


def test_export():
    lmf = LightMicroscopyHDF('/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/check_roi_bigger_cache.h5')

    # list_of_transforms = [TransformTuple(0, 'auto_to_template', 'affine', True)]
    list_of_transforms = []
    e_cmd = export.ExportCmd(channel_name='cfos', output_path='../results/check_roi_0.nii.gz',
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=0,
                      list_of_transforms=list_of_transforms,
                      phys_origin=None, phys_size=None,
                      segmentation_name=None, region_id=None,
                      grid_size=[0, 0, 0], overlap_mm=None)

    lmf.export_image(e_cmd)


def test_region_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)
    e_cmd = export.ExportCmd(channel_name='fos', output_path='../results/exp_4/extracted_reg_382_level_0.nii.gz',
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=0,
                      list_of_transforms=[],
                      phys_origin=None, phys_size=None,
                      segmentation_name='right_signal_segmentation', region_id=382)

    lmf.export_image(e_cmd)


def test_chunk_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)
    '''e_cmd = ExportCmd(channel_name='fos', output_path='../results/exp_4/chunk_level_0.nii.gz',
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=0,
                      list_of_transforms=[],
                      phys_origin=[2.84, 1.81, -4.524], phys_size=[.5, .5, .5],
                      segmentation_name=None, region_id=None)'''

    e_cmd = export.ExportCmd(channel_name='cfos', output_path='../results/13_D_chunk_level_0_cfos.nii.gz',
                      output_resolution=None,
                      input_orientation='LAI', input_resolution_level=0,
                      list_of_transforms=[],
                      phys_origin=[2.895, 4.675, -4.1], phys_size=[0.75, .75, .75],
                      segmentation_name=None, region_id=None)

    lmf.export_image(e_cmd)


def test_grid_of_chunks_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)
    e_cmd = export.ExportCmd(channel_name='cfos', output_path='../results/grid_chunks/whatever',
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=1,
                      list_of_transforms=[],
                      phys_origin=[2.5, 1.0, -4.0], phys_size=[0.75, .75, .75],
                      segmentation_name=None, region_id=None, grid_size=[3, 3, 3], overlap_mm=0.2)

    lmf.export_image(e_cmd)


def test_whole_grid_of_chunks_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)
    e_cmd = export.ExportCmd(channel_name='fos', output_path='../results/whole_grid/overlap',
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=2,
                      list_of_transforms=[],
                      phys_origin=[2.5, 1.0, -4.0], phys_size=[2.2, 2.2, 2.2],
                      segmentation_name=None, region_id=None, grid_size=[0, 0, 0], overlap_mm=0.4)

    lmf.export_image(e_cmd)


def test_whole_grid_of_chunks_transformed_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)

    list_of_transforms = [transforms.TransformTuple(0, 'cfos_to_auto', 'affine', True),
                          transforms.TransformTuple(1, 'inverse_warp_template', 'df', True),
                          transforms.TransformTuple(2, 'cfos_to_template', 'affine', True)]

    e_cmd = export.ExportCmd(channel_name='fos', output_path='../results/whole_grid/transform',
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=2,
                      list_of_transforms=list_of_transforms,
                      # phys_origin=[2.5, 1.0, -4.0], phys_size=[2.2, 2.2, 2.2],
                      phys_origin=[7.475, 5.975, -4.475], phys_size=[1.25, 1.25, 1.25],
                      segmentation_name=None, region_id=None, grid_size=[2, 2, 2], overlap_mm=0.0)

    lmf.export_image(e_cmd)


def test_chunk_transform_export(hdf_path, channel_name='cfos', output_path='../results/check_margis_2.nii.gz'):
    lmf = LightMicroscopyHDF(hdf_path)
    list_of_transforms = [transforms.TransformTuple(0, 'cfos_to_auto', 'affine', True),
                          transforms.TransformTuple(1, 'inverse_warp_template', 'df', True),
                          transforms.TransformTuple(2, 'cfos_to_template', 'affine', True)]

    list_of_transforms = []
    e_cmd = export.ExportCmd(channel_name=channel_name,
                      output_path=output_path,
                      output_resolution=[0.006, 0.006, 0.006],
                      input_orientation='PSR', input_resolution_level=None,
                      list_of_transforms=list_of_transforms,
                      # phys_origin=[7.475, 5.975, -4.475], phys_size=[2.25, 1.25, 3.25],
                      # phys_origin=[7.475, 5.975, -4.475], phys_size=[1., 1., 1.],
                      phys_origin=[2.5, 1.0, -4.0], phys_size=[5., 5., 5.],
                      # phys_origin=[5.325, 6.6, -1.975], phys_size=[1.7, 1.5, 1.5],
                      segmentation_name=None, region_id=None)

    lmf.export_image(e_cmd)


def test_auto_chunk_transform_export(hdf_path):
    lmf = LightMicroscopyHDF(hdf_path)
    list_of_transforms = [transforms.TransformTuple(0, 'inverse_warp_template', 'df', True),
                          transforms.TransformTuple(1, 'auto_to_template', 'affine', True)]
    e_cmd = export.ExportCmd(channel_name='autofluo', output_path='../results/exp_4/auto_2.nii.gz',
                      output_resolution=None,
                      input_orientation='PSR', input_resolution_level=0,
                      list_of_transforms=list_of_transforms,
                      phys_origin=[7.475, 5.975, -4.475], phys_size=[2.25, 1.25, 3.25],
                      segmentation_name=None, region_id=None)

    lmf.export_image(e_cmd)


def test_slice_export(hdf_path, channel, output_path):
    lmf = LightMicroscopyHDF(hdf_path)
    e_cmd = export.ExportSlicesCmd(channel_name=channel, output_path=output_path, input_orientation="PSR",
                            input_resolution_level=1, slicing_range=[20, 400, 2], axis=1,
                            extract_roi=[1, 1, 500, 500],
                            ref_channel=None, ref_level=None)

    lmf.export_slices(e_cmd)


if __name__ == '__main__':
    logger.info("Testing....")


    input_path = '/media/sbednarek/4BCFEE837AD4D9DD/Diana/new_autofluo_2018_01/fos_8/Z_planes/auto_8_Z000.ome.tif'
    json_path = '/home/sbednarek/DEV/lsfm/results/fos_8_metadata.json'
    input_path = '/data/sbednarek/pnas23/200302-1-1_N09_N58181_mGRE_M4D.nhdr'
    json_path = '/data/sbednarek/pnas23/test.json'
    #input_path = '/data/sbednarek/neuroinfB_slash_data/data/lsfmpy_tutorial/lsfm_image_server/example_data/cfos/Z000000.tif'
    #json_path = '/data/sbednarek/neuroinfB_slash_data/data/lsfmpy_tutorial/lsfm_image_server/cfos.json'
    input_path = '/data/sbednarek/pnas23/200302-1-1_N09_N58181_mGRE_M4D.nii.gz'
    json_path = '/data/sbednarek/pnas23/test_nii.json'
    #input_path = '/home/sbednarek/neuroinfC/sbednarek/200302-1-1_N11_N58211NLSAM_autof_M4D.nhdr'
    #json_path = '/home/sbednarek/neuroinfC/sbednarek/200302-1-1_N11_N58211NLSAM_autof_M4D.json'
    

    #input_path = '/home/sbednarek/neuroinfC/sbednarek/200302-1-1_N07_N58211NLSAM_tdi3-color_blue_M4D.nhdr'
    #json_path = '/home/sbednarek/neuroinfC/sbednarek/200302-1-1_N07_N58211NLSAM_tdi3-color_blue_M4D.json'
    #input_path = '/home/sbednarek/neuroinfC/sbednarek/200302-1-1_N08_N58211NLSAM_md_M4D.nhdr'
    #json_path = '/home/sbednarek/neuroinfC/sbednarek/200302-1-1_N08_N58211NLSAM_md_M4D.json'
    #input_path = '/home/sbednarek/neuroinfC/sbednarek/200302-1-1_N05_N58211NLSAM_rd_M4D.nhdr'
    #json_path = '/home/sbednarek/neuroinfC/sbednarek/200302-1-1_N05_N58211NLSAM_rd_M4D.json'
    #input_path = '/home/sbednarek/neuroinfC/sbednarek/200302-1-1_N09_N58181_mGRE_M4D.nhdr'
    #json_path = '/home/sbednarek/neuroinfC/sbednarek/200302-1-1_N09_N58181_mGRE_M4D.json'
    
   

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
    auto_seg_json = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/metadata/segementation_in_autofluo_11.json'
    json_seg_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/metadata/segmentation_in_signal_11.json'
    df_forward_path = '/home/sbednarek/DEV/lsfm/resources/deformation_field_test/30_transforms/forward_warp_structural_001.nii.gz'
    df_forward_json_path = '../results/metadata/forward_warp.json'

    exp_4_json = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/resources/647-exp4.json'
    exp_4_path = '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/resources/exp4_cfos/Z000000.tif'

    #test_proxy(input_path, json_path)
    test_lmdhf(input_path, json_path, '/data/sbednarek/pnas23/tmpN71.h5')

    # test_proxy(ome_input_path, ome_json_path)
    # test_proxy(single_input_path, single_json_path)
    # test_proxy(nifti_path, nifti_json_path)

    # test_lmdhf(single_input_path, single_json_path, '../results/single_test_DIRECTION.h5')
    # test_lmdhf(ome_input_path, ome_json_path, '../results/testing_ome2.h5')
    # test_lmdhf(nifti_path, nifti_json_path, '../results/testing_nifti.h5')
    # test_lmdhf(df_path, df_json_path, '../results/df_test2.h5', True)
    # test_lmdhf(mysz_path, mysz_json_path, '../results/myszdf_RAS.h5')

    '''test_lmdhf('/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/resources/Ctrl_4_nowyROI_647_1_MMStack_Pos08.ome.tif',
               '/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/resources/check_roi.json',
               '../results/check_roi_big_chunks.h5')'''

    # test_export()

    # test_export('../results/myszdf_RAS.h5')

    # test_export('../results/single_test.h5')

    # test_cfos_export('/home/sbednarek/DEV/lsfm/results/cfos.h5')

    # test_cfos_add_displacement_field(df_path, df_json_path, '/home/sbednarek/DEV/lsfm/results/cfos.h5')
    # test_cfos_apply_displacement_field('/home/sbednarek/DEV/lsfm/results/cfos.h5')

    # write_template_to_hdf(template_path, json_template_path, '../results/template.h5')

    # test_template_export('../results/template.h5')

    # test_lmdhf(df_forward_path, df_forward_json_path, '../results/template.h5', True)

    # write_cfos_affine('../results/template.h5')

    # test_template_export('../results/template.h5')

    # test_lmdhf(exp_4_path, exp_4_json, '../results/exp_4_bigger_close.h5', False)

    # write_exp4_affine('../results/exp_4_new.h5')

    # write_exp4_affine('/home/sbednarek/DEV/lsfm/results/experimental_4.h5')

    # add_displacement_field(df_path, df_json_path, '../results/exp_4_new.h5')

    # add_displacement_field(df_path, df_json_path, '/home/sbednarek/DEV/lsfm/results/experimental_4.h5')

    # test_exp4_export('../results/exp_4_new.h5')

    # test_exp4_cfos_export('../results/exp_4_new.h5')

    # write_seg_to_hdf(seg_path, json_seg_path, '../results/exp_4_new.h5')

    #test_region_export('../results/exp_4_new.h5')

    # test_chunk_export('../results/exp_4_new.h5')

    #test_chunk_transform_export('../results/exp_4_new.h5')
    # test_auto_chunk_transform_export('../results/exp_4_new.h5')

    # test_chunk_transform_export('/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4_new.h5',
    #                            'fos', output_path='/home/sbednarek/DEV/lsfm_schema/lsfm_image_server/results/exp_4/check_margis_1.nii.gz')

    # test_chunk_transform_export('/mnt/nicl/home/pmajka/results/experimental_5.h5',
    #                            output_path='../results/experimental_5/in_cfos_whole_amygdala_native_cfos.nii.gz')

    # write_transforms_signal('/mnt/nicl/home/pmajka/results/experimental_5.h5', 'cfos',
    #                        '/home/sbednarek/DEV/lsfm_pipeline/000012/lsfm_image_server/30_transforms/')

    # test_chunk_transform_export('/mnt/nicl/home/pmajka/results/experimental_3.h5',
    #                            output_path='../results/experimental_3/in_cfos_whole_amygdala_native_cfos.nii.gz')

    # test_chunk_transform_export('/mnt/nicl/home/pmajka/results/experimental_2.h5',
    #                           output_path='../results/experimental_2/testing_chunk_4.nii.gz')

    # test_chunk_transform_export('/mnt/nicl/home/pmajka/results/control_2.h5',
    #                            output_path='../results/control_2/in_cfos_whole_amygdala_native_cfos.nii.gz')


    # test_slice_export('../results/exp_4_new.h5', 'autofluo', '../results/slices/img_%04d.tif')

    # write_seg_to_hdf(auto_seg_path, auto_seg_json, '../results/exp_4_new.h5')


    # test_chunk_export('/mnt/nicl/home/sbednarek/lsfm/000013.h5')

    # test_grid_of_chunks_export('/mnt/nicl/home/pmajka/results/experimental_2.h5')
    #test_whole_grid_of_chunks_export('../results/exp_4_new.h5')
    # test_whole_grid_of_chunks_transformed_export('../results/exp_4_new.h5')
