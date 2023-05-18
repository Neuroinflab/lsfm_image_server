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
from . import processing

import numpy as np
import nibabel as nib
from . import constants as const

from .transforms import TransposeTuple, CompositeTransform

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)


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
        self.affine_vp, self.affine_pv = processing.get_affines(self.voxel_size,
                                                                    self.origin)
        self.ct = CompositeTransform(self.export_cmd.list_of_transforms)

    def get_axes_transpose(self, ndim=3):
        transpose, flip = np.arange(ndim), np.ones(ndim)
        inv_transpose, inv_flip = np.arange(ndim), np.ones(ndim)
        transpose[:3], flip[:3] = processing.get_transpose(self.export_cmd.input_orientation)
        inv_transpose[:3] = np.argsort(transpose[:3])
        inv_flip[:3] = flip[inv_transpose[:3]] < 0

        return transpose[:3], flip[:3], inv_transpose[:3], inv_flip[:3], transpose, flip, inv_transpose, inv_flip

    def _compute_indices(self):

        self.start_os_mm = self.export_cmd.phys_origin
        self.chunk_size_in_voxels = (self.export_cmd.phys_size // self.voxel_size).astype(int)
        self.end_os_mm = np.array(self.export_cmd.phys_origin) + \
                         np.array(self.export_cmd.phys_size) * np.array([1, 1, -1.])

        self.start_is_mm = np.abs(self.ct.composite.TransformPoint(self.start_os_mm)) * np.array([1, 1, -1.])
        self.start_index = np.abs(np.around(self.affine_pv.TransformPoint(self.start_is_mm)).astype(int))

        self.end_index = self.start_index + self.chunk_size_in_voxels
        self.start_index = np.array(self.start_index)[self.transpose_tuple.inv_transpose]
        self.end_index = np.array(self.end_index)[self.transpose_tuple.inv_transpose]

    def _prepare_chunk(self):
        self.chunk = self.pyramid_level.get_level_chunk_data(self.start_index, self.end_index,
                                                             self.transpose_tuple.inv_flip)

    def _transpose_chunk(self):

        self.chunk = self.chunk.transpose(self.transpose_tuple.ltranspose)
        self.chunk = processing.reverse_axes(self.chunk, self.transpose_tuple.flip < 0)

    def _resample_chunk(self):

        logger.debug("_resample_chunk START_IS_MM {}".format(self.start_is_mm))

        self.chunk_affine = nib.affines.from_matvec(np.diag(self.voxel_size) * const.SIGN_RAS_A, self.start_is_mm)
        self.meta_image = processing.MetaImage(self.chunk, self.chunk_affine,
                                               self.pyramid_level.pixel_type,
                                               const.DIRECTION_RAS,
                                               self.pyramid_level.is_segmentation,
                                               self.pyramid_level.is_nifti)

        if self.export_cmd.resample_image:
            self.scale_factors = self.voxel_size / self.export_cmd.output_resolution
        else:
            self.scale_factors = np.array([1., 1., 1.])

        self.meta_image = processing.resample_image(self.meta_image, self.scale_factors)

    def _transform_chunk(self):

        self.chunk = processing.apply_transformations(self.meta_image, self.ct)

    def _resize_chunk(self):
        pass

    def process(self):
        self._compute_indices()
        self._prepare_chunk()
        self._transpose_chunk()
        self._resample_chunk()
        self._transform_chunk()
        self._resize_chunk()

        return self.chunk


class ImageRegionExporter(ImageExporter):

    def _compute_indices(self):
        self.start, self.end = processing.get_bounding_box(self.export_cmd.segmentation,
                                                           self.export_cmd.region_id)

        self.start_is_mm = self.export_cmd.segmentation.TransformIndexToPhysicalPoint(self.start)
        self.end_is_mm = self.export_cmd.segmentation.TransformIndexToPhysicalPoint(self.end)
        if self.export_cmd.region_size:
            logger.debug("start_is before margin: {}".format(self.start_is_mm))
            logger.debug("end_is_mm before margin: {}".format(self.end_is_mm))
            seg_spacing = np.array(self.export_cmd.segmentation.GetSpacing())
            region_size = np.array(self.end) - np.array(self.start)
            region_size_mm = region_size * seg_spacing
            output_size_mm = self.export_cmd.output_resolution * self.export_cmd.region_size
            margin_mm = (output_size_mm - region_size_mm) / 2.
            margin_mm += self.export_cmd.output_resolution
            if np.any(margin_mm < 0):
                logger.info("WARNING, ROI TO SMALL")
            self.start_is_mm -= (margin_mm * np.array([1., 1., -1.]))
            self.end_is_mm += (margin_mm * np.array([1., 1., -1.]))
            logger.debug("start_is after margin: {}".format(self.start_is_mm))
            logger.debug("end_is_mm after margin: {}".format(self.end_is_mm))
        self.start_index = np.abs(np.around(self.affine_pv.TransformPoint(self.start_is_mm)).astype(int))
        self.end_index = np.abs(np.around(self.affine_pv.TransformPoint(self.end_is_mm)).astype(int))
        self.start_index = np.array(self.start_index)[self.transpose_tuple.inv_transpose]
        self.end_index = np.array(self.end_index)[self.transpose_tuple.inv_transpose]

    def _resize_chunk(self):
        if self.export_cmd.region_size:
            if not np.all(np.array(self.chunk.GetSize()) == self.export_cmd.region_size):
                logger.debug("Fixing size..")
                self.chunk = processing.trim_image_to_size(self.chunk,
                                                           self.export_cmd.region_size)


class ImageWholeExporter(ImageExporter):

    def _compute_indices(self):
        if self.export_cmd.phys_origin is not None:
            self.start_os_mm = self.export_cmd.phys_origin
            self.start_is_mm = self.ct.composite.TransformPoint(self.start_os_mm)
        else:
            self.start_is_mm = self.origin
            self.start_os_mm = self.ct.affine_composite.GetInverse().TransformPoint(self.start_is_mm)

        self.start_index = np.abs(np.around(self.affine_pv.TransformPoint(self.start_is_mm)).astype(int))
        self.chunk_size_in_voxels = self.shape
        self.end_index = self.start_index + self.chunk_size_in_voxels
        self.start_index = np.array(self.start_index)[self.transpose_tuple.inv_transpose]
        self.end_index = np.array(self.end_index)[self.transpose_tuple.inv_transpose]


class ExportCmd(object):
    def __init__(self, channel_name, output_path, output_resolution, input_orientation,
                 input_resolution_level, list_of_transforms, phys_origin, phys_size,
                 segmentation_name, region_id, grid_size, overlap_mm, region_size):
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
        self.region_size = region_size
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
