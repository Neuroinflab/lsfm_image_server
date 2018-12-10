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
import fire
import json
import logging
import dump_metadata as dm

from lmdf import LightMicroscopyHDF
from transforms import TransformTuple
from image import ImageProxy
from export import ExportSlicesCmd, ExportCmd


class LmdfIO_CLI(object):
    def __init__(self):
        self.hdf_path = None
        self.logger = logging.getLogger("LmdfIO_CLI")
        self.list_of_transforms = []

    def write(self, hdf_path, channel_name, image_path, bdv_xml, metadata_path=None,
              slab_memory_size=2., file_name_format=None, is_multichannel=False, is_segmentation=False):
        """
        Write image data, affine transformations, displacement fields and segmentation
        into HDF5 file.

        :param hdf_path: path to the HDF5 file
        :param channel_name: name of the new channel
        :param image_path: path to the source file or first file in case of a Z-stack
        :param metadata_path: path to the json metadata file
        :param bdv_xml: name of the accompanying xml file for BigDataViewer compatibility (not a path)
        :param slab_memory_size: how much physical memory can be used (default 2GB)
        :param file_name_format: define file name format for image series (e.g. Z%06d.tif)
        :param is_multichannel: flag for mutlichannel files (e.g. displacement fields)
        :param is_segmentation: flag for segmentation files
        :return: None
        """
        meta = dm.ImageMetaData(image_path)
        self.hdf_path = hdf_path
        if not metadata_path:
            meta = dm.getImageMetaDataClass(image_path)(image_path)
        else:
            try:

                with open(metadata_path) as fp:
                    json_meta = json.load(fp)
                meta.update(json_meta)
            except dm.InvalidImageSourceException:
                self.logger.error("Image type not recognized or file does not exist", exc_info=False)
            except IOError:
                self.logger.error("File could not be opened", exc_info=True)

        self.logger.debug(meta)

        ip = ImageProxy.get_image_proxy_class(meta)(channel_name, 'whatever', file_name_format, meta, bdv_xml,
                                                    slab_memory_size, 'RAS', is_multichannel, is_segmentation)

        lmf = LightMicroscopyHDF(self.hdf_path)
        lmf.write_channel(ip)
        lmf.close()

    def write_affine(self, hdf_path, channel_name, affine_name, affine_path):
        """
        Reads ITK produced affine file and writes it in nifti and itk compatible formats to a given channel.

        :param hdf_path: path to the HDF5 file
        :param channel_name: Identifies the channel to which affine will be written.
        :param affine_name: Name of the affine, must be unique for channel.
        :param affine_path: Path to a .txt file containing the affine parameters.
        :return: None
        """
        self.hdf_path = hdf_path
        lmf = LightMicroscopyHDF(self.hdf_path)
        channel = lmf.get_channel(channel_name)
        channel.write_affine(affine_name, affine_path)
        lmf.close()

    def export_slices(self, hdf_path, channel_name, input_orientation,
                      input_resolution_level, slicing_range, axis, output_path=None, extract_roi=None):
        """

        :param hdf_path: path to the source HDF5 file
        :param channel_name: channel from which to export slices
        :param output_path: formatted output path
        :param input_orientation: the internal (anatomical) orientation of the image (def. RAS)
        :param input_resolution_level: export data from a given, pre-resampled level of resolution pyramid
        :param slicing_range: start,stop,step or None (None - will export everything)
        :param axis: the z axis
        :param extract_roi: ox,oy,sx,sy or None (None - will export with maximum height and width)
        :return:
        """

        if not os.path.isfile(hdf_path):
            sys.exit("File not found!")

        self.hdf_path = hdf_path
        lmf = LightMicroscopyHDF(self.hdf_path, access_mode='a')

        export_cmd = ExportSlicesCmd(channel_name=channel_name, output_path=output_path,
                                     input_orientation=input_orientation,
                                     input_resolution_level=input_resolution_level,
                                     slicing_range=slicing_range, axis=axis,
                                     extract_roi=extract_roi,
                                     ref_channel=None, ref_level=None)

        return lmf.export_slices(export_cmd)
        # lmf.close()

    def export(self, hdf_path, channel_name, output_path=None, output_resolution=None, input_orientation='RAS',
               input_resolution_level=None, phys_origin=None, phys_size=None,
               segmentation_name=None, region_id=None, grid_size=None, overlap_mm=None):
        """
        Export image data at desired resolution, volume, and physical coordinates, with optional transformation
        to a reference space. Save to a file or to memory for further processing.

        :param hdf_path: path to source HDF5 file
        :param channel_name: channel from which to export data
        :param output_path: where to write output data
        :param output_resolution: resolution to which requested data subset will be resampled (in mm!)
        :param input_orientation: the internal (anatomical) orientation of the image (def. RAS)
        :param input_resolution_level: export data from a given, pre-resampled level of resolution pyramid
        :param phys_origin: physical origin (in RAS) of the image subset in output physical space (in mm!)
        :param phys_size:  physical size of the image subset (in mm!)
        :param segmentation_name: segmentation based on which a subset region of image volume will be exported
        :param region_id: id of the region in segmentation specified in segmentation_name parameter
        :param grid_size: x,y,z -> number of chunks in each direction(RAS), if 0,0,0 whole chunked image will be exported
        :param overlap_mm: overlap for chunk export (in mm!)
        :return:
        """

        if not os.path.isfile(hdf_path):
            sys.exit("File not found!")

        self.hdf_path = hdf_path
        lmf = LightMicroscopyHDF(self.hdf_path, 'a')
        self.logger.debug("resolution: {}".format(output_resolution))

        export_cmd = ExportCmd(channel_name=channel_name, output_path=output_path,
                               output_resolution=output_resolution,
                               input_orientation=input_orientation, input_resolution_level=input_resolution_level,
                               list_of_transforms=self.list_of_transforms,
                               phys_origin=phys_origin, phys_size=phys_size,
                               segmentation_name=segmentation_name, region_id=region_id,
                               grid_size=grid_size, overlap_mm=overlap_mm)

        return lmf.export_image(export_cmd)

    def add_transform(self, ordn, name, transform_type, invert=True):
        """
        For export purposes, define the transforms to be used on source data.
         Transforms include affines and displacement
        fields, which must be already present in source HDF5 file. Use this command multiple times to chain a series
        of transforms.

        :param ordn: int, defines the order in which transforms will be executed
        :param name: str, name of the transform
        :param transform_type: str, affine or df
        :param invert: boolean, should the inverse of the transform be used
        :return: None
        """
        transform = TransformTuple(ordn, name, transform_type, invert)
        self.list_of_transforms.append(transform)
        self.logger.debug(self.list_of_transforms)
        return self

    def info(self, hdf_path, channel_name=None):
        """
        Print available information about HDF5 file, or specified channel.

        :param hdf_path: HDF5 file of interest
        :param channel_name: Print information only concerning specified channel.
        :return: None
        """

        if not os.path.isfile(hdf_path):
            sys.exit("File not found!")

        self.hdf_path = hdf_path
        lmf = LightMicroscopyHDF(self.hdf_path)

        if channel_name is not None:
            print lmf.get_channel(channel_name)
        else:
            print lmf

        lmf.close()

    def show_log(self, hdf_path):

        self.hdf_path = hdf_path
        lmf = LightMicroscopyHDF(self.hdf_path)
        logs = lmf.get_logs()
        lmf.close()

        for k, v in logs:
            print("------------ {} ------------\n".format(k))
            print(v)
            print("\n########################################################\n\n\n")

    def add_metadata(self, hdf_path, m_key, m_value):
        """
        Add arbitrary metadata as key value pair.

        :param hdf_path: path to hdf5 file
        :param m_key: metadata key
        :param m_value: metadata value
        :return: None
        """
        self.hdf_path = hdf_path
        lmf = LightMicroscopyHDF(self.hdf_path, 'a')
        lmf.add_metadata(m_key, m_value)
        lmf.close()

    def __repr__(self):
        return "Command line interface for the Light Microscopy HDF5 Software"


if __name__ == '__main__':
    fire.Fire(LmdfIO_CLI)
