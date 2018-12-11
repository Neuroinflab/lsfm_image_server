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

from lxml import etree
from utils import PathUtil

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)


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

        '''LOWEST_RES = 128
        subdiv_16_16_16 = np.array([16, 16, 16])
        subdiv_32_32_4 = np.array([[4, 32, 32],
                                   [32, 4, 32],
                                   [32, 32, 4]])

        subdiv_128_128_128 = np.array([128, 128, 128])
        subdiv_256_256_32 = np.array([[32, 256, 256],
                                      [256, 32, 256],
                                      [256, 256, 32]])'''

        LOWEST_RES = 256
        subdiv_16_16_16 = np.array([16, 16, 16])
        subdiv_32_32_4 = np.array([[4, 32, 32],
                                   [32, 4, 32],
                                   [32, 32, 4]])

        subdiv_128_128_128 = np.array([256, 256, 256])
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
                if shape_min > LOWEST_RES:
                    subdivisions.append(subdiv_256_256_32[voxel_scale.argmax()])
                else:
                    subdivisions.append(subdiv_32_32_4[voxel_scale.argmax()])
            else:
                if shape_min > LOWEST_RES:
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
        """
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
        """
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