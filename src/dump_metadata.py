#!/usr/bin/env python
# -*- coding: utf-8 -*

"""
This is a framework for parsing metadata from image files,
printing them on screen and dumping to json file for further processing / usage.

Framework can be extended by adding classes for particular file types.
Each class is completely responsible for parsing its files appropriately.

Example usage:
dump_metadata.py --input-file=example.nii.gz --output-file=metadata.json
"""

import sys
import json
import argparse

import numpy as np
import nibabel as nib
import tifffile as tf

from os import path
from lxml import objectify
from utils import InputImageType

VALID_EXTENSIONS = ["tif",
                    "tiff",
                    "nii",
                    "nii.gz",
                    "ome.tif",
                    "ome.tiff"]

FILE_TYPE = {'tif': InputImageType.Tiff,
             'tiff': InputImageType.Tiff,
             'ome.tif': InputImageType.OmeTiff,
             'ome.tiff': InputImageType.OmeTiff,
             'nii': InputImageType.Nifti,
             'nii.gz': InputImageType.Nifti}


def is_image_source_valid(file_path):
    """
    Checks if file path provided by user is a valid image file from
    which metadata can be harvested

    :param file_path: str
                path to image file
    :return: bool
    """
    if not path.isfile(file_path):
        return False
    _, fname = path.split(file_path)
    ext = fname.split('.', 1)[1]
    return ext.lower() in VALID_EXTENSIONS


class InvalidImageSourceException(Exception):
    pass


class ImageMetaData(dict):
    """
    Stores metadata parsed from image files
    """

    def __init__(self, file_path):
        dict.__init__(self)
        if not is_image_source_valid(file_path):
            raise InvalidImageSourceException
        self['file_path'] = file_path

    def __repr__(self):
        return "\n".join(["%s = %s" % (k, v) for k, v in self.items()])

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)


class TiffImageMetaData(ImageMetaData):
    """
    Helper class delegating tiff parsing to specialized tiff metadata classes
    """

    def __setitem__(self, key, value):
        if key == "file_path" and value:
            tiff_img = tf.TiffFile(value)
            tiff_page = tiff_img.pages[0]

            if tiff_page.is_imagej:
                print('is imageJ')
                print value
                # print ImageJTiffImageMetaData(value)
                self = ImageJTiffImageMetaData(value)

            elif tiff_page.is_ome:
                return OmeTiffImageMetaData

            else:
                return ImageMetaData


class OmeTiffImageMetaData(ImageMetaData):
    """
    Dict based class for parsing and storing Ome tiff file metadata
    """

    tagInfoMapfloat = {"voxel_size_x": "PhysicalSizeX",
                       "voxel_size_y": "PhysicalSizeY",
                       "voxel_size_z": "PhysicalSizeZ"}
    tagInfoMapint = {"image_size_x": "SizeX",
                     "image_size_y": "SizeY",
                     "image_size_z": "SizeZ"}

    def __parse(self, file_path):
        tiff_file = tf.TiffFile(file_path)
        page = tiff_file.pages[0]
        ome_xml = objectify.XML(page.tags['image_description'].value)
        pixel_info = ome_xml['Image']['Pixels']  # AttributeError
        for tag, ome_key in self.tagInfoMapfloat.items():
            self[tag] = np.float64(pixel_info.attrib[ome_key])
        for tag, ome_key in self.tagInfoMapint.items():
            self[tag] = np.int(pixel_info.attrib[ome_key])
        self["bits_per_sample"] = page.bits_per_sample
        self["file_type"] = self.__class__.__name__[:-13]

    def __setitem__(self, key, value):
        if key == "file_path" and value:
            self.__parse(value)
        ImageMetaData.__setitem__(self, key, value)


class ImageJTiffImageMetaData(ImageMetaData):
    """
    Dict based class for parsing and storing ImageJ tiff file metadata
    """

    def __parse(self, file_path):
        tiff_file = tf.TiffFile(file_path)
        page = tiff_file.pages[0]

        shape = page.shape
        if len(shape) == 3:

            Z, Y, X = page.shape
            self["image_size_x"] = np.int(X)
            self["image_size_y"] = np.int(Y)
            self["image_size_z"] = np.int(Z)

            x_denominator, x_nominator = page.tags['x_resolution'].value
            voxel_size_x = x_nominator / x_denominator

            y_denominator, y_nominator = page.tags['y_resolution'].value
            voxel_size_y = y_nominator / y_denominator

            self["voxel_size_x"] = np.float64(voxel_size_x)
            self["voxel_size_y"] = np.float64(voxel_size_y)
            self["voxel_size_z"] = np.float64(page.imagej_tags['spacing'])
            self["bits_per_sample"] = page.bits_per_sample
            self["file_type"] = self.__class__.__name__[:-13]

        elif len(shape) == 2:
            print 'cannot read all necessary information, please fill in missing values in .json file'
            Y, X = page.shape
            self["image_size_x"] = np.int(X)
            self["image_size_y"] = np.int(Y)
            self["image_size_z"] = 0
            self["voxel_size_x"] = 0
            self["voxel_size_y"] = 0
            self["voxel_size_z"] = 0
            self["bits_per_sample"] = page.bits_per_sample
            self["file_type"] = "MultipleTiff"

    def __setitem__(self, key, value):
        if key == "file_path" and value:
            self.__parse(value)
        ImageMetaData.__setitem__(self, key, value)


class NiftiImageMetaData(ImageMetaData):
    """
    Dict based class for parsing and storing nifti file metadata
    """

    def __parse(self, file_path):
        nii_file = nib.load(file_path)
        self["image_size_x"] = np.int(nii_file.header['dim'][1])
        self["image_size_y"] = np.int(nii_file.header['dim'][2])
        self["image_size_z"] = np.int(nii_file.header['dim'][3])
        self["voxel_size_x"] = np.float64(nii_file.header['pixdim'][1])
        self["voxel_size_y"] = np.float64(nii_file.header['pixdim'][2])
        self["voxel_size_z"] = np.float64(nii_file.header['pixdim'][3])
        self["bits_per_sample"] = np.int(nii_file.header['bitpix'])
        self["file_type"] = self.__class__.__name__[:-13]

    def __setitem__(self, key, value):
        if key == "file_path" and value:
            self.__parse(value)
        ImageMetaData.__setitem__(self, key, value)


def getImageMetaDataClass(file_path, module=sys.modules[ImageMetaData.__module__]):
    """
    Returns appropriate metadata class based on file extension.

    :param file_path: str
            path to the image file
    :param module: module
            where to look for class definition
    :return: class

    """
    if not is_image_source_valid(file_path):
        raise InvalidImageSourceException

    _, filename = path.split(file_path)
    ext = filename.split('.', 1)[1]
    file_type = FILE_TYPE[ext]

    # various tiff file types cannot be differentiated
    # based solely on file extension
    if file_type == InputImageType.Tiff:
        pages = tf.TiffFile(file_path).pages
        page = pages[0]
        if page.is_ome:
            if len(pages) > 1: # it's streamable
                file_type = InputImageType.OmeTiff
            else:
                file_type = InputImageType.ImageJTiff
        elif page.is_imagej:
            file_type = InputImageType.ImageJTiff
        else:
            # best guess we can make - missing info will be provided by the user
            file_type = InputImageType.ImageJTiff

    subclass = "%sImageMetaData" % file_type.name
    return hasattr(module, subclass) and getattr(module, subclass) or ImageMetaData


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read relevant metadata from image \
                                                  and dump them in json format',
                                     prog='dump_metadata')

    parser.add_argument('--input-file',
                        required=True,
                        help="Path to the image file",
                        type=str)

    parser.add_argument('--output-file',
                        required=False,
                        default=None,
                        help="Path to output json file",
                        type=str)

    args = parser.parse_args()
    file_path = args.input_file
    meta_data = getImageMetaDataClass(file_path)(file_path)

    print meta_data
    if args.output_file is not None:
        with open(args.output_file, 'w') as fp:
            json.dump(meta_data, fp)
