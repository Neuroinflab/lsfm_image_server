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

import logging

import numpy as np
import SimpleITK as sitk
import nibabel as nib
import constants as const

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def reverse_axes(data, reverse=[False, False, False]):
    if reverse[0]:
        data = data[::-1, :, :]
    if reverse[1]:
        data = data[:, ::-1, :]
    if reverse[2]:
        data = data[:, :, ::-1]
    return data


def get_bounding_box(segmentation, reg_id):
    label_stats_filter = sitk.LabelStatisticsImageFilter()
    label_stats_filter.Execute(segmentation, segmentation)
    x1, x2, y1, y2, z1, z2 = label_stats_filter.GetBoundingBox(reg_id)
    start_point = (x1, y1, z1)
    end_point = (x2, y2, z2)

    return start_point, end_point


def get_affines(voxel_size, origin):
    affine_voxel_to_physical = sitk.AffineTransform(3)
    params = tuple(np.hstack([np.diag(voxel_size).ravel(), origin]))
    affine_voxel_to_physical.SetParameters(params)
    affine_physical_to_voxel = affine_voxel_to_physical.GetInverse()

    return affine_voxel_to_physical, affine_physical_to_voxel


def create_dummy_image(size, direction, spacing, origin, pixel_id):
    logger.debug("size: {}".format(size))
    logger.debug("direction: {}".format(direction))
    logger.debug("spacing: {}".format(spacing))
    logger.debug("origin: {}".format(origin))
    logger.debug("pixel_id: {}".format(pixel_id))
    img = sitk.Image(int(size[0]), int(size[1]), int(size[2]), int(pixel_id))
    img.SetDirection(direction)
    img.SetSpacing(spacing)
    img.SetOrigin(origin)

    return img


def trim_image_to_size(img, size):
    roi_filter = sitk.RegionOfInterestImageFilter()
    roi_filter.SetIndex((0, 0, 0))
    roi_filter.SetSize(size)
    out_img = roi_filter.Execute(img)
    return out_img


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
        logger.debug("margin for segmentation reduced at index, YOU ARE ON THE EDGE")
        m_index[np.argwhere(m_index < 0)] = 0
        reduce_size = np.array(index) - m_index
        m_size -= reduce_size

    if np.any((m_index + m_size) > img_size):
        logger.debug("margin for segmentation reduced at upper bound, YOU ARE ON THE EDGE")
        out_of_bounds = np.argwhere((m_index + m_size) > img_size)
        m_size[out_of_bounds] = img_size[out_of_bounds]

    index = list(m_index)
    size = list(m_size)

    roi_filter.SetIndex(index)
    roi_filter.SetSize(size)
    out_img = roi_filter.Execute(img)

    return out_img


def extract_labeled_region(reg_id, img, labels, margin=0):
    """
    Extracts region of interest from img identified by region_id, which corresponds to
    label intensity value in labels'
    :param reg_id: unit16 indicating designated region
    :param img: sitk.Image
    :param labels: sitk.Image of segmentation
    :param margin:
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
    img = trim_image_to_label_region(img, labels, reg_id, margin=margin)
    # trim segmentation to some reasonable bounding box
    labels = trim_image_to_label_region(labels, labels, reg_id, margin=margin)

    threshold_filter.SetLower(reg_id)
    threshold_filter.SetUpper(reg_id)
    labels = threshold_filter.Execute(labels)

    rescale_filter.SetOutputMaximum(1)
    labels = rescale_filter.Execute(labels)
    labels = sitk.Cast(labels, img.GetPixelID())

    region = multiply_filter.Execute(img, labels)

    return region


def compute_offset(sitk_image, spacing):
    offset = np.array(sitk_image.GetDirection()).reshape(3, 3) * spacing

    for i in range(offset.shape[0]):
        offset[i] *= 0.5

    return offset


def get_size_origin(img, transform, spacing):
    """
    Recompute the size and origin for image transformation.
    For RAS orientation only.

    :type img: sitk.Image
    :type transform: sitk.Transform
    :type spacing: np.array
    :param img: input image for transformation
    :param transform: itk (affine) transform
    :param spacing: spacing of the output image
    :return:
    """
    size = img.GetSize()
    # find all corners of an image in coordinate space
    index = []
    for i in xrange(8):
        for j in xrange(3):
            if (i >> j) % 2:
                index.append(size[j] - 1)
            else:
                index.append(0)

    index = np.array(index)
    index = index.reshape(-1, 3)

    # find all corners of the image in physical output space
    points = []
    corners = []
    for p in index:
        point = img.TransformIndexToPhysicalPoint(p)
        points.append(point)
        corners.append(transform.TransformPoint(point))

    _max = np.ones(3) * np.finfo(np.float32).min
    _min = np.ones(3) * np.finfo(np.float32).max

    for c in corners:
        for j in xrange(3):
            _min[j] = c[j] if c[j] < _min[j] else _min[j]
            _max[j] = c[j] if c[j] > _max[j] else _max[j]

    new_size = np.round(((_max - _min) / spacing) + 1)
    _min[2] = _max[2]  # RAS direction
    new_origin = _min

    return new_size, new_origin


def apply_transformations(meta_image, composite_transform):

    input_image = meta_image.get_sitk_image()
    output_size, output_origin = get_size_origin(input_image,
                                                 composite_transform.affine_composite.GetInverse(),
                                                 np.array(input_image.GetSpacing()))

    logger.debug("output_size: {}".format(output_size))
    logger.debug("output_origin: {}".format(output_origin))
    logger.debug("input spacing: {}".format(input_image.GetSpacing()))

    resampler = sitk.ResampleImageFilter()
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(0)

    resampler.SetOutputOrigin(output_origin)
    resampler.SetOutputSpacing(input_image.GetSpacing())
    resampler.SetSize(map(int, output_size))
    resampler.SetOutputDirection(input_image.GetDirection())

    resampler.SetTransform(composite_transform.composite)

    out = resampler.Execute(input_image)
    return out


def resample_sitk(sitk_image, scale_factors, interpolation_type):

    scale_factors = np.array(scale_factors, dtype=np.float64)

    input_size = np.array(sitk_image.GetSize(), dtype=np.float64)
    output_size = np.int16(np.floor(input_size * scale_factors))
    output_size[0] = int(round(input_size[0] * scale_factors[0]))

    voxel_size = np.array(sitk_image.GetSpacing())

    output_spacing = voxel_size * input_size
    output_spacing = output_spacing / output_size

    input_offset = compute_offset(sitk_image, voxel_size)
    output_offset = compute_offset(sitk_image, output_spacing)

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

    input_offset = compute_offset(sitk_image, meta_image.voxel_size)
    output_offset = compute_offset(sitk_image, output_spacing)

    output_origin = meta_image.origin - input_offset + output_offset
    output_origin = np.diag(output_origin)

    output_affine = nib.affines.from_matvec(np.diag(output_spacing) * const.SIGN_RAS_A,
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

    logger.debug("spacing after resampling: {}".format(output_image.GetSpacing()))

    output_image = sitk.GetArrayFromImage(output_image)
    output_image = output_image.transpose(meta_image.transpose)

    out_meta_image = MetaImage(output_image, output_affine, meta_image.pixel_type,
                               meta_image.direction, meta_image.is_segmentation)

    return out_meta_image


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
        logger.debug("copying data")
        img_data = self.data.copy()
        transpose_data = [2, 1, 0]
        if self.sitk_data.component_type == 'vector':
            if len(img_data.shape) == 4:
                # img_data = img_data[:, :, :, 0, :]
                transpose_data = np.arange(len(img_data.shape))
                transpose_data[:3] = [2, 1, 0]
        logger.debug("transposing data")
        img_data = img_data.transpose(self.transpose)
        logger.debug("converting to sitk image array")
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
        transpose, flip = get_transpose(source_orientation)
        logger.debug("transpose: {}\n flip: {}".format(transpose, flip))
        _transpose = np.arange(len(self.data_shape))
        _transpose[:3] = transpose
        self.data = self.data.transpose(_transpose)
        self.reverse_axes(flip < 0)

        self.voxel_size = self.voxel_size[transpose]
        self.origin = np.abs(self.origin[transpose]) * const.SIGN_LPI_T
        affine = np.abs(np.diag(self.voxel_size)) * const.SIGN_LPI_A
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
        self.direction = const.DIRECTION_LPI
        self.affine[:3, :3] = np.abs(self.affine[:3, :3]) * const.SIGN_LPI_A
        self.affine[:3, 3] = np.abs(self.affine[:3, 3]) * const.SIGN_LPI_T

    def set_direction_RAS(self):

        self.direction = const.DIRECTION_RAS
        self.affine[:3, :3] = np.abs(self.affine[:3, :3]) * const.SIGN_RAS_A
        if not self.is_nifti:
            self.affine[:3, 3] = np.abs(self.affine[:3, 3]) * const.SIGN_RAS_T
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