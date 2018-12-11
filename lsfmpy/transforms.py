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
import nibabel as nib
import SimpleITK as sitk

import constants as const

from collections import namedtuple

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)

TransformTuple = namedtuple('TransformTuple', 'order name type invert')
TransposeTuple = namedtuple('TransposeTuple', 'transpose, flip, inv_transpose, inv_flip,'
                            'ltranspose, lflip, linv_transpose, linv_flip')


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
        dir_matrix *= const.RAS_INVERSE_DIRECTION[:, np.newaxis]
        return dir_matrix

    @staticmethod
    def compute_offset(matrix, center, translation):
        """
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
        """

        offset = translation + center
        offset -= matrix.dot(center)
        return offset

    @staticmethod
    def get_affine_matrix(fname):
        """
        Reads from a file an itk matrix defining affine transformation and converts
        it to a nifti-compatible format

        Parameters
        ----------
        param fname: file name (txt)
        type fname: str

        Returns
        -------
        rtype: ndarray of float (4, 4)
        """
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
        if len(self.transforms) == 0:
            logger.debug("creating identity transform")
            return sitk.Transform(sitk.AffineTransform(3))
        self.sort_transforms()
        ct = sitk.Transform(self.transforms[0].transform)
        for transform in self.transforms[1:]:
            ct.AddTransform(transform.transform)

        return ct

    @property
    def affine_composite(self):
        if len(self.transforms) == 0:
            logger.debug("creating identity transform")
            return sitk.Transform(sitk.AffineTransform(3))
        act = None
        self.sort_transforms()
        for transform in self.transforms:
            if transform.transform_type == 'affine':
                if act is None:
                    act = sitk.Transform(transform.transform)
                else:
                    act.AddTransform(transform.transform)

        return act
