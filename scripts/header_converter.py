#!/usr/bin/env python
# -*- coding: utf-8 -*

import argparse
import numpy as np
import nibabel as nib
import SimpleITK as sitk

# What does this guy do?
# What are assumptions wrt. to input data
# How to use it?
# What to use it for?
# How not to use it?
#
#
# Essentially: we need an abstract here.
# An user-end summary. The comments within
# individual functions act as documentation
# for developers and for troubleshooting.


# original ITK_DIRECTION = np.array([-1, -1, 1])
ITK_DIRECTION = np.array([-1, -1, 1])
INPUT_DIRECTION = np.array([1, 1, 1])
#INPUT_DIRECTION = np.array([-1,-1,-1])
#ITK_DIRECTION = np.array([-1, -1, -1])


def voxel_sizes_from_affine(affine):
    """
    Returns voxel size in each direction given affine (mapping between
    voxel coordinates and world (mm) coordinates.
    Parameters
    ----------
    affine : 2D array-like
        Affine transformation array.  Usually shape (4, 4), but can be any 2D
        array.
    Returns
    -------
    vox_sizes : 1D array
        Voxel sizes for each input axis of affine.  Usually 1D array length 3,
        but in general has length (N-1) where input `affine` is shape (M, N).
    """
    top_left = affine[:-1, :-1]
    return np.sqrt(np.sum(top_left ** 2, axis=0))


def get_unity_transform(phys_aff):
    """
    Returns affine transform mapping from physical space to coordinate space of
    uniformly sized voxels (1x1x1) in arbitrary units with origin at (0,0,0).
    The orientation of output space is assumed to be RAS by itk convention.
    Parameters
    ----------
    phys_aff : 2D array-like
        Affine transformation array.
    Returns
    -------
    affine :  2D array-like
    """
    scale_vector = np.ones(3) / voxel_sizes_from_affine(phys_aff) * INPUT_DIRECTION
    Tvector = phys_aff[:3, 3] * scale_vector
    Tvector_dir = ITK_DIRECTION / np.sign(np.diag(phys_aff)[:3])
    Tvector *= Tvector_dir

    unit_affine = nib.affines.from_matvec(np.diag(scale_vector), Tvector)

    return unit_affine


def get_translation_transform(shift):
    """
    Returns translation transform shifting the origin to given voxel
    Parameters
    ----------
    shift : 1D array-like
        coordinates of the new origin in itk RAS oriented space
    Returns
    -------
    affine :  2D array-like
    """
    itk_params = np.hstack([np.eye(3, 3).flatten(), shift * ITK_DIRECTION * -1])

    new_transform = sitk.AffineTransform(3)
    new_transform.SetParameters(itk_params)
    new_transform.SetFixedParameters((0., 0., 0.))

    return new_transform


def get_affine_from_file(fname):
    """
    Reads affine from nifty file and returns it
    Parameters
    ----------
    fname : str
        path to nifti file
    Returns
    -------
    affine :  2D array-like
    """
    img = nib.load(fname)
    return img.affine


def get_affine_as_itk(affine):
    """
    Creates ITK compatible affine transform from numpy 2D array
    Parameters
    ----------
    phys_aff : 2D array-like
        Affine transformation array.
    Returns
    -------
    new_transform :  sitk.AffineTransform
    """
    itk_params = np.hstack([affine[:3, :3].flatten(), affine[:3, 3]])

    new_transform = sitk.AffineTransform(3)
    new_transform.SetParameters(itk_params)
    new_transform.SetFixedParameters((0., 0., 0.))

    return new_transform


def write_affine_as_itk(affine, fname):
    """
    Writes numpy 2D array as itk transform to txt file.
    Parameters
    ----------
    affine : 2D array-like
    fname: str
        path to output txt file
    Returns
    -------
    new_transform :  sitk.AffineTransform
    """
    itk_params = np.hstack([affine[:3, :3].flatten(), affine[:3, 3]])

    new_transform = sitk.AffineTransform(3)
    new_transform.SetParameters(itk_params)
    new_transform.SetFixedParameters((0., 0., 0.))

    sitk.WriteTransform(new_transform, fname)

    return new_transform


def compute_transform(input_fn, output_fn, shift):
    """
    #XXX ??
    """
    transform = get_unity_transform(get_affine_from_file(input_fn))
    translation = get_translation_transform(shift)
    translation = nib.affines.from_matvec(
        np.array(translation.GetMatrix()).reshape(3, 3),
        translation.GetTranslation())

    combined = translation.dot(transform)
    write_affine_as_itk(combined, output_fn)


parser = \
    argparse.ArgumentParser(description='Compute ITK transform to unit voxel space')

parser.add_argument('--input-file-mm',
                    required=True,
                    help="path to input file",
                    type=str)

parser.add_argument('--output-file',
                    required=True,
                    help='path to output file',
                    type=str)

parser.add_argument('--voxel-shift',
                    required=False,
                    help='shift in unit voxel space to match position of the template',
                    nargs=3,
                    default=[0, 0, 0],
                    type=int)

args = parser.parse_args()

if __name__ == '__main__':

    compute_transform(args.input_file_mm,
                      args.output_file,
                      args.voxel_shift)
