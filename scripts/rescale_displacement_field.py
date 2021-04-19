#!/usr/bin/env python
# encoding: utf-8



"""
This is a very simple script which purpose is to, simply, rescale the values of
the displacement field so that the warp are applicable in the physical
coordinates, not only in the unit coordinates. The rescaling is performed based
on the provided affine transformation. In such case the individual components
can be scaled with a different factor.

Note that currently only the scaling factors of the affine transformation are
taken into account. No shearing is utilized and, obviously, the rotation has no
influence.

Currently this script uses an ungly combination of SimpleITK / Niftii
combination of modules to carry out the scaling. This is bad and ugly and, in
the future should be reimplemented using pure SimpleITK or, even, better just
ITK.

TODO: XXX: Add assertions 1) Header, dimensionality.
"""

import os
import sys
import argparse

import nibabel as nib
import SimpleITK as sitk


def rescale_displacement_field(args):
    """
    Rescales the displacement field according to the provided reference affine
    transformation.
    """

    input_df = args.input_displacement
    output_df = args.output_displacement
    transform_file = args.reference_transform


    displacement_field = nib.load(input_df)
    displacement = displacement_field.get_data()

    transform = sitk.ReadTransform(transform_file)
    scaling_factors = [transform.GetParameters()[x] for x in [0, 4, 8]]
    print(("SCALING FACTORS: {}".format(scaling_factors)))
    print(("DISPLACEMENT SHAPE: {}".format(displacement.shape)))


    displacement[:, :, :, :, 0] /= scaling_factors[0]
    displacement[:, :, :, :, 1] /= scaling_factors[1]
    displacement[:, :, :, :, 2] /= scaling_factors[2]


    img = nib.Nifti1Image(displacement, displacement_field.affine, header=displacement_field.header)
    nib.save(img, args.output_displacement)


def parse_arguments():
    """
    """

    parser = argparse.ArgumentParser(
        description='Rescales the values of the displacement field based on the provided affine transformation.')

    parser.add_argument("--input-filename", "-i", dest="input_displacement",
        type=str, help="Filename of the input displacement field.",
        required=True)

    parser.add_argument("--output-filename", "-o", dest="output_displacement",
        required=True, help="Output displacement field filename.", type=str)

    parser.add_argument("--affine-transformation", "-t", dest="reference_transform",
        required=True, help="Reference affine transformation (txt file).", type=str)

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_arguments()
    rescale_displacement_field(args)
