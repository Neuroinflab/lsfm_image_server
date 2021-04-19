#!/usr/bin/env python
# encoding: utf-8

import sys

import nibabel as nib
from skimage import data, exposure, img_as_float

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    clip_limit_ = float(sys.argv[3])
except:
    clip_limit_ = 0.013

image = nib.load(input_file)
imd = image.get_data()

im1 = exposure.rescale_intensity(imd, (imd.min(), imd.max()), (0, 1))
for i in range(imd.shape[2]):
    print(i)
    im1[:,:,i] = exposure.equalize_adapthist(im1[:,:,i], clip_limit=clip_limit_)

im2 = exposure.rescale_intensity(imd, (imd.min(), imd.max()), (0, 1))
for i in range(imd.shape[1]):
    print(i)
    im2[:,i,:] = exposure.equalize_adapthist(im2[:,i,:], clip_limit=clip_limit_)

im3 = exposure.rescale_intensity(imd, (imd.min(), imd.max()), (0, 1))
for i in range(imd.shape[0]):
    print(i)
    im3[i,:,:] = exposure.equalize_adapthist(im3[i,:,:], clip_limit=clip_limit_)


img = (im1 + im2 + im3) / 3.
img = nib.Nifti1Image(img, image.affine, header=image.header)
nib.save(img, output_file)
