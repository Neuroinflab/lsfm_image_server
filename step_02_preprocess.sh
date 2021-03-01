#!/bin/bash -xe

source defaults.sh

# This script is for preprocessing the input images.
# Both 'signal' and 'structural'

# The assumptions on the input data are that
# - the data is spatially consistent, as is has been extracted from the HDF5 container.
# - the niftii files are in proper RAS orientation.

MASK_STRUCT_PHYSICAL001_LO=${DIR_INPUT_DATA}"/01_structural/001_25_mm_mask.nii.gz"
MASK_STRUCT_UNIT001_LO=${DIR_IMG_PREPROCESSED}"/structural_001_25_unit_mask.nii.gz"

MASK_SIGNAL_PHYSICAL001_LO=${DIR_INPUT_DATA}"/00_signal/001_25_mm_mask.nii.gz"
MASK_SIGNAL_UNIT001_LO=${DIR_IMG_PREPROCESSED}"/signal_001_25_unit_mask.nii.gz"

debug "Transforming 25um images to the unit space."

TransformToNativeSpace \
    ${IMG_SIGNAL_PHYSICAL001_LO} \
    ${IMG_SIGNAL_UNIT001_LO} \
    ${TR_SIGNAL_UNIT_001_LO}

TransformToNativeSpace \
    ${IMG_STRUCT_PHYSICAL001_LO} \
    ${IMG_STRUCT_UNIT001_LO} \
    ${TR_STRUCT_UNIT_001_LO}


debug "Preprocessing structural images."

N4BiasFieldCorrection \
    -d 3 \
    -s 4 \
    -i ${IMG_STRUCT_UNIT001_LO} \
    --output ${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01.nii.gz

python scripts/clahe.py \
    ${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01.nii.gz \
    ${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01_clahe.nii.gz \
      0.01 

# requires c3d1.1
c3d \
    ${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01_clahe.nii.gz \
    -median 1x1x1 \
    -o ${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01_clahe_median.nii.gz \

ImageMath 3 \
    ${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01_clahe_median_grad.nii.gz \
    Grad \
    ${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01_clahe_median.nii.gz \
    1.5 1



debug "Preprocessing atlas image."

ImageMath 3 \
    ${DIR_TEMPLATE}/template_grad.nii.gz \
    Grad \
    ${FILE_TEMPLATE} \
    1.5 1

# This script is for preprocessing the mask for structure image by transforming it
# to the unit coordinate system.

# The assumptions on the input data are that
# - the data is spatially consistent, as is has been extracted from the HDF5 container.
# - the niftii files are in proper RAS orientation.



debug "Transforming 25um mask to the unit space."

TransformToNativeSpace \
    ${MASK_STRUCT_PHYSICAL001_LO} \
    ${MASK_STRUCT_UNIT001_LO} \
    ${TR_STRUCT_UNIT_001_LO}


debug "Transforming signal 25um mask to the unit space."


TransformToNativeSpace \
    ${MASK_SIGNAL_PHYSICAL001_LO} \
    ${MASK_SIGNAL_UNIT001_LO} \
    ${TR_STRUCT_UNIT_001_LO}

