#!/bin/bash -x 

source defaults.sh

CheckExecutables

RESLICE_BACKGROUND=0

FIXED_IMAGE=${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01_clahe_median.nii.gz
FIXED_IMAGE_STRUCT_GRAD=${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01_clahe_median_grad.nii.gz

MOVING_IMAGE=${DIR_TEMPLATE}/template.nii.gz
MOVING_IMAGE_GRAD=${DIR_TEMPLATE}/template_grad.nii.gz

# If there is masking infolved in the coregistration process,
# calculate the masked images.
if [ -e ${DIR_IMG_PREPROCESSED}/structural_001_25_unit_mask.nii.gz ]
then
    c3d ${DIR_IMG_PREPROCESSED}/structural_001_25_unit_mask.nii.gz \
        ${FIXED_IMAGE} \
        -times \
        -o ${DIR_IMG_PREPROCESSED}/fixed_image_for_registration.nii.gz
    FIXED_IMAGE=${DIR_IMG_PREPROCESSED}/fixed_image_for_registration.nii.gz

    c3d ${DIR_IMG_PREPROCESSED}/structural_001_25_unit_mask.nii.gz \
        ${FIXED_IMAGE_STRUCT_GRAD} \
        -times \
        -o ${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01_clahe_median_grad_masked.nii.gz
    FIXED_IMAGE_STRUCT_GRAD=${DIR_IMG_PREPROCESSED}/structural_001_25um_N4i01_clahe_median_grad_masked.nii.gz

fi

if [ -e ${DIR_TEMPLATE}/template_mask.nii.gz ]
then
    c3d ${DIR_TEMPLATE}/template_mask.nii.gz \
        ${MOVING_IMAGE} \
        -times \
        -o ${DIR_TEMPLATE}/template_masked_for_registration.nii.gz
    MOVING_IMAGE=${DIR_TEMPLATE}/template_masked_for_registration.nii.gz

    c3d ${DIR_TEMPLATE}/template_mask.nii.gz \
        ${MOVING_IMAGE_GRAD} \
        -times \
        -o ${DIR_TEMPLATE}/template_gradient_for_registration_masked.nii.gz
    MOVING_IMAGE_GRAD=${DIR_TEMPLATE}/template_gradient_for_registration_masked.nii.gz

fi

REG_PARAMETER='1.75,0.75'
SYN_PARAMETER='4'

bash -x step_05_core_registration.sh \
    ${FIXED_IMAGE} ${MOVING_IMAGE} \
    ${FIXED_IMAGE_STRUCT_GRAD} ${MOVING_IMAGE_GRAD} \
    ${REG_PARAMETER} ${SYN_PARAMETER}

