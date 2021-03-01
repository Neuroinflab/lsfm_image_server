#!/bin/bash -x

source defaults.sh

FIXED=${IMG_SIGNAL_UNIT001_LO}
MOVING=${IMG_STRUCT_UNIT001_LO}

FIXED_MASK=${DIR_IMG_PREPROCESSED}/signal_001_25_unit_mask.nii.gz
MOVING_MASK=${DIR_IMG_PREPROCESSED}/structural_001_25_unit_mask.nii.gz

FIXED_PHYSICAL=${IMG_SIGNAL_PHYSICAL001_LO}
MOVING_PHYSICAL=${IMG_STRUCT_PHYSICAL001_LO}

# This is tricky and was implemented for the purpose of using unit
# images of different resolution.
# If the registration has been calculating using LO images, the LO
# unit transforms have to be selected, otherwise the HI unit transformation
# should be selected.
#C_STRUCT_TO_UNIT_TRANSFORM=${TR_STRUCT_UNIT_001_HI}
C_STRUCT_TO_UNIT_TRANSFORM=${TR_STRUCT_UNIT_001_LO}
mkdir -p ${DIR_COREGISTRATION}/01_signal_and_structural/

USE_RIGID_AFFINE='false'
DO_CALCULATE_TRANSFORM='true'
DO_RESLICE_FIXED_TO_MOVING='true'
DO_RESLICE_MOVING_TO_FIXED='true'


OUTPUTNAME=${DIR_COREGISTRATION}/01_signal_and_structural/f-signal001_m-structural001

if [ ${DO_CALCULATE_TRANSFORM} = 'true' ]
then

    if [ ${USE_RIGID_AFFINE} = 'true' ]
    then
        RIGID=" --rigid-affine true --affine-gradient-descent-option 0.5x0.95x1.e-4x1.e-4"
    else
        RIGID=" --rigid-affine false "
    fi

    if [ -f ${OUTPUTNAME}Affine.txt ];
    then
        warning "The ${OUTPUTNAME}Affine.txt already exists. "
        warning "The existing calculations will be renamed to: ${DIR_COREGISTRATION}/01_signal_and_structural_previous_${CURRENT_TIMESTAMP}"
        mv -v ${DIR_COREGISTRATION}/01_signal_and_structural/ \
           ${DIR_COREGISTRATION}/01_signal_and_structural_previous_${CURRENT_TIMESTAMP}
    fi

    mkdir -p ${DIR_COREGISTRATION}/01_signal_and_structural/
    
    ANTS 3 \
        -m MI[${FIXED_MASK},${MOVING_MASK},1,32] \
        -o ${OUTPUTNAME} \
        -i 0 \
        --use-Histogram-Matching \
        --number-of-affine-iterations 10000x10000x10000x10000x10000 \
        ${RIGID}

    # Copying the masks to the working directory of the coregistration to
    # preserve them for the purpose of potential debugging
    cp -v ${FIXED_MASK} ${MOVING_MASK} ${DIR_COREGISTRATION}/01_signal_and_structural/

    # Generate affine transformation in physical space:
    ComposeMultiTransform 3 \
        ${DIR_COREGISTRATION}/01_signal_and_structural/f-signal001_m-structural001_physical_affine.txt \
        ${C_STRUCT_TO_UNIT_TRANSFORM} \
        ${OUTPUTNAME}Affine.txt \
        -i ${C_STRUCT_TO_UNIT_TRANSFORM}
fi


if [ ${DO_RESLICE_MOVING_TO_FIXED} = 'true' ]
then
    # Reslicing in the unit space:
    WarpImageMultiTransform 3 \
        ${MOVING} \
        ${OUTPUTNAME}AffineDeformedUnitSpace.nii.gz \
        -R ${FIXED} \
        -i ${OUTPUTNAME}Affine.txt

    # Reslicing in the physical space:
    WarpImageMultiTransform 3 \
        ${MOVING_PHYSICAL} \
        ${OUTPUTNAME}AffineDeformedPhysicalSpace.nii.gz \
        -R ${FIXED_PHYSICAL} \
        ${C_STRUCT_TO_UNIT_TRANSFORM} \
        ${OUTPUTNAME}Affine.txt \
        -i ${C_STRUCT_TO_UNIT_TRANSFORM}
fi

if [ ${DO_RESLICE_FIXED_TO_MOVING} = 'true' ]
then
    WarpImageMultiTransform 3 \
        ${FIXED_PHYSICAL} \
        ${DIR_IMG_PREPROCESSED}/signal_001_in_structural_001.nii.gz \
        -R ${IMG_STRUCT_UNIT001_LO} \
        -i ${TR_STRUCT_UNIT_001_LO} \
        -i ${DIR_COREGISTRATION}/01_signal_and_structural/f-signal001_m-structural001_physical_affine.txt \


fi
