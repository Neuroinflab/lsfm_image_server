#!/bin/bash -x
CASE_ID=`cat case_id`
DIR_RAMDISK=/dev/shm/
CURRENT_TIMESTAMP=`date +%Y-%m-%d-%H-%M-%S`

DIR_TEMPLATE_SRC=../allen_ccf3_template/output_25/

DIR_INPUT_DATA=01_input_data/
DIR_INPUT_DATA_STRUCT=${DIR_INPUT_DATA}/01_structural/
DIR_INPUT_DATA_SIGNAL=${DIR_INPUT_DATA}/00_signal/
 
DIR_TEMPLATE=00_reference/
DIR_IMG_PREPROCESSED=02_preprocessed_input_images/
DIR_COREGISTRATION=03_coregistration/
DIR_TRANSFORMS=04_transforms/
DIR_RESULTS=05_results/

mkdir -p ${DIR_IMG_PREPROCESSED}
mkdir -p ${DIR_COREGISTRATION}
mkdir -p ${DIR_COREGISTRATION}/01_signal_and_structural/
mkdir -p ${DIR_COREGISTRATION}/10_structural_and_template/

mkdir -p ${DIR_INPUT_DATA}
mkdir -p ${DIR_INPUT_DATA_STRUCT} ${DIR_INPUT_DATA_SIGNAL}
mkdir -p ${DIR_TEMPLATE}
mkdir -p ${DIR_RESULTS}

mkdir -p ${DIR_TRANSFORMS}
mkdir -p ${DIR_COREGISTRATION} 


C_TEMPLATE_HEMISPHERE="left"

FILE_DEFAULT_SEGMENTATION=${DIR_TEMPLATE}/segmentation.nii.gz
FILE_SEGMENTATION=${FILE_DEFAULT_SEGMENTATION}

FILE_DEFAULT_TEMPLATE=${DIR_TEMPLATE}/template.nii.gz
FILE_TEMPLATE=${FILE_DEFAULT_TEMPLATE}

FILE_DEFAULT_TEMPLATE_PHYSICAL=${DIR_TEMPLATE}/template_25um_both_hemispheres.nii.gz
FILE_TEMPLATE_PHYSICAL=${FILE_DEFAULT_TEMPLATE_PHYSICAL}

FILE_DEFAULT_SEGMENTATION_PHYSICAL=${DIR_TEMPLATE}/segmentation_25um_both_hemispheres.nii.gz
FILE_SEGMENTATION_PHYSICAL=${FILE_DEFAULT_SEGMENTATION_PHYSICAL}

AFFINE_ATLAS_LEFT_TO_UNIT_SPACE=${DIR_TEMPLATE}/template_${C_TEMPLATE_HEMISPHERE}_hemisphere_to_unit_space_Affine.txt
AFFINE_ATLAS_LEFT_TO_UNIT_SPACE_TRANSLATION=${DIR_TEMPLATE}/template_${C_TEMPLATE_HEMISPHERE}_hemisphere_to_unit_space_translation_Affine.txt

FILE_MANUAL_AFFINE=${DIR_COREGISTRATION}/10_structural_and_template/manual_initial_affine.txt

function log {
    echo "[$(date --rfc-3339=seconds)]: $*"
}

function debug {
    echo -e "\e[90m`log $*` \e[22m\e[0m"
}

function error {
    echo -e "\e[31m`log $*` \e[22m\e[0m"
}

function warning {
    echo -e "\e[33m`log $*` \e[22m\e[0m"
}

function positive {
    echo -e "\e[32m`log $*` \e[22m\e[0m"
}

function timestamp {
      date +"%s"
}

function CheckExecutables {
    local EXECUTABLES='mkdir rm cp tar WarpImageMultiTransform c3d antsRegistration'

    for executable in ${EXECUTABLES}
    do
        prog_path=`which ${executable}`
        if [ $? -ne 0 ];
        then
            error "(!) The ${executable} is not available."
            error "(!) The environment cannot be set. Exiting."
            exit 1
        fi
    done
}

function _GetGenericInputImage {
    local channel=$1
    local resolution=$2
    local space=$3
    local dir_type=$4

    echo ${DIR_INPUT_DATA}"/${dir_type}/${channel}_${resolution}_${space}.nii.gz"
}

function GetStructuralImage {
    local channel=$1
    local resolution=$2
    local space=$3

    echo `_GetGenericInputImage ${channel} ${resolution} ${space} 01_structural`
}

function GetSignalImage {
    local channel=$1
    local resolution=$2
    local space=$3

    echo `_GetGenericInputImage ${channel} ${resolution} ${space} 00_signal`
}

function GetUnitTransformFilename {
	local dir_type=$1
    local channel=$2
    local resolution=$3

	echo ${DIR_TRANSFORMS}/${dir_type}_${channel}_affine_unit_to_physical_${resolution}.txt	
}


function GetDisplacementFieldFilename {
	local dir_type=$1
	local channel=$2
	local direction=$3

	echo ${DIR_TRANSFORMS}/${direction}_warp_${dir_type}_${channel}.nii.gz
}

function GetResultImageFilename {
	local moving=$1
	local fixed=$2
	local transform_type=$3

	echo ${DIR_RESULTS}/${moving}_in_${fixed}_${transform_type}.nii.gz
}


function TransformToNativeSpace {
    local FN_PHYSICAL=$1
    local FN_UNIT=$2
    local FN_TRANSFORM=$3

    python scripts/header_converter.py \
        --input-file-mm ${FN_PHYSICAL} \
        --output-file ${FN_TRANSFORM}

    c3d ${FN_PHYSICAL} \
        -spacing 1x1x1mm \
        -origin 0x0x0mm \
        -scale 0 \
        -orient RAS \
        -type ushort \
        -o ${FN_UNIT}

    # One would think that it would be sufficient to use the convert3d here.
    # However, this is not that easy as we do not have a possibility of inverting
    # the transform in the convert3d.
    WarpImageMultiTransform 3 \
        ${FN_PHYSICAL} \
        ${FN_UNIT} \
        -R ${FN_UNIT} \
        -i ${FN_TRANSFORM}
}

C_HI_RES=10
C_LO_RES=25
C_HI_RES_MM=`python -c "print(${C_HI_RES}./1000)" | tr -d '\n'`
C_LO_RES_MM=`python -c "print(${C_LO_RES}./1000)" | tr -d '\n'`

IMG_STRUCT_PHYSICAL001_LO=`GetStructuralImage 001 ${C_LO_RES} mm`
IMG_STRUCT_PHYSICAL001_HI=`GetStructuralImage 001 ${C_HI_RES} mm`
IMG_STRUCT_UNIT001_LO=`GetStructuralImage 001 ${C_LO_RES} unit`
IMG_STRUCT_UNIT001_HI=`GetStructuralImage 001 ${C_HI_RES} unit`

IMG_SIGNAL_PHYSICAL001_LO=`GetSignalImage 001 ${C_LO_RES} mm`
IMG_SIGNAL_PHYSICAL001_HI=`GetSignalImage 001 ${C_HI_RES} mm`
IMG_SIGNAL_UNIT001_LO=`GetSignalImage 001 ${C_LO_RES} unit`
IMG_SIGNAL_UNIT001_HI=`GetSignalImage 001 ${C_HI_RES} unit`

TR_STRUCT_UNIT_001_LO=`GetUnitTransformFilename structural 001 ${C_LO_RES}`
TR_STRUCT_UNIT_001_HI=`GetUnitTransformFilename structural 001 ${C_HI_RES}`
TR_SIGNAL_UNIT_001_LO=`GetUnitTransformFilename signal 001 ${C_LO_RES}`
TR_SIGNAL_UNIT_001_HI=`GetUnitTransformFilename signal 001 ${C_HI_RES}`
