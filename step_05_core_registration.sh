#!/bin/bash -x

source defaults.sh

FIXED_IMAGE=$1
MOVING_IMAGE=$2
FIXED_IMAGE_STRUCT_GRAD=$3
MOVING_IMAGE_GRAD=$4

REG_PARAMETER=$5
SYN_PARAMETER=$6


FLAG_AFFINE_ALIGNMENT="true"
FLAG_DEFORMABLE_COREGISTRATION="true"
FLAG_AFFINE_CORRECTION="true"
FLAG_CONVERT_DISPLACEMENT_FIELD="true"
FLAG_RESLICE_IMAGES="true"

WORKSPACE=`pwd`
FIXED_MD5=`md5sum ${FIXED_IMAGE} | cut -f1 -d" "`
MOVING_MD5=`md5sum ${MOVING_IMAGE} | cut -f1 -d" "`
WORKDIR=${WORKSPACE}/workdir_fixed_${FIXED_MD5}_moving_${MOVING_MD5}/

mkdir -p ${WORKDIR}

AFFINE_REG_FIXED=${FIXED_IMAGE}
AFFINE_REG_MOVING=${MOVING_IMAGE}
REG_OUTPUT=${WORKDIR}/registration_${REG_PARAMETER}_${SYN_PARAMETER}/

AFFINE_REG_FIXED=${DIR_IMG_PREPROCESSED}/structural_001_25_unit_mask.nii.gz
AFFINE_REG_MOVING=${DIR_TEMPLATE}/template_mask.nii.gz

SCRIPT_RESCALE_DF=scripts/rescale_displacement_field.py

mkdir -p ${REG_OUTPUT}


if [ -e ${REG_OUTPUT}/Deformed.nii.gz ];
then
   warning "Deformable transformations seem to be already calculated."
   warning "Press Ctrl + C to stop if you don't want to continue."
   sleep 5
fi


function AffineCoregistration {
    # This function runs the affine coregistration. There is nothing particularly strange here.
    # However the process itself is a little bit complex. 
    #

    if [ -e ${FILE_MANUAL_AFFINE} ]
    then
        cp -v ${FILE_MANUAL_AFFINE} ${REG_OUTPUT}/0GenericAffine.txt
    else
        if [ ! -f ${REG_OUTPUT}/initial_affine.txt ];
        then
            debug "No initial affine transformation found. Calculating one."
            antsRegistration \
		--dimensionality 3 \
		--float 0 \
		--output ${REG_OUTPUT}/translation_ \
		--interpolation Linear \
		--use-histogram-matching 0 \
		--initial-moving-transform [${AFFINE_REG_FIXED},${AFFINE_REG_MOVING},1] \
		--transform Rigid[0.1] \
		--metric MI[${AFFINE_REG_FIXED},${AFFINE_REG_MOVING},1,32,Regular,0.25] \
		--convergence [1x1x1x1,1e-6,10] \
		--shrink-factors 4x3x2x1 \
		--smoothing-sigmas 0x0x0x0vox

        ConvertTransformFile 3 \
                 ${REG_OUTPUT}/translation_0GenericAffine.mat \
                 ${REG_OUTPUT}/initial_affine.txt

        else
            debug "Initial affine transformation found: "
            debug "${REG_OUTPUT}/initial_affine.txt"
        fi

        if [ ! -f ${REG_OUTPUT}/0GenericAffine.txt ];
        then
            antsRegistration \
                --dimensionality 3 \
                --initial-moving-transform ${REG_OUTPUT}/initial_affine.txt \
                --metric "Mattes[${AFFINE_REG_FIXED},${AFFINE_REG_MOVING},1,32,regular]" \
                --transform "Affine[0.5]" \
                --convergence 10000x10000x10000x00000 \
                --shrink-factors 4x3x2x1 \
                --smoothing-sigmas 0x0x0x0 \
                --output ${REG_OUTPUT}/

            ConvertTransformFile 3 \
                ${REG_OUTPUT}/0GenericAffine.mat \
                ${REG_OUTPUT}/0GenericAffine.txt

            WarpImageMultiTransform 3 \
                ${AFFINE_REG_FIXED} \
                ${REG_OUTPUT}/fixed_to_atlas_affine_check.nii.gz \
                -R ${AFFINE_REG_MOVING} \
                --use-NN \
                -i ${REG_OUTPUT}/0GenericAffine.txt
            
        fi
    fi

}


function DeformableCoregistration {


    if [ ! -f ${REG_OUTPUT}WarpX.nii.gz ];
    then
        ANTS 3 -v \
            -m CC[${FIXED_IMAGE},${MOVING_IMAGE},1.0,${SYN_PARAMETER}] \
            -m CC[${FIXED_IMAGE_STRUCT_GRAD},${MOVING_IMAGE_GRAD},1.0,${SYN_PARAMETER}] \
            -m MSQ[${DIR_IMG_PREPROCESSED}/structural_001_25_unit_mask.nii.gz,${DIR_TEMPLATE}/template_mask.nii.gz,1,5] \
            -t SyN[0.25] \
            -r Gauss[${REG_PARAMETER}] \
            -a ${REG_OUTPUT}/0GenericAffine.txt \
            --continue-affine false \
            --number-of-affine-iterations 10000x10000x10000x10000x10000 \
            --rigid-affine false \
            --use-Histogram-Matching 1 \
            --affine-gradient-descent-option 0.1x0.5x5.e-6x5.e-6 \
            -i ${LSFM_ANTS_ITERATIONS} \
            -o ${REG_OUTPUT} \
            --use-all-metrics-for-convergence

            # This part is just to make sure that the coregistration went fine
            # in the unit space. Essentially, this is just For debugging purposes
            # and for tuning the registration parameters.
            WarpImageMultiTransform 3 \
                ${MOVING_IMAGE} ${REG_OUTPUT}/DeformedJustAffine.nii.gz \
                -R ${FIXED_IMAGE} \
                ${REG_OUTPUT}Affine.txt

            WarpImageMultiTransform 3 \
                ${MOVING_IMAGE} ${REG_OUTPUT}/DeformedPureFixed.nii.gz \
                -R ${FIXED_IMAGE} \
                ${REG_OUTPUT}Warp.nii.gz ${REG_OUTPUT}Affine.txt
            
            WarpImageMultiTransform 3 \
                ${MOVING_IMAGE} ${WORKDIR}/DeformedJustAffineStruct.nii.gz \
                -R ${FIXED_IMAGE} \
                ${REG_OUTPUT}Affine.txt

            WarpImageMultiTransform 3 \
                ${MOVING_IMAGE} ${REG_OUTPUT}/Deformed.nii.gz \
                -R ${FIXED_IMAGE} \
                ${REG_OUTPUT}Warp.nii.gz ${REG_OUTPUT}Affine.txt

            WarpImageMultiTransform 3 \
                ${FILE_SEGMENTATION} ${REG_OUTPUT}/segmentation.nii.gz \
                -R ${FIXED_IMAGE} \
                --use-NN \
                ${REG_OUTPUT}Warp.nii.gz ${REG_OUTPUT}Affine.txt

    else
        debug "Deformable coregistration seems to be already done:"
        debug "${REG_OUTPUT}/initial_affine.txt"
    fi
}


function AffineCoregistrationCorrection {
    # This part is just to make sure that the coregistration went fine
    # Mapping the template / segmentation and structural 001 between one another.
    # This is a bit tricky. Here we define an 'uncorrected' physical affine transformation
    # between the template and the physical space.

    ComposeMultiTransform 3 \
        ${DIR_TRANSFORMS}/template_left_physical_to_structural_001_Affine_uncorrected.txt \
        ${AFFINE_ATLAS_LEFT_TO_UNIT_SPACE} \
        -i ${AFFINE_ATLAS_LEFT_TO_UNIT_SPACE_TRANSLATION} \
        ${REG_OUTPUT}Affine.txt \
        -i ${AFFINE_ATLAS_LEFT_TO_UNIT_SPACE} \

    antsApplyTransforms -d 3 \
        --input-image-type 0 \
        --input ${FILE_TEMPLATE_PHYSICAL} \
        --reference-image ${FIXED_IMAGE} \
        --output ${WORKDIR}/template_in_structural_001_Affine_unc.nii.gz \
        --transform [${TR_STRUCT_UNIT_001_LO},1] \
        --transform [${DIR_TRANSFORMS}/template_left_physical_to_structural_001_Affine_uncorrected.txt,0] \
        --default-value 0

    ANTS 3 \
        -m MI[${WORKDIR}/DeformedJustAffineStruct.nii.gz,${WORKDIR}/template_in_structural_001_Affine_unc.nii.gz,1,32] \
        -o ${DIR_TRANSFORMS}/template_to_structural_001_subpixel_correction_ \
        -i 0 \
        --use-Histogram-Matching \
        --number-of-affine-iterations 10000 \
        --rigid-affine false \
        --affine-gradient-descent-option 0.05x0.5x1.e-4x1.e-4

    antsApplyTransforms -d 3 \
        --input-image-type 0 \
        --input ${WORKDIR}/template_in_structural_001_Affine_unc.nii.gz \
        --reference-image ${FIXED_IMAGE} \
        --output ${WORKDIR}/subpixel_correction_verification.nii.gz \
        --transform [${DIR_TRANSFORMS}/template_to_structural_001_subpixel_correction_Affine.txt,0] \
        --default-value 0

    ComposeMultiTransform 3 \
        ${DIR_TRANSFORMS}/template_left_physical_to_structural_001_Affine.txt \
        ${AFFINE_ATLAS_LEFT_TO_UNIT_SPACE} \
        -i ${AFFINE_ATLAS_LEFT_TO_UNIT_SPACE_TRANSLATION} \
        ${DIR_TRANSFORMS}/template_to_structural_001_subpixel_correction_Affine.txt \
        ${REG_OUTPUT}Affine.txt \
        -i ${AFFINE_ATLAS_LEFT_TO_UNIT_SPACE} &> tmptmp.txt


    antsApplyTransforms -d 3 \
        --input-image-type 0 \
        --input ${FILE_TEMPLATE_PHYSICAL} \
        --reference-image ${IMG_STRUCT_PHYSICAL001_LO} \
        --output ${WORKDIR}/template_in_structural_001_Affine.nii.gz \
        --transform [${DIR_TRANSFORMS}/template_left_physical_to_structural_001_Affine.txt,0] \
        --default-value 0

    # Unfortunately, this part has to be removed. We cannot merge together
    # partial affine transformations into a sigle one, unfortunately. The problem is that,
    # In between of applying the affine transformations, we need to apply a deformable
    # one.
    # Instead, we just copy the auto -> to template transformation to 
    ComposeMultiTransform 3 \
        ${DIR_TRANSFORMS}/template_left_physical_to_signal_001_Affine.txt \
        ${DIR_TRANSFORMS}/template_left_physical_to_structural_001_Affine.txt \

    # Unfortunately, ants cannot read even symbolic links. We have to create a hardcopy.
    cp -v \
        ${DIR_COREGISTRATION}/01_signal_and_structural/f-signal001_m-structural001_physical_affine.txt \
        ${DIR_TRANSFORMS}/structure_001_to_signal_001_physical_affine.txt


}


function MapDisplacementToPhysical {
    # Now, one of the mos difficult part: transforming the displacement
    # field affinely and then rescaling it so it can act
    # in the physical coordinates.

    # The displacement fields are processed in the ramdrive
    # and just then they are transfered to their final location:
    temp_imgname="$(mktemp ${DIR_RAMDISK}/temp_imgname.XXXXXX.nii)"

    WarpImageMultiTransform 3 \
        ${REG_OUTPUT}Warp.nii.gz \
        ${temp_imgname} \
        -R ${IMG_STRUCT_PHYSICAL001_LO} \
        ${TR_STRUCT_UNIT_001_LO}

    python ${SCRIPT_RESCALE_DF} \
        --input-filename ${temp_imgname} \
        --output-filename `GetDisplacementFieldFilename structural 001 forward` \
        --affine-transformation ${TR_STRUCT_UNIT_001_LO}

    WarpImageMultiTransform 3 \
        ${REG_OUTPUT}InverseWarp.nii.gz \
        ${temp_imgname} \
        -R ${IMG_STRUCT_PHYSICAL001_LO} \
        ${TR_STRUCT_UNIT_001_LO}

    python ${SCRIPT_RESCALE_DF} \
        --input-filename ${temp_imgname} \
        --output-filename `GetDisplacementFieldFilename structural 001 inverse` \
        --affine-transformation ${TR_STRUCT_UNIT_001_LO}

    rm -rfv ${temp_imgname}
}


function ResliceImages {
    #
    #
    # Copy the affine transformation in unit space to the
    # output transformations_directory. This transformation
    # will not be used for any other purpose than just the debugging.
    cp -v \
        ${REG_OUTPUT}/0GenericAffine.txt \
        ${DIR_TRANSFORMS}/template_to_structural_unit_space_Affine.txt

    antsApplyTransforms -d 3 \
        --input-image-type 0 \
        --input ${FILE_TEMPLATE_PHYSICAL} \
        --reference-image ${IMG_STRUCT_PHYSICAL001_LO} \
        --output `GetResultImageFilename template structural-001-25 deformable` \
        --transform `GetDisplacementFieldFilename structural 001 forward` \
        --transform ${DIR_TRANSFORMS}/template_left_physical_to_structural_001_Affine.txt \
        --default-value 0

    antsApplyTransforms -d 3 \
        --input-image-type 0 \
        --input ${FILE_TEMPLATE_PHYSICAL} \
        --output `GetResultImageFilename template signal-001-25 deformable` \
        --reference-image ${IMG_SIGNAL_PHYSICAL001_LO} \
        --transform ${DIR_TRANSFORMS}/structure_001_to_signal_001_physical_affine.txt \
        --transform `GetDisplacementFieldFilename structural 001 forward` \
        --transform ${DIR_TRANSFORMS}/template_left_physical_to_signal_001_Affine.txt \
        --default-value 0

    antsApplyTransforms -d 3 \
        --interpolation NearestNeighbor \
        --input-image-type 0 \
        --input ${FILE_SEGMENTATION_PHYSICAL} \
        --reference-image ${IMG_STRUCT_PHYSICAL001_LO} \
        --output `GetResultImageFilename segmentation structural-001-25 deformable` \
        --transform `GetDisplacementFieldFilename structural 001 forward` \
        --transform ${DIR_TRANSFORMS}/template_left_physical_to_structural_001_Affine.txt \
        --default-value 0

    antsApplyTransforms -d 3 \
        --interpolation NearestNeighbor \
        --input-image-type 0 \
        --input ${FILE_SEGMENTATION_PHYSICAL} \
        --output `GetResultImageFilename segmentation signal-001-25 deformable` \
        --reference-image ${IMG_SIGNAL_PHYSICAL001_LO} \
        --transform ${DIR_TRANSFORMS}/structure_001_to_signal_001_physical_affine.txt \
        --transform `GetDisplacementFieldFilename structural 001 forward` \
        --transform ${DIR_TRANSFORMS}/template_left_physical_to_signal_001_Affine.txt \
        --default-value 0

    # And now: mapping the LSFM data to the template.

    antsApplyTransforms -d 3 \
        --input-image-type 0 \
        --reference-image ${FILE_TEMPLATE_PHYSICAL} \
        --input ${IMG_STRUCT_PHYSICAL001_LO} \
        --output `GetResultImageFilename structural-001-25 template deformable` \
        --transform [${DIR_TRANSFORMS}/template_left_physical_to_structural_001_Affine.txt,1] \
        --transform `GetDisplacementFieldFilename structural 001 inverse` \
        --default-value 0

    antsApplyTransforms -d 3 \
        --input-image-type 0 \
        --reference-image ${FILE_TEMPLATE_PHYSICAL} \
        --input ${IMG_SIGNAL_PHYSICAL001_LO} \
        --output `GetResultImageFilename signal-001-25 template deformable` \
        --transform [${DIR_TRANSFORMS}/template_left_physical_to_signal_001_Affine.txt,1] \
        --transform `GetDisplacementFieldFilename structural 001 inverse` \
        --transform [${DIR_TRANSFORMS}/structure_001_to_signal_001_physical_affine.txt,1] \
        --default-value 0




}



if [ ${FLAG_AFFINE_ALIGNMENT} = 'true' ]
then
    AffineCoregistration
fi

if [ ${FLAG_DEFORMABLE_COREGISTRATION} = 'true' ]
then
    DeformableCoregistration
fi

if [ ${FLAG_AFFINE_CORRECTION} = 'true' ]
then
   AffineCoregistrationCorrection 
fi

if [ ${FLAG_CONVERT_DISPLACEMENT_FIELD} = 'true' ]
then
   MapDisplacementToPhysical
fi

if [ ${FLAG_RESLICE_IMAGES} = 'true' ]
then
   ResliceImages
fi
