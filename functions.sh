#!/bin/sh

source defaults.sh

LSFM_COREGISTRATION_AFFINE=1x0x0x0x0x0
LSFM_ANTS_ITERATIONS_MED=100x100x00
LSFM_ANTS_ITERATIONS_HIGH=200x200x100

export LSFM_ANTS_ITERATIONS="${LSFM_ANTS_ITERATIONS:-$LSFM_COREGISTRATION_AFFINE}"


function lsfm_set_mapping_affine
{
    log "Coregistration with atlas set to affine."
    export LSFM_ANTS_ITERATIONS=${LSFM_COREGISTRATION_AFFINE}
}

function lsfm_set_mapping_overnight
{
    log "Coregistration with atlas set to low quality deformable."
    export LSFM_ANTS_ITERATIONS=${LSFM_ANTS_ITERATIONS_MED}
}

function lsfm_set_mapping_high_quality
{
    log "Coregistration with atlas set to high quality deformable."
    export LSFM_ANTS_ITERATIONS=${LSFM_ANTS_ITERATIONS_HIGH}
}

