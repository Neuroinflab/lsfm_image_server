#!/bin/bash -xe

# The purpose of this step is to prepare the experimental imaging dataset
# for processing. The images to coregister originate from the HDF5 format
# file. This script extracts the imging data in the desired resolution,
# suitable for the purpose of the coregistration with the atlas.

source defaults.sh

HDF5_SRC_PATH=./${CASE_ID}.h5
HDF5_CFOS_CHANNEL_NAME=cfos
HDF5_AUTO_CHANNEL_NAME=autofluo


lsfmpy export \
    --hdf-path ${HDF5_SRC_PATH} \
    --channel-name ${HDF5_CFOS_CHANNEL_NAME} \
    --output-path ${IMG_SIGNAL_PHYSICAL001_LO} \
    --output-resolution ${C_LO_RES_MM},${C_LO_RES_MM},${C_LO_RES_MM} \
    --input-orientation RPI 

lsfmpy export \
    --hdf-path ${HDF5_SRC_PATH} \
    --channel-name ${HDF5_AUTO_CHANNEL_NAME} \
    --output-path ${IMG_STRUCT_PHYSICAL001_LO} \
    --output-resolution ${C_LO_RES_MM},${C_LO_RES_MM},${C_LO_RES_MM} \
    --input-orientation RPI


