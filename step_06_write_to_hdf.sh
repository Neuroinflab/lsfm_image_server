#!/bin/bash -xe

# The purpose of this step is to prepare the experimental imaging dataset
# for processing. The images to coregister originate from the HDF5 format
# file. This script extracts the imging data in the desired resolution,
# suitable for the purpose of the coregistration with the atlas.

source defaults.sh

HDF5_SRC_PATH=./${CASE_ID}.h5
HDF5_CFOS_CHANNEL_NAME=cfos
HDF5_AUTO_CHANNEL_NAME=autofluo


# to save displacement fields and segmentations we first need their metadata
dump_metadata --input-file `GetDisplacementFieldFilename structural 001 forward` \
	--output-file forward_warp.json << EOF
n
EOF

dump_metadata --input-file `GetDisplacementFieldFilename structural 001 inverse` \
	--output-file inverse_warp.json << EOF
n
EOF

dump_metadata --input-file `GetResultImageFilename segmentation signal-001-25 deformable` \
	--output-file cfos_segmentation.json << EOF
n
EOF

dump_metadata --input-file `GetResultImageFilename segmentation structural-001-25 deformable` \
	--output-file auto_segmentation.json << EOF
n
EOF

# write the affine for mapping between signal and structural channel
lsfmpy write-affine \
    --hdf-path ${HDF5_SRC_PATH} \
    --channel-name ${HDF5_CFOS_CHANNEL_NAME} \
    --affine-name  signal_to_structural \
    --affine-path ${DIR_TRANSFORMS}/structure_001_to_signal_001_physical_affine.txt

sleep 1

lsfmpy write-affine \
    --hdf-path ${HDF5_SRC_PATH} \
    --channel-name ${HDF5_CFOS_CHANNEL_NAME} \
    --affine-name  structural_to_template \
    --affine-path ${DIR_TRANSFORMS}/template_left_physical_to_structural_001_Affine.txt

sleep 1

lsfmpy write-affine \
    --hdf-path ${HDF5_SRC_PATH} \
    --channel-name ${HDF5_AUTO_CHANNEL_NAME} \
    --affine-name  structural_to_template \
    --affine-path ${DIR_TRANSFORMS}/template_left_physical_to_structural_001_Affine.txt

sleep 1

### This section saves the displacement fields ###

lsfmpy write \
    --hdf-path ${HDF5_SRC_PATH} \
    --channel-name forward_warp \
    --image-path `GetDisplacementFieldFilename structural 001 forward` \
    --metadata-path forward_warp.json \
    --is-multichannel True \
    --bdv-xml None 

lsfmpy write \
    --hdf-path ${HDF5_SRC_PATH} \
    --channel-name inverse_warp \
    --image-path `GetDisplacementFieldFilename structural 001 inverse` \
    --metadata-path inverse_warp.json \
    --is-multichannel True \
    --bdv-xml None 

### We can also save atlas segmentation mapped to each channel ###

lsfmpy write \
    --hdf-path ${HDF5_SRC_PATH} \
    --channel-name cfos_segmentation \
    --image-path `GetResultImageFilename segmentation signal-001-25 deformable` \
    --metadata-path cfos_segmentation.json \
    --is-segmentation True \
    --bdv-xml None 

lsfmpy write \
    --hdf-path ${HDF5_SRC_PATH} \
    --channel-name auto_segmentation \
    --image-path `GetResultImageFilename segmentation structural-001-25 deformable` \
    --metadata-path auto_segmentation.json \
    --is-segmentation True \
    --bdv-xml None
 
rm -rfv {cfos_segmentation,auto_segmentation,forward_warp,inverse_warp}.json
