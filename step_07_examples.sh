#!/bin/bash -xe

source defaults.sh

HDF5_SRC_PATH=../${CASE_ID}.h5
EXAMPLE_DIR=07_examples

mkdir -p ${EXAMPLE_DIR}


#Autofluorescence channel mapped to template at 10 µm:
lsfmpy \
    add-transform \
	--name inverse_warp \
	--ordn 0 \
	--transform-type df \
    - add-transform \
	--name structural_to_template \
	--ordn 1 \
	--transform-type affine \
	--invert True \
    - export \
	--channel-name autofluo \
	--output-resolution 0.01,0.01,0.01 \
	--output-path ${EXAMPLE_DIR}/${CASE_ID}_auto_in_template_10um.nii.gz \
	--hdf-path ${HDF5_SRC_PATH} \
	--input-orientation RPI

#Cfos channel mapped to template at 10 µm:
lsfmpy \
    add-transform \
	--transform-type affine \
	--name signal_to_structural \
	--ordn 0 \
	--invert True \
    - add-transform \
	--name inverse_warp \
	--ordn 1 \
	--transform-type df \
    -  add-transform \
	--name structural_to_template \
	--ordn 2 \
	--transform-type affine \
	--invert True \
    - export \
	--channel-name cfos \
	--output-resolution 0.01,0.01,0.01 \
	--output-path ${EXAMPLE_DIR}/${CASE_ID}_cfos_in_template_10um.nii.gz \
	--hdf-path ${HDF5_SRC_PATH} \
	--input-orientation RPI

#Basolateral Amygdala, anterior part in signal channel at 10 µm:
lsfmpy \
    export \
	--hdf-path ${HDF5_SRC_PATH} \
	--channel-name cfos \
	--output-path ${EXAMPLE_DIR}/${CASE_ID}_cfos_BLAa_10um.nii.gz \
	--output-resolution 0.01,0.01,0.01 \
	--segmentation-name cfos_segmentation \
	--region_id 303 \
	--input-orientation RPI

#Anterodorsal thalamic nucleus in signal channel mapped to template at 2 µm:
lsfmpy \
    add-transform \
	--transform-type affine \
	--name signal_to_structural \
	--ordn 0 \
	--invert True \
    - add-transform \
	--name inverse_warp \
	--ordn 1 \
	--transform-type df \
    -  add-transform \
	--name structural_to_template \
	--ordn 2 \
	--transform-type affine \
	--invert True \
    - export \
	--channel-name cfos \
	--output-resolution 0.002,0.002,0.002 \
	--output-path ${EXAMPLE_DIR}/${CASE_ID}_cfos_in_template_AD_2um.nii.gz \
	--hdf-path ${HDF5_SRC_PATH} \
	--input-orientation RPI \
	--segmentation-name cfos_segmentation \
	--region_id 64

# Native resolution chunk at arbitrary physical location
lsfmpy \
    export \
	--hdf-path ${HDF_SRC_PATH} \
	--channel-name cfos \
	--input-orientation RPI \
	--phys-origin 1.76,6.91,2.01 \
	--phys-size 0.1,0.1,0.1 \
	--input-resolution-level 0 \
	--output-path ${EXAMPLE_DIR}/${CASE_ID}_cfos_chunk_native.nii.gz

