&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
![Python](https://img.shields.io/badge/python-v3.8+-blue.svg)
[![pypy](https://badge.fury.io/py/lsfmpy.svg)](https://badge.fury.io/py/lsfmpy)
[![DOI](https://img.shields.io/badge/DOI-10.18150%2FNIDUBWC-informational)](https://doi.org/10.18150/NIDUBW)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) 


<img src="https://raw.githubusercontent.com/Neuroinflab/lsfm_image_server/lsfmpy_dev/media/thumbnail.png" align="left" alt="logo" width="90" height="90"/>

LSFMPy - python library for processing 3D light sheet fluorescence microscopy images
=====================

<p align="center"><img width=65% src=https://github.com/Neuroinflab/lsfm_image_server/blob/lsfmpy_dev/media/lsfmpy-github-fig-contents.png?raw="true"></img></p>

### Contents

 * [Basic Overview](#basic-overview)
 * [Setup](#setup)
 * [Sample data](#sample-data)
 * [Image conversion](#image-conversion)
 * [Image registration](#image-registration)
 * [Usage examples](#usage-examples)

### Basic Overview
The primary motivation for creating this package is to enable efficient manipulation of large, three-dimensional image datasets whilst preserving precise spatial information and facilitating spatial transformations. This is particularly useful in neuroscience, where the additional layer of information about localization in specific brain regions is crucial for proper interpretation and integration of experimental data.
 

### Setup
```bash
pip install lsfmpy
```
- Three new commands should be now available:
```text
lsfmpy
lsfmpy_create_autocomplete
dump_metadata
```
- To enable the command line autocompletion, run:
```bash
lsfmpy_create_autocomplete --output-file autocomplete_lsfmpy.sh
source autocomplete_lsfmpy.sh
```
Additionally, you will need the [Advanced Normalization Tools](https://picsl.upenn.edu/software/ants/) (ANTS) and [Convert 3D](http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D) software packages available in your PATH environmental variable.

### Sample data
- You can download sample data from [Repository for Open Data](https://doi.org/10.18150/NIDUBW) and follow our tutorial to process it:
```bash
git clone -b tutorial https://github.com/Neuroinflab/lsfm_image_server.git
cd lsfm_image_server
export TUTORIAL=${PWD}
```
- Download data with curl (around 60 GB):
```bash
bash -xe step_00_download_data.sh
```
- Move data and unpack:
```bash
mkdir ${TUTORIAL}/example_data/autofluorescence
mv auto_channel_part_*.zip ${TUTORIAL}/example_data/autofluorescence
mkdir ${TUTORIAL}/example_data/cfos
mv cfos_channel_part_*.zip ${TUTORIAL}/example_data/cfos
cd ${TUTORIAL}/example_data/autofluroescence
unzip *. zip
cd ${TUTORIAL}/example_data/cfos
unzip *. zip
```
- Clean up:
```bash
# rm -f ${TUTORIAL}/example_data/autofluorescence/auto_channel_part_*.zip 
# rm -f ${TUTORIAL}/example_data/cfos/cfos_channel_part_*.zip
```
- The structure of the tutorial package:

```bash
lsfmpy_tutorial/
└── defaults.sh		# set up constants and logging
└── case_id 		   # unique specimen ID
└── step_00_download_data.sh 
└── step_01_set_up_images.sh 
└── step_02_preprocess.sh
└── step_03_signal_and_structural_coregistration.sh
└── step_04_set_up_template_coregistration.sh
└── step_05_core_registration.sh
└── step_06_write_to_hdf.sh
└── step_07_examples.sh
└── functions.sh 	# convenience functions for setting registration fidelity
├── scripts       	# python scripts (should work with python 2.7 and 3.8)
│   └── header_converter.py
│   └── clahe.py
│   └── rescale_displacement_field.py
│   └── tqdm.py
├── example_data/
│   └── autofluorescence/
│   │   └── Z000000.tif
│   │   └── ...
│   └── cfos/
│   │   └── Z000000.tif
│   │   └── ...
└── 00_reference/
│   └── template.nii.gz 
│   └── segmentation.nii.gz 
│   └── template_25um_both_hemispheres.nii.gz 
│   └── segmentation_25um_both_hemispheres.nii.gz 
│   └── labels.txt # labels for segmentation
│   └── template_mask.nii.gz  	# template mask clipped to match the extent of the target brain image
└── 01_input_data/
    └── 00_signal/
    │   └── 001_25_mm_mask.nii.gz # mask outlining the brain matter in signal channel
    └── 01_structural/
        └── 001_25_mm_mask.nii.gz # mask outlining the brain matter in autofluorescence channel
	
```
### Image conversion
**Prepare the metadata json file**

You will be asked to fill the following information for the autofluorescence image:
- image_size_z 490
- voxel_size_y 0.00145
- voxel_size_x 0.00145
- voxel_size_z 0.01
```bash
cd ${TUTORIAL}
dump_metadata \
	--input-file ./example_data/autofluorescence/Z000000.tif \
	--output-file ./autofluo.json
```
You will be asked to fill the following information for the cfos image:
- image_size_z 1265
- voxel_size_y 0.00145
- voxel_size_x 0.00145
- voxel_size_z 0.004
```bash
cd ${TUTORIAL}
dump_metadata \
	--input-file ./example_data/cfos/Z000000.tif \
	--output-file ./cfos.json
```
**Image conversion**

Let's start with the autofluorescence channel:
```bash
lsfmpy \
    # use the write command
    write \
        # path to the target hdf file,
	  # which will be created if not found
        --hdf-path ./000001.h5 \
	  # path to image metadata
        --metadata-path ./autofluo.json \
	  #  file name pattern 
        --file-name-format Z%06d.tif \
	  # channel name for this modality
        --channel-name autofluo \
	  # path to first file in the sequence
        --image-path ./example_data/autofluorescence/Z000000.tif \
   	  # file name for BigDataViewer-compatible xml
        --bdv-xml 000001.xml
```

```bash
lsfmpy write \
        --hdf-path ./000001.h5 \
        --metadata-path ./cfos.json \
        --file-name-format Z%06d.tif \
        --channel-name cfos \
        --image-path ./example_data/cfos/Z000000.tif \
        --bdv-xml 000001.xml
```
### Image registration

1. Export images at 25 µm isotropic resolution:
```bash
bash step_01_set_up_images.sh 
```

2. Run image preprocessing before the registration:
```bash
bash step_02_preprocess.sh
```

3. Co-register the two channels to make sure that we know the spatial transformations between the structural (autofluorescence) and signal (cfos) channels:

```bash
bash step_03_signal_and_structural_coregistration.sh
```
4. Finally, register the structural channel to the template.
You can choose from three different default settings, determining registration fidelity and speed: 
- lsfm_set_mapping_affine       - coregistration with atlas set to affine only
- lsfm_set_mapping_overnight    - coregistration with atlas set to low quality deformable
- lsfm_set_mapping_high_quality - coregistration with atlas set to high quality deformable

```bash
source functions.sh

# load setting for high quality mapping:
lsfm_set_mapping_high_quality
bash step_04_set_up_template_coregistration.sh
# before proceeding with registration, the script will check if following executables are available:
# mkdir rm cp tar WarpImageMultiTransform c3d antsRegistration
# this will also execute step_05_core_registration.sh
```
```shell
lsfmpy_tutorial/
└── example_data/	
└── 00_reference/	# preprocessed allen ccfv3	
└── 01_input_data/	# images exported for registration
└── workdir_fixed_{md5}_moving_{md5}/ 
└── 02_preprocessed_input_images/ # preprocessed images
└── 03_coregistration/ 		# signal channel mapped to structural channel
└── 04_transforms/		# affine and deformable transformations
└── 05_results/			# both channels mapped to template, 
				# template and segmentation mapped to the physical space of 
				# signal and structural image
```

5. After reviewing the result, save the computed transforms and the segmentation to the LSFMPy hdf file:
```bash
bash step_06_write_to_hdf.sh
```
### Usage examples

1.  Export the autofluorescence channel at 10 µm resolution, mapped to template

<p align="center"><img width=65% src=https://raw.githubusercontent.com/Neuroinflab/lsfmpy_image_server/lsfm_dev/media/auto_in_temp_10um.gif></img></p>

```bash
lsfmpy \
    # you can iteratively add multiple transformations
    # both affine and deformable, indicating the order
    # in which they should be executed 
    add-transform \
	    # name of the transformation as it was written to the hdf file
          --name inverse_warp \
          # indicate the position of the transformation in the chain:
          --ordn 0 \
          # indicate type of transformation (deformable)
    --transform-type df \


    # add another transformation
    - add-transform \
    --name structural_to_template \
    # this is second transformation in chain
    --ordn 1 \
    # this is an affine transformation
    --transform-type affine \
    # and it should be inverted
    --invert True \
	


    # now we can export the data
    - export \
	    # from the autofluorescence channel
    --channel-name autofluo \
    # at 10 µm isotropic resolution
    --output-resolution 0.01,0.01,0.01 \
    # indicate output file
    --output-path ./000001_auto_in_template_10um.nii.gz \
    # indicate source hdf
    --hdf-path ./000001.h5 \
    # indicate anatomical orientation of the image
    # (if you are not sure, set the orientation to RAS
          # and check results in ITKSnap)
    --input-orientation RPI

```
2. Export Anterodorsal nucleus region resampled to 2 µm resolution in signal channel, mapped to template

<p align="center"><img width=65% src=https://raw.githubusercontent.com/Neuroinflab/lsfm_image_server/lsfmpy_dev/media/AD_cfos_in_template_2um.gif></img></p>

```bash
lsfmpy \
    # specify the mapping from signal channel to the template
    --add-transform \
        --transform-type affine \
        --name signal_to_structural \
        --ordn 0 \
        --invert True \
    –-add-transform \
        --name inverse_warp \
        --ordn 1 \
        --transform-type df \
    –-add-transform \
        --name structural_to_template \
        --ordn 2 \
        --transform-type affine \
        --invert True \
    # export Anterodorsal nucleus from the signal channel at 2 µm voxel size
    -export \
        --channel-name cfos \
        --output-resolution 0.002,0.002,0.002 \
        --output-path ./000001_cfos_in_template_AD_2um.nii.gz \
        --hdf-path ./000001.h5 \
        --input-orientation RPI \
        --segmentation-name cfos_segmentation \
        --region_id 64

```
