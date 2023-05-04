&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
![Python](https://img.shields.io/badge/python-v3.8+-blue.svg)
[![pypy](https://badge.fury.io/py/lsfmpy.svg)](https://badge.fury.io/py/lsfmpy)
[![DOI](https://img.shields.io/badge/DOI-10.18150%2FNIDUBWC-informational)](https://doi.org/10.18150/NIDUBW)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) 


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://raw.githubusercontent.com/Neuroinflab/lsfm_image_server/lsfm_schema/media/thumbnail.png" align="left" alt="logo" width="90" height="90"/>
Tutorial for the LSFMPy library
=====================

<p> &nbsp;</p>

### Contents

 * [Overview](#overview)
 * [Technical notes](#technical-notes)
 * [Setup & Installation](#setup)
 * [Sample data download](#sample-data)
 * [The structure of the tutorial package](#the-structure-of-the-tutorial-package)
 * [Image conversion to the LSFMPy format](#image-conversion) 
 * [Image registration to the Allen Mouse Common Coordinate Framework v3 (CCF v3) template](#image-registration)
 * [Usage examples](#usage-examples)

### Overview
The following tutorial introduces the LSFMPy package and the accompanying registration pipeline.
Data files used in this tutorial can be downloaded from the [Repository for Open Data](https://doi.org/10.18150/NIDUBW)
(mind the
substantial size of the dataset: 63 GB). Images of a mouse brain hemisphere were acquired with
a light sheet fluorescence microscope at two different excitation wavelengths, resulting in two com-
plementary imaging channels. The first one represents the brain’s autofluorescence (488 nm), which
reflects its anatomical structure. The second image shows the fluorescence from immunostained c-Fos
transcription factor (638 nm), a known neuroplasticity marker.
 
### Technical notes

* Basic familiarity with command line interfaces and with the Linux operating system is advised
* For viewing nifti files, we recommend the [ITK-SNAP](http://www.itksnap.org/) application
* For viewing hdf5 files created with the LSFMPy, we recommend standard tools such as [HD-
FView](https://www.hdfgroup.org/downloads/hdfview/) or [ViTables](https://vitables.org/)
* The default anatomical orientation for LSFMPy datasets and associated software follows the RAS
(right-anterior-superior) convention
* Removing large groups or datasets from an hdf5 file will not immediately free up disk space
occupied by the file; only after using the [h5repack](https://manpages.ubuntu.com/manpages/lunar/man1/h5repack.1.html) tool will the disk space usage be reduced
to reflect the updated hdf5 file contents

### Setup & Installation
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

### The structure of the tutorial package:

The `example_data` folder contains two subfolders with image series of autofluorescence and cfos signal
of an entire mouse brain hemisphere imaged in the sagittal plane [Fig. 1A-B](). The TIFF files have
been already stitched with the freely available `BigStitcher` (Hörl et al., 2019) plugin for Fiji ImageJ
and denoised with the `VSNR algorithm` (Fehrenbach et al., 2012). For image registration, the
00_reference folder contains preprocessed template (`template.nii.gz`) and segmentation of the
Allen Mouse CCF v3 at 25 µm voxel resolution (`segmentation.nii.gz`), together with a label
description file (`labels.txt`). For your convenience, the template mask has been already manually
adjusted to match the brain being registered [Fig. 1C](). In the `01_input_data` directory, where
00_signal corresponds to the cfos channel, and 01_structural corresponds to the autofluorescence
channel, you will find masks for each imaging channel.
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
