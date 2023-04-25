
![Python](https://img.shields.io/badge/python-v3.8+-blue.svg)
[![pypy](https://badge.fury.io/py/lsfmpy.svg)](https://badge.fury.io/py/lsfmpy)
[![DOI](https://img.shields.io/badge/DOI-10.18150%2FNIDUBWC-informational)](https://doi.org/10.18150/NIDUBW)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


<img src="https://raw.githubusercontent.com/Neuroinflab/lsfm_image_server/lsfm_schema/media/thumbnail.png" align="left" alt="logo" width="70" height="70"/>

LSFMPy - python library for processing 3D light sheet fluorescence microscopy images
=====================

<p align="center"><img width=55% src=https://github.com/Neuroinflab/lsfm_image_server/blob/lsfm_schema/media/lsfmpy-github-fig-contents.png?raw="true"></img></p>

### Contents

 * [Basic Overview](#basic-overview)
 * [Setup](#setup)
 * [Sample data](#sample-data)
 * [Image conversion](#image-conversion)
 * [Image registration](#image-registration)
 * [Usage examples](#usage-examples)

### Basic Overview
The primary motivation for creating this package is to enable efficient manipulation of large, three-dimensional image datasets whilst preserving precise spatial information and facilitating spatial transformations. This is particularly useful in neuroscience, where the additional layer of information about localization in specific brain regions is crucial for proper data interpretation and analyzis.
 

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

### Image conversion


### Image registration


### Usage examples
