
![Python](https://img.shields.io/badge/python-v3.8+-blue.svg)
[![pypy](https://badge.fury.io/py/lsfmpy.svg)](https://badge.fury.io/py/lsfmpy)
[![DOI](https://img.shields.io/badge/DOI-10.18150%2FNIDUBWC-informational)](https://doi.org/10.18150/NIDUBW)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


<img src="https://raw.githubusercontent.com/Neuroinflab/lsfm_image_server/lsfm_schema/media/thumbnail.png" align="left" alt="logo" width="70" height="70"/>

LSFMPy - python library for processing 3D light sheet fluorescence microscopy images
=====================

<p align="center"><img width=85% src=https://github.com/Neuroinflab/lsfm_image_server/blob/lsfm_schema/media/lsfmpy-github-fig-contents.png?raw="true"></img></p>

### Contents

 * [Basic Overview](#basic-overview)
 * [Setup](#setup)
 * [Image conversion](#image-conversion)
 * [Image registration](#image-registration)
 * [Usage examples](#usage-examples)

### Basic Overview

 

### Setup
```bash
pip install lsfmpy
```
Three new commands should be now available:
```text
lsfmpy
lsfmpy_create_autocomplete
dump_metadata
```
To enable the command line autocompletion, run:
```bash
lsfmpy_create_autocomplete --output-file autocomplete_lsfmpy.sh
source autocomplete_lsfmpy.sh
```
Additionally, you will need the [Advanced Normalization Tools](https://picsl.upenn.edu/software/ants/) (ANTS) and [Convert 3D](http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D) software packages available in your PATH environmental variable.

### Image conversion


### Image registration


### Usage examples
