
![Python](https://img.shields.io/badge/python-v3.8+-blue.svg)
[![pypy](https://badge.fury.io/py/lsfmpy.svg)](https://badge.fury.io/py/lsfmpy)
[![DOI](https://img.shields.io/badge/DOI-10.18150%2FNIDUBWC-informational)](https://doi.org/10.18150/NIDUBW)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


<img src="https://raw.githubusercontent.com/Neuroinflab/lsfm_image_server/lsfm_schema/media/thumbnail.png" align="left" alt="logo" width="70" height="70"/>

LSFMPy - python library for processing 3D light sheet fluorescence microscopy images
=====================
### Basic Overview

<p align="center"><img width=75% src=https://github.com/Neuroinflab/lsfm_image_server/blob/lsfm_schema/media/lsfmpy-github-fig-contents.png?raw="true"></img></p>
 

### Last Stable Release
```bash
pip install lsfmpy
```


Functionality:
 
--------------

1. Writing:

* a series of tiffs, a  stack of images in a single tiff or a nifti file (with RAS orientation) can be compressed to hdf5

  * channel name is required
  * a json file with metadata is required
  * an xml file for big data viewer compatibility will be created (and recreated every time a new channel is added)
  * user can specify how much memory should be used (for file reading, in GB) at once
  * deafult 2 GB
  * does not take into account memory for processing this amount of data
  * an affine from a file (itk compatible) can be written for a specified channel


2. Exporting:

* 3D images

  * reoriented from original orientation to RAS
  * resampled (down- or upsampled) to desired resolution
  * from desired level of internal pyramid of resolutions
  * works for single compononet images, multicomponent images, segmentation
  * affine and displacement fields can be specified for transformation of image during export. the order, type and inversion of transforms must be given. many transforms can be chained (a composite transform will be created and executed)

  * if channel segmentation (from registration) is present in the hdf file, it can be used to specify which region (anatomical structure) should be exported (a bounding box), this also can be transformed to a reference space during export

  * a subset of image volume, in the form of a cuboid, with origin specified in the output space (i.e. reference space if transformations are required) and of given physical size and resolution (in mm) can also be exported
  * if a grid (number of chunks in each direction in output space) and overlap are also given when exporting a subset of image volume, appropriate grid of chunks will be exported (doesn't work with transformations)
  * if grid = 0,0,0 whole image volume will be exported as a grid of chunks (origin parameter is ignored)


slices

* slices along a given (0,1,2) axis can be exported from pre-resampled levels, with	given range

* if range is None, all slices will be exported

* roi can be specified via ox,oy,sx,sy (origin and size) in pixel coordinates


Export functions can either write data to specified path name, or if output_path isNone, will return the data (or generator for multiple chunks/slices)


viewing:

* hdf5 files can be viewed with BigDataViewer plugin for Fiji via accompanying xml file


info:

* informations about the data written to file (or for specified channel) can be printed


metadata:

* arbitrary metadata (as a key-value pair) can be added to file

provenance:

* all operations on hdf5 file, along with command line, will be logged (unix systems)

interactive:

* console commands with -- --interactive will return data in IPython notebook

dump_metadata:

* some metadata can be recovered from some of the tiff/nifti files, and a json file will be created, it needs to be inspected, and missing values must be filled in

h5_external:
* hdf5 files can be repacked, compressed or decompressed with this utility (repacking is useful for recovering space after deletion of large amounts of data that 	otherwise will not be freed)
