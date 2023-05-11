#!/usr/bin/env python
# -*- coding: utf-8 -*

import os
import h5py
import numpy as np

from .utils import PathUtil


class LSFMInfo(object):
    """

    """

    def __init__(self, file_path):

        assert os.path.isfile(file_path)

        self._file_path = file_path
        self._f = h5py.File(file_path, 'r')
        self._lsfm = self._f['/LSFM']

    @property
    def case_ids(self):
        return list(self._lsfm.keys())

    @staticmethod
    def _print_attrs(name, obj):
        print(name)
        for key, val in list(obj.attrs.items()):
            print(("    %s: %s" % (key, val)))

    def print_pyramid(self, case_id):
        case = self._f[PathUtil.get_lsfm_image_path(case_id)]
        print("              ")
        print("     /\.      ")
        print("    /  \`.    ")
        print("   /    \ `.  ")
        print("  /______\/ \n")

        print(("Resolution pyramid structure for: ", case.name))
        origin = np.array(case.attrs['origin'],
                          dtype=np.float32)
        shape = np.array(case.attrs['shape'],
                         dtype=np.int)
        spacing = np.array(case.attrs['spacing'],
                            dtype=np.float32)
        phys_size = shape[:,:3] * spacing
        n_levels = origin.shape[0]

        print(("number of levels: ", n_levels))
        for i in range(n_levels):
            print("\n")
            print(("level no. {}".format(i)))
            print(("shape:         {: <7} {: <7} {: <7}".format(shape[i][0], shape[i][1], shape[i][2])))
            print(("spacing:       {:.5f} {:.5f} {:.5f}".format(spacing[i][0], spacing[i][1], spacing[i][2])))
            print(("origin:        {:+.4f} {:+.4f} {:+.4f}".format(origin[i][0], origin[i][1], origin[i][2])))
            print(("physical size: {:.5f} {:.5f} {:.5f}".format(phys_size[i][0], phys_size[i][1], phys_size[i][2])))
            print(("*" * 50))


    def show_structure(self):
        self._f.visititems(LSFMInfo._print_attrs)

