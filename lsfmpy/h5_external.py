#!/usr/bin/env python
# -*- coding: utf-8 -*


import argparse
import subprocess

from shutil import copyfile
from shutil import rmtree
from distutils import spawn


class H5_Repacker(object):

    _CMD_REPACK = 'h5repack'
    _CMD_COMPRESS = 'GZIP=6'
    _CMD_DECOMPRESS = 'NONE'

    @staticmethod
    def tools_available():
        if spawn.find_executable('h5repack') is None:
            return False
        return True

    @classmethod
    def repack(cls, fname_in, fname_out):
        """

        :param fname_in:
        :return:
        """

        if H5_Repacker.tools_available():
            subprocess.Popen([cls._CMD_REPACK, fname_in, fname_out]).wait()

    @classmethod
    def compress(cls, fname_in, fname_out):
        """

        :param fname_in:
        :return:
        """
        if H5_Repacker.tools_available():
            subprocess.Popen([cls._CMD_REPACK, '-f', cls._CMD_COMPRESS, fname_in, fname_out]).wait()

    @classmethod
    def decompress(cls, fname_in, fname_out):
        """

        :param fname_in:
        :return:
        """
        if H5_Repacker.tools_available():
            subprocess.Popen([cls._CMD_REPACK, '-f', cls._CMD_DECOMPRESS, fname_in, fname_out]).wait()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Repack HDF5 file to delete data from it, or to compress it',
                                     prog='h5_external')

    parser.add_argument('--input-file',
                        required=True,
                        help="Path to the hdf file",
                        type=str)

    parser.add_argument('--output-file',
                        required=True,
                        default=None,
                        help="Path to output hdf5 file",
                        type=str)

    compression_parser = parser.add_mutually_exclusive_group()

    compression_parser.add_argument('--compress',
                        action='store_true')

    compression_parser.add_argument('--decompress',
                        action='store_true')

    compression_parser.set_defaults(compress=False,
                        decompress=False)

    args = parser.parse_args()

    if args.compress:
        H5_Repacker.compress(args.input_file,
                             args.output_file)
    elif args.decompress:
        H5_Repacker.decompress(args.input_file,
                               args.output_file)
    else:
        H5_Repacker.repack(args.input_file,
                           args.output_file)
