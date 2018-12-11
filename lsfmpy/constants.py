#!/usr/bin/env python
# -*- coding: utf-8 -*
###############################################################################
#                                                                             #
#    LSFM image server                                                        #
#                                                                             #
#    Copyright (C) 2017-2018 Sylwia Bednarek, Piotr Majka                     #
#    (Laboratory of Neuroinformatics; Nencki Institute of Experimental        #
#    Biology of Polish Academy of Sciences)                                   #
#                                                                             #
#    This software is free software: you can redistribute it and/or modify    #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This software is distributed in the hope that it will be useful,         #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this software.  If not, see http://www.gnu.org/licenses/.     #
#                                                                             #
###############################################################################

import numpy as np

M_ORIGIN = 'origin'
M_SHAPE = 'shape'
M_VOXEL_SIZE = 'spacing'

RAS_INVERSE_DIRECTION = np.array([-1, -1, 1])
RAS_DIRECTION = np.array([1, 1, -1])

INVERSE_AFFINE_SIGN = np.array([[1, 1, -1, -1],
                                [1, 1, -1, -1],
                                [-1, -1, 1, 1],
                                [1, 1, 1, 1]])

SIGN_LPI_RAS_T = np.array([1, 1, -1])
SIGN_LPI_RAS_A = np.array([[-1, 1, -1], [1, -1, -1], [1, 1, -1]])

SIGN_RAS_T = np.array([-1, -1, -1])
SIGN_RAS_A = np.array([[-1, 1, -1],
                       [1, -1, -1],
                       [1, 1, -1]])
DIRECTION_RAS = (1.0, 0.0, 0.0,
                 0.0, 1.0, 0.0,
                 0.0, 0.0, -1.0)

SIGN_LPI_T = np.array([-1, -1, 1])
SIGN_LPI_A = np.array([[1, 1, -1], [1, 1, -1], [-1, -1, 1]])
DIRECTION_LPI = (-1.0, 0.0, 0.0,
                 0.0, -1.0, 0.0,
                 0.0, 0.0, 1.0)