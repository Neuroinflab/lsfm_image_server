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