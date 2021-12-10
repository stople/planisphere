#!/usr/bin/python3
# solarsystem_log.py
# -*- coding: utf-8 -*-
#


# ----------------------------------------------------------------------------

"""
Vector transform functions
"""

import cairocffi as cairo


# import math
from math import sin, cos
import numpy as np

# Cairo matrices is 2x3
# Extend this with a 3x4 matrix

import cairocffi.matrix

# convert between a full 3d matrix and the 2d matrix used by cairo


class Matrix3d:
    #def __init__(self):

    matrix = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

    # Transform by applying transformation matrix
    # Transformation matrix must be 4x4 elements
    def transform(self, new_matrix):
        self.matrix = np.dot(new_matrix, self.matrix)

    def rotate_x(self, radians):
        rot_x_matrix = np.array([
            [1,            0,             0, 0],
            [0, cos(radians), -sin(radians), 0],
            [0, sin(radians),  cos(radians), 0],
            [0,            0,             0, 1]
        ])
        self.transform(rot_x_matrix)

    def rotate_y(self, radians):
        rot_y_matrix = np.array([
            [ cos(radians), 0, sin(radians), 0],
            [            0, 1,            0, 0],
            [-sin(radians), 0, cos(radians), 0],
            [            0, 0,            0, 1]
        ])
        self.transform(rot_y_matrix)

    def rotate_z(self, radians):
        rot_z_matrix = np.array([
            [cos(radians), -sin(radians), 0, 0],
            [sin(radians),  cos(radians), 0, 0],
            [           0,             0, 1, 0],
            [           0,             0, 0, 1]
        ])
        self.transform(rot_z_matrix)

    def flip_y(self):
        flip_y_matrix = np.array([
            [1,  0, 0, 0],
            [0, -1, 0, 0],
            [0,  0, 1, 0],
            [0,  0, 0, 1]
        ])
        self.transform(flip_y_matrix)

    def scale(self, x, y, z):
        scale_matrix = np.array([
            [x, 0, 0, 0],
            [0, y, 0, 0],
            [0, 0, z, 0],
            [0, 0, 0, 1]
        ])
        self.transform(scale_matrix)

    def translate(self, x, y, z):
        translate_matrix = np.array([
            [1, 0, 0, x],
            [0, 1, 0, y],
            [0, 0, 1, z],
            [0, 0, 0, 1]
        ])
        print(self.matrix)
        print(translate_matrix)
        self.transform(translate_matrix)
        print(self.matrix)

    def translate_old(self, x, y, z):
        translate_matrix = np.array([
            [0, 0, 0, x],
            [0, 0, 0, y],
            [0, 0, 0, z],
            [0, 0, 0, 0]
        ])
        print(self.matrix)
        print(translate_matrix)
        self.matrix = self.matrix + translate_matrix
        print(self.matrix)

    #  xx | xy | xz | x0
    #  yx | yy | yz | y0
    #  zx | zy | xz | z0
    #   0 |  0 |  0 |  1

    def get_cairo_matrix(self):
        m = self.matrix
        xx = m[0, 0]
        xy = m[0, 1]
        x0 = m[0, 3]
        yx = m[1, 0]
        yy = m[1, 1]
        y0 = m[1, 3]

        matrix = cairocffi.matrix.Matrix(xx=xx, yx=yx, xy=xy, yy=yy, x0=x0, y0=y0)
        return matrix
