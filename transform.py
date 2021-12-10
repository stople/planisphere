#!/usr/bin/python3
# solarsystem_log.py
# -*- coding: utf-8 -*-
#


# ----------------------------------------------------------------------------

"""
Vector transform functions
"""

import math


# Rotate 2d
#
# Rotate the 2d point [x, y] around origo with the angle (a) specified in radians.
#
# The point can be expressed with polar coordinates, radius r and angle t.
# x = r * cos(t), y = r * sin(t).
#
# The rotated point have the same r, but angle (t) have been incremented with a.
#
# Rotated point:
# x' = r * cos(t + a), y' = r * sin(t + a)
# Using trigonometric identities, this becomes:
# x' = r * [cos(t) * cos(a) - sin(t) * sin(a)], y' = r * [sin(t) * cos(a) + cos(t) * sin(a)]
# Using the expressions for x and y, these identities are found:
# cos(t) = x / r, sin(t) = y / r
# This results in:
# x' = r * [x / r * cos(a) - y / r * sin(a)], y' = r * [y / r * cos(a) + x / r * sin(a)]
# x' = x * cos(a) - y * sin(a), y' = y * cos(a) + x * sin(a)
def rotate_2d(pos, angle):
    x = pos[0]
    y = pos[1]
    x_rot = x * math.cos(angle) - y * math.sin(angle)
    y_rot = y * math.cos(angle) + x * math.sin(angle)
    return [x_rot, y_rot]


def rotate_3d_x(pos, angle):
    yz = [pos[1], pos[2]]
    xy_rot = rotate_2d(yz, angle)
    return [pos[0], xy_rot[0], xy_rot[1]]


def rotate_3d_y(pos, angle):
    zx = [pos[2], pos[0]]
    zx_rot = rotate_2d(zx, angle)
    return [xy_rot[2], pos[1], xy_rot[0]]


def rotate_3d_z(pos, angle):
    xy = [pos[0], pos[1]]
    xy_rot = rotate_2d(xy, angle)
    return [xy_rot[0], xy_rot[1], pos[2]]


def add_2d(pos1, pos2):
    xy = [pos1[0] + pos2[0], pos1[1] + pos2[1]]
    return xy


def sub_2d(pos1, pos2):
    xy = [pos1[0] - pos2[0], pos1[1] - pos2[1]]
    return xy


def add_3d(pos1, pos2):
    xyz = [pos1[0] + pos2[0], pos1[1] + pos2[1], pos1[2] + pos2[2]]
    return xyz


def sub_3d(pos1, pos2):
    xyz = [pos1[0] - pos2[0], pos1[1] - pos2[1], pos1[2] - pos2[2]]
    return xyz


# Cartesian coordinates to polar coordinates
#
# Angle: Radians, 0 <= a < 2 * pi
def cartesian_to_polar_2d(pos):
    x = pos[0]
    y = pos[1]

    angle = math.atan2(y, x)  # between -PI and PI
    if angle < 0:
        angle += math.pi * 2  # between 0 and 2 PI

    radius = math.sqrt(x * x + y * y)

    return [angle, radius]


# Cartesian coordinates to polar coordinates
#
# Angle: Radians, 0 <= a < 2 * pi
def polar_to_cartesian_2d(pos):
    angle = pos[0]
    radius = pos[1]

    x = math.cos(angle) * radius
    y = math.sin(angle) * radius

    return [x, y]


# Astronomical functions below

# Input:  Longitude (0 - 360)
#         Latitude  (-90 - +90)
#         Radius    (any unit)
# Output: xyz       (same unit as radius)

def longitude_latitude_distance_to_cartesian(pos):
    longitude_rad = math.radians(pos[0])
    latitude_rad = math.rad(pos[1])
    distance = pos[2]

    # 1: Get Z
    z = distance * math.sin(latitude_rad)

    # 2: Get xy
    pos_xy = polar_to_cartesian_2d([longitude_rad, distance])
    x = pos_xy[0]
    y = pos_xy[1]

    return [x, y, z]


# Latitude/longitude, in degrees
def cartesian_to_latitude_longitude_dist(pos):
    x = pos[0]
    y = pos[1]
    z = pos[2]

    # 1: Calculate longitude (xy plane) (0 <= long < 360)
    xy_polar = cartesian_to_polar_2d([x, y])
    longitude = math.degrees(xy_polar[0])
    xy_length = xy_polar[1]

    # 2: Calculate latitude (-90 <= lat <= 90)
    xy_z_polar = cartesian_to_polar_2d([xy_length, z])
    latitude = math.degrees(xy_z_polar[0])
    if latitude > 180:
        latitude -= 360
    radius = xy_z_polar[1]

    return [longitude, latitude, radius]

# Input:  Equatorial position, cartesian
# Output: Ecliptic position, cartesian

def equatorial_to_ecliptic_heliocentric(equatorial):
    # Points toward vernal equinox

    ecliptic = rotate_3d_x(equatorial, )
