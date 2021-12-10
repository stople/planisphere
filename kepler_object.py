#!/usr/bin/python3
# solarsystem_log.py
# -*- coding: utf-8 -*-
#


# ----------------------------------------------------------------------------

"""
File containing the definition of all planets in the solar system

Dataset: "Keplerian elements" (https://ssd.jpl.nasa.gov/?planet_pos).
Currently only table 1 (1800 AD to 2050 AD) is implemented.
Data uses the ecliptic plane J2000.0 as reference.

"""



"""
Render the solar system with logarithmic scales
"""
import math
import transform

import re
from math import pi, sin, cos, atan2, hypot

from numpy import arange

import calendar
from bright_stars_process import fetch_bright_star_list
from constants import unit_deg, unit_rev, unit_mm, unit_cm, r_1, r_gap, central_hole_size, radius
from graphics_context import BaseComponent
from settings import fetch_command_line_arguments
from text import text
from themes import themes


J2000_JULIAN_DAY = 2451545
OBLIQUITY = 23.43928  # Earth's axial tilt


class KeplerianElement():
    def __init__(self, init_state, delta_pr_century, wrap):
        self.init_state = float(init_state)
        self.delta_pr_century = float(delta_pr_century)
        self.wrap = float(wrap)

    def get_at_delta_time(self, delta_time):
        ret = self.init_state + delta_time * self.delta_pr_century
        if self.wrap != 0:
            ret %= self.wrap
        return ret

    def get_at_J2000(self, julian_time):
        return self.get_at_delta_time((julian_time - J2000_JULIAN_DAY) / (100 * 365.25))


# Inverse kepler
#
# Kepler's equation:  M = E - e sin E
# E: Eccentric anomaly. M: Mean anomaly.
# This function calculates E if M is known, which cannot be done algebraicly.
# Using a numerical method, as found on Wikipedia (Kepler's equation)
def inverse_kepler(eccentricity, mean_anomaly):
    MAX_ITERATIONS = 100
    MIN_ACCURACY = 0.0001

    m_rad = math.radians(mean_anomaly)
    e_rad = m_rad

    for i in range(1, MAX_ITERATIONS):
        old_e_rad = e_rad
        e_rad = m_rad + eccentricity * math.sin(e_rad)
        if abs(e_rad - old_e_rad) < MIN_ACCURACY:
            break
    return math.degrees(e_rad)




class KeplerianObject():
    """
    Object following a Kepler orbit (major planets in the solar system)
    """

    def __init__(self, parameter_string):
        lines = parameter_string.splitlines()
        self.name = lines[0][0:9]

        line1data = lines[0][9:].split(sep=None, maxsplit=-1)  # Split string on spaces, ignore consecutive spaces
        line2data = lines[1].split(sep=None, maxsplit=-1)  # Split string on spaces, ignore consecutive spaces

        # Semi major axis (a). Unit: AU.
        #
        # Distance from ellipse center to border at maximum.
        #
        # (Step 1: Define the max radius of the orbit)
        self.semi_major_axis             = KeplerianElement(line1data[0], line2data[0], 0)

        # Eccentricity (e). Unit: None.
        #
        # Shape of orbit. 0: Perfect circle. Less than 1: Ellipse. (1: Parabolic. >1: Hyperbolic. Not applicable here)
        # Barycenter (one of the focus points) is e * a away from ellipse center (assuming ellipse, 0 <= e < 1)
        #
        # (Step 2: Define the shape of the orbit (as an ellipse, focal points 0-1 diameter apart)
        self.eccentricity                = KeplerianElement(line1data[1], line2data[1], 0)

        # Inclination (i). Unit: deg.
        #
        # Angle of orbit plane vs ecliptic plane.
        self.inclination                 = KeplerianElement(line1data[2], line2data[2], 360)

        # Mean longitude (L). Unit: deg.
        #
        # Longitude of the object (position) along reference plane (ecliptic plane),
        # relative to the reference direction (vernal equinox),
        # assuming no inclination and no eccentricity.
        #
        # EDIT: This seem to be "mean angle along orbit plane since mean angle pointed at ascending node"
        # EDIT2: Not so fast... there are 2 periapsis angles
        self.mean_longitude              = KeplerianElement(line1data[3], line2data[3], 360)

        # Longitude of periapsis (ω with strike over). Unit: deg.
        #
        # Longitude of periapsis (position where the object is closest to the barycenter)
        # along reference plane (ecliptic plane),
        # relative to the reference direction (vernal equinox),
        # assuming no inclination and no eccentricity.
        #
        # (Step 3:
        self.longitude_of_periapsis      = KeplerianElement(line1data[4], line2data[4], 360)

        # Longitude of the ascending node (Ω). Unit: deg.
        #
        # Longitude along reference plane (ecliptic plane),
        # relative to the reference direction (vernal equinox),
        # where the object is ascending through the reference plane.
        #
        # (Step 6: Define an axis in this direction from the sun. Tilt the planet's plane along this axis.)
        self.longitude_of_ascending_node = KeplerianElement(line1data[5], line2data[5], 360)


        # How to compute position from Keplerian elements:
        #
        # Find mean anomaly (M) (angle since last periapsis, with constant (average) angular velocity) using L and
        # Find eccentric anomaly (E) (angle, relating to orbit arc distance?) using M and e
        # Compute the cartesian coordinates in orbit plane, with x pointing to periapsis, using a, e and E
        # 1: Rotate plane around z with [argument_of_periapsis], point x toward Ascending node
        # 2: Tilt plane around x with [inclination], to the reference plane (ecliptic)
        # 3: Rotate plane around z with [longitude of ascending node], point x toward reference direction (vernal equinox)

    # Argument of periapsis (ω). Unit: deg.
    #
    # Angle along orbit plane, from ascending node to periapsis.
    def get_argument_of_periapsis(self, julian_time):
        longitude_of_periapsis = self.longitude_of_periapsis.get_at_J2000(julian_time)
        longitude_of_ascending_node = self.longitude_of_ascending_node.get_at_J2000(julian_time)
        periapsis_argument = longitude_of_periapsis - longitude_of_ascending_node
        return periapsis_argument % 360


    # Mean anomaly (M with strike over). Unit: deg
    #
    # Assuming a virtual (mean) orbit along the orbit plane, fixed at periapsis,
    # with no eccentricity (circular orbit).
    # Angle (anomaly) from periapsis to the object.
    def get_mean_anomaly(self, julian_time):
        mean_longitude = self.mean_longitude.get_at_J2000(julian_time)
        longitude_of_periapsis = self.longitude_of_periapsis.get_at_J2000(julian_time)
        mean_anomaly = mean_longitude - longitude_of_periapsis
        return ((mean_anomaly + 180) % 360) - 180


    # Eccentric anomaly (E). Unit: deg
    #
    # Angle (anomaly) from periapsis to the object along the reference plane, measured as following:
    # Draw a circle with center at the center of the ellipse, with radius equal to semi major axis (a).
    # Draw a line perpendicular to the major axis, crossing the object, and hitting the circle further out.
    # Measure the angle from periapsis, via center, to where the perpendicular line hits the circle.
    def get_eccentric_anomaly(self, julian_time):
        eccentricity = self.eccentricity.get_at_J2000(julian_time)
        mean_anomaly = self.get_mean_anomaly(julian_time)
        eccentric_anomaly = inverse_kepler(eccentricity, mean_anomaly)
        return eccentric_anomaly

    # Position (x' y' z') along own orbit plane. Unit: AU
    #
    # Center on barycenter.
    # x in direction of perihelion.
    def get_position_own_orbit(self, julian_time):
        semi_major_axis = self.semi_major_axis.get_at_J2000(julian_time)
        eccentric_anomaly = self.get_eccentric_anomaly(julian_time)
        eccentricity = self.eccentricity.get_at_J2000(julian_time)
        x = semi_major_axis * (math.cos(math.radians(eccentric_anomaly)) - eccentricity)
        y = semi_major_axis * math.sqrt(1 - eccentricity * eccentricity) * math.sin(math.radians(eccentric_anomaly))
        z = 0
        return [x, y, z]


    # Position (x y z) relative to reference plane (ecliptic), heliocentric. Unit: AU
    #
    # Position on orbit plane.
    # Center on barycenter.
    # x in reference direction (vernal equinox).
    # z is distance above reference plane.
    #
    # Transformation:
    # 1: Rotate plane around z with [argument_of_periapsis], point x toward Ascending node
    # 2: Tilt plane around x with [inclination], to the reference plane (ecliptic)
    # 3: Rotate plane around z with [longitude of ascending node], point x toward reference direction (vernal equinox)
    def get_position_ecliptic(self, julian_time):
        pos_orbit_plane_periapsis = self.get_position_own_orbit(julian_time)

        p_rad = math.radians(self.get_argument_of_periapsis(julian_time))
        i_rad = math.radians(self.inclination.get_at_J2000(julian_time))
        l_rad = math.radians(self.longitude_of_ascending_node.get_at_J2000(julian_time))

        # 1: Rotate plane around z with [argument_of_periapsis], point x toward Ascending node
        pos_orbit_plane_asc_node = transform.rotate_3d_z(pos_orbit_plane_periapsis, p_rad)

        # 2: Tilt plane around x with [inclination], to the reference plane (ecliptic)
        pos_reference_plane_asc_node = transform.rotate_3d_x(pos_orbit_plane_asc_node, i_rad)

        # 3: Rotate plane around z with [longitude of ascending node], point x toward reference direction (vernal equinox)
        pos_ecliptic = transform.rotate_3d_z(pos_reference_plane_asc_node, l_rad)

        return pos_ecliptic

    def get_orbit_time_in_days(self):
        mean_longitude_delta_pr_century = self.mean_longitude.delta_pr_century
        mean_longitude_delta_pr_day = mean_longitude_delta_pr_century / (100 * 365.25)
        orbits_pr_day = mean_longitude_delta_pr_day / 360
        days_pr_orbit = 1 / orbits_pr_day
        return days_pr_orbit

    def get_days_since_periapsis(self, julian_time):
        mean_anomaly = self.get_mean_anomaly(julian_time)
        days_pr_orbit = self.get_orbit_time_in_days()
        days_since_periapsis = ((mean_anomaly % 360) / 360) * days_pr_orbit
        return days_since_periapsis


    # Position (x y z) relative to reference plane, geocentric. Unit: AU
    #
    # Position on equatorial plane.
    # Center on earth.
    # x in reference direction (vernal equinox).
    # z is distance above reference plane.
    #
    # Transformation:
    # 1: Subtract the position of the Earth (heliocentric to geocentric)
    def get_position_ecliptic_geocentric(self, julian_time, earth_position):
        pos_ecliptic = self.get_position_ecliptic(julian_time)

        # 1: Subtract the position of the Earth (heliocentric to geocentric)
        pos_ecliptic_geocentric = transform.sub_3d(pos_ecliptic, earth_position)

        return pos_ecliptic_geocentric


    # Position (x y z) relative to equatorial plane. Unit: AU
    #
    # Position on equatorial plane.
    # Center on earth.
    # x in reference direction (vernal equinox).
    # z in direcion of celestial north pole (Polaris)
    #
    # Transformation:
    # 1: Tilt plane around x with [obliquity] (axial tilt), to the equatorial plane
    def get_position_equatorial(self, julian_time, earth_position):
        pos_ecliptic_geocentric = self.get_position_ecliptic_geocentric(julian_time, earth_position)

        o_rad = math.radians(OBLIQUITY)

        # 1: Tilt plane around x with [obliquity] (axial tilt), to the equatorial plane
        pos_equatorial = transform.rotate_3d_x(pos_ecliptic_geocentric, o_rad)

        return pos_equatorial


    def get_celestial_position(self, julian_time, earth_position):
        position_equatorial = self.get_position_equatorial(julian_time, earth_position)
        celestial_pos = transform.cartesian_to_latitude_longitude_dist(position_equatorial)
        return celestial_pos


    def debug(self):
        print(self.mean_longitude.get_at_delta_time(0))
        print(self.mean_longitude.get_at_delta_time(1/100/365))
        print(self.mean_longitude.get_at_delta_time(2/100/365))


'''
    @property
    def julian_time(self):
        return self.__julian_time

    @julian_time.setter
    def julian_time(self, value):
        self.__
        return self.__julian_time
'''






class SolarSystem():

    def __init__(self, julian_time):

        # Data copied from table1
        '''
               a              e               I                L            long.peri.      long.node.
           AU, AU/Cy     rad, rad/Cy     deg, deg/Cy      deg, deg/Cy      deg, deg/Cy     deg, deg/Cy
-----------------------------------------------------------------------------------------------------------
        '''
        table1_raw = '''
Mercury   0.38709927      0.20563593      7.00497902      252.25032350     77.45779628     48.33076593
          0.00000037      0.00001906     -0.00594749   149472.67411175      0.16047689     -0.12534081
Venus     0.72333566      0.00677672      3.39467605      181.97909950    131.60246718     76.67984255
          0.00000390     -0.00004107     -0.00078890    58517.81538729      0.00268329     -0.27769418
EM Bary   1.00000261      0.01671123     -0.00001531      100.46457166    102.93768193      0.0
          0.00000562     -0.00004392     -0.01294668    35999.37244981      0.32327364      0.0
Mars      1.52371034      0.09339410      1.84969142       -4.55343205    -23.94362959     49.55953891
          0.00001847      0.00007882     -0.00813131    19140.30268499      0.44441088     -0.29257343
Jupiter   5.20288700      0.04838624      1.30439695       34.39644051     14.72847983    100.47390909
         -0.00011607     -0.00013253     -0.00183714     3034.74612775      0.21252668      0.20469106
Saturn    9.53667594      0.05386179      2.48599187       49.95424423     92.59887831    113.66242448
         -0.00125060     -0.00050991      0.00193609     1222.49362201     -0.41897216     -0.28867794
Uranus   19.18916464      0.04725744      0.77263783      313.23810451    170.95427630     74.01692503
         -0.00196176     -0.00004397     -0.00242939      428.48202785      0.40805281      0.04240589
Neptune  30.06992276      0.00859048      1.77004347      -55.12002969     44.96476227    131.78422574
          0.00026291      0.00005105      0.00035372      218.45945325     -0.32241464     -0.00508664
Pluto    39.48211675      0.24882730     17.14001206      238.92903833    224.06891629    110.30393684
         -0.00031596      0.00005170      0.00004818      145.20780515     -0.04062942     -0.01183482
        '''

        '''

Pluto    39.48211675      0.24882730     17.14001206      238.92903833    224.06891629    110.30393684
         -0.00031596      0.00005170      0.00004818      145.20780515     -0.04062942     -0.01183482
        '''

        table1_list = table1_raw.splitlines()

        test = KeplerianObject(table1_list[1] + "\n" + table1_list[2])

        self.planet = []

        for i in range(0, 9):
            self.planet.append(KeplerianObject(table1_list[i * 2 + 1] + "\n" + table1_list[i * 2 + 2]))

    def get_planet_by_index(self, index):
        return self.planet[index]




# Test


def debugPlanets():
    venus = "Venus     0.72333566      0.00677672      3.39467605      181.97909950    131.60246718     76.67984255\n"
    venus += "          0.00000390     -0.00004107     -0.00078890    58517.81538729      0.00268329     -0.27769418"

    venusObj = KeplerianObject(venus)
    venusObj.debug()

    day = 2459170

    solar_system = SolarSystem(day)

    earth = solar_system.get_planet_by_index(2)
    earth_position = earth.get_position_ecliptic(day)

    # solar_system.planet[0].debug()

    for i in range(0, 9):
        planet = solar_system.get_planet_by_index(i)
        print("Planet: " + planet.name)

        print("degrees pr day: ", (planet.mean_longitude.delta_pr_century / 36525))
        print("days pr degree: ", (36525 / planet.mean_longitude.delta_pr_century))

        print("delta50y,m0f ", planet.mean_longitude.delta_pr_century / 2)
        print("delta50y,m1f ", 50 * 365.25 / (36525 / planet.mean_longitude.delta_pr_century))
        daysPrDegreeInt = round(36525 * 32 / planet.mean_longitude.delta_pr_century) / 32
        if (daysPrDegreeInt == 0):
            daysPrDegreeInt = 1
        print("delta50y,m1i ", 50 * 365.25 / daysPrDegreeInt)


        '''

        print("semi_major_axis")
        print(planet.semi_major_axis.get_at_J2000(day))
        print("eccentricity")
        print(planet.eccentricity.get_at_J2000(day))
        print("inclination")
        print(planet.inclination.get_at_J2000(day))
        print("mean_longitude")
        print(planet.mean_longitude.get_at_J2000(day))
        print("longitude_of_periapsis")
        print(planet.longitude_of_periapsis.get_at_J2000(day))
        print("longitude_of_ascending_node")
        print(planet.longitude_of_ascending_node.get_at_J2000(day))

        print("get_argument_of_periapsis")
        print(planet.get_argument_of_periapsis(day))
        print("get_mean_anomaly")
        print(planet.get_mean_anomaly(day))
        print("get_eccentric_anomaly")
        print(planet.get_eccentric_anomaly(day))
        print("get_position_own_orbit")
        print(planet.get_position_own_orbit(day))

        print("get_position_ecliptic")
        print(planet.get_position_ecliptic(day))

        print("get_position_ecliptic_geocentric")
        print(planet.get_position_ecliptic_geocentric(day, earth_position))
        print("get_position_equatorial")
        print(planet.get_position_equatorial(day, earth_position))
        print("get_celestial_position")
        print(planet.get_celestial_position(day, earth_position))
        '''
        print("")

# Debug
def computeJackTable():
    day = 2459170
    solar_system = SolarSystem(day)
    for i in range(0, 9):
        planet = solar_system.get_planet_by_index(i)
        print("\"" + planet.name + "\", ")
        print(planet.semi_major_axis.init_state)
        print(round(planet.semi_major_axis.init_state * 512))
        print(round(planet.semi_major_axis.init_state * math.sqrt(1 - planet.eccentricity.init_state * planet.eccentricity.init_state) * 512))
        print(math.sqrt(1 - planet.eccentricity.init_state * planet.eccentricity.init_state) * 128)


    for i in range(0,91):
        print("let sinTable[" + str(i) + "] = " + str(int(round(math.sin(math.radians(i)) * 256))) + ";")

    for i in range(0,129):
        print("let atanTable[" + str(i) + "] = " + str(int(round(math.degrees(math.atan(i / 128))))) + ";")

    for i in range(0, 9):
        planet = solar_system.get_planet_by_index(i)

        # bindeg: Binary degrees (0-32767), where 32768 is a full orbit

        # compute "period", a number of days where delta bindeg is as close to 32767 as possible
        deltaDegPrDay = planet.mean_longitude.delta_pr_century / 36525
        deltaBinDegPrDay = deltaDegPrDay * 32768 / 360
        period = int(32768 / deltaBinDegPrDay)
        if (deltaBinDegPrDay * period > 32767):
            period = period - 1

        if (period > 32767):
            period = 32767


        print("let planets[" + str(i) + "] = Planet.new(\"" + planet.name + "\", "
              + str(round(planet.semi_major_axis.init_state * 512)) + ", "
              + str(round(planet.eccentricity.init_state * 32768)) + ", "
              + str(round(planet.inclination.init_state) % 360) + ", "
              + str(round(planet.longitude_of_periapsis.init_state) % 360) + ", "
              + str(round(planet.longitude_of_ascending_node.init_state) % 360) + ", "
              + str(round((planet.mean_longitude.init_state % 360) * 32768 / 360)) + ", "
              + str(round(period * deltaBinDegPrDay)) + ", "
              + str(period) + ", "
              + str(round(planet.semi_major_axis.init_state * math.sqrt(1 - planet.eccentricity.init_state * planet.eccentricity.init_state) * 512))
              + ");")





if __name__ == "__main__":
    # debugPlanets()
    computeJackTable()
