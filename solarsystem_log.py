#!/usr/bin/python3
# solarsystem_log.py
# -*- coding: utf-8 -*-
#


# ----------------------------------------------------------------------------

"""
Render the solar system with logarithmic scales
"""

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

import kepler_object
import math
import transform
import datetime

AU_PR_LIGHT_YEAR = 63241
KM_PR_AU = 149597870.700

class SolarSystemLog(BaseComponent):
    """
    Render the star wheel for the planisphere.
    """


    def default_filename(self):
        """
        Return the default filename to use when saving this component.
        """
        return "star_wheel_log"

    def bounding_box(self, settings):
        """
        Return the bounding box of the canvas area used by this component.

        :param settings:
            A dictionary of settings required by the renderer.
        :return:
         Dictionary with the elements 'x_min', 'x_max', 'y_min' and 'y_max' set
        """
        return {
            'x_min': -r_1 - 4 * unit_mm,
            'x_max': r_1 + 4 * unit_mm,
            'y_min': -r_1 - 4 * unit_mm,
            'y_max': r_1 + 4 * unit_mm
        }

    def do_rendering(self, settings, context):
        """
        This method is required to actually render this item.

        :param settings:
            A dictionary of settings required by the renderer.
        :param context:
            A GraphicsContext object to use for drawing
        :return:
            None
        """

        is_southern = settings['latitude'] < 0
        language = settings['language']
        latitude = abs(settings['latitude'])
        theme = themes[settings['theme']]

        context.set_font_size(1.2)

        # Radius of outer edge of star chart
        r_2 = r_1 - r_gap

        # Radius of day-of-month ticks from centre of star chart
        r_3 = r_1 * 0.1 + r_2 * 0.9

        # Radius of every fifth day-of-month tick from centre of star chart
        r_4 = r_1 * 0.2 + r_2 * 0.8

        # Radius of lines between months on date scale
        r_5 = r_1

        # Radius for writing numeric labels for days of the month
        r_6 = r_1 * 0.4 + r_2 * 0.6

        # Shade background to month scale
        shading_inner_radius = r_1 * 0.55 + r_2 * 0.45
        context.begin_path()
        context.circle(centre_x=0, centre_y=0, radius=r_1)
        context.circle(centre_x=0, centre_y=0, radius=shading_inner_radius)
        context.fill(color=theme['shading'])

        # Draw the outer edge of planisphere
        context.begin_path()
        context.circle(centre_x=0, centre_y=0, radius=r_1)
        context.fill(color=theme['background'])

        # Draw the central hole in the middle of the planisphere
        context.begin_sub_path()
        context.circle(centre_x=0, centre_y=0, radius=central_hole_size)
        context.stroke(color=theme['edge'])

        # Combine these two paths to make a clipping path for drawing the star wheel
        context.clip()

        # Draw lines of constant declination at 15 degree intervals.
        '''
        for dec in arange(-80, 85, 15):
            # Convert declination into radius from the centre of the planisphere
            r = radius(dec=dec, latitude=latitude)
            if r > r_2:
                continue
            context.begin_path()
            context.circle(centre_x=0, centre_y=0, radius=r)
            context.stroke(color=theme['grid'])
        '''

        #day = 2489156
        #day = 2459170

        #day_datetime = datetime.datetime(2020, 5, 17)
        day_datetime = datetime.datetime.now()
        day_number_in_year = day_datetime.date().timetuple().tm_yday
        day = kepler_object.J2000_JULIAN_DAY + (day_datetime.date().timetuple().tm_year - 2000) * 365.25 + day_number_in_year

        # Scale of the system:
        #max_render_radius = 40  # max radius to render, in AU
        #max_render_radius = 300  # max radius to render, in AU
        #max_render_radius = 100000  # max radius to render, in AU
        max_render_radius = AU_PR_LIGHT_YEAR * 25000
        max_render_radius_log = math.log(max_render_radius, 10) + 1  # max radius to render, in log(AU). Add 1 to multiply with 10
        length_of_au_log = r_2 / max_render_radius_log
        tick_length = r_2 / 10  # Length of 10 AU ticks
        min_tick_distance = tick_length / 2  # Minimum distance between minor ticks, skip category if smaller

        solar_system = kepler_object.SolarSystem(day)




        tick_length = 10
        half_tick_length = tick_length / 2

        # Draw logarithmic circles. 0.1 to 40
        for log_line in arange(2, 11, 1):  # 2,3,4,5,6,7,8,9,10
            len_in_au_log = math.log(log_line, 10)

            while len_in_au_log < max_render_radius_log * 10:  # * 10 to draw all the way out

                r = len_in_au_log * length_of_au_log
                '''
                context.begin_path()
                context.circle(centre_x=0, centre_y=0, radius=r)
                context.stroke(line_width=1, color=theme['grid'])
                '''
                # Tick marks for each AU
                # Positive X
                context.begin_path()
                context.move_to(x=r, y=-half_tick_length)
                context.line_to(x=r, y=half_tick_length)
                context.stroke(line_width=1, dotted=False, color=theme['grid'])

                # Negative X
                context.begin_path()
                context.move_to(x=-r, y=-half_tick_length)
                context.line_to(x=-r, y=half_tick_length)
                context.stroke(line_width=1, dotted=False)

                # Positive Y
                context.begin_path()
                context.move_to(y=r, x=-half_tick_length)
                context.line_to(y=r, x=half_tick_length)
                context.stroke(line_width=1, dotted=False)

                # Negative Y
                context.begin_path()
                context.move_to(y=-r, x=-half_tick_length)
                context.line_to(y=-r, x=half_tick_length)
                context.stroke(line_width=1, dotted=False)

                len_in_au_log += 1




        # 1: Draw x/y axis for ecliptic coordinate system

        # Draw x axis
        context.begin_path()
        context.move_to(x=0, y=-r_1)
        context.line_to(x=0, y=r_1)
        context.stroke(line_width=1, dotted=False)

        # Draw y axis
        context.begin_path()
        context.move_to(x=-r_1, y=0)
        context.line_to(x=r_1, y=0)
        context.stroke(line_width=1, dotted=False)


        # Draw orbit for each planet, logarithmic scale

        # Draw orbits

        planet_pos_tick_length = 0.002
        planet_pos_half_tick_length = planet_pos_tick_length / 2

        current_year = 2022
        #current_month = 1

        for i in range(0, 9):
            planet = solar_system.get_planet_by_index(i)
            print("Planet: " + planet.name)


            semi_major_axis = planet.semi_major_axis.get_at_J2000(day)

            print("semi_major_axis")
            print(semi_major_axis)

            # Log. Add 1 to multiply with 10
            log_semi_major_axis = math.log(semi_major_axis, 10) + 1
            print("log_semi_major_axis")
            print(log_semi_major_axis)

            r = log_semi_major_axis * length_of_au_log
            print("radius")
            print(r)

            # Draw the orbit
            context.begin_path()
            context.circle(centre_x=0, centre_y=0, radius=r)
            context.stroke(color=theme['edge'])

            # Draw location at current day

            #stroke_planet_orbit(planet, day)
            stroke_planet_orbit_get(planet, day, r, planet_pos_tick_length, 4, context, theme['edge'])

            # Draw strokes, 1/4 turn. Resolution: Day, 5 days, 10 days, month, half year, 5 year, 10 year, max 3 types

            steps = [0, 0, 1, 1, 3, 4, 5, 5, 5]

            stroke_planet_orbit_multiple(planet, day_datetime.date(), r, planet_pos_tick_length / 2, steps[i], context, theme['edge'])

        # Data for 6.12.21, https://ssd.jpl.nasa.gov/horizons/app.html
        draw_object_km("Voyager 1", -4.498040770062099E+09, -1.839704434327980E+10, 1.331230823914344E+10, length_of_au_log, context, theme['edge'])
        draw_object_km("Voyager 2", 5.325151335924605E+09, -1.435487832307063E+10, -1.169558464606284E+10, length_of_au_log, context, theme['edge'])

        draw_circle("Termination shock", 80, length_of_au_log, context, theme['edge'])

        draw_circle("Light year", 63241, length_of_au_log, context, theme['edge'])

        #draw_object_ly("Voyager 2", 5.325151335924605E+09, -1.435487832307063E+10, -1.169558464606284E+10, length_of_au_log, context, theme['edge'])

        # Equatorial_dir: 3d polar orientation. Right Ascension in hours (0-23), Declination in degrees (-90 to +90)
        # Equatorial_polar: 3d polar position, with radius in any unit
        # Equatorial_cartesian_au: 3d polar position, with radius in AU.

        sagittarius_equatorial = [17.76, -29.01]
        sagittarius_ecliptic = equatorial_polar_to_ecliptic_cartesian(sagittarius_equatorial)

        sagittarius_geo_polar_h_deg_ly = [17.76, -29.01, 26.4]
        sagittarius_j2000_ly = [17.76, -29.01, 26.4]


# Input:  Right Ascension in hours (0-23).
#         Declination in degrees (-90 to +90)
#         Radius in any unit
# Output: x,y,z in the same unit as radius
def equatorial_polar_to_ecliptic_cartesian(equatorial_polar):
    equatorial_right_ascension_rad = equatorial_polar[0] / 24 * (math.pi * 2)  # Hours to radians
    equatorial_declination_rad = equatorial_polar[1] / 360 * (math.pi * 2)  # Degrees to radians
    equatorial_radius = equatorial_polar[2]



def draw_circle(name, radius_au, length_of_au_log, context, color):
    # Log. Add 1 to multiply with 10
    radius_in_au_log = math.log(radius_au, 10) + 1

    r = radius_in_au_log * length_of_au_log

    # Draw the orbit
    context.begin_path()
    context.circle(centre_x=0, centre_y=0, radius=r)
    context.stroke(color=color)



def km_to_au(distance):
    return distance / KM_PR_AU

def ly_to_au(distance):
    return distance * AU_PR_LIGHT_YEAR

def au_to_km(distance):
    return distance * KM_PR_AU


def draw_object_ly(name, xpos_ly, ypos_ly, zpos_ly, length_of_au_log, context, color):
    xpos_au = ly_to_au(xpos_ly)
    ypos_au = ly_to_au(ypos_ly)
    zpos_au = ly_to_au(zpos_ly)
    draw_object_km(name, xpos_au, ypos_au, zpos_au, length_of_au_log, context, color)

def draw_object_km(name, xpos_km, ypos_km, zpos_km, length_of_au_log, context, color):
    xpos_au = km_to_au(xpos_km)
    ypos_au = km_to_au(ypos_km)
    zpos_au = km_to_au(zpos_km)
    draw_object_au(name, xpos_au, ypos_au, zpos_au, length_of_au_log, context, color)


# Input is in J2000.0 ecliptic frame
def draw_object_au(name, xpos, ypos, zpos, length_of_au_log, context, color):
    pos_au = [xpos, ypos, zpos]
    polar_pos = transform.cartesian_to_polar_2d(pos_au)  # ignore Z
    radius_in_au = polar_pos[1]
    radius_in_au_log = math.log(radius_in_au, 10) + 1
    polar_pos[1] = radius_in_au_log * length_of_au_log

    pos_log = transform.polar_to_cartesian_2d(polar_pos)
    pos_log[1] *= -1  # Invert so that positive Y direction is upwards

    context.begin_path()
    context.circle(centre_x=pos_log[0], centre_y=pos_log[1], radius=0.002)
    context.stroke(color=color)

    return

def stroke_planet_orbit_multiple(planet, day, radius, stroke_len, steps, context, color):
    # Day, 5 day, 10 day, month, 6 month, year, 5 year, 10 year
    #days_pr_step = [1, 5, 10, 30, 180, 360, 360 * 5, 360 * 10]
    #years_pr_step = [(1 / 360), (1 / 72), (1 / 36), (1 / 12), (1 / 2), 1, 5, 10]



    # temp_day = day.date()


    # Step back in time one quarter orbit time

    orbit_time = datetime.timedelta(days=planet.get_orbit_time_in_days())
    one_eighth_orbit_time = orbit_time / 4

    # first_day = day - one_eighth_orbit_time
    last_day = day + one_eighth_orbit_time

    # temp_day = first_day
    #temp_day = datetime.datetime.today()

    # step back in time until highest day significance

    temp_day = day

    while 1:
        if get_day_signification(temp_day) >= steps + 2:
            temp_day = temp_day - datetime.timedelta(days=1)
            break
        temp_day = temp_day - datetime.timedelta(days=1)


    planet_pos_ecliptic_start = planet.get_position_ecliptic(datetime_to_julian_time(temp_day))
    polar_pos_start = transform.cartesian_to_polar_2d(planet_pos_ecliptic_start)
    angle_start = polar_pos_start[0]


    # Step forward until last day
    while 1:
    #while temp_day < last_day:
        temp_day = temp_day + datetime.timedelta(days=1)

        current_day_significance = get_day_signification(temp_day)
        '''
        current_day_significance = 0 # day
        if temp_day.day % 5 == 1 and temp_day.day != 31:
            current_day_significance = 1  # 5 day (day 1, 6, 11, 16, 21, 26 inside a month)
        if temp_day.day % 10 == 1 and temp_day.day != 31:
            current_day_significance = 2  # 10 day (day 1, 11, 21 inside a month)
        if temp_day.day == 1:
            current_day_significance = 3  # month
            if temp_day.month % 6 == 1:
                current_day_significance = 4  # half year (month 1, 7)
            if temp_day.month == 1:
                current_day_significance = 5  # year
                if temp_day.year % 5 == 0:
                    current_day_significance = 6  # 5 year (2020, 2025, 2030, ...)
                if temp_day.year % 10 == 0:
                    current_day_significance = 7  # 10 year (2020, 2030, 2040, ...)
        '''
        if current_day_significance < steps:
            continue

        if current_day_significance > steps + 2:
            current_day_significance = steps + 2

        julian_day = datetime_to_julian_time(temp_day)

        '''
        day_number_in_year = temp_day.timetuple().tm_yday
        julian_day = kepler_object.J2000_JULIAN_DAY + (
                    temp_day.timetuple().tm_year - 2000) * 365.25 + day_number_in_year
        '''
        factor = [1, 2, 4]

        # troke_len_now = stroke_len * (current_day_significance - steps + 1)
        stroke_len_now = stroke_len * factor[current_day_significance - steps]

        stroke_planet_orbit_get(planet, julian_day, radius, stroke_len_now, 1, context, color)

        if current_day_significance >= steps + 2:
            planet_pos_ecliptic = planet.get_position_ecliptic(julian_day)
            polar_pos = transform.cartesian_to_polar_2d(planet_pos_ecliptic)
            angle = polar_pos[0]

            delta_angle = angle - angle_start
            if delta_angle < 0:
                delta_angle += math.pi * 2

            if delta_angle > math.pi: # Half orbit
                break

    return


def get_day_signification(day):
    current_day_significance = 0  # day
    if day.day % 5 == 1 and day.day != 31:
        current_day_significance = 1  # 5 day (day 1, 6, 11, 16, 21, 26 inside a month)
    if day.day % 10 == 1 and day.day != 31:
        current_day_significance = 2  # 10 day (day 1, 11, 21 inside a month)
    if day.day == 1:
        current_day_significance = 3  # month
        if day.month % 6 == 1:
            current_day_significance = 4  # half year (month 1, 7)
        if day.month == 1:
            current_day_significance = 5  # year
            if day.year % 5 == 0:
                current_day_significance = 6  # 5 year (2020, 2025, 2030, ...)
            if day.year % 10 == 0:
                current_day_significance = 7  # 10 year (2020, 2030, 2040, ...)
    return current_day_significance

def stroke_planet_orbit_get(planet, day, radius, width, stroke_width, context, color):
    planet_pos_ecliptic = planet.get_position_ecliptic(day)
    # x = planet_pos_ecliptic[0]
    # y = 0 - planet_pos_ecliptic[1]  # Invert so that positive direction is upwards

    # Convert to polar coordinates
    polar_pos = transform.cartesian_to_polar_2d(planet_pos_ecliptic)  # ignore Z

    # Use the same radius as the orbit
    polar_pos[1] = radius

    # Make radius logarithmic
    # polar_pos[1] = (math.log(polar_pos[1], 10) + 1) * length_of_au_log

    polar_pos_start = polar_pos.copy()
    # polar_pos_start[1] -= planet_pos_half_tick_length

    polar_pos_end = polar_pos.copy()
    # polar_pos_end[1] += planet_pos_half_tick_length

    polar_pos_start[1] = polar_pos_start[1] - width / 2
    polar_pos_end[1] = polar_pos_end[1] + width / 2

    pos_start = transform.polar_to_cartesian_2d(polar_pos_start)
    pos_end = transform.polar_to_cartesian_2d(polar_pos_end)

    pos_start[1] *= -1  # Invert so that positive Y direction is upwards
    pos_end[1] *= -1  # Invert so that positive Y direction is upwards

    context.begin_path()
    context.move_to(x=pos_start[0], y=pos_start[1])
    context.line_to(x=pos_end[0], y=pos_end[1])

    context.stroke(line_width=stroke_width, dotted=False, color=color)
    return

def datetime_to_julian_time(day_datetime):
    day_number_in_year = day_datetime.timetuple().tm_yday
    day_julian = kepler_object.J2000_JULIAN_DAY + (
                day_datetime.timetuple().tm_year - 2000) * 365.25 + day_number_in_year
    return day_julian


# Do it right away if we're run as a script
if __name__ == "__main__":
    # Fetch command line arguments passed to us
    arguments = fetch_command_line_arguments(default_filename=SolarSystemLog().default_filename())

    # Render the star wheel for the planisphere
    SolarSystemLog(settings={
        'latitude': arguments['latitude'],
        'language': 'en',
        'theme': arguments['theme'],
    }).render_to_file(
        filename=arguments['filename'],
        img_format=arguments['img_format'],

    )
