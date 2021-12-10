#!/usr/bin/python3
# solarsystem.py
# -*- coding: utf-8 -*-
#


# ----------------------------------------------------------------------------

"""
File containing the definition of all planets in the solar system

Dataset: "Keplerian elements, table 1 (1800 AD to 2050 AD)"
https://ssd.jpl.nasa.gov/?planet_pos

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

import cairocffi as cairo

import cairo_matrix_3d as matrix3d



class SolarSystemLog(BaseComponent):
    """
    Render the star wheel for the planisphere.
    """

    def default_filename(self):
        """
        Return the default filename to use when saving this component.
        """
        return "star_wheel2"

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

        '''
        # Draw lines of constant declination at 15 degree intervals.
        for dec in arange(-80, 85, 1):
            # Convert declination into radius from the centre of the planisphere
            r = radius(dec=dec, latitude=latitude)
            if r > r_2:
                continue
            context.begin_path()
            context.circle(centre_x=0, centre_y=0, radius=r)
            context.stroke(color=theme['grid'])
        '''

        #day = 2489156
        day = 2459170

        # Scale of the system:
        #max_render_radius = 45  # max radius to render, in AU
        max_render_radius = 6  # max radius to render, in AU
        length_of_au = r_2 / max_render_radius
        tick_length = r_2 / 10  # Length of 10 AU ticks
        min_tick_distance = tick_length / 2  # Minimum distance between minor ticks, skip category if smaller

        solar_system = kepler_object.SolarSystem(day)

        # Orientation of drawing:
        # Drawing similar to the unit circle.
        # The reference direction (x, vernal equinox) is to the right.
        # The other direction of the ecliptic plane (y), is up.

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

        # Draw ticks along axis
        # Strategy: Draw 10 - 5 - 1 - 0.5 - 0.1 AU lines.
        # Skip drawing lines if distance is closer than

        tick_values = [10, 5, 1, 0.5, 0.1]

        tick_length_temp = tick_length / 0.4

        value = 1000
        types_rendered = 0

        for i in range(0, 5):
            prev_value = value
            value = tick_values[i]

            if value > max_render_radius:
                continue

            types_rendered += 1
            if types_rendered > 3:
                break

            tick_distance = value * length_of_au
            #if tick_distance < min_tick_distance:
            #    break

            tick_length_temp *= 0.4
            half_tick_length = tick_length_temp / 2

            pos = 0
            k = 0
            RENDER_LIMIT = 50
            #delta_pos = RENDER_LIMIT / value
            while pos < RENDER_LIMIT:
                pos += value
                k += 1

                if pos % prev_value == 0:
                    continue

                '''
                for pos in range(0, int(50 / value)):

                if pos % prev_value == 0:
                    continue

                '''

                # Tick marks for each AU
                # Positive X
                context.begin_path()
                context.move_to(x=pos * length_of_au, y=-half_tick_length)
                context.line_to(x=pos * length_of_au, y=half_tick_length)
                context.stroke(line_width=1, dotted=False)

                # Negative X
                context.begin_path()
                context.move_to(x=-pos * length_of_au, y=-half_tick_length)
                context.line_to(x=-pos * length_of_au, y=half_tick_length)
                context.stroke(line_width=1, dotted=False)

                # Positive Y
                context.begin_path()
                context.move_to(y=pos * length_of_au, x=-half_tick_length)
                context.line_to(y=pos * length_of_au, x=half_tick_length)
                context.stroke(line_width=1, dotted=False)

                # Negative Y
                context.begin_path()
                context.move_to(y=-pos * length_of_au, x=-half_tick_length)
                context.line_to(y=-pos * length_of_au, x=half_tick_length)
                context.stroke(line_width=1, dotted=False)

        # Draw orbits
        for i in range(0, 9):
            planet = solar_system.get_planet_by_index(i)
            print("Planet: " + planet.name)

            print("longitude_of_periapsis")
            print(planet.longitude_of_periapsis.get_at_J2000(day))
            print("get_mean_anomaly")
            print(planet.get_mean_anomaly(day))
            print("get_eccentric_anomaly")
            print(planet.get_eccentric_anomaly(day))
            print("get_position_own_orbit")
            print(planet.get_position_own_orbit(day))

            print("get_position_ecliptic")
            print(planet.get_position_ecliptic(day))

            # Simplified: Assume orbit in ecliptic plane
            '''
            planet_pos_ecliptic = planet.get_position_ecliptic(day)
            x = planet_pos_ecliptic[0]
            y = 0 - planet_pos_ecliptic[1]  # Invert so that positive direction is upwards

            semi_major_axis = planet.semi_major_axis.get_at_J2000(day)
            inclination = planet.inclination.get_at_J2000(day)
            semi_major_axis_xy = semi_major_axis * math.cos(math.radians(inclination))

            #context.begin_sub_path()
            context.begin_path()
            context.eccentric_ellipse(centre_x=0, centre_y=-0, semi_major_axis=semi_major_axis_xy * length_of_au,
                                      eccentricity=planet.eccentricity.get_at_J2000(day), radians=-math.radians(planet.longitude_of_periapsis.get_at_J2000(day)))
            '''
            semi_major_axis = planet.semi_major_axis.get_at_J2000(day)
            eccentricity = planet.eccentricity.get_at_J2000(day)
            argument_of_periapsis = planet.get_argument_of_periapsis(day)
            inclination = planet.inclination.get_at_J2000(day)
            longitude_of_ascending_node = planet.longitude_of_ascending_node.get_at_J2000(day)

            mat3d = matrix3d.Matrix3d()
            mat3d.scale(1, math.sqrt(1 - eccentricity * eccentricity), 1)
            print("raw")
            print(mat3d.matrix)
            mat3d.translate(-eccentricity * length_of_au * semi_major_axis, 0, 0)
            print("trans")
            print(mat3d.matrix)
            mat3d.rotate_z(math.radians(argument_of_periapsis))
            print("rotz")
            print(mat3d.matrix)
            mat3d.rotate_x(math.radians(inclination))
            print("rotx")
            print(mat3d.matrix)
            mat3d.rotate_z(math.radians(longitude_of_ascending_node))
            print("rotz")
            print(mat3d.matrix)

            #mat3d.rotate_z(0.5)

            mat3d.flip_y()  # Flip Y because how this is rendered

            mat2d = mat3d.get_cairo_matrix()
            #mat2d = mat3d.get_cairo_matrix_flipped_y()
            print(mat2d)
            mat2d.scale(length_of_au * semi_major_axis, length_of_au * semi_major_axis)
            print(mat2d)

            context.begin_path()
            context.ellipse_with_matrix(mat2d)
            context.stroke(color=theme['edge'])


        # Draw planets
        for i in range(0, 9):
            planet = solar_system.get_planet_by_index(i)
            print("Planet: " + planet.name)

            planet_pos_ecliptic = planet.get_position_ecliptic(day)
            x = planet_pos_ecliptic[0]
            y = 0 - planet_pos_ecliptic[1]  # Invert so that positive direction is upwards

            #context.begin_sub_path()
            context.begin_path()
            context.circle(centre_x=x * length_of_au, centre_y=y * length_of_au, radius=central_hole_size)
            context.stroke(color=theme['edge'])

        # Draw periapsis
        for i in range(0, 9):
            planet = solar_system.get_planet_by_index(i)
            print("Planet: " + planet.name)

            print("Orbit time: ")
            print(planet.get_orbit_time_in_days())

            days_since_periapsis = planet.get_days_since_periapsis(day)
            print("Days since periapsis:")
            print(days_since_periapsis)

            periapsis_day = day - days_since_periapsis
            planet_pos_ecliptic_periapsis = planet.get_position_ecliptic(periapsis_day)
            x = planet_pos_ecliptic_periapsis[0]
            y = 0 - planet_pos_ecliptic_periapsis[1]  # Invert so that positive direction is upwards

            polar_periapsis = transform.cartesian_to_polar_2d([x, y])
            pp_a = polar_periapsis[0]  # Angle
            pp_r = polar_periapsis[1]  # Radius

            pp_r_start = pp_r * length_of_au - tick_length / 4
            pp_r_end   = pp_r * length_of_au + tick_length / 4

            context.begin_path()
            context.move_to(x=pp_r_start * math.cos(pp_a), y=pp_r_start * math.sin(pp_a))
            context.line_to(x=pp_r_end   * math.cos(pp_a), y=pp_r_end   * math.sin(pp_a))
            context.stroke(line_width=1, dotted=False)

            # Draw ticks along axis
            # Strategy: Draw 365.25 - (365.25/12) - 10 - 5 - 1 days.
            # Skip drawing lines if distance is closer than
            '''
            tick_values = [365.15, (365.25/12), 10, 5, 1]

            tick_length_temp = tick_length / 0.4

            value = 1000
            types_rendered = 0

            for i in range(0, 5):
                prev_value = value
                value = tick_values[i]

                if value > max_render_radius:
                    continue

                types_rendered += 1
                if types_rendered > 3:
                    break

                tick_distance = value * length_of_au
                # if tick_distance < min_tick_distance:
                #    break

                tick_length_temp *= 0.4
                half_tick_length = tick_length_temp / 2

                pos = 0
                k = 0
                RENDER_LIMIT = 50
                # delta_pos = RENDER_LIMIT / value
                while pos < RENDER_LIMIT:
                    pos += value
                    k += 1

                    if pos % prev_value == 0:
                        continue
            '''

        '''
        #context.begin_sub_path()
        context.begin_path()
        #context.eccentric_ellipse(centre_x=r_2/2, centre_y=-r_2/2, semi_major_axis=r_2, eccentricity=0.8, radians=0.9)


        mat3d = matrix3d.Matrix3d()
        mat3d.scale(1, 0.5, 1)
        print("raw")
        print(mat3d.matrix)
        mat3d.translate(-1 * length_of_au, 0, 0)
        print("trans")
        print(mat3d.matrix)
        mat3d.rotate_x(pi/4)
        print("rotx")
        print(mat3d.matrix)
        mat3d.rotate_z(pi/4)
        print("rotz")
        print(mat3d.matrix)
        mat2d = mat3d.get_cairo_matrix()
        print(mat2d)
        mat2d.scale(length_of_au, length_of_au)
        print(mat2d)
        context.ellipse_with_matrix(mat2d)

        context.stroke(line_width=5, color=theme['edge'])
        '''

        '''
        # Draw constellation stick figures
        for line in open("raw_data/constellation_stick_figures.dat", "rt"):
            line = line.strip()

            # Ignore blank lines and comment lines
            if (len(line) == 0) or (line[0] == '#'):
                continue

            # Split line into words.
            # These are the names of the constellations, and the start and end points for each stroke.
            [name, ra1, dec1, ra2, dec2] = line.split()

            # If we're making a southern hemisphere planisphere, we flip the sky upside down
            if is_southern:
                ra1 = -float(ra1)
                ra2 = -float(ra2)
                dec1 = -float(dec1)
                dec2 = -float(dec2)

            # Project RA and Dec into radius and azimuth in the planispheric projection
            r_point_1 = radius(dec=float(dec1), latitude=latitude)
            if r_point_1 > r_2:
                continue

            r_point_2 = radius(dec=float(dec2), latitude=latitude)
            if r_point_2 > r_2:
                continue

            p1 = (-r_point_1 * cos(float(ra1) * unit_deg), -r_point_1 * sin(float(ra1) * unit_deg))
            p2 = (-r_point_2 * cos(float(ra2) * unit_deg), -r_point_2 * sin(float(ra2) * unit_deg))

            # Impose a maximum length of 4 cm on constellation stick figures; they get quite distorted at the edge
            if hypot(p2[0] - p1[0], p2[1] - p1[1]) > 4 * unit_cm:
                continue

            # Stroke a line
            context.begin_path()
            context.move_to(x=p1[0], y=p1[1])
            context.line_to(x=p2[0], y=p2[1])
            context.stroke(color=theme['stick'], line_width=1, dotted=True)

        # Draw stars from Yale Bright Star Catalogue
        for star_descriptor in fetch_bright_star_list()['stars'].values():
            [ra, dec, mag] = star_descriptor[:3]

            # Discard stars fainter than mag 4
            if mag == "-" or float(mag) > 4.0:
                continue

            ra = float(ra)
            dec = float(dec)

            # If we're making a southern hemisphere planisphere, we flip the sky upside down
            if is_southern:
                ra *= -1
                dec *= -1

            r = radius(dec=dec, latitude=latitude)
            if r > r_2:
                continue

            # Represent each star with a small circle
            context.begin_path()
            context.circle(centre_x=-r * cos(ra * unit_deg), centre_y=-r * sin(ra * unit_deg),
                           radius=0.18 * unit_mm * (5 - mag))
            context.fill(color=theme['star'])

        # Write constellation names
        context.set_font_size(0.7)
        context.set_color(theme['constellation'])

        # Open a list of the coordinates where we place the names of the constellations
        for line in open("raw_data/constellation_names.dat"):
            line = line.strip()

            # Ignore blank lines and comment lines
            if (len(line) == 0) or (line[0] == '#'):
                continue

            # Split line into words
            [name, ra, dec] = line.split()[:3]

            # Translate constellation name into the requested language, if required
            if name in text[language]['constellation_translations']:
                name = text[language]['constellation_translations'][name]

            ra = float(ra) * 360. / 24
            dec = float(dec)

            # If we're making a southern hemisphere planisphere, we flip the sky upside down
            if is_southern:
                ra = -ra
                dec = -dec

            # Render name of constellation, with _s turned into spaces
            name2 = re.sub("_", " ", name)
            r = radius(dec=dec, latitude=latitude)
            if r > r_2:
                continue
            p = (-r * cos(ra * unit_deg), -r * sin(ra * unit_deg))
            a = atan2(p[0], p[1])
            context.text(text=name2, x=p[0], y=p[1], h_align=0, v_align=0, gap=0, rotation=unit_rev / 2 - a)

        # Calendar ring counts clockwise in northern hemisphere; anticlockwise in southern hemisphere
        s = -1 if not is_southern else 1

        def theta2014(d):
            """
            Convert Julian Day into a rotation angle of the sky about the north celestial pole at midnight,
            relative to spring equinox.

            :param d:
                Julian day
            :return:
                Rotation angle, radians
            """
            return (d - calendar.julian_day(year=2014, month=3, day=20, hour=16, minute=55, sec=0)) / 365.25 * unit_rev

        # Write month names around the date scale
        context.set_font_size(2.3)
        context.set_color(theme['date'])
        for mn, (mlen, name) in enumerate(text[language]['months']):
            theta = s * theta2014(calendar.julian_day(year=2014, month=mn + 1, day=mlen // 2, hour=12, minute=0, sec=0))

            # We supply circular_text with a negative radius here, as a fudge to orientate the text with bottom-inwards
            context.circular_text(text=name, centre_x=0, centre_y=0, radius=-(r_1 * 0.65 + r_2 * 0.35),
                                  azimuth=theta / unit_deg + 180,
                                  spacing=1, size=1)

        # Draw ticks for the days of the month
        for mn, (mlen, name) in enumerate(text[language]['months']):
            # Tick marks for each day
            for d in range(1, mlen + 1):
                theta = s * theta2014(calendar.julian_day(year=2014, month=mn + 1, day=d, hour=0, minute=0, sec=0))

                # Days of the month which are multiples of 5 get longer ticks
                R = r_3 if (d % 5) else r_4

                # The last day of each month is drawn as a dividing line between months
                if d == mlen:
                    R = r_5

                # Draw line
                context.begin_path()
                context.move_to(x=r_2 * cos(theta), y=-r_2 * sin(theta))
                context.line_to(x=R * cos(theta), y=-R * sin(theta))
                context.stroke(line_width=1, dotted=False)

            # Write numeric labels for the 10th, 20th and last day of each month
            for d in [10, 20, mlen]:
                theta = s * theta2014(calendar.julian_day(year=2014, month=mn + 1, day=d, hour=0, minute=0, sec=0))
                context.set_font_size(1.2)

                # First digit
                theta2 = theta + 0.15 * unit_deg
                context.text(text="%d" % (d / 10), x=r_6 * cos(theta2), y=-r_6 * sin(theta2),
                             h_align=1, v_align=0,
                             gap=0,
                             rotation=-theta + pi / 2)

                # Second digit
                theta2 = theta - 0.15 * unit_deg
                context.text(text="%d" % (d % 10), x=r_6 * cos(theta2), y=-r_6 * sin(theta2),
                             h_align=-1, v_align=0,
                             gap=0,
                             rotation=-theta + pi / 2)

        # Draw the dividing line between the date scale and the star chart
        context.begin_path()
        context.circle(centre_x=0, centre_y=0, radius=r_2)
        context.stroke(color=theme['date'], line_width=1, dotted=False)
'''




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
