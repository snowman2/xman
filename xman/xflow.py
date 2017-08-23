# -*- coding: utf-8 -*-
"""xman

    Module for solving for flow using Manning's
    equation and plotting solutions.

    Author:  Alan D Snow, 2017.
    License: BSD-3-Clause
"""
from descartes.patch import PolygonPatch
from shapely.geometry import LineString, Polygon
import matplotlib.pyplot as plt
import numpy as np


def get_xpolygon_bottom(cross_section_polygon, depth_from_bottom):
    """
    Returns a subset of the polygon based on the distance from the bottom.

    Parameters:
    -----------
    cross_section_polygon: :obj:`shapely.geometry.polygon.Polygon`
        Polygon of cross section.
    depth_from_bottom: float
        Distance from the bottom of the polygon to retain.

    Returns:
    --------
    :obj:`shapely.geometry.polygon.Polygon`
        Subset of the polygon.
    """

    # determine extraction polygon
    x_coords, y_coords = cross_section_polygon.exterior.coords.xy
    max_y = np.amax(y_coords)
    min_y = np.amin(y_coords)
    max_x = np.amax(x_coords)
    min_x = np.amin(x_coords)

    y_extract_bottom = min_y + depth_from_bottom
    y_extract_top = max_y + 1

    extraction_poly_xy_list = [(min_x - 1, y_extract_top),
                               (max_x + 1, y_extract_top),
                               (max_x + 1, y_extract_bottom),
                               (min_x - 1, y_extract_bottom)]

    extraction_poly = Polygon(extraction_poly_xy_list)

    # get bottom polygon to depth
    return cross_section_polygon.difference(extraction_poly)


def flow_from_xsection(cross_section_polygon,
                       water_depth,
                       channel_slope,
                       mannings_roughness):
    """
    Returns the flow in m^3/s based on the cross section polygon,
    the water depth, and the slope.

    Parameters:
    -----------
    cross_section_polygon: :obj:`shapely.geometry.polygon.Polygon`
        A polygon of the cross section with the top left coordinate
        at the beginning of the array.
    water_depth: float
        Depth from bottom of the cross section.
    channel_slope: float,
        Slope of the channel in m/m.
    mannings_roughness: float
        Value for Manning's roughness (N).

    Returns:
    --------
    float
        Flow in m^3/s.
    """
    bottom_poly = get_xpolygon_bottom(cross_section_polygon,
                                      water_depth)

    # get wetted perimeter line
    bot_x, bot_y = bottom_poly.exterior.coords.xy
    bottom_line = LineString(zip(bot_x[:-1], bot_y[:-1]))

    # return flow
    wetted_perimiter = bottom_line.length
    cross_section_area = bottom_poly.area
    return ((1 / mannings_roughness) * cross_section_area *
            (cross_section_area / wetted_perimiter) ** (2 / 3) *
            channel_slope ** 0.5)


def find_depth_velocity_area(cross_section_polygon,
                             target_flow,
                             max_delta=0.01,
                             channel_slope=0.001,
                             mannings_roughness=0.013):
    """
    Iterative solver that calculates the depth, velocity,
    and area of the cross section filled with water.

    Parameters:
    -----------
    cross_section_polygon: :obj:`shapely.geometry.polygon.Polygon`
        A polygon of the cross section with the top left coordinate
        at the beginning of the array.
    target_flow: float
        This is the desired flow to match in m^3/s.
    max_delta: float, optional
        This is the maximum difference between the target_flow and
        the calculated flow. Default ia 0.01.
    channel_slope: float, optional
        This is the slope of the channel (m/m). Default is 0.001.
    mannings_roughness: float, optional
        This is the manning's n roughness value of the channel.
        Default os 0.013.

    Returns:
    --------
    depth: float
        The depth of the water in the channel.
    area: float
        The area of the cross section filled with water.
    velocity: float
        The velocity (m/s) of water moving through the channel.
    bottom_poly: :obj:`shapely.geometry.polygon.Polygon`
        A polygon where the channel is filled with water.


    Example::

        from shapely.geometry import Polygon
        from xman.xflow import find_depth_velocity_area


        target_flow = 35.
        max_delta = 1.
        xy_list = [(0, 20), (1, 13), (2, 12),
                   (3, 11), (4, 10), (5, 4),
                   (6, 7), (7, 15), (8, 16), (9, 20)]
        cross_section_polygon = Polygon(xy_list)
        info = find_depth_velocity_area(cross_section_polygon,
                                        target_flow,
                                        max_delta)

    """
    calc_flow = -9999.

    y_coords = cross_section_polygon.exterior.coords.xy[1]
    previous_depth_max = np.amax(y_coords)
    previous_depth_min = np.amin(y_coords)
    depth_guess = (previous_depth_max + previous_depth_min) / 2.

    # check flow bounds
    max_flow = flow_from_xsection(cross_section_polygon,
                                  previous_depth_max,
                                  channel_slope,
                                  mannings_roughness)

    if target_flow > max_flow:
        raise ValueError("The target_flow {0} is greater than "
                         "the bankful discharge {1}. Overflow ..."
                         .format(target_flow, max_flow))

    # iterate to find cross section matching discharge
    # with mannings equation
    while np.abs(calc_flow - target_flow) > max_delta:
        calc_flow = flow_from_xsection(cross_section_polygon,
                                       depth_guess,
                                       channel_slope,
                                       mannings_roughness)
        if calc_flow > target_flow:
            previous_depth_max = depth_guess
            depth_guess = (depth_guess + previous_depth_min) / 2.
        elif calc_flow < target_flow:
            previous_depth_min = depth_guess
            depth_guess = (depth_guess + previous_depth_max) / 2.

    bottom_poly = get_xpolygon_bottom(cross_section_polygon, depth_guess)
    bot_y = bottom_poly.exterior.coords.xy[1]
    depth = np.amax(bot_y) - np.amin(bot_y)
    area = bottom_poly.area

    return depth, area, (target_flow / area), bottom_poly


def plot_filled_cross_section(cross_section_polygon, water_filled_polygon):
    """
    Plot the original cross section with the region filled with water.

    Parameters:
    -----------
    cross_section_polygon: :obj:`shapely.geometry.polygon.Polygon`
        Polygon of cross section.
    water_filled_polygon: :obj:`shapely.geometry.polygon.Polygon`
        Polygon of water filled region.
    """
    plt.figure()
    plot_ax = plt.axes()
    plot_ax.set_aspect('equal')

    # build the polygon from exterior points
    patch = PolygonPatch(cross_section_polygon,
                         facecolor=[1, 1, 1],
                         edgecolor=[0, 0, 0],
                         alpha=0.7,
                         zorder=2)
    plot_ax.add_patch(patch)
    patch = PolygonPatch(water_filled_polygon,
                         facecolor=[0, 0, 0.7],
                         edgecolor=[0, 0, 0],
                         alpha=0.7,
                         zorder=2)
    plot_ax.add_patch(patch)

    # use bounding box to set plot limits
    minx, miny, maxx, maxy = cross_section_polygon.bounds
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)
    plt.show()


def find_depth_velocity_area_plot(cross_section_polygon,
                                  target_flow,
                                  max_delta=0.01,
                                  channel_slope=0.001,
                                  mannings_roughness=0.013):
    """
    This workflow calculates the depth, velocity,
    and area of the cross section filled with water.
    Then, it plots the final cross section and prints
    the results to the console.

    Parameters:
    -----------
    cross_section_polygon: :obj:`shapely.geometry.polygon.Polygon`
        A polygon of the cross section with the top left coordinate
        at the beginning of the array.
    target_flow: float
        This is the desired flow to match in m^3/s.
    max_delta: float, optional
        This is the maximum difference between the target_flow and
        the calculated flow. Default ia 0.01.
    channel_slope: float, optional
        This is the slope of the channel (m/m). Default is 0.001.
    mannings_roughness: float, optional
        This is the manning's n roughness value of the channel.
        Default os 0.013.

    Example::

        from shapely.geometry import Polygon
        from xman.xflow import find_depth_velocity_area_plot


        target_flow = 35.
        max_delta = 1.
        xy_list = [(0, 20), (1, 13), (2, 12),
                   (3, 11), (4, 10), (5, 4),
                   (6, 7), (7, 15), (8, 16), (9, 20)]
        cross_section_polygon = Polygon(xy_list)
        find_depth_velocity_area_plot(cross_section_polygon,
                                      target_flow,
                                      max_delta)

    """
    depth, area, velocity, bottom_poly = find_depth_velocity_area(
        cross_section_polygon,
        target_flow,
        max_delta,
        channel_slope,
        mannings_roughness,
    )
    print("DEPTH    {:.2f}".format(depth))
    print("AREA     {:.2f}".format(area))
    print("VELOCITY {:.2f}".format(velocity))
    plot_filled_cross_section(cross_section_polygon, bottom_poly)
