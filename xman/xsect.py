# -*- coding: utf-8 -*-
"""xman

    Module for extracting a cross seciton
    from a DEM.

    Author:  Alan D Snow, 2017.
    License: BSD-3-Clause
"""
from gazar.grid import GDALGrid, utm_proj_from_latlon
import numpy as np
from shapely.geometry import LineString, Polygon
from osgeo import osr

def x_section_from_latlon(elevation_file,
                          x_section_lat0,
                          x_section_lon0,
                          x_section_lat1,
                          x_section_lon1,
                          as_polygon=False,
                          auto_clean=False):

    """
    This workflow extracts a cross section from a DEM
    based on the input latitude and longitude point pairs.

    Parameters:
    -----------
    elevation_file: str
        Path to the elevation DEM.
    x_section_lat0: float
        THe first coordinate latitude.
    x_section_lon0: float
        THe first coordinate longitude.
    x_section_lat1: float
        THe second coordinate latitude.
    x_section_lon1: float
        THe second coordinate longitude.
    as_polygon: bool, optional
        If True, will return cross section as a
        :obj:`shapely.geometry.Polygon`. Default is False.
    auto_clean: bool, optional
        If True, will attempt to clean any issues from the polygon.
        Default is False.

    Returns:
    --------
    list or :obj:`shapely.geometry.Polygon`
        Cross section information.
        The list will be xy coordinate pairs.


    Example::

        from shapely.geometry import Polygon
        from xman.xsect import x_section_from_latlon


        elevation_file = '/path/to/elevation.tif'
        lat1 = 34.105265417341442
        lon1 = 38.993958690587505
        lat2 = 34.107264451129197
        lon2 = 38.99355588515526)
        x_sect_list = x_section_from_latlon(elevation_file,
                                            lat1,
                                            lon1,
                                            lat2,
                                            lon2)

    """
    utm_proj = utm_proj_from_latlon(x_section_lat0, x_section_lon0,
                                    as_osr=True)
    sp_ref = osr.SpatialReference()
    sp_ref.ImportFromEPSG(4326)
    geo_to_utm_trans = osr.CoordinateTransformation(sp_ref, utm_proj)

    x_line_m = LineString((
        geo_to_utm_trans.TransformPoint(x_section_lon0, x_section_lat0)[:2],
        geo_to_utm_trans.TransformPoint(x_section_lon1, x_section_lat1)[:2]
    ))

    elevation_utm_ggrid = GDALGrid(elevation_file).to_projection(utm_proj)

    x_sect_list = []

    for x_step in np.linspace(0, x_line_m.length, num=20):
        x_point = x_line_m.interpolate(x_step)
        x_sect_list.append((
            x_step, elevation_utm_ggrid.get_val_coord(x_point.x, x_point.y)
        ))

    if as_polygon or auto_clean:
        x_sect_poly = Polygon(x_sect_list)
        if not x_sect_poly.is_valid and auto_clean:
            x_sect_poly = x_sect_poly.buffer(0)
            print("WARNING: Cross section cleaned up.")
            if hasattr(x_sect_poly, 'geoms'):
                if len(x_sect_poly.geoms) > 1:
                    largest_poly = x_sect_poly.geoms[0]
                    for geom_poly in x_sect_poly.geoms[1:]:
                        if geom_poly.area > largest_poly.area:
                            largest_poly = geom_poly
                    x_sect_poly = largest_poly

        if as_polygon:
            return x_sect_poly

        x_coords, y_coords = x_sect_poly.exterior.coords.xy
        return list(zip(x_coords, y_coords))

    return x_sect_list