import math

from shapely.geometry import Polygon
from shapely.prepared import prep
import numpy as np

# from https://www.matecdev.com/posts/shapely-polygon-gridding.html


def partition(geom, delta, mode='grid'):
    prepared_geom = prep(geom)
    if mode == 'grid':
        grid = list(filter(prepared_geom.intersects,
                        _grid_bounds(geom, delta)))
    elif mode == 'center':
        grid = list(filter(prepared_geom.intersects,
                        _center_bounds(geom, delta)))
        
    return grid

def _center_bounds(geom, delta):
    delta = int(delta/2)
    minx, miny, maxx, maxy = geom.bounds
    x_delta = (maxx - minx)/delta
    y_delta = (maxy - miny)/delta
    x_delta = math.ceil(x_delta) if 0.9 < x_delta%1 < 1 else x_delta
    y_delta = math.ceil(y_delta) if 0.9 < y_delta%1 < 1 else y_delta
    nx = int(x_delta) + 1
    ny = int(y_delta) + 1
    gx, gy = np.linspace(minx, maxx, nx), np.linspace(miny, maxy, ny)
    grid = []
    poly = Polygon([[gx[1], gy[1]],
                    [gx[1], gy[-2]],
                    [gx[-2], gy[-2]],
                    [gx[-2], gy[1]]])
    grid.append(poly)
    return grid

def _grid_bounds(geom, delta):
    minx, miny, maxx, maxy = geom.bounds
    x_delta = (maxx - minx)/delta
    y_delta = (maxy - miny)/delta
    x_delta = math.ceil(x_delta) if 0.9 < x_delta%1 < 1 else x_delta
    y_delta = math.ceil(y_delta) if 0.9 < y_delta%1 < 1 else y_delta
    nx = int(x_delta) + 1
    ny = int(y_delta) + 1
    gx, gy = np.linspace(minx, maxx, nx), np.linspace(miny, maxy, ny)
    grid = []
    for i in range(len(gx)-1):
        for j in range(len(gy)-1):
            polyij = Polygon([[gx[i],gy[j]],
                              [gx[i],gy[j+1]],
                              [gx[i+1],gy[j+1]],
                              [gx[i+1],gy[j]]])
            grid.append(polyij)
    return grid
