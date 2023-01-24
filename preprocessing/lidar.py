import os
import pdal
import sys
from osgeo import gdal, osr
from pathlib import Path
from tqdm import tqdm
import geopandas as gpd
import json

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = conda_env_path/'share'/'proj'

def _get_polygon_str(x_cord, y_cord):
    polygon_str = 'POLYGON(('
    for x, y in zip(list(x_cord), list(y_cord)):
        polygon_str += f'{x} {y}, '
    polygon_str = polygon_str[:-2]
    polygon_str += '))'
    return polygon_str

def clip_laz_by_plots(laz_path, site_plots_path, output_laz_path):
    laz_path = Path(laz_path)
    output_laz_path = Path(output_laz_path)
    laz_file_paths = [f for f in laz_path.glob('*colorized.laz')]
    shp_file = [i for i in site_plots_path.glob('*.shp')][0]
    polygons_utm = gpd.read_file(shp_file)
    polygons_str_list = [_get_polygon_str(polygons_utm.geometry.iloc[i].exterior.coords.xy[0].tolist(),
                                          polygons_utm.geometry.iloc[i].exterior.coords.xy[1].tolist())
                         for i in range(polygons_utm.shape[0])
                        ]

    for laz_file_path in tqdm(laz_file_paths):
        pdal_json = {
            "pipeline": [
                {
                    "type": "readers.las",
                    "filename": str(laz_file_path)
                },
                {
                    "type": "filters.crop",
                    "polygon": polygons_str_list,
                },
                {
                    "type": "writers.las",
                    "filename": str(output_laz_path/f"{laz_file_path.stem}.laz"),
                    "extra_dims": "all"
                }
            ]
        }
        pdal_json_str = json.dumps(pdal_json)
        pipeline = pdal.Pipeline(pdal_json_str)
        count = pipeline.execute()
        if count == 0:
            os.remove(str(output_laz_path/f"{laz_file_path.stem}.laz"))



