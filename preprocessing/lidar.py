import geopandas as gpd
import json
import os
import pdal
import sys
import whitebox

from pathlib import Path
from tqdm import tqdm

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = str(conda_env_path/'share'/'proj')


def _get_polygon_str(x_cord, y_cord):
    polygon_str = 'POLYGON(('
    for x, y in zip(list(x_cord), list(y_cord)):
        polygon_str += f'{x} {y}, '
    polygon_str = polygon_str[:-2]
    polygon_str += '))'
    return polygon_str


def clip_laz_by_plots(laz_path,
                      site_plots_path,
                      site, year,
                      output_laz_path,
                      end_result=False):
    laz_path = Path(laz_path)
    site_plots_path = Path(site_plots_path)
    year = str(year)
    output_folder = 'clipped_lidar' if not end_result else 'output'
    pp_laz_path = Path(output_laz_path)/site/year/'clipped_lidar'
    pp_laz_path.mkdir(parents=True, exist_ok=True)
    output_laz_path = Path(output_laz_path)/site/year/output_folder
    output_laz_path.mkdir(parents=True, exist_ok=True)

    laz_file_paths = [f for f in laz_path.glob('*colorized.laz')]
    shp_file = [i for i in site_plots_path.glob('*.shp')][0]
    polygons_utm = gpd.read_file(shp_file)
    polygons_str_list = [_get_polygon_str(polygons_utm.geometry.
                                          iloc[i].exterior.coords.
                                          xy[0].tolist(),
                                          polygons_utm.geometry.
                                          iloc[i].exterior.coords.
                                          xy[1].tolist())
                         for i in range(polygons_utm.shape[0])]

    # crop each lidar file with all plots
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
                    "filename": str(pp_laz_path/laz_file_path.name),
                    "extra_dims": "all"
                }
            ]
        }
        pdal_json_str = json.dumps(pdal_json)
        pipeline = pdal.Pipeline(pdal_json_str)
        count = pipeline.execute()
        if count == 0:
            os.remove(str(pp_laz_path/laz_file_path.name))

    # merge clipped lidar files into merge.laz
    merged_laz_file = str(pp_laz_path/'merge.laz')
    laz_files = [str(i) for i in pp_laz_path.glob('*.laz')]
    pdal_json = {
        "pipeline": laz_files + [
            {
                "type": "filters.merge"
            },
            {
                "type": "writers.las",
                "filename": merged_laz_file,
                "extra_dims": "all"
            }
        ]
    }
    pdal_json_str = json.dumps(pdal_json)
    pipeline = pdal.Pipeline(pdal_json_str)
    pipeline.execute()

    # clip merged lidar file again into plot level lidar files
    for row in tqdm(polygons_utm.itertuples()):
        pdal_json = {
            "pipeline": [
                {
                    "type": "readers.las",
                    "filename": merged_laz_file
                },
                {
                    "type": "filters.crop",
                    "polygon": _get_polygon_str(row.geometry.exterior.
                                                coords.xy[0].tolist(),
                                                row.geometry.exterior.
                                                coords.xy[1].tolist()),
                },
                {
                    "type": "writers.las",
                    "filename": f"{str(output_laz_path)}/{row.plotID}.laz",
                    "extra_dims": "all"
                }
            ]
        }
        pdal_json_str = json.dumps(pdal_json)
        pipeline = pdal.Pipeline(pdal_json_str)
        pipeline.execute()
    return str(output_laz_path)


def normalize_laz(laz_path,
                  site,
                  year,
                  output_path,
                  end_result=False):
    laz_path = Path(laz_path)
    year = str(year)
    output_folder = 'subtracted_ground_lidar' if not end_result else 'output'
    output_path = Path(output_path)/site/year/output_folder
    output_path.mkdir(parents=True, exist_ok=True)
    laz_file_paths = [i for i in laz_path.glob('*.laz')]

    wht = whitebox.WhiteboxTools()
    for laz_path in laz_file_paths:
        wht.height_above_ground(
            i=str(laz_path),
            output=str(output_path/f'{laz_path.stem}.laz')
        )
    return output_path
