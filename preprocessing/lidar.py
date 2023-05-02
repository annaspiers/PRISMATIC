import json
import logging
import geopandas as gpd
import os
import pdal
import sys
import whitebox

from pathlib import Path
from tqdm import tqdm
from utils.download_functions import download_aop_files

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = str(conda_env_path/'share'/'proj')

log = logging.getLogger(__name__)
wht = whitebox.WhiteboxTools()
wht.set_verbose_mode(False)


def download_lidar(site, date, lidar_path):
    lidar_path = Path(lidar_path)
    product_code = 'DP1.30003.001'
    path = lidar_path/site/date
    file_types = ['prj', 'kml', 'shx', 'shp', 'dbf', 'laz']
    for file_type in file_types:
        if file_type == 'laz':
            p = path/'laz'
        else:
            p = path/'shape'
        download_aop_files(product_code,
                           site,
                           date,
                           str(p),
                           match_string=file_type,
                           check_size=False)

    product_codes = ['DP3.30024.001', 'DP3.30025.001']
    file_type = 'tif'
    for product_code in product_codes:
        p = path/file_type
        download_aop_files(product_code,
                           site,
                           date,
                           str(p),
                           match_string=file_type,
                           check_size=False)
    return str(path/'laz'), str(path/'tif')


def clip_lidar_by_plots(laz_path,
                        tif_path,
                        site_plots_path,
                        site,
                        year,
                        output_laz_path,
                        end_result=False):
    log.info(f'Processing LiDAR data for site: {site} / year: {year}')
    laz_path = Path(laz_path)
    tif_path = Path(tif_path)
    site_plots_path = Path(site_plots_path)
    year = str(year)
    output_folder = 'clipped_lidar' if not end_result else 'output'
    pp_laz_path = Path(output_laz_path)/site/year/'clipped_lidar'
    pp_laz_path.mkdir(parents=True, exist_ok=True)
    output_laz_path = Path(output_laz_path)/site/year/output_folder
    output_laz_path.mkdir(parents=True, exist_ok=True)

    laz_file_paths = [f for f in laz_path.glob('*colorized.laz')]
    shp_file = [i for i in
                site_plots_path.glob('*.shp')
                if 'plots' in str(i)][0]
    polygons_utm = gpd.read_file(shp_file)

    log.info('Cropping lidar files given all plots...')
    for laz_file_path in tqdm(laz_file_paths):
        wht.clip_lidar_to_polygon(
            str(laz_file_path),
            str(shp_file),
            str(pp_laz_path/laz_file_path.name)
        )

    log.info('Merging clipped lidar files...')
    merged_laz_file = str(pp_laz_path/'merge.laz')
    laz_files = [str(i) for i in
                 pp_laz_path.glob('*.laz')
                 if 'merge' not in str(i)]
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

    log.info('Clipping lidar into plot level lidar files...')
    for row in tqdm(polygons_utm.itertuples()):
        pdal_json = {
            "pipeline": [
                {
                    "type": "readers.las",
                    "filename": merged_laz_file
                },
                {
                    "type": "filters.crop",
                    "polygon": _get_polygon_str(row.geometry.
                                                exterior.coords
                                                .xy[0].tolist(),
                                                row.geometry.
                                                exterior.coords
                                                .xy[1].tolist()),
                },
                {
                    "type": "writers.las",
                    "filename": (f'{str(output_laz_path)}'
                                 f'/{row.plotID}_{row.subplotID}.laz'),
                    "extra_dims": "all"
                }
            ]
        }
        pdal_json_str = json.dumps(pdal_json)
        pipeline = pdal.Pipeline(pdal_json_str)
        pipeline.execute()

    shp_paths = []
    for row in tqdm(polygons_utm.itertuples()):
        p = polygons_utm.query(f"plotID == '{row.plotID}' "
                               f"and subplotID == '{row.subplotID}'")
        saved_path = shp_file.parent/f'{row.plotID}_{row.subplotID}.shp'
        p.to_file(saved_path)
        shp_paths.append(saved_path)

    laz_files = [str(i) for i in
                 pp_laz_path.glob('*.laz')
                 if 'merge' not in str(i)]
    for laz_file in tqdm(laz_files):
        tif_filename = '_'.join(Path(laz_file)
                                .stem
                                .replace('DP1', 'DP3')
                                .split('_')[:-4])
        tif_files = [i for i in tif_path.glob(f'{tif_filename}*')]
        for tif_file in tqdm(tif_files):
            for shp_path in tqdm(shp_paths):
                output_file_name = \
                    f"{shp_path.stem}_{tif_file.stem.split('_')[-1]}.tif"
                wht.clip_raster_to_polygon(
                    str(tif_file),
                    str(shp_path),
                    str(output_laz_path/output_file_name)
                )
    log.info(f'Processed LiDAR data for site: {site} / year: {year}')
    return str(output_laz_path)


def normalize_laz(laz_path,
                  site,
                  year,
                  output_path,
                  end_result=False):
    log.info(f'Normalizing LiDAR data for site: {site} / year: {year}')
    laz_path = Path(laz_path)
    year = str(year)
    output_folder = 'normalized_lidar' if not end_result else 'output'
    output_path = Path(output_path)/site/year/output_folder
    output_path.mkdir(parents=True, exist_ok=True)
    laz_file_paths = [i for i in laz_path.glob('*.laz')]

    for laz_path in tqdm(laz_file_paths):
        wht.height_above_ground(
            i=str(laz_path),
            output=str(output_path/f'{laz_path.stem}.laz')
        )
    log.info(f'Normalized LiDAR data for site: {site} / year: {year}')
    return output_path


def _get_polygon_str(x_cord, y_cord):
    polygon_str = 'POLYGON(('
    for x, y in zip(list(x_cord), list(y_cord)):
        polygon_str += f'{x} {y}, '
    polygon_str = polygon_str[:-2]
    polygon_str += '))'
    return polygon_str
