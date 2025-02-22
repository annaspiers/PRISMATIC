import json
import logging
import geopandas as gpd
import os
import pdal
import sys
import numpy as np
# import whitebox
import rpy2.robjects as ro
import shutil
import zipfile

from pathlib import Path
from tqdm import tqdm
from utils.download_functions import download_aop_files

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = str(conda_env_path/'share'/'proj')

log = logging.getLogger(__name__)
# wht = whitebox.WhiteboxTools()
# wht.set_verbose_mode(False)


def download_aop_bbox(site, year, path, hs_type, coords_bbox):
    """Download lidar data (raw and raster) for site-year

    Parameters
    ----------
    site : str
        Site name
    year : str
        Lidar year
    path : str
        Path to store the downloaded data
    coords_bbox : str
        UTM coordinates for bounding box of data subset to download 

    Returns
    -------
    (str, str)
        Path to the result lidar folder ('*/laz'), hyperspectral folder 
        ('*/hs_tile' or '*/hs_flightline'), and the result raster folder ('*/tif')
    """
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'leaf_area_density_helper.R'))
    download_by_tile = ro.r('download_by_tile')
    move_downloaded_files = ro.r('move_downloaded_files')

    # expand coords_bbox 
    # Generate easting and northing coordinates, 1000m apart
    easting = np.arange(coords_bbox[0], coords_bbox[1] , 1000)
    northing = np.arange(coords_bbox[2], coords_bbox[3] , 1000)

    # Create a meshgrid of coordinates
    easting_grid, northing_grid = np.meshgrid(easting, northing)

    # Flatten the grids into 1D arrays
    easting_flat = easting_grid.flatten()
    northing_flat = northing_grid.flatten()

    path = Path(path)/site/year
    
    # lidar and mosaic/tile shp
    product_code = 'DP1.30003.001' #lidar
    # p = path/'laz'
    download_by_tile(dpID=product_code, 
                        site=site, year=year, 
                        eastings=ro.IntVector(easting_flat),
                        northings=ro.IntVector(northing_flat), 
                        savepath=str(path) )
    # Move shp stuff to appropriate folders
    file_types = ['prj', 'shx', 'shp', 'dbf', 'merged_tiles', 'laz']
    for file_type in file_types:
        if file_type == 'laz':
            file_type = '_classified_point_cloud_colorized.laz'
            # p = path/'laz'
            move_downloaded_files(dir_out = str(path), dp_id = product_code
                     ,dp_name = "laz", file_pattern = file_type, delete_orig=True)
        else:
            # p = path/'shape'
            move_downloaded_files(dir_out = str(path), dp_id = product_code
                     ,dp_name = "shape", file_pattern = file_type, delete_orig=False)

    # rasters
    product_codes = ['DP3.30024.001', 'DP3.30025.001', #DTM, slope/aspect
                     'DP3.30015.001', 'DP3.30010.001'] #chm, rgb
    for product_code in product_codes:
        download_by_tile(dpID=product_code, 
                        site=site, year=year, 
                        eastings=ro.IntVector(easting_flat),
                        northings=ro.IntVector(northing_flat), 
                        savepath=str(path) )
        move_downloaded_files(dir_out = str(path), dp_id = product_code
                     ,dp_name = "tif", file_pattern =  ".tif", delete_orig=True)
        
    product_code = 'DP3.30026.001' #veg indices
    download_by_tile(dpID=product_code, 
                        site=site, year=year, 
                        eastings=ro.IntVector(easting_flat),
                        northings=ro.IntVector(northing_flat), 
                        savepath=str(path) )
    move_downloaded_files(dir_out = str(path), dp_id = product_code
                     ,dp_name = "tif", file_pattern = ".zip", delete_orig=True, unzip=True)
    
    # zip_files = [file for file in os.listdir(p) if file.endswith('.zip')] # get the list of files
    # for zip_file in zip_files:  #for each zipfile
    #     with zipfile.ZipFile(Path(p/zip_file)) as item: # treat the file as a zip
    #         item.extractall(p)  # extract it in the working directory
    #         item.close()
    vi_error_files = [file for file in os.listdir(Path(path/"tif")) if file.endswith('_error.tif')]
    for vi_error_file in vi_error_files: 
        os.remove(Path(path/"tif"/vi_error_file)) # remove error files

    # hyperspectral
    if hs_type=="tile":
        product_code = 'DP3.30006.001'
    elif hs_type=="flightline": 
        product_code = 'DP1.30006.001'
    else:
        print("must specify hs_type argument")
    download_by_tile(dpID=product_code, 
                        site=site, year=year, 
                        eastings=ro.IntVector(easting_flat),
                        northings=ro.IntVector(northing_flat), 
                        savepath=str(path) )
    move_downloaded_files(dir_out = str(path), dp_id = product_code
                     ,dp_name = "hs_"+hs_type, file_pattern =  ".h5", delete_orig=True)
    
    return str(path/'laz'), str(path/hs_type), str(path/'tif')


def download_lidar(site, year, lidar_path, use_tiles_w_veg):
    """Download lidar data (raw and raster) for site-year

    Parameters
    ----------
    site : str
        Site name
    year : str
        Lidar year
    lidar_path : str
        Path to store the downloaded data

    Returns
    -------
    (str, str)
        Path to the result lidar folder ('*/laz')
        and the result raster folder ('*/tif')
    """

    lidar_path = Path(lidar_path)
    path = lidar_path/site/year
        
    product_code = 'DP1.30003.001' #lidar
    file_types = ['prj', 'shx', 'shp', 'dbf', 'merged_tiles', 'laz']
    for file_type in file_types:
        if file_type == 'laz':
            file_type = '_classified_point_cloud_colorized.laz'
            p = path/'laz'
        else:
            p = path/'shape'
        download_aop_files(product_code,
                        site,
                        year,
                        str(p),
                        match_string=file_type,
                        check_size=False)

    product_codes = ['DP3.30024.001', 'DP3.30025.001'] #DTM, slope/aspect
    file_type = 'tif'
    for product_code in product_codes:
        p = path/file_type
        download_aop_files(product_code,
                        site,
                        year,
                        str(p),
                        match_string=file_type,
                        check_size=False,
                        use_tiles_w_veg=use_tiles_w_veg)
            
    return str(path/'laz'), str(path/'tif')

# def download_lidar(site, year, lidar_path, use_tiles_w_veg, coords_bbox=[]):
#     """Download lidar data (raw and raster) for site-year

#     Parameters
#     ----------
#     site : str
#         Site name
#     year : str
#         Lidar year
#     lidar_path : str
#         Path to store the downloaded data

#     Returns
#     -------
#     (str, str)
#         Path to the result lidar folder ('*/laz')
#         and the result raster folder ('*/tif')
#     """

#     if coords_bbox:
#         # Generate easting and northing coordinates, 1000m apart
#         easting = np.arange(coords_bbox[0], coords_bbox[1] + 1, 1000)
#         northing = np.arange(coords_bbox[2], coords_bbox[3] + 1, 1000)

#         # Create a meshgrid of coordinates
#         easting_grid, northing_grid = np.meshgrid(easting, northing)

#         # Flatten the grids into 1D arrays
#         easting_flat = easting_grid.flatten()
#         northing_flat = northing_grid.flatten()

#         # Combine the coordinates into a single array of strings
#         coords = np.array([f"{easting}_{northing}" for easting, northing in zip(easting_flat, northing_flat)])

#     lidar_path = Path(lidar_path)
#     path = lidar_path/site/year
    
#     for coord in coords:
#         product_code = 'DP1.30003.001' #lidar
#         file_types = ['prj', 'shx', 'shp', 'dbf', 'laz']
#         for file_type in file_types:
#             if file_type == 'laz':
#                 file_type = '_classified_point_cloud_colorized.laz'
#                 p = path/'laz'
#             else:
#                 p = path/'shape'
#             download_aop_files(product_code,
#                             site,
#                             year,
#                             str(p),
#                             match_string=[coord,file_type],
#                             check_size=False)
            
#         file_types = ['merged_tiles']
#         for file_type in file_types:
#             p = path/'shape'
#             download_aop_files(product_code,
#                             site,
#                             year,
#                             str(p),
#                             match_string=[coord,file_type],
#                             check_size=False)

#         product_codes = ['DP3.30024.001', 'DP3.30025.001'] #DTM, slope/aspect
#         file_type = 'tif'
#         for product_code in product_codes:
#             p = path/file_type
#             download_aop_files(product_code,
#                             site,
#                             year,
#                             str(p),
#                             match_string=[coord,file_type],
#                             check_size=False,
#                             use_tiles_w_veg=use_tiles_w_veg)
            
#     return str(path/'laz'), str(path/'tif')


def clip_lidar_by_plots(laz_path,
                        tif_path,
                        site_plots_path,
                        site,
                        year,
                        output_laz_path,
                        end_result=False):
    """_summary_

    Parameters
    ----------
    laz_path : str
        Path to the lidar folder
        Format '*/laz'
    tif_path : str
        Path to the raster folder
        Format '*/tif'
    site_plots_path : str
        Path to the shp folder
        Format '*/shape'
    site : str
        Site name
    year : str
        Inventory year
    output_laz_path : str
         Default output data root
    end_result : bool, optional
        If this is the end result,
        redirect the result into '*/output' folder.

    Returns
    -------
    str
        Path to the folder that the result of this function is saved to
    """
    log.info(f'Clipping laz and tif data to plots for site: {site} / year: {year}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'leaf_area_density_helper.R'))
    clip_lidar_to_polygon_lidR = ro.r('clip_lidar_to_polygon_lidR')
    clip_raster_to_polygon_lidR = ro.r("clip_raster_to_polygon_lidR")
    
    laz_path = Path(laz_path)
    tif_path = Path(tif_path)
    site_plots_path = Path(site_plots_path)
    year = str(year)
    output_folder = 'clipped_lidar_tiles' if not end_result else 'clipped_to_plots'
    pp_laz_path = Path(output_laz_path)/site/year/'clipped_lidar_tiles'
    pp_laz_path.mkdir(parents=True, exist_ok=True)
    output_laz_path = Path(output_laz_path)/site/year/output_folder
    output_laz_path.mkdir(parents=True, exist_ok=True)

    laz_file_paths = [f for f in laz_path.glob('*colorized.laz')]
    shp_file = [i for i in site_plots_path.glob('plots.shp')][0] 
    #otherwise only the first plot was being saved
                # [i for i in
                # site_plots_path.glob('*.shp')
                # if 'plots' in str(i)][0]
    polygons_utm = gpd.read_file(shp_file)
    # if ic_type=="rs_inv_plots":
    # if 'plotID' in polygons_utm.columns and 'subplotID' in polygons_utm.columns: 
    #     polygons_utm['plotID'] = polygons_utm['plotID'] + '_' + polygons_utm['subplotID']
    #ais ^ fix this so that plotID and subplotID are connected earlier on and so I don't ahve to do this ad hoc


    log.info('Cropping lidar files given all plots...')
        
    #for laz_file_path in tqdm(laz_file_paths):
    clip_lidar_to_polygon_lidR(str(laz_path),
            str(shp_file),
            str(output_laz_path)) 
        #ais ^ need to make this ore efficient - gets the job done for now though

        # wht.clip_lidar_to_polygon(
        #     str(laz_file_path),
        #     str(shp_file),
        #     str(pp_laz_path/laz_file_path.name)
        # )

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

    log.info('Clipping lidar into plot level files...')
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
                } #,
                # {
                #     "type": "writers.las",
                #     "filename": (f'{str(output_laz_path)}'
                #                  f'/{row.plotID}_{row.subplotID}.laz'),
                #     "extra_dims": "all"
                # }
            ]
        }
        pdal_json_str = json.dumps(pdal_json)
        pipeline = pdal.Pipeline(pdal_json_str)
        pipeline.execute()

    shp_paths = []
    for row in tqdm(polygons_utm.itertuples()):
        p = polygons_utm.query(f"plotID == '{row.plotID}' "
                               f"and subplotID == '{row.subplotID}'")
        saved_path = shp_file.parent/f'{row.plotID}.shp'
        # saved_path = shp_file.parent/f'{row.plotID}_{row.subplotID}.shp'
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
                clip_raster_to_polygon_lidR( str(tif_file),
                    str(shp_path),
                    str(output_laz_path/output_file_name))
                # wht.clip_raster_to_polygon(
                #     str(tif_file),
                #     str(shp_path),
                #     str(output_laz_path/output_file_name)
                # )
    log.info(f'Processed LiDAR data for site: {site} / year: {year}')
    return str(output_laz_path)


def normalize_laz(laz_path,
                  site,
                  year,
                  output_path,
                  end_result=False):
    """Normalize laz files

    Parameters
    ----------
    laz_path : str
        Path to the lidar folder
        Format '*/laz'
    site : str
        Site name
    year : str
        Inventory year
    output_path : str
        Default output data root
    end_result : bool, optional
        If this is the end result,
        redirect the result into '*/output' folder.

    Returns
    -------
    str
        Path to the folder that the result of this function is saved to
    """
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'leaf_area_density_helper.R'))
    normalize_lidR = ro.r('normalize_lidR')

    log.info(f'Normalizing LiDAR data for site: {site} / year: {year}')
    laz_path = Path(laz_path)
    year = str(year)
    output_folder = 'normalized_lidar_tiles' if not end_result else 'output'
    output_path = Path(output_path)/site/year/output_folder
    output_path.mkdir(parents=True, exist_ok=True)
    laz_file_paths = [i for i in laz_path.glob('*.laz')]
    
    for laz_path in tqdm(laz_file_paths):
        if not os.path.isfile(str(output_path/f'{laz_path.stem}.laz')):
            normalize_lidR(str(laz_path),str(output_path/f'{laz_path.stem}'))
        # wht.height_above_ground(
        #     i=str(laz_path),
        #     output=str(output_path/f'{laz_path.stem}.laz')
        # )
        
        # since whitebox isn't working, use lidR
        #copy laz files to noramlized folder
        #then normalize them
    log.info(f'Normalized LiDAR data for site: {site} / year: {year}')
    return output_path


def _get_polygon_str(x_cord, y_cord):
    """
    Prepare the polygon string for clipping with PDAL
    """
    polygon_str = 'POLYGON(('
    for x, y in zip(list(x_cord), list(y_cord)):
        polygon_str += f'{x} {y}, '
    polygon_str = polygon_str[:-2]
    polygon_str += '))'
    return polygon_str
