import logging
import geopandas as gpd
import pandas as pd
import requests
from pathlib import Path
from zipfile import ZipFile
from utils.plot_partition import partition

NEON_POLYGONS_LINK = ('https://www.neonscience.org/'
                      'sites/default/files/All_NEON_TOS_Plots_V9_0.zip')
OUTPUT_FOLDERNAME = 'All_NEON_TOS_Plots_V9'
EPSG = 'epsg:32611'
PLOTS_FOLDER = 'plots'
INVENTORY_PLOTS_FOLDER = 'inventory_plots'
INVENTORY_PARTITIONED_PLOTS_FOLDER = 'inventory_partitioned_plots'
PLOT_PARTITION_SIZE = 20

log = logging.getLogger(__name__)


def download_polygons(data_path):
    data_path = Path(data_path)
    zip_filename = Path(NEON_POLYGONS_LINK).name
    zip_data_path = data_path/zip_filename
    data = requests.get(NEON_POLYGONS_LINK)
    with open(zip_data_path, 'wb') as f:
        f.write(data.content)

    with ZipFile(zip_data_path, 'r') as zf:
        zf.extractall(path=data_path)

    output_folder_path = str(data_path/OUTPUT_FOLDERNAME)
    log.info(f'Downloaded NEON plots saved at: {output_folder_path}')
    return output_folder_path


def preprocess_polygons(input_data_path, output_data_path):
    input_data_path = Path(input_data_path)
    output_data_path = Path(output_data_path)
    shp_file = [i for i in input_data_path.glob('*Polygons*.shp')][0]
    polygons = gpd.read_file(shp_file)
    polygons_utm = polygons.to_crs(EPSG)
    output_data_path = output_data_path/PLOTS_FOLDER
    output_data_path.mkdir(parents=True, exist_ok=True)
    polygons_utm.to_file(output_data_path/'plots.shp')
    log.info(f'Processed NEON plots saved at: {output_data_path}')
    return str(output_data_path)


def preprocess_polygons_site_inventory(input_data_path,
                                       site, year,
                                       output_data_path,
                                       inventory_file_path=None):
    input_data_path = Path(input_data_path)
    output_data_path = Path(output_data_path)
    year = str(year)
    polygons_utm = gpd.read_file(input_data_path/'plots.shp')
    polygons = polygons_utm[polygons_utm.siteID == site]
    polygons = polygons[polygons.pointID.isin(['31', '32',
                                               '40', '41'])]
    if inventory_file_path:
        inventory = pd.read_csv(inventory_file_path)
        inventory_avail_plots = inventory.plotID.unique()
        polygons = polygons[polygons.plotID.isin(inventory_avail_plots)]

    output_folder_path = output_data_path/site/year/INVENTORY_PLOTS_FOLDER
    output_folder_path.mkdir(parents=True, exist_ok=True)
    polygons.to_file(output_folder_path/'plots.shp')
    return str(output_folder_path)


def preprocess_polygons_site_inventory_partition(input_data_path,
                                                 sampling_effort_path,
                                                 inventory_path,
                                                 site, year,
                                                 output_data_path):
    input_data_path = Path(input_data_path)
    output_data_path = Path(output_data_path)
    year = str(year)
    # sampling_effort = pd.read_csv(sampling_effort_path)
    # inventory = pd.read_csv(inventory_path)
    polygons = gpd.read_file(input_data_path/'plots.shp')

    name = []
    ps = []
    for polygon in polygons.itertuples():
        plot_id = polygon.plotID
        plot_polygon = polygon.geometry
        plot_partitions = partition(plot_polygon,
                                    PLOT_PARTITION_SIZE,
                                    mode='center')
        p = plot_partitions[0]
        name.append(plot_id)
        ps.append(p)
        # idxs = ['31', '40', '32', '41']
        # plot_partitions = partition(plot_polygon,
        #                             PLOT_PARTITION_SIZE)
        # for i, p in zip(idxs, plot_partitions):
        #     name.append(f'{plot_id}_{i}')
        #     ps.append(p)

    df = gpd.GeoDataFrame(data=zip(name, ps), columns=['plotID', 'geometry'],
                          crs=polygons.crs)
    output_folder_path = \
        output_data_path/site/year/INVENTORY_PARTITIONED_PLOTS_FOLDER
    output_folder_path.mkdir(parents=True, exist_ok=True)
    df.to_file(output_folder_path/'plots.shp')
    return str(output_folder_path)
