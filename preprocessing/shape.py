import geopandas as gpd
import pandas as pd
import requests
from pathlib import Path
from zipfile import ZipFile

NEON_POLYGONS_LINK = 'https://www.neonscience.org/sites/default/files/All_NEON_TOS_Plots_V9_0.zip'
OUTPUT_FOLDERNAME = 'All_NEON_TOS_Plots_V9'
EPSG = 'epsg:32611'
PLOTS_FOLDER = 'plots'
INVENTORY_PLOTS_FOLDER = 'inventory_plots'

def download_neon_polygons(data_path):
    data_path = Path(data_path)
    zip_filename = Path(NEON_POLYGONS_LINK).name
    zip_data_path = data_path/zip_filename
    data = requests.get(NEON_POLYGONS_LINK)
    with open(zip_data_path, 'wb') as f:
        f.write(data.content)

    with ZipFile(zip_data_path, 'r') as zf:
        zf.extractall(path=data_path)

    output_folder_path = data_path/OUTPUT_FOLDERNAME
    return output_folder_path


def preprocessing_neon_polygons_site(data_path, site):
    data_path = Path(data_path)
    downloaded_data_path = data_path/OUTPUT_FOLDERNAME
    shp_file = [i for i in downloaded_data_path.glob('*.shp')][0]
    polygons = gpd.read_file(shp_file)
    polygons_site = polygons[polygons.siteID == site]
    polygons_site_utm = polygons_site.to_crs(EPSG)
    polygons_site_utm.to_file(data_path/site/PLOTS_FOLDER/'plots.shp')

def preprocessing_neon_polygons_site_inventory(data_path, site, inventory_file_path):
    data_path = Path(data_path)
    polygons_site_utm = gpd.read_file(data_path/site/PLOTS_FOLDER/f'{PLOTS_FOLDER}.shp')
    inventory = pd.read_csv(inventory_file_path)
    inventory_avail_plots = inventory.plotID.unique()
    filtered_polygons_site_utm = polygons_site_utm[polygons_site_utm.plotID.isin(inventory_avail_plots)]
    filtered_polygons_site_utm.to_file(data_path/site/INVENTORY_PLOTS_FOLDER/'plots.shp')

download_neon_polygons('/home/toanngo/Documents/GitHub/prisma/preprocessing/data')