import numpy as np
import geopandas as gpd
import pandas as pd
import requests
from pathlib import Path
from shapely import Point
from zipfile import ZipFile
from utils.plot_partition import partition
import matplotlib.pyplot as plt

NEON_POLYGONS_LINK = ('https://www.neonscience.org/'
                      'sites/default/files/All_NEON_TOS_Plots_V9_0.zip')
OUTPUT_FOLDERNAME = 'All_NEON_TOS_Plots_V9'
EPSG = 'epsg:32611'
PLOTS_FOLDER = 'plots'
INVENTORY_PLOTS_FOLDER = 'inventory_plots'
PLOT_PARTITION_SIZE = 20


def download_polygons(data_path):
    data_path = Path(data_path)
    zip_filename = Path(NEON_POLYGONS_LINK).name
    zip_data_path = data_path/zip_filename
    data = requests.get(NEON_POLYGONS_LINK)
    with open(zip_data_path, 'wb') as f:
        f.write(data.content)

    with ZipFile(zip_data_path, 'r') as zf:
        zf.extractall(path=data_path)

    output_folder_path = data_path/OUTPUT_FOLDERNAME
    return str(output_folder_path)


def preprocess_polygons(input_data_path,
                        sampling_effort_path,
                        inventory_path,
                        site, year,
                        output_data_path):
    input_data_path = Path(input_data_path)
    output_data_path = Path(output_data_path)
    year = str(year)

    shp_file = [i for i in input_data_path.glob('*Polygons*.shp')][0]
    polygons = gpd.read_file(shp_file)
    sampling_effort = pd.read_csv(sampling_effort_path)
    veg_df = pd.read_csv(inventory_path)
    avail_veg_df = veg_df[(~pd.isna(veg_df.basalStemDiameter) |
                           ~pd.isna(veg_df.stemDiameter))]
    geometry = [Point(xy) for xy in zip(avail_veg_df['adjDecimalLongitude'],
                                        avail_veg_df['adjDecimalLatitude'])]
    polygons_site = polygons[polygons.siteID == site]
    avail_veg_gdf = gpd.GeoDataFrame(avail_veg_df,
                                     crs=polygons_site.crs,
                                     geometry=geometry)
    avail_veg_gdf = avail_veg_gdf.to_crs(EPSG)
    polygons_site_utm = polygons_site.to_crs(EPSG)
    polygons_site_utm = polygons_site_utm[polygons_site_utm
                                          .pointID
                                          .isin(['31', '32', '40', '41'])]
    polygons_site_utm = polygons_site_utm[polygons_site_utm
                                          .plotID
                                          .isin(avail_veg_df.plotID.unique())]
    names = []
    ps = []
    processed_plots = {}
    for plot_id, group in polygons_site_utm.groupby('plotID'):
        veg_gdf = avail_veg_gdf[avail_veg_gdf.plotID == plot_id]
        query = f'plotID == "{plot_id}" and plotType == "distributed"'
        sampling_area = sampling_effort.query(query)\
            .totalSampledAreaTrees.values[0]
        sampling_side = int(np.sqrt(sampling_area))

        for row in group.itertuples():
            p = plot_polygon = row.geometry
            tree_inside_plot = sum(plot_polygon.contains(veg_gdf.geometry))
            if plot_id not in processed_plots or \
               tree_inside_plot > processed_plots[plot_id]:

                processed_plots[plot_id] = tree_inside_plot
                if row.geometry.area > sampling_area:
                    # perform clipping
                    pplot = partition(plot_polygon,
                                      sampling_side,
                                      mode='center')
                    p = pplot[0]
                    tree_inside_subplot = sum(p.contains(veg_gdf.geometry))

                    idxs = ['31', '40', '32', '41']
                    pplots = partition(plot_polygon,
                                       sampling_side)
                    for _, plot in zip(idxs, pplots):
                        n = sum(plot.contains(veg_gdf.geometry))
                        if tree_inside_subplot < n:
                            tree_inside_subplot = n
                            p = plot
                names.append(plot_id)
                ps.append(p)
                fig, ax = plt.subplots(figsize=(5, 5))
                gpd.GeoSeries(p).boundary.plot(ax=ax)
                gpd.GeoSeries([plot_polygon]).boundary.plot(ax=ax, color="red")
                veg_gdf.plot(ax=ax, color="red")
                plt.title(label=f"Site: {plot_id}")
                output_folder_path = \
                    output_data_path/site \
                    / year/INVENTORY_PLOTS_FOLDER
                fig.savefig(output_folder_path/f'{plot_id}.png')
                plt.close()

    df = gpd.GeoDataFrame(data=zip(names, ps), columns=['plotID', 'geometry'],
                          crs=polygons.crs)
    output_folder_path = \
        output_data_path/site/year/INVENTORY_PLOTS_FOLDER
    output_folder_path.mkdir(parents=True, exist_ok=True)
    df.to_file(output_folder_path/'plots.shp')
    return str(output_folder_path)
