import logging
import json
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
INVENTORY_PLOTS_FOLDER = 'inventory_plots'
METADATA = 'metadata'
SUBPLOT_AREA_UNIT = 400

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


def preprocess_polygons(input_data_path,
                        sampling_effort_path,
                        inventory_path,
                        site, year,
                        output_data_path):
    log.info(f'Processing polygons for site: {site} / '
             f'year: {year} given inventory')
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
    subplot_ids = []
    ps = []
    processed_plots = {}
    for plot_id, group in polygons_site_utm.groupby('plotID'):
        veg_plot_metadata = {
            'plot_id': plot_id,
            'total_ind': None,
            'percent_ind_have_location': None,
            'total_ind_stem_gt_10': None,
            'total_ind_stem_lt_10': None,
            'sampling_effort_trees': None,
            'sampling_effort_sapling': None,
            'clipped_subplot_position': None
        }
        veg_gdf = avail_veg_gdf[avail_veg_gdf.plotID == plot_id]
        veg_plot_metadata['total_ind'] = veg_gdf.shape[0]
        veg_plot_metadata['percent_ind_have_location'] = \
            round(sum(~veg_gdf.geometry.is_empty)
                  / veg_gdf.shape[0]*100,
                  2)
        veg_plot_metadata['total_ind_stem_gt_10'] = \
            sum(veg_gdf.stemDiameter >= 10)
        veg_plot_metadata['total_ind_stem_lt_10'] = \
            sum(veg_gdf.stemDiameter < 10)

        sampling_area_trees = group.area.values[0]
        sampling_area_sapling = group.area.values[0]
        try:
            query = f'plotID == "{plot_id}"'
            sampling_area = sampling_effort.query(query)
            v = (sampling_area
                 .totalSampledAreaTrees
                 .values[0])
            if np.isnan(v):
                log.warning(f'Site: {site}, plot: {plot_id} '
                            'has empty sampling_area_trees')
            else:
                sampling_area_trees = v
            v = (sampling_area
                 .totalSampledAreaShrubSapling
                 .values[0])
            if np.isnan(v):
                log.warning(f'Site: {site}, plot: {plot_id} '
                            'has empty sampling_area_sapling')
            else:
                sampling_area_sapling = v
        except IndexError:
            log.warning(f'Site: {site}, plot: {plot_id} '
                         'does not have in: pp_sampling_effort.csv')
        sampled_subplots = sampling_area.subplotsSampled.values[0]
        subplots = [v.strip() for v in sampled_subplots.split('|')]
        sampling_area_trees = SUBPLOT_AREA_UNIT
        sampling_area_sapling = SUBPLOT_AREA_UNIT
        if np.isnan(sampling_area_trees) and len(subplots) < 4:
            sampling_area_trees = SUBPLOT_AREA_UNIT * len(subplots)
        if np.isnan(sampling_area_sapling) and len(subplots) < 4:
            sampling_area_sapling = SUBPLOT_AREA_UNIT * len(subplots)

        veg_plot_metadata['sampling_effort_trees'] = sampling_area_trees
        veg_plot_metadata['sampling_effort_sapling'] = sampling_area_sapling

        for row in group.itertuples():
            p = plot_polygon = row.geometry
            veg_gdf_position_list = veg_gdf.geometry[~veg_gdf
                                                     .geometry
                                                     .is_empty]
            tree_in_plot = sum(plot_polygon.contains(veg_gdf_position_list))
            
            if plot_id not in processed_plots or \
               tree_in_plot > processed_plots[plot_id]:
                processed_plots[plot_id] = tree_in_plot
                if row.geometry.area > sampling_area_trees:
                    pplots = partition(plot_polygon,
                                      np.sqrt(SUBPLOT_AREA_UNIT),
                                      subplots)
                if len(pplots) == 1:
                    p = pplots[0]
                    names.append(plot_id)
                    ps.append(p)
                    subplot_ids.append('central')
                else:
                    for subplot, p in zip(subplots, pplots):
                        names.append(plot_id)
                        ps.append(p)
                        subplot_ids.append(subplot)
                veg_plot_metadata['clipped_subplot_position'] = subplots
                if tree_in_plot == 0:
                    log.warning(f'{plot_id} does not have any tree location. '
                                'Selecting the first plot encountered as '
                                'default. Clipping to the central area given '
                                'the totalSampledAreaTrees in sampling_effort.'
                                )

                # save result for diagnostics
                fig, ax = plt.subplots(figsize=(5, 5))
                gpd.GeoSeries([plot_polygon]).boundary.plot(ax=ax, color="red")
                for p in pplots:
                    gpd.GeoSeries(p).boundary.plot(ax=ax)
                veg_gdf.plot(ax=ax, color="red")
                plt.title(label=f"Site: {plot_id} {sampled_subplots}")
                output_folder_path = \
                    (output_data_path/'diagnostics'/site
                     / year/INVENTORY_PLOTS_FOLDER)
                output_folder_path.mkdir(parents=True, exist_ok=True)
                fig.savefig(output_folder_path/f'{plot_id}.png')
                plt.close()
                output_folder_metadata_path = \
                    (output_data_path/'diagnostics'/site
                     / year/METADATA)
                output_folder_metadata_path.mkdir(parents=True, exist_ok=True)
                f_name = output_folder_metadata_path/f'{plot_id}.json'
                with open(f_name, 'w') as f:
                    json.dump(veg_plot_metadata, f, indent=4)

    df = gpd.GeoDataFrame(data=zip(names, subplot_ids, ps), columns=['plotID', 'subplotID', 'geometry'],
                          crs=polygons_site_utm.crs)
    output_folder_path = \
        output_data_path/site/year/INVENTORY_PLOTS_FOLDER
    output_folder_path.mkdir(parents=True, exist_ok=True)
    df.to_file(output_folder_path/'plots.shp')
    log.info(f'Processed polygons saved at: {output_folder_path}\n'
             'Diagnostics and metadata saved at: '
             f'{output_folder_metadata_path}')
    return str(output_folder_path)
