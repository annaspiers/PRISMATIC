from preprocessing.inventory import download_veg_structure_data, \
                                    preprocess_veg_structure_data
from preprocessing.plots import download_neon_polygons, \
                                preprocess_neon_polygons, \
                                preprocess_neon_polygons_site_inventory, \
                                preprocess_neon_polygons_site_inventory_partition
from preprocessing.lidar import clip_laz_by_plots, \
                                clip_laz_by_inventory_plots, \
                                subtract_ground_plots
from preprocessing.biomass import preprocess_biomass

from constants import site, year, data_path, laz_path


# download and preprocessing inventory data
download_veg_structure_data(site)
inventory_file_path = preprocess_veg_structure_data(site, year, data_path)
# # download NEON plots/watch
neon_plots_path = download_neon_polygons(data_path)
# # preprocessing NEON plots
pp_neon_plots_path = preprocess_neon_polygons(neon_plots_path, data_path)
# preprocessing NEON plots given site ID and inventory
pp_neon_plots_site_inventory_path = \
    preprocess_neon_polygons_site_inventory(pp_neon_plots_path,
                                               site,
                                               year,
                                               data_path,
                                               inventory_file_path)
sampling_effort_path = f'{data_path}/{site}/plot_sampling_effort.csv'
# partition plots to subplots
pp_neon_partition_plots_site_inventory_path = \
preprocess_neon_polygons_site_inventory_partition(pp_neon_plots_site_inventory_path,
                                                     sampling_effort_path,
                                                     inventory_file_path,
                                                     site,
                                                     year,
                                                     data_path)

# # clip lidar data
normalized_laz_path = subtract_ground_plots(laz_path,
                                            site,
                                            year,
                                            data_path,
                                            end_result=False)
merged_laz_file = clip_laz_by_plots(normalized_laz_path,
                                    pp_neon_plots_site_inventory_path,
                                    site,
                                    year,
                                    data_path)
clipped_laz_inventory_path = \
    clip_laz_by_inventory_plots(merged_laz_file,
                                pp_neon_partition_plots_site_inventory_path,
                                site,
                                year,
                                data_path,
                                end_result=True)

# # biomass
preprocess_biomass(inventory_file_path,
                      pp_neon_partition_plots_site_inventory_path,
                      sampling_effort_path,
                      site, year,
                      data_path,
                      end_result=True)
