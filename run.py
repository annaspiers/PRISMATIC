from preprocessing.inventory import download_veg_structure_data, \
                                    preprocessing_veg_structure_data
from preprocessing.plots import download_neon_polygons, \
                                preprocessing_neon_polygons, \
                                preprocessing_neon_polygons_site_inventory
from preprocessing.lidar import clip_laz_by_plots, \
                                clip_laz_by_inventory_plots
from preprocessing.biomass import preprocessing_biomass

from constants import site, year, data_path, laz_path, output_laz_path


# download and preprocessing inventory data
download_veg_structure_data(site)
inventory_file_path = preprocessing_veg_structure_data(site, year, data_path)
# inventory_file_path = '/home/toanngo/Documents/GitHub/prisma/preprocessing/data/SOAP/2019/pp_veg_structure.csv'
# download NEON plots/watch
neon_plots_path = download_neon_polygons(data_path)
# neon_plots_path = '/home/toanngo/Documents/GitHub/prisma/preprocessing/data/All_NEON_TOS_Plots_V9'
# preprocessing NEON plots
pp_neon_plots_path = preprocessing_neon_polygons(neon_plots_path, data_path)
# pp_neon_plots_path = '/home/toanngo/Documents/GitHub/prisma/preprocessing/data/plots'
# preprocessing NEON plots given site ID and inventory
pp_neon_plots_site_inventory_path = preprocessing_neon_polygons_site_inventory(pp_neon_plots_path,
                                                                               site,
                                                                               year,
                                                                               data_path,
                                                                               inventory_file_path)
# pp_neon_plots_site_inventory_path = '/home/toanngo/Documents/GitHub/prisma/preprocessing/data/SOAP/2019/inventory_plots'
# clip lidar data
merged_laz_file = clip_laz_by_plots(laz_path,
                                    pp_neon_plots_site_inventory_path,
                                    site,
                                    year,
                                    data_path)
# merged_laz_file = '/home/toanngo/Documents/GitHub/prisma/preprocessing/data/SOAP/2019/clipped_lidar/merge.laz'
clipped_laz_inventory_path = clip_laz_by_inventory_plots(merged_laz_file,
                                                         pp_neon_plots_site_inventory_path,
                                                         site,
                                                         year,
                                                         data_path)

# biomass
preprocessing_biomass(inventory_file_path,
                      pp_neon_plots_site_inventory_path,
                      site, year,
                      data_path)
