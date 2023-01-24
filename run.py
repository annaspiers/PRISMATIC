from preprocessing.inventory import download_veg_structure_data, \
                                    preprocessing_veg_structure_data
from preprocessing.plots import download_neon_polygons, \
                                preprocessing_neon_polygons, \
                                preprocessing_neon_polygons_site_inventory
from preprocessing.lidar import clip_laz_by_plots

site = 'SOAP'
year = 2019
data_path = '/home/toanngo/Documents/GitHub/prisma/preprocessing/data'
laz_path = ''
output_laz_path = ''

# download and preprocessing inventory data
download_veg_structure_data(site)
inventory_file_path = preprocessing_veg_structure_data(site, year)

# download NEON plots
neon_plots_path = download_neon_polygons(data_path)

# preprocessing NEON plots
pp_neon_plots_path = preprocessing_neon_polygons(neon_plots_path)

# preprocessing NEON plots given site ID and inventory
pp_neon_plots_site_inventory_path = preprocessing_neon_polygons_site_inventory(pp_neon_plots_path,
                                                                               site,
                                                                               inventory_file_path)

# clip lidar data
clip_laz_by_plots(laz_path,
                  pp_neon_plots_site_inventory_path,
                  output_laz_path)
