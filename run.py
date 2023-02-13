from preprocessing.inventory import download_veg_structure_data, \
                                    preprocess_veg_structure_data
from preprocessing.plots import download_polygons, \
                                preprocess_polygons
from preprocessing.lidar import clip_laz_by_plots, \
                                normalize_laz
from preprocessing.biomass import preprocess_biomass

from constants import site, year, data_path, laz_path


# download and preprocessing inventory data
download_veg_structure_data(site)
inventory_file_path = preprocess_veg_structure_data(site, year, data_path)

# download NEON plots/watch
neon_plots_path = download_polygons(data_path)
# preprocessing NEON plots
sampling_effort_path = f'{data_path}/{site}/plot_sampling_effort.csv'
# partition plots to subplots
partitioned_plots_path = \
    preprocess_polygons(neon_plots_path,
                        sampling_effort_path,
                        inventory_file_path,
                        site,
                        year,
                        data_path)

# clip lidar data
normalized_laz_path = normalize_laz(laz_path,
                                    site,
                                    year,
                                    data_path,
                                    end_result=False)
output_laz_path = clip_laz_by_plots(normalized_laz_path,
                                    partitioned_plots_path,
                                    site,
                                    year,
                                    data_path,
                                    end_result=True)

# biomass
preprocess_biomass(inventory_file_path,
                   partitioned_plots_path,
                   sampling_effort_path,
                   site, year,
                   data_path,
                   end_result=True)
