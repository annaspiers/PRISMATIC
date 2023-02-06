import logging
import hydra

from preprocessing.inventory import download_veg_structure_data, \
                                    preprocessing_veg_structure_data
from preprocessing.plots import download_neon_polygons, \
                                preprocessing_neon_polygons, \
                                preprocessing_neon_polygons_site_inventory
from preprocessing.lidar import clip_laz_by_plots, \
                                clip_laz_by_inventory_plots
from preprocessing.biomass import preprocessing_biomass

log = logging.getLogger(__name__)


@hydra.main(version_base='1.2', config_path='conf', config_name='config')
def main(cfg):
    log.info(f'Run with configuration: {cfg}')

    site = cfg.sites.site
    year = cfg.sites.year
    lidar_path = cfg.paths.lidar_path
    output_data_path = cfg.paths.data_path

    download_veg_structure_data(site)
    inventory_file_path = preprocessing_veg_structure_data(site,
                                                           year,
                                                           output_data_path)

    neon_plots_path = download_neon_polygons(output_data_path)
    pp_neon_plots_path = preprocessing_neon_polygons(neon_plots_path,
                                                     output_data_path)

    pp_neon_plots_site_inventory_path = \
        preprocessing_neon_polygons_site_inventory(pp_neon_plots_path,
                                                   site,
                                                   year,
                                                   output_data_path,
                                                   inventory_file_path)

    merged_lidar_file = clip_laz_by_plots(lidar_path,
                                          pp_neon_plots_site_inventory_path,
                                          site,
                                          year,
                                          output_data_path)

    clipped_laz_inventory_path = \
        clip_laz_by_inventory_plots(merged_lidar_file,
                                    pp_neon_plots_site_inventory_path,
                                    site,
                                    year,
                                    output_data_path)


    preprocessing_biomass(inventory_file_path,
                          pp_neon_plots_site_inventory_path,
                          site,
                          year,
                          output_data_path)


if __name__ == '__main__':
    main()
