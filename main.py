import logging
import hydra

from preprocessing.inventory import download_veg_structure_data, \
                                    preprocess_veg_structure_data
from preprocessing.plots import download_polygons, \
                                preprocess_polygons
from preprocessing.lidar import clip_laz_by_plots, \
                                normalize_laz
from preprocessing.biomass import preprocess_biomass

log = logging.getLogger(__name__)


@hydra.main(version_base='1.2', config_path='conf', config_name='config')
def main(cfg):
    log.info(f'Run with configuration: {cfg}')

    site = cfg.sites.site
    year = cfg.sites.year
    lidar_path = cfg.paths.lidar_path
    data_path = cfg.paths.data_path

    # process inventory
    download_veg_structure_data(site, data_path)
    inventory_file_path, \
        sampling_effort_path = preprocess_veg_structure_data(site,
                                                             year,
                                                             data_path)

    # process plots
    neon_plots_path = download_polygons(data_path)
    partitioned_plots_path = \
        preprocess_polygons(neon_plots_path,
                            sampling_effort_path,
                            inventory_file_path,
                            site,
                            year,
                            data_path)

    # clip lidar data
    normalized_laz_path = normalize_laz(lidar_path,
                                        site,
                                        year,
                                        data_path,
                                        end_result=False)
    clip_laz_by_plots(normalized_laz_path,
                      partitioned_plots_path,
                      site,
                      year,
                      data_path,
                      end_result=True)

    # biomass
    preprocess_biomass(inventory_file_path,
                       partitioned_plots_path,
                       sampling_effort_path,
                       site,
                       year,
                       data_path,
                       end_result=True)

    log.info('DONE')


if __name__ == '__main__':
    main()
