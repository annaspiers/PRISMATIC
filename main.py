import logging
import hydra

from preprocessing.inventory import download_veg_structure_data, \
                                    preprocess_veg_structure_data, \
                                    download_trait_table
from preprocessing.plots import download_polygons, \
                                preprocess_polygons
from preprocessing.lidar import clip_lidar_by_plots, \
                                download_lidar, \
                                normalize_laz
from preprocessing.biomass import preprocess_biomass
from utils.utils import build_cache, force_rerun

log = logging.getLogger(__name__)


@hydra.main(version_base='1.2', config_path='conf', config_name='config')
def main(cfg):
    log.info(f'Run with configuration: {cfg}')
    root_lidar_path = cfg.paths.lidar_path
    data_path = cfg.paths.data_path

    global_force_rerun = cfg.sites.global_run_params.force_rerun
    global_run = cfg.sites.global_run_params.run
    if isinstance(global_run, str):
        global_run = [global_run]

    for site, v in cfg.sites.run.items():
        if not global_run or site in global_run:
            for year_inventory, p in v.items():
                rerun_status = {}
                for k, v in p.force_rerun.items():
                    rerun_status[k] = v or global_force_rerun.get(k, False)
                year_lidar = p.year_lidar
                log.info(f'Run process for site: {site}, '
                         f'year: {year_inventory}, year_lidar: {year_lidar}, '
                         f'with rerun status: {rerun_status}')
                cache = build_cache(site, year_inventory, year_lidar,
                                    data_path, root_lidar_path)

                # download neon_trait_table
                url = cfg.others.neon_trait_table.neon_trait_link
                trait_table_path = (force_rerun(cache,
                                                force={
                                                    'download_trait_table':
                                                    (cfg
                                                     .others
                                                     .neon_trait_table
                                                     .force_rerun)})
                                    (download_trait_table)
                                    (url, data_path))

                # download lidar
                laz_path, tif_path = (force_rerun(cache,
                                                  force=rerun_status)
                                      (download_lidar)
                                      (site,
                                       year_lidar,
                                       root_lidar_path))

                # process inventory
                _, _ = (force_rerun(cache,
                                    force=rerun_status)
                        (download_veg_structure_data)
                        (site, data_path))
                inventory_file_path, \
                    sampling_effort_path = (force_rerun(cache,
                                                        force=rerun_status)
                                            (preprocess_veg_structure_data)
                                            (site,
                                             year_inventory,
                                             data_path))

                # process plots
                neon_plots_path = (force_rerun(cache,
                                               force=rerun_status)
                                   (download_polygons)
                                   (data_path))
                partitioned_plots_path = \
                    (force_rerun(cache,
                                 force=rerun_status)
                     (preprocess_polygons)
                     (neon_plots_path,
                      sampling_effort_path,
                      inventory_file_path,
                      site,
                      year_inventory,
                      data_path))

                # clip lidar data
                normalized_laz_path = (force_rerun(cache,
                                                   force=rerun_status)
                                       (normalize_laz)(laz_path,
                                                       site,
                                                       year_inventory,
                                                       data_path))
                (force_rerun(cache,
                             force=rerun_status)
                 (clip_lidar_by_plots)
                 (normalized_laz_path,
                  tif_path,
                  partitioned_plots_path,
                  site,
                  year_inventory,
                  data_path,
                  end_result=True))

                # biomass
                (force_rerun(cache,
                             force=rerun_status)
                 (preprocess_biomass)(inventory_file_path,
                                      partitioned_plots_path,
                                      sampling_effort_path,
                                      site,
                                      year_inventory,
                                      data_path,
                                      neon_trait_table_path=trait_table_path,
                                      end_result=True))

    log.info('DONE')


if __name__ == '__main__':
    main()
