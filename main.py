import logging
import hydra

import os 
os.environ['R_HOME'] = os.path.join(os.environ['CONDA_PREFIX'], 'lib/R') 
# so that this line works in inventory.py: import rpy2.robjects as robjects
# solution from https://github.com/rpy2/rpy2/issues/882
# generalized with help from https://stackoverflow.com/questions/36539623/how-do-i-find-the-name-of-the-conda-environment-in-which-my-code-is-running

from preprocessing.inventory import download_veg_structure_data, \
                                    preprocess_veg_structure_data, \
                                    download_trait_table
from preprocessing.plots import download_polygons, \
                                preprocess_polygons
from preprocessing.lidar import clip_lidar_by_plots, \
                                download_lidar, \
                                normalize_laz
from preprocessing.biomass import preprocess_biomass
from preprocessing.lad import preprocess_lad
from preprocessing.hyperspectral import download_hs_L3_tiles, \
                                        prep_aop_imagery, \
                                        create_training_data, \
                                        train_pft_classifier

from utils.utils import build_cache, force_rerun

log = logging.getLogger(__name__)


@hydra.main(version_base='1.2', config_path='conf', config_name='config')
def main(cfg):
    log.info(f'Run with configuration: {cfg}')
    # data_raw_path = cfg.paths.data_raw_path
    data_raw_aop_path = cfg.paths.data_raw_aop_path
    data_raw_inv_path = cfg.paths.data_raw_inv_path
    data_int_path = cfg.paths.data_int_path
    data_final_path = cfg.paths.data_final_path

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
                year_aop = p.year_aop
                log.info(f'Run process for site: {site}, '
                         f'year: {year_inventory}, year_aop: {year_aop}, '
                         f'with rerun status: {rerun_status}')
                cache = build_cache(site, year_inventory, year_aop, data_raw_aop_path, 
                                    data_raw_inv_path, data_int_path, data_final_path)
                
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
                                    (url, data_raw_inv_path))

                # process inventory
                _, _ = (force_rerun(cache,
                                    force=rerun_status)
                        (download_veg_structure_data) #ais how to automate this to enter 72?
                        (site, data_raw_inv_path))
                inventory_file_path, \
                    sampling_effort_path = (force_rerun(cache,
                                                        force=rerun_status)
                                            (preprocess_veg_structure_data)
                                            (site,
                                             year_inventory,
                                             data_raw_inv_path))
                neon_plots_path = (force_rerun(cache,
                                               force=rerun_status)
                                   (download_polygons)
                                   (data_raw_inv_path))
                
        # LIDAR
                # download lidar
                laz_path, tif_path = (force_rerun(cache,
                                                  force=rerun_status)
                                      (download_lidar)
                                      (site,
                                       year_aop,
                                       data_raw_aop_path))
                
                # process plots                
                partitioned_plots_path = \
                    (force_rerun(cache,
                                 force=rerun_status)
                     (preprocess_polygons)
                     (neon_plots_path,
                      sampling_effort_path,
                      inventory_file_path,
                      site,
                      year_inventory,
                      data_int_path))

                # clip lidar data
                normalized_laz_path = (force_rerun(cache,
                                                   force=rerun_status)
                                       (normalize_laz)(laz_path,
                                                       site,
                                                       year_inventory,
                                                       data_int_path))
                clipped_laz_path = (force_rerun(cache,
                                                force=rerun_status)
                                    (clip_lidar_by_plots)
                                    (normalized_laz_path,
                                     tif_path,
                                     partitioned_plots_path,
                                     site,
                                     year_inventory,
                                     data_int_path,
                                     end_result=True))

                # leaf area density
                (force_rerun(cache,
                             force=rerun_status)
                 (preprocess_lad)
                 (clipped_laz_path,
                  inventory_file_path,
                  site,
                  year_inventory,
                  data_int_path,
                  end_result=False))

                # biomass
                biomass_path = (force_rerun(cache,
                             force=rerun_status)
                                (preprocess_biomass)
                                (inventory_file_path,
                                      partitioned_plots_path,
                                      sampling_effort_path,
                                      site,
                                      year_inventory,
                                      data_int_path,
                                      neon_trait_table_path=trait_table_path,
                                      end_result=False)) 
                #ais check warnings for utils.allometry - lots of species not detected in neon_trait_table
                
        # HYPERSPECTRAL
                # download hs data
                hs_L3_path, tif_path = (force_rerun(cache,
                                                  force=rerun_status)
                                      (download_hs_L3_tiles)
                                      (site,
                                       year_aop,
                                       data_raw_aop_path))
                
                # prep NEON AOP data for classifier
                stacked_aop_path = (force_rerun(cache,
                                                  force=rerun_status)
                                      (prep_aop_imagery)
                                      (site,
                                       year_inventory,
                                       hs_L3_path, 
                                       tif_path,
                                       data_int_path))
                
                training_crown_shp_path, \
                    training_spectra_csv_path = (force_rerun(cache,
                                                  force=rerun_status)
                                      (create_training_data)
                                      (site,
                                       year_inventory,
                                       biomass_path,
                                       data_int_path,
                                       stacked_aop_path,
                                       px_thresh=2,
                                       use_case="train",
                                       aggregate_from_1m_to_2m_res=False)) 
                
                # Train model
                rf_model_path = (force_rerun(cache,
                                                  force=rerun_status)
                                      (train_pft_classifier)
                                      (site,
                                       year_inventory,
                                       stacked_aop_path,
                                       training_crown_shp_path, 
                                       training_spectra_csv_path, 
                                       data_int_path,
                                       pcaInsteadOfWavelengths=True, 
                                       ntree=5000, 
                                       randomMinSamples=False, 
                                       independentValidationSet=True)) # ais specify in yaml whether using corrected/uncorrected hs data

                


    # when adding new content to main.py
        # import functions at top of main.py
        # add functions to run in body of main,py
        # track in utils.py
        # track rerun in sites.yaml
        # add to function_workflow.drawio and data structure diagram
        # push to github
            
    log.info('DONE')

if __name__ == '__main__':
    main()
