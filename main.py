import logging
import hydra
import os 
os.environ['R_HOME'] = os.path.join(os.environ['CONDA_PREFIX'], 'lib/R') 
# so that this line works in inventory.py: import rpy2.robjects as robjects
# solution from https://github.com/rpy2/rpy2/issues/882
# generalized with help from https://stackoverflow.com/questions/36539623/how-do-i-
# find-the-name-of-the-conda-environment-in-which-my-code-is-running

from utils.utils import build_cache, force_rerun
from initialize.inventory import download_veg_structure_data, \
                                    preprocess_veg_structure_data, \
                                    download_trait_table
from initialize.plots import download_polygons, \
                                preprocess_polygons
from initialize.lidar import clip_lidar_by_plots, \
                                download_lidar, \
                                normalize_laz
from initialize.biomass import preprocess_biomass
from initialize.lad import preprocess_lad
from initialize.hyperspectral import download_hs_L3_tiles, \
                                        prep_aop_imagery, \
                                        create_training_data, \
                                        train_pft_classifier
from initialize.generate_initial_conditions import generate_initial_conditions

log = logging.getLogger(__name__)

@hydra.main(version_base='1.2', config_path='conf', config_name='config')
def main(cfg):
    log.info(f'Run with configuration: {cfg}')
    data_raw_aop_path = cfg.paths.data_raw_aop_path
    data_raw_inv_path = cfg.paths.data_raw_inv_path
    data_int_path = cfg.paths.data_int_path
    data_final_path = cfg.paths.data_final_path

    # Parameters to change manually
    use_case = "train" # train, predict
    ic_type = "rs_inv_plots" # field_inv_plots, rs_inv_plots, rs_random_plots
    aggregate_from_1m_to_2m_res = False
    px_thresh = 2
    ntree = 5000
    randomMinSamples = False
    independentValidationSet = True
    pcaInsteadOfWavelengths = True
    n_plots = 1000
    min_distance = 20
    plot_length = 20

    global_force_rerun = cfg.sites.global_run_params.force_rerun
    global_run = cfg.sites.global_run_params.run
    if isinstance(global_run, str):
        global_run = [global_run]

    # First download/process raw data for all site-years to get stacked AOP product
    for site, v in cfg.sites.run.items():
        if not global_run or site in global_run:
            for year_inventory, p in v.items():
                rerun_status = {}
                for k, v in p.force_rerun.items():
                    rerun_status[k] = v or global_force_rerun.get(k, False)
                year_aop = p.year_aop
                log.info(f'Download raw data for site: {site}, '
                         f'year: {year_inventory}, year_aop: {year_aop}, '
                         f'with rerun status: {rerun_status}')
                cache = build_cache(site=site, 
                                    year_inventory=year_inventory, 
                                    year_aop=year_aop, 
                                    data_raw_aop_path=data_raw_aop_path, 
                                    data_raw_inv_path=data_raw_inv_path, 
                                    data_int_path=data_int_path, 
                                    data_final_path=data_final_path, 
                                    use_case=use_case, 
                                    ic_type=ic_type)
                
                # download neon_trait_table
                url = cfg.others.neon_trait_table.neon_trait_link
                trait_table_path = (force_rerun(cache, force={
                                                    'download_trait_table':
                                                    (cfg
                                                     .others
                                                     .neon_trait_table
                                                     .force_rerun)})
                                                     (download_trait_table)
                                                     (download_link=url, 
                                                      data_path=data_raw_inv_path))

                
                
                # download lidar
                laz_path, tif_path = (force_rerun(cache, force=rerun_status)
                                      (download_lidar)
                                      (site=site,
                                       year=year_aop,
                                       lidar_path=data_raw_aop_path))
                
                # download hs data
                hs_L3_path, tif_path = (force_rerun(cache, force=rerun_status)
                                        (download_hs_L3_tiles)
                                        (site=site,
                                         year=year_aop,
                                         data_raw_aop_path=data_raw_aop_path))
                
                _, _ = (force_rerun(cache, force=rerun_status)
                        (download_veg_structure_data) #ais how to automate this to seloect US(OR)?
                        (site=site, 
                         data_path=data_raw_inv_path))
                
                neon_plots_path = (force_rerun(cache, force=rerun_status)
                                   (download_polygons)
                                   (data_path=data_raw_inv_path))
                
        # Process intermediate and final data
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
                    cache = build_cache(site=site, 
                                    year_inventory=year_inventory, 
                                    year_aop=year_aop, 
                                    data_raw_aop_path=data_raw_aop_path, 
                                    data_raw_inv_path=data_raw_inv_path, 
                                    data_int_path=data_int_path, 
                                    data_final_path=data_final_path, 
                                    use_case=use_case, 
                                    ic_type=ic_type)
                    
                    # process plots   
                inventory_file_path, \
                    sampling_effort_path = (force_rerun(cache, force=rerun_status)
                                            (preprocess_veg_structure_data)
                                            (site=site,
                                             year_inv=year_inventory,
                                             year_aop=year_aop,
                                             data_path=data_raw_inv_path))
                
                partitioned_plots_path = (force_rerun(cache, force=rerun_status)
                                          (preprocess_polygons)
                                          (input_data_path=neon_plots_path,
                                           sampling_effort_path=sampling_effort_path,
                                           inventory_path=inventory_file_path,
                                           site=site,
                                           year=year_inventory,
                                           output_data_path=data_int_path))

                # clip lidar data
                normalized_laz_path = (force_rerun(cache, force=rerun_status)
                                       (normalize_laz)
                                       (laz_path=laz_path,
                                        site=site,
                                        year=year_inventory,
                                        output_path=data_int_path))
                clipped_laz_path = (force_rerun(cache, force=rerun_status)
                                    (clip_lidar_by_plots)
                                    (laz_path=normalized_laz_path,
                                     tif_path=tif_path,
                                     site_plots_path=partitioned_plots_path,
                                     site=site,
                                     year=year_inventory,
                                     output_laz_path=data_int_path,
                                     end_result=True))
                
                # prep NEON AOP data for classifier
                stacked_aop_path = (force_rerun(cache, force=rerun_status)
                                    (prep_aop_imagery)
                                    (site=site,
                                     year=year_inventory,
                                     hs_L3_path=hs_L3_path, 
                                     tif_path=tif_path,
                                     data_int_path=data_int_path))
                    
                # leaf area density
                (force_rerun(cache, force=rerun_status)
                         (preprocess_lad)
                         (laz_path=clipped_laz_path,
                          inventory_path=inventory_file_path,
                          site=site,
                          year=year_inventory,
                          output_data_path=data_int_path,
                          use_case=use_case,
                          end_result=False))

                # biomass
                biomass_path = (force_rerun(cache, force=rerun_status)
                                (preprocess_biomass)
                                (data_path=inventory_file_path,
                                 site_plots_path=partitioned_plots_path,
                                 sampling_effort_path=sampling_effort_path,
                                 site=site,
                                 year=year_inventory,
                                 output_data_path=data_int_path,
                                 neon_trait_table_path=trait_table_path,
                                 end_result=False)) 
                #ais check warnings for utils.allometry - lots of species not detected in neon_trait_table
                
                training_crown_shp_path, \
                training_spectra_csv_path = (force_rerun(cache, 
                                                         force=rerun_status)
                                             (create_training_data)
                                             (site=site,
                                              year=year_inventory,
                                              biomass_path=biomass_path,
                                              data_int_path=data_int_path,
                                              data_final_path=data_final_path,
                                              stacked_aop_path=stacked_aop_path,
                                              px_thresh=px_thresh,
                                              use_case="train",
                                              ic_type=ic_type,
                                              aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res))            
                
                # Train model
                rf_model_path = (force_rerun(cache,  force=rerun_status)
                                     (train_pft_classifier)
                                     (site=site,
                                      year=year_inventory,
                                      stacked_aop_path=stacked_aop_path,
                                      training_shp_path=training_crown_shp_path,  
                                      training_spectra_path=training_spectra_csv_path, 
                                      data_int_path=data_int_path,
                                      pcaInsteadOfWavelengths=pcaInsteadOfWavelengths, 
                                      ntree=ntree, 
                                      randomMinSamples=randomMinSamples, 
                                      independentValidationSet=independentValidationSet)) 
                    # ais ^ specify in yaml whether using corrected/uncorrected hs data

                    # INITIAL CONDITIONS
                if use_case=="predict":
                # Generate initial conditions
                    cohort_path, patch_path = (force_rerun(cache,  force=rerun_status)
                                                 (generate_initial_conditions)
                                                 (site=site,
                                                  year_inv=year_inventory,
                                                  year_aop = year_aop,
                                                  data_raw_aop_path = data_raw_aop_path,
                                                  data_int_path=data_int_path, 
                                                  data_final_path=data_final_path,
                                                  rf_model_path = rf_model_path, 
                                                  stacked_aop_path = stacked_aop_path,
                                                  biomass_path = biomass_path,
                                                  use_case=use_case,
                                                  ic_type=ic_type,
                                                  n_plots = n_plots,
                                                  min_distance = min_distance, 
                                                  plot_length = plot_length, 
                                                  aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res, 
                                                  pcaInsteadOfWavelengths=pcaInsteadOfWavelengths))

    # when adding new content to main.py
        # import functions at top of main.py
        # add functions to run in body of main.py
        # track in utils.py
        # track rerun in sites.yaml
        # add to function_workflow.drawio and data structure diagram
        # push to github
            
    log.info('DONE')

if __name__ == '__main__':
    main()

# PROJECT FOLLOWUPS
# What function generates the individual plot shp files? We know where plots.shp comes from
# E.g., /Users/AISpiers/dev/RS-PRISMATIC/inititalize/data/intermediate/SOAP/2019/inventory_plots/SOAP_003_central.dbf
# Solution: generated in clip_lidar_by_plots()
# Don’t save intermediate files saved in “/Users/AISpiers/dev/RS-PRISMATIC/inititalize/data/intermediate/SOAP/2019/clipped_lidar_tiles”
# Should the y-axis be biomass instead of number of individuals? /Users/AISpiers/dev/RS-PRISMATIC/inititalize/data/intermediate/diagnostics/SOAP/2019/biomass_plots/SOAP.png
# Calculate_ila is a dummy function, returning 0.1
# Diagnostics lad profiles are plotted with legend cut off
# The function download_urls in utils/download_functions is not called anywhere
# No function is defined in utils/rank.py - remove it?
# Do I ever use “vst_training.csv” - could delete where I save it in hyperspectral_helper.R
# Where to save and use pft_list?
# pft_list <-  c("pine_PFT","cedar_PFT","fir_PFT","oak_PFT","shrub_PFT")
