import logging
import hydra
import os 
os.environ['R_HOME'] = os.path.join(os.environ['CONDA_PREFIX'], 'lib/R') 

# # when adding new content to main.py
#     # import functions at top of main.py
#     # add functions to run in body of main.py
#     # track in utils.py
#     # track rerun in sites.yaml
#     # add to function_workflow.drawio and data structure diagram
#     # push to github

from utils.utils import build_cache_site, build_cache_all, force_rerun
from initialize.inventory import download_veg_structure_data, \
                                    download_trait_table, \
                                    preprocess_veg_structure_data
from initialize.plots import download_polygons, \
                                preprocess_polygons
from initialize.lidar import download_lidar, \
                                clip_lidar_by_plots, \
                                normalize_laz
from initialize.biomass import preprocess_biomass
from initialize.lad import preprocess_lad
from initialize.hyperspectral import download_hyperspectral, \
                                        generate_pft_reference, \
                                        create_tree_crown_polygons, \
                                        prep_aop_imagery, \
                                        extract_spectra_from_polygon, \
                                        correct_flightlines, \
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
    use_case    = cfg.others.use_case 
    ic_type     = cfg.others.ic_type 
    hs_type     = cfg.others.hs_type 
    month_window= cfg.others.month_window 
    n_plots     = cfg.others.n_plots 
    plot_length = cfg.others.plot_length 
    px_thresh   = cfg.others.px_thresh 
    ntree       = cfg.others.ntree 
    min_distance    = cfg.others.min_distance 
    use_tiles_w_veg = cfg.others.use_tiles_w_veg  
    randomMinSamples = cfg.others.randomMinSamples 
    aggregate_from_1m_to_2m_res = cfg.others.aggregate_from_1m_to_2m_res 
    independentValidationSet    = cfg.others.independentValidationSet
    pcaInsteadOfWavelengths     = cfg.others.pcaInsteadOfWavelengths 
    balance_training_to_min_PFT = cfg.others.balance_training_to_min_PFT
    
    global_force_rerun = cfg.sites.global_run_params.force_rerun
    global_run = cfg.sites.global_run_params.run
    if isinstance(global_run, str):
        global_run = [global_run]

    # First download/process raw data for all site-years to get stacked AOP product
    for site, v in cfg.sites.run.items():
        print(cfg.sites.run.items())
        if not global_run or site in global_run:
            for year_inventory, p in v.items():
                rerun_status = {}
                for k, v in p.force_rerun.items():
                    rerun_status[k] = v or global_force_rerun.get(k, False)
                year_aop = p.year_aop
                log.info(f'Download raw data for site: {site}, '
                         f'year: {year_inventory}, year_aop: {year_aop}, '
                         f'with rerun status: {rerun_status}')
                cache = build_cache_site(site=site, 
                                    year_inventory=year_inventory, 
                                    year_aop=year_aop, 
                                    data_raw_aop_path=data_raw_aop_path, 
                                    data_raw_inv_path=data_raw_inv_path, 
                                    data_int_path=data_int_path, 
                                    data_final_path=data_final_path, 
                                    use_case=use_case, 
                                    ic_type=ic_type,
                                    hs_type=hs_type)
                
                # download lidar
                laz_path, tif_path = (force_rerun(cache, force=rerun_status)
                                      (download_lidar)
                                      (site=site,
                                        year=year_aop,
                                        lidar_path=data_raw_aop_path,
                                        use_tiles_w_veg=use_tiles_w_veg))
                
                # download hs data
                hs_path = (force_rerun(cache, force=rerun_status)
                                        (download_hyperspectral)
                                        (site=site,
                                         year=year_aop,
                                         data_raw_aop_path=data_raw_aop_path,
                                         hs_type=hs_type))
                
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
                                
                _, _ = (force_rerun(cache, force=rerun_status)
                        (download_veg_structure_data)
                        (site=site, 
                         data_path=data_raw_inv_path))
                
                neon_plots_path = (force_rerun(cache, force=rerun_status)
                                   (download_polygons)
                                   (data_path=data_raw_inv_path))



    sites = []
    for site, v in cfg.sites.run.items():
        sites.append(site)
        pft_reference_path = (force_rerun(cache, force=rerun_status)
                                            (generate_pft_reference)
                                            (sites=sites,
                                            data_raw_inv_path=data_raw_inv_path, 
                                            data_int_path=data_int_path,
                                            trait_table_path=trait_table_path))



    # Process intermediate data per site
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
                    cache = build_cache_site(site=site, 
                                    year_inventory=year_inventory, 
                                    year_aop=year_aop, 
                                    data_raw_aop_path=data_raw_aop_path, 
                                    data_raw_inv_path=data_raw_inv_path, 
                                    data_int_path=data_int_path, 
                                    data_final_path=data_final_path, 
                                    use_case=use_case, 
                                    ic_type=ic_type,
                                    hs_type = hs_type )
                    
    #                 # Reset variable names to proper site/year
    #                 # ais is there a cleaner way to do this? ask sy-toan
    #                 laz_path=os.path.join(data_raw_aop_path,site,year_aop,"laz")
    #                 if hs_type=="tile":
    #                     hs_path=os.path.join(data_raw_aop_path,site,year_aop,"hs_tile_h5")
    #                 else:
    #                     hs_path=os.path.join(data_raw_aop_path,site,year_aop,"hs_flightline_h5")
    #                 tif_path=os.path.join(data_raw_aop_path,site,year_aop,"tif")
    #                 neon_plots_path=os.path.join(data_raw_inv_path,"All_NEON_TOS_Plots_V9")
                    
    #                 # process plots   
    #             inventory_file_path, \
    #                 sampling_effort_path = (force_rerun(cache, force=rerun_status)
    #                                         (preprocess_veg_structure_data)
    #                                         (site=site,
    #                                          year_inv=year_inventory,
    #                                          year_aop=year_aop,
    #                                          data_path=data_raw_inv_path,
    #                                          month_window=month_window))
                
    #             partitioned_plots_path = (force_rerun(cache, force=rerun_status)
    #                                       (preprocess_polygons)
    #                                       (input_data_path=neon_plots_path,
    #                                        sampling_effort_path=sampling_effort_path,
    #                                        inventory_path=inventory_file_path,
    #                                        site=site,
    #                                        year=year_inventory,
    #                                        output_data_path=data_int_path))

    #             # clip lidar data
    #             normalized_laz_path = (force_rerun(cache, force=rerun_status)
    #                                    (normalize_laz)
    #                                    (laz_path=laz_path,
    #                                     site=site,
    #                                     year=year_inventory,
    #                                     output_path=data_int_path))
                
    #             # ais walk through workflow below with a fine-toothed comb 
    #             clipped_laz_path = (force_rerun(cache, force=rerun_status)
    #                                 (clip_lidar_by_plots)
    #                                 (laz_path=normalized_laz_path,
    #                                  tif_path=tif_path,
    #                                  site_plots_path=partitioned_plots_path,
    #                                  site=site,
    #                                  year=year_inventory,
    #                                  output_laz_path=data_int_path,
    #                                  end_result=True))
    #             # leaf area density
    #             (force_rerun(cache, force=rerun_status)
    #                      (preprocess_lad)
    #                      (laz_path=clipped_laz_path,
    #                       inventory_path=inventory_file_path,
    #                       site=site,
    #                       year=year_inventory,
    #                       output_data_path=data_int_path,
    #                       use_case=use_case,
    #                       end_result=False))
    #             # biomass
    #             biomass_path = (force_rerun(cache, force=rerun_status)
    #                             (preprocess_biomass)
    #                             (data_path=inventory_file_path,
    #                              site_plots_path=partitioned_plots_path,
    #                              sampling_effort_path=sampling_effort_path,
    #                              site=site,
    #                              year=year_inventory,
    #                              output_data_path=data_int_path,
    #                              neon_trait_table_path=trait_table_path,
    #                              end_result=False)) 
    #             # ais check warnings for utils.allometry - lots of species not detected in neon_trait_table
                
    #             if hs_type=="flightline":
    #                 (force_rerun(cache, force=rerun_status)
    #                      (correct_flightlines)
    #                      (site=site,
    #                      year_inv=year_inventory, 
    #                      year_aop=year_aop, 
    #                      data_raw_aop_path=data_raw_aop_path,
    #                      data_int_path=data_int_path))
                
                
    #             training_crown_shp_path = (force_rerun(cache, 
    #                                                      force=rerun_status)
    #                                          (create_tree_crown_polygons)
    #                                          (site=site,
    #                                           year=year_inventory,
    #                                           data_raw_inv_path=data_raw_inv_path, 
    #                                           data_int_path=data_int_path, 
    #                                           biomass_path=biomass_path,
    #                                           pft_reference_path=pft_reference_path,
    #                                           px_thresh=px_thresh))  
                
    #             # prep NEON AOP data for classifier
    #             stacked_aop_path = (force_rerun(cache, force=rerun_status)
    #                                 (prep_aop_imagery)
    #                                 (site=site,
    #                                  year=year_inventory,
    #                                  hs_type=hs_type,
    #                                  hs_path=hs_path,
    #                                  tif_path=tif_path,
    #                                  data_int_path=data_int_path,
    #                                 use_tiles_w_veg=use_tiles_w_veg))
                
    #             training_spectra_csv_path = (force_rerun(cache, 
    #                                    force=rerun_status)
    #                                    (extract_spectra_from_polygon)
    #                                    (site=site,
    #                                     year=year_inventory,
    #                                     shp_path=training_crown_shp_path,
    #                                     data_int_path=data_int_path,
    #                                     data_final_path=data_final_path,
    #                                     stacked_aop_path=stacked_aop_path,
    #                                     use_case=use_case,
    #                                     ic_type=ic_type,
    #                                     aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res)) 

    # # ^ intermediate data finished processing for all site/years. Next, train RF 

    for k, v in p.force_rerun.items():
        rerun_status[k] = v or global_force_rerun.get(k, False)
    # log.info(f'Run process for all sites, '
    #     f'with rerun status: {rerun_status}')
        #ais ask sy-toan how to code this to train RF only once,, and not for every function in sites.yaml (stepped through k)
    
        from pathlib import Path #ais temporary
        trait_table_path=Path("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/NEON_trait_table.csv") #ais temporary 
        
        cache = build_cache_all(
                            trait_table_path=trait_table_path,
                            data_raw_inv_path=data_raw_inv_path, 
                            data_int_path=data_int_path, 
                            data_final_path=data_final_path, 
                            use_case=use_case, 
                            ic_type=ic_type)                
            
        # Train model
        rf_model_path = (force_rerun(cache,  force=rerun_status)
                        (train_pft_classifier)
                        (sites=sites,
                        data_int_path=data_int_path,
                        pcaInsteadOfWavelengths=pcaInsteadOfWavelengths, 
                        ntree=ntree, 
                        randomMinSamples=randomMinSamples, 
                        independentValidationSet=independentValidationSet,
                        balance_training_to_min_PFT=balance_training_to_min_PFT)) 
        # ais why does this run 5 times, for each site



        # INITIAL CONDITIONS
        # ais set this up so that I can specify site/year or generate initial conditions for ALL sites
        rf_model_path="/pscratch/sd/a/aspiers/data/intermediate/rf_dir/rf_model_tree_crowns_training.joblib"
        site='SJER' #ais how to specify this automatically? - when I specify 'usecase=predict'
        year_inventory='2021'  #ais how to specify this automatically?
        year_aop = '2021-03'  #ais how to specify this automatically?
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
                                                stacked_aop_path = os.path.join(data_int_path,site,year_inventory,'stacked_aop'),
                                                biomass_path = os.path.join(data_int_path,site,year_inventory,"biomass/pp_veg_structure_IND_IBA_IAGB_live.csv") ,
                                                use_case=use_case,
                                                ic_type=ic_type,
                                                n_plots = n_plots,
                                                min_distance = min_distance, 
                                                plot_length = plot_length, 
                                                aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res, 
                                                pcaInsteadOfWavelengths=pcaInsteadOfWavelengths))
                    # ais do I actually need to specify year_inv and year_aop?
            
    log.info('DONE')

if __name__ == '__main__':
    main()