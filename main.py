import logging
import hydra
import os 
from pathlib import Path
os.environ['R_HOME'] = os.path.join(os.environ['CONDA_PREFIX'], 'lib/R') 

# # when adding new content to main.py
#     # import functions at top of main.py
#     # add functions to run in body of main.py
#     # track in utils.py
#     # track rerun in sites.yaml
#     # add to function_workflow.drawio and data structure diagram
#     # push to github

from utils.utils import build_cache_site, force_rerun 
from initialize.inventory import download_veg_structure_data, \
                                    download_trait_table, \
                                    prep_veg_structure
from initialize.plots import download_polygons, \
                                prep_polygons
from initialize.lidar import download_aop_bbox, \
                                download_lidar, \
                                clip_lidar_by_plots, \
                                normalize_laz
from initialize.biomass import prep_biomass
from initialize.lad import prep_lad
from initialize.hyperspectral import download_hyperspectral, \
                                        prep_manual_training_data, \
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
    data_int_path     = cfg.paths.data_int_path
    data_final_path   = cfg.paths.data_final_path

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
    multisite   = cfg.others.multisite
    coords_bbox = cfg.others.coords_bbox 
    # balance_training_to_min_PFT = cfg.others.balance_training_to_min_PFT
    
    global_force_rerun = cfg.sites.global_run_params.force_rerun
    global_run = cfg.sites.global_run_params.run
    if isinstance(global_run, str):
        global_run = [global_run]

    if use_case=="train": 
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
                                        hs_type=hs_type,
                                        coords_bbox=coords_bbox)
                    
                    if not coords_bbox:
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
                    
                    else:
                        # download lidar and hyperspectral data 
                        hs_path, laz_path, tif_path  = (force_rerun(cache, force=rerun_status)
                                                (download_aop_bbox)
                                                (site=site,
                                                year=year_aop,
                                                path=data_raw_aop_path,
                                                hs_type=hs_type,
                                                coords_bbox=coords_bbox))
                                                        
                    _, _ = (force_rerun(cache, force=rerun_status)
                            (download_veg_structure_data)
                            (site=site, 
                            data_path=data_raw_inv_path))
                    
                    # process plots   
                    inventory_file_path, \
                        sampling_effort_path = (force_rerun(cache, force=rerun_status)
                                                (prep_veg_structure)
                                                (site=site,
                                                year_inv=year_inventory,
                                                year_aop=year_aop,
                                                data_path=data_raw_inv_path,
                                                month_window=month_window))
                    
                    neon_plots_path = (force_rerun(cache, force=rerun_status)
                                    (download_polygons)
                                    (data_path=data_raw_inv_path))

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
                                        hs_type = hs_type,
                                         coords_bbox=coords_bbox )
                        
                        # Reset variable names to proper site/year
                        # ais is there a cleaner way to do this? ask sy-toan
                        laz_path=os.path.join(data_raw_aop_path,site,year_aop,"laz")
                        if hs_type=="tile":
                            hs_path=os.path.join(data_raw_aop_path,site,year_aop,"hs_tile")
                        else:
                            hs_path=os.path.join(data_raw_aop_path,site,year_aop,"hs_flightline")
                        tif_path=os.path.join(data_raw_aop_path,site,year_aop,"tif")
                        neon_plots_path=os.path.join(data_raw_inv_path,"All_NEON_TOS_Plots_V9")
                       
                                           
                    partitioned_plots_path = (force_rerun(cache, force=rerun_status)
                                            (prep_polygons)
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
                    
                    # ais walk through workflow below with a fine-toothed comb 
                    # ais is this for both laz and tif files? or just laz?
                    # if for both laz and tif, then the output, clipped_laz_path, should be used more than just 1 time in following function...
                    clipped_laz_path = (force_rerun(cache, force=rerun_status)
                                        (clip_lidar_by_plots)
                                        (laz_path=normalized_laz_path,
                                        tif_path=tif_path,
                                        site_plots_path=partitioned_plots_path,
                                        site=site,
                                        year=year_inventory,
                                        output_laz_path=data_int_path,
                                        end_result=True))
                    
                    # leaf area density
                    (force_rerun(cache, force=rerun_status)
                            (prep_lad)
                            (laz_path=clipped_laz_path,
                            inventory_path=inventory_file_path,
                            site=site,
                            year=year_inventory,
                            output_path=data_int_path,
                            use_case=use_case)) #ais or should I just use "train" here - whats best practice

                    # biomass
                    biomass_path = (force_rerun(cache, force=rerun_status)
                                    (prep_biomass)
                                    (data_path=inventory_file_path,
                                    site_plots_path=partitioned_plots_path,
                                    sampling_effort_path=sampling_effort_path,
                                    site=site,
                                    year=year_inventory,
                                    data_int_path=data_int_path,
                                    neon_trait_table_path=trait_table_path,
                                    end_result=False)) 
                    # ais check warnings for utils.allometry - lots of species not detected in neon_trait_table
                    
                    if ic_type!="field_inv_plots": #only run the following functions if RS data is being used
                        if hs_type=="flightline":
                            (force_rerun(cache, force=rerun_status)
                                (correct_flightlines)
                                (site=site,
                                year_inv=year_inventory, 
                                year_aop=year_aop, 
                                data_raw_aop_path=data_raw_aop_path,
                                data_int_path=data_int_path))
                        
                        
                        training_crown_shp_path = (force_rerun(cache, 
                                                                force=rerun_status)
                                                    (prep_manual_training_data)
                                                    (site=site,
                                                    year=year_inventory,
                                                    data_raw_inv_path=data_raw_inv_path, 
                                                    data_int_path=data_int_path, 
                                                    biomass_path=biomass_path))  
                        
                        # prep NEON AOP data for classifier
                        stacked_aop_path = (force_rerun(cache, force=rerun_status)
                                            (prep_aop_imagery)
                                            (site=site,
                                            year=year_inventory,
                                            hs_type=hs_type,
                                            hs_path=hs_path,
                                            tif_path=tif_path,
                                            data_int_path=data_int_path,
                                            use_tiles_w_veg=use_tiles_w_veg))
                        
                        training_spectra_csv_path = (force_rerun(cache, 
                                            force=rerun_status)
                                            (extract_spectra_from_polygon)
                                            (site=site,
                                                year=year_inventory,
                                                shp_path=training_crown_shp_path,
                                                data_int_path=data_int_path,
                                                data_final_path=data_final_path,
                                                stacked_aop_path=stacked_aop_path,
                                                use_case=use_case,
                                                aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res,
                                                ic_type=ic_type)) 
        # ^ intermediate data finished processing for all site/years. Next, train RF 

    
        sites = []
        for site, v in cfg.sites.run.items():
            sites.append(site)
            
        for k, v in p.force_rerun.items():
            rerun_status[k] = v or global_force_rerun.get(k, False)
            log.info(f'Run process for all sites, '
            f'with rerun status: {rerun_status}')
            
        cache = build_cache_site(site=site, 
                                        year_inventory=year_inventory, 
                                        year_aop=year_aop, 
                                        data_raw_aop_path=data_raw_aop_path, 
                                        data_raw_inv_path=data_raw_inv_path, 
                                        data_int_path=data_int_path, 
                                        hs_type=hs_type,
                                        coords_bbox=coords_bbox)
                
        if ic_type!="field_inv_plots": #only run if using RS data
            rf_model_path = (force_rerun(cache,  force=rerun_status)
                            (train_pft_classifier)
                            (sites=sites,
                            data_int_path=data_int_path,
                            pcaInsteadOfWavelengths=pcaInsteadOfWavelengths, 
                            ntree=ntree, 
                            randomMinSamples=randomMinSamples, 
                            independentValidationSet=independentValidationSet)) 
        # ais why does this run 5 times, for each site



    if use_case=="predict": 
    
        #ais think about how to loop through ic_type - ic_type should be in build_cache_predict     - [ ] 
            #ais - fix code so ‘generate initial conditions’ steps through all sites, years, and ic_types’
        for site, v in cfg.sites.run.items():
            print(cfg.sites.run.items())
            if not global_run or site in global_run:
                for year_inventory, p in v.items():
                    rerun_status = {}
                    for k, v in p.force_rerun.items():
                        rerun_status[k] = v or global_force_rerun.get(k, False)
                    year_aop = p.year_aop
                    log.info(f'Download raw data for site: {site}, ' #ais change words
                            f'year: {year_inventory}, year_aop: {year_aop}, '
                            f'with rerun status: {rerun_status}')
                    cache = build_cache_site(site=site, 
                                        year_inventory=year_inventory, 
                                        year_aop=year_aop, 
                                        data_raw_aop_path=data_raw_aop_path, 
                                        data_raw_inv_path=data_raw_inv_path, 
                                        data_int_path=data_int_path, 
                                        hs_type=hs_type,
                                        coords_bbox=coords_bbox)
            # cache = build_cache_predict(site=site, #ais work on swapping this out for build_cache_site
            #                             year_inventory=year_inventory, 
            #                             use_case=use_case,
            #                             data_raw_aop_path=data_raw_aop_path, 
            #                             data_int_path=data_int_path, 
            #                             data_final_path=data_final_path,
            #                             ic_type=ic_type)
            # #ais ask sy-toan for help so I don't have to manually assign all these arguments

            # site='TEAK' #SOAP SJER
            # year_inventory='2021'
            # year_aop = '2021-07' #2021-07 2021-03

                    rf_model_path=os.path.join(data_int_path,'rf_dir/rf_model.joblib')
                    ic_type_path = str(data_final_path+'/'+site+'/'+year_inventory+'/'+ic_type) #needs to be str to be input
                    Path(ic_type_path).mkdir(parents=True, exist_ok=True)

                                                    
                    # Generate initial conditions
                    cohort_path, patch_path = (force_rerun(cache,  force=rerun_status)
                                            (generate_initial_conditions)
                                            (site=site,
                                                year_inv= year_inventory,
                                                year_aop = year_aop,
                                                data_raw_aop_path = data_raw_aop_path,
                                                data_int_path=data_int_path, 
                                                data_final_path=data_final_path,
                                                rf_model_path = rf_model_path, 
                                                stacked_aop_path = os.path.join(data_int_path,site,year_inventory,'stacked_aop'),
                                                biomass_path = os.path.join(data_int_path,site,year_inventory,"biomass/pp_veg_structure_IND_IBA_IAGB_live.csv"),
                                                use_case=use_case,
                                                ic_type=ic_type,
                                                ic_type_path=ic_type_path,
                                                n_plots = n_plots,
                                                min_distance = min_distance, 
                                                plot_length = plot_length, 
                                                aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res, 
                                                pcaInsteadOfWavelengths=pcaInsteadOfWavelengths,
                                                multisite=multisite))
            
    log.info('DONE')

if __name__ == '__main__':
    main()
    