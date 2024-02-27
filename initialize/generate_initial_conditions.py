import logging
import os
import sys
import glob

import rpy2.robjects as ro

from pathlib import Path
from preprocessing.lad import preprocess_lad

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = str(conda_env_path/'share'/'proj')

log = logging.getLogger(__name__)

def generate_initial_conditions(site, year_inv, year_aop, data_raw_aop_path, data_int_path,
                                data_final_path, rf_model_path, stacked_aop_path, biomass_path,
                                use_case, ic_type, n_plots, min_distance, plot_length, 
                                aggregate_from_1m_to_2m_res, pcaInsteadOfWavelengths):                                            
    """Generate FATES initial conditions from remote sensing data in three scenarios:
            1) inventory-based IC's over inventory plots (ic_type = field_inv_plots)
            2) remote-sensing-based IC's from plots used for forest inventory 
                    (ic_type = rs_inv_plots)
            3) remote-sensing-based IC's from plots randomly generated across 
                the AOP spatial extent (ic_type = rs_random_plots)

    We generate files formatted to initialize FATES: cohort (.css) and patch (.pss) files
    """

    if not os.path.exists(Path(data_final_path,site,year_inv)):
        os.makedirs(Path(data_final_path,site,year_inv))    
    
    log.info(f'Generating initial conditions for: {site} {year_inv}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'initial_conditions_helper.R'))
    r_source(str(Path(os.getcwd()+'/preprocessing/hyperspectral_helper.R')))
    r_source(str(Path(os.getcwd()+'/utils/hyperspectral_tools.R')))

    ic_type_path = str(data_final_path+'/'+site+'/'+year_inv+'/'+ic_type) #needs to be str to be input
    Path(ic_type_path).mkdir(parents=True, exist_ok=True)

    generate_cohort_patch_files = ro.r('generate_cohort_patch_files') 
    if (ic_type=="field_inv_plots"):
        cohort_path, patch_path = generate_cohort_patch_files(site=site, 
                                    year=year_inv, 
                                    data_int_path=data_int_path, 
                                    biomass_path=biomass_path,
                                    ic_type=ic_type, 
                                    ic_type_path=ic_type_path)

    if (ic_type=="rs_inv_plots" or ic_type=="rs_random_plots"):        
        ### 1) Specify spatial extent of initial conditions and prepare spatial data by
        # clipping point clouds to and extracting raster data from plots
        prepare_spatial_data_for_ic = ro.r('prepare_spatial_data_for_ic') 
        plots_shp_path, lad_laz_clipped_path, extracted_features_csv = prepare_spatial_data_for_ic(
                                            site, year_inv, year_aop, ic_type, data_raw_aop_path, 
                                            data_int_path, stacked_aop_path, data_final_path, 
                                            ic_type_path, use_case, n_plots, min_distance, 
                                            plot_length, aggregate_from_1m_to_2m_res) 
        log.info('Plots, clipped laz, and clipped stacked rasters prepared for generating initial conditions')

        ### 2) Predict pfts on new plots
        # proxy to neon-veg-SOAPpfts/11-predict-pft.R
        # if ic_type==rs_inv_plots, use known PFTs from NEON inventory data 
        # if ic_type==rs_random_plots, predict PFTs with RF classifier
        predict_plot_level_PFTs = ro.r('predict_plot_level_PFTs') 
        classified_PFTs_path = predict_plot_level_PFTs(site, 
                                                    year_inv, 
                                                    data_int_path, 
                                                    data_final_path, 
                                                    stacked_aop_path, 
                                                    ic_type_path, 
                                                    rf_model_path,
                                                    extracted_features_csv,
                                                    ic_type, 
                                                    pcaInsteadOfWavelengths)
        log.info('Classified PFTs for FATES initializaiton stored in: '
                f'{classified_PFTs_path}')

        ### 3) Generate lad profiles on new plots 
        lad_laz_clipped_path = lad_laz_clipped_path
        # if ic_type=="rs_inv_plots", use existing LAD profiles
        # or if ic_type=="rs_random_plots", generate lad profiles...
        if ic_type=="rs_random_plots":
            if not glob.glob(str(ic_type_path+'/clipped_to_plots/'+'*lad.json') ):
                # Generate new LAD profiles
                preprocess_lad(laz_path=lad_laz_clipped_path, 
                            inventory_path=None, 
                            site=site, 
                            year=year_inv, 
                            output_data_path=ic_type_path, 
                            use_case="predict",
                            end_result=True)
                #ais why can LAD not be processed for some sites?
    
        ### 4) Generate cohort and patch files
            # second half of step (1)
        cohort_path, patch_path = generate_cohort_patch_files(site=site, 
                                    year=year_inv, 
                                    data_int_path=data_int_path, 
                                    biomass_path=biomass_path,
                                    ic_type=ic_type, 
                                    ic_type_path=ic_type_path, 
                                    plots_shp_path=plots_shp_path,
                                    classified_PFTs_path=classified_PFTs_path,
                                    lad_laz_clipped_path=lad_laz_clipped_path)
    
    log.info('Cohort and patch files for FATES initializaiton generated')

    return cohort_path, patch_path

# create generate_fates_init_conditions script in new subfolder (not preprocessing/)
    # GENERATE FATES INITIAL CONDITIONS
        # with RS data
        # 1) generate random plots: first half of preprocessing/generate_ic_remotesensing.R
        # 2) predict pfts on new plots: neon-veg-SOAPpfts/11-predict-pft.R
        # 3) generate lad profiles on new plots: preprocessing/main.py with preprocess_lad==T for force run in sites.yaml
        # 4) generate cohort and patch files: second half of step (1)

# # TEST CASE
# site='SOAP',
# year_inv='2019',
# year_aop = '2019-06',
# use_case = "predict"
# ic_type = "rs_inv_plots" # field_inv_plots, rs_inv_plots, rs_random_plots
# aggregate_from_1m_to_2m_res = False
# px_thresh = 2
# ntree = 5000
# pcaInsteadOfWavelengths = True
# n_plots = 1000
# min_distance = 20
# plot_length = 20
# data_raw_aop_path = "/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/aop",
# data_int_path="/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/intermediate", 
# data_final_path='/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/final',
# rf_model_path = "/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/intermediate/SOAP/2019/training/rf_tree_crowns_training/rf_model_tree_crowns-training.RData",
# stacked_aop_path = "/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/intermediate/SOAP/2019/stacked_aop",
# biomass_path='/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/intermediate/SOAP/2019/biomass',
                       
# test = generate_initial_conditions(site=site[0],
#                             year_inv=year_inv[0],
#                             year_aop = year_aop[0],
#                             data_raw_aop_path = data_raw_aop_path[0],
#                             data_int_path=data_int_path[0], 
#                             data_final_path=data_final_path[0],
#                             rf_model_path = rf_model_path[0],
#                             stacked_aop_path = stacked_aop_path[0],
#                             biomass_path=biomass_path[0],
#                                         use_case=use_case,
#                                         ic_type=ic_type,
#                                         n_plots = n_plots,
#                                         min_distance = min_distance, 
#                                         plot_length = plot_length, 
#                                         aggregate_from_1m_to_2m_res = aggregate_from_1m_to_2m_res, 
#                                         pcaInsteadOfWavelengths = pcaInsteadOfWavelengths)

# 2+2
# # to do later ais (preprocessing/hyperspectral.py)
# add argument name before input in every function call in main.py 
# save tifs into their own data product subfolders e.g. raw/aop/tif/chm
# combine lidar and hs aop downloads into one step? think on this
# add 06-plot_aop_imagery function for example tile - is this helpful -maybe combine with script 10 function?
    # even since we have the full tile saved in stacked aop?
# create function for script 10 - save plots in diagnostics folder?
# check through functions in 00 script