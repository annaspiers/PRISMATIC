import logging
import os
import sys
import glob

import pandas as pd
import geopandas as gpd
import rasterio
import rasterio.mask
from rasterio.plot import show
import numpy as np
from shapely.geometry import Point
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
import joblib
from joblib import dump, load

import rpy2.robjects as ro

from pathlib import Path

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = str(conda_env_path/'share'/'proj')

# custom functions in another script
from initialize.lad import prep_lad
from initialize.hyperspectral import filter_out_wavelengths, \
                                        extract_spectra_from_polygon

log = logging.getLogger(__name__)

def generate_initial_conditions(site, year_inv, year_aop, data_raw_aop_path, data_int_path,

                                data_final_path, rf_model_path, stacked_aop_path, biomass_path,
                                use_case, ic_type, ic_type_path, n_plots, min_distance, plot_length, 
                                aggregate_from_1m_to_2m_res, pcaInsteadOfWavelengths, multisite):                                            
    """Generate FATES initial conditions from remote sensing data in three scenarios:
            1) inventory-based IC's over inventory plots (ic_type = field_inv_plots)
            2) remote-sensing-based IC's from plots used for forest inventory 
                    (ic_type = rs_inv_plots)
            3) remote-sensing-based IC's from plots gridded wall to wall across 
                    the eddy covariance tower spatial extent (ic_type = rs_tower_ftpt)
            4) remote-sensing-based IC's from plots randomly generated across 
                the AOP spatial extent (ic_type = rs_random_plots)
            5) remote-sensing-based IC's from plots gridded wall to wall across 
                the AOP spatial extent (ic_type = rs_wall2wall)

    We generate files formatted to initialize FATES: cohort (.css) and patch (.pss) files
    """

    if not os.path.exists(Path(data_final_path,site,year_inv)):
        os.makedirs(Path(data_final_path,site,year_inv))   

    log.info(f'Generating initial conditions for: {site} {year_inv}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'initial_conditions_helper.R'))
    generate_pcss = ro.r('generate_pcss') 
    
    if (ic_type=="field_inv_plots"):
        cohort_path, patch_path = generate_pcss(site=site, 
                                    year=year_inv, 
                                    data_int_path=data_int_path, 
                                    biomass_path=biomass_path,
                                    ic_type=ic_type, 
                                    ic_type_path=ic_type_path)

    else :  #any of the rs_* ic_types      
        ### 1) Specify spatial extent of initial conditions and prepare spatial data by
        # clipping point clouds to and extracting raster data from plots
        plots_shp_path, plots_laz_path, extracted_features_csv = prepare_spatial_data_for_ic(
                                            site, year_inv, year_aop, ic_type, data_raw_aop_path, 
                                            data_int_path, stacked_aop_path, data_final_path, 
                                            ic_type_path, use_case, n_plots, min_distance, 
                                            plot_length, aggregate_from_1m_to_2m_res) 
        log.info('Plots, clipped laz, and clipped stacked rasters prepared for generating initial conditions')

        ### 2) Predict pfts on new plots
        # proxy to neon-veg-SOAPpfts/11-predict-pft.R
        # if ic_type==rs_inv_plots, use known PFTs from NEON inventory data 
        # if ic_type==rs_random_plots, predict PFTs with RF classifier
        classified_PFTs_path = predict_plot_level_PFTs(stacked_aop_path=stacked_aop_path, 
                                                    ic_type=ic_type, 
                                                    ic_type_path=ic_type_path, 
                                                    rf_model_path=rf_model_path,
                                                    extracted_features_filename=extracted_features_csv[0],
                                                    pcaInsteadOfWavelengths=pcaInsteadOfWavelengths)
        log.info('Classified PFTs for FATES initializaiton stored in: '
                f'{classified_PFTs_path}')

        ### 3) Generate lad profiles on new plots 
        # if ic_type=="rs_inv_plots", use existing LAD profiles
        # or if ic_type=="rs_tower_ftpt","rs_random_plots","rs_wall2wall" generate lad profiles...
        if ic_type!="rs_inv_plots":
            if not glob.glob(str(ic_type_path+'/clipped_to_plots/'+'*lad.json') ):
                # Generate new LAD profiles
                prep_lad(laz_path=plots_laz_path, 
                            inventory_path=None, 
                            site=site, 
                            year=year_inv, 
                            output_path=ic_type_path, 
                            use_case="predict")
                #ais why can LAD not be processed for some sites?

        ### 4) Generate cohort and patch files
            # second half of step (1)
        cohort_path, patch_path = generate_pcss(site=site, 
                                        year=year_inv, 
                                        data_int_path=data_int_path, 
                                        biomass_path=biomass_path,
                                        ic_type=ic_type, 
                                        ic_type_path=ic_type_path, 
                                        plots_shp_path=plots_shp_path,
                                        classified_PFTs_path=classified_PFTs_path,
                                        plots_laz_path=plots_laz_path,
                                        multisite=multisite)
    
    log.info('Cohort and patch files for FATES initialization generated')

    return cohort_path, patch_path


def prepare_spatial_data_for_ic(site, year_inv, year_aop, ic_type, data_raw_aop_path, 
                                data_int_path, stacked_aop_path, data_final_path, ic_type_path,
                                use_case, n_plots, min_distance, plot_length, 
                                aggregate_from_1m_to_2m_res):
    # This function generates random plots across the remote sensing 
    # spatial extent and clips laz and rasters 
    generate_random_plots = ro.r('generate_random_plots') 
    generate_wall2wall_plots = ro.r('generate_wall2wall_plots') 
    clip_norm_laz_to_shp = ro.r('clip_norm_laz_to_shp') 
    # clip_lidar_to_polygon_lidR = ro.r('clip_lidar_to_polygon_lidR') 
    
    ### Specify spatial extent of initial conditions
    if ic_type == "rs_inv_plots": 
        # use existing inventory plot extents
        aoi_shp_dir = os.path.join(data_raw_aop_path, site, year_aop, "shape")
        for filename in os.listdir(aoi_shp_dir):
            if filename.endswith("merged_tiles.shp"):
                aoi_shp_path = os.path.join(aoi_shp_dir, filename)
        plots_shp_path = os.path.join(data_int_path, site, str(year_inv), "inventory_plots", "plots.shp")
        plots_laz_path = os.path.join(data_int_path, site, str(year_inv), "clipped_to_plots")

        plots_shp = gpd.read_file(plots_shp_path)
        aoi_shp   = gpd.read_file(aoi_shp_path)
        # Save figure with plots and AOI
        fig, ax = plt.subplots()
        aoi_shp.plot(ax=ax, color='black', edgecolor='black', linewidth=1)
        plots_shp.plot(ax=ax, color='red', edgecolor='red', linewidth=.75)
        plt.title('AOI and Plots for initial condition type ' + ic_type)
        plt.savefig(os.path.join(ic_type_path,'plot_coverage_across_AOI.png'), dpi=300, bbox_inches='tight')
                
    else:         
        plots_shp_path = os.path.join(ic_type_path, "ic_type_plots.shp")
        plots_laz_path = os.path.join(ic_type_path, "clipped_to_plots")

        if not os.path.exists(plots_shp_path):
            if not os.path.exists(plots_laz_path): os.makedirs(plots_laz_path)
            # Generate plots across AOI

            if ic_type == "rs_random_plots":
                aoi_shp_dir = os.path.join(data_raw_aop_path, site, year_aop, "shape")
                for filename in os.listdir(aoi_shp_dir):
                    if filename.endswith("merged_tiles.shp"):
                        aoi_shp_path = os.path.join(aoi_shp_dir, filename)
                generate_random_plots(aoi_shp_path, ic_type,ic_type_path, n_plots, min_distance, plot_length)
                    
            if ic_type == "rs_tower_ftpt":
                aoi_shp_dir = os.path.join(data_raw_aop_path, site)
                for filename in os.listdir(aoi_shp_dir):
                    if filename.endswith("90percent_flux.shp"):
                        aoi_shp_path = os.path.join(aoi_shp_dir, filename)                
                generate_wall2wall_plots(aoi_shp_path,plot_length,ic_type,ic_type_path)

            elif ic_type == "rs_wall2wall":
                aoi_shp_dir = os.path.join(data_raw_aop_path, site, year_aop, "shape")
                for filename in os.listdir(aoi_shp_dir):
                    if filename.endswith("merged_tiles.shp"):
                        aoi_shp_path = os.path.join(aoi_shp_dir, filename)
                generate_wall2wall_plots(aoi_shp_path,plot_length,ic_type,ic_type_path) 
            
        # Clip normalized point clouds to randomized plots    
        #if the number of clipped plots does not exqual the number of plot shapes, continue clipping laz
        if len(os.listdir(plots_laz_path)) != len(gpd.read_file(plots_shp_path)): 
            clip_norm_laz_to_shp(site, year_inv, data_int_path, ic_type_path, plots_shp_path)
        #ais sometimes thereare going to be some patches dropped - what to do?
       
    # ais tried changing this function to python, but it ran incredibly slowly and didn't work 
    # put too much time into trying to make it work (see scrap below)
    extracted_features_csv = extract_spectra_from_polygon(site=site,
                                                           year=year_inv,
                                                           data_int_path=data_int_path,
                                                           data_final_path=data_final_path,
                                                           stacked_aop_path=stacked_aop_path,
                                                           shp_path=plots_shp_path,
                                                           use_case=use_case,
                                                           aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res,                                                           
                                                           ic_type=ic_type)

    return (plots_shp_path, plots_laz_path, extracted_features_csv)



def predict_plot_level_PFTs(stacked_aop_path, ic_type_path, rf_model_path, 
                            extracted_features_filename,
                            ic_type, pcaInsteadOfWavelengths):

    classified_PFTs_path = os.path.join(ic_type_path,"pfts_per_plot.csv")

    if not os.path.exists(classified_PFTs_path):

        wavelengths = pd.read_csv(os.path.join(stacked_aop_path, 
                                               "wavelengths.txt"))['wavelengths'].tolist()
        stacked_aop_layer_names = pd.read_csv(os.path.join(stacked_aop_path, 
                                                           "stacked_aop_layer_names.txt"))['stacked_aop_layer_names'].tolist()
        #ais remove the following rows earlier in the workflow so I don't need to do it here
        # stacked_aop_layer_names.drop(columns=['pixelNumber','eastingIDs','northingIDs'], inplace=True)
        stacked_aop_layer_names = [x for x in stacked_aop_layer_names if x not in ['pixelNumber','eastingIDs','northingIDs']]
            
        # Filter out unwanted wavelengths
        wavelength_lut = filter_out_wavelengths(wavelengths=wavelengths, layer_names=stacked_aop_layer_names)
        
        # features to use in the RF models
        featureNames = ["shapeID"] + wavelength_lut['xwavelength'].tolist() + [name for name in stacked_aop_layer_names if not name.isdigit()]
        
        # Prep extracted features csv for RF model
        extracted_features_df = pd.read_csv(extracted_features_filename)
        # filter the data to contain only the features of interest 
        features_df = extracted_features_df[featureNames] 

        # Remove any rows with NA   
        features_df.dropna(inplace=True)
        #ais why are there rows with NAs?
        
        ### Predict PFT ID for FATES patches 
        print("Predicting plot-level PFTs")

         # Load RF model
        rf_model = joblib.load(rf_model_path) 
                
        if pcaInsteadOfWavelengths:
            # remove the individual spectral reflectance bands from the training data
            features_noWavelengths = features_df.drop([col for col in features_df.columns if col.startswith("X")], axis=1)
            
            # Define evaluation procedure for CV
            # cv = StratifiedKFold(n_splits = 3) #ais make k_fold a global param #, shuffle=True, random_state=10210)
            
            # PCA
            scaler = StandardScaler()
            features_scaled = scaler.fit_transform(features_df[[col for col in features_df.columns if col.startswith("X")]])
            nPCs = sum(1 for item in rf_model.feature_names_in_ if item.startswith('PC'))
            pca = PCA(n_components=nPCs) 
            #ais ^ I may have 9 PCs used in training the RF model, but to get to 99% variance for only rs_inv plots, 
            # one PC may be sufficient, so I'm forcing the same number of PCs as model
            features_pca = pca.fit_transform(features_scaled) 
            # nPCs = features_pca.shape[1] 
            features = pd.concat([features_noWavelengths.reset_index(drop=True), pd.DataFrame(features_pca, columns=[f"PC{i+1}" for i in range(nPCs)])],axis=1)
        else:
            features = features_df
            
        # Predict PFT ID
        features['pred_PFT'] = rf_model.predict(features.drop('shapeID',axis=1)) 
        
        # Summarize as percentages per plot
        features_pft_by_pct = features.groupby(['shapeID', 'pred_PFT']).size().reset_index(name='count')
        features_pft_by_pct['pct'] = features_pft_by_pct['count'] / features_pft_by_pct.groupby('shapeID')['count'].transform('sum')
        
        # Write to CSV
        features_pft_by_pct.to_csv(classified_PFTs_path, index=False)

        # Notes w Nicola
        # save geometry of each pixel in addition to the coordinates of the plot (for plotting RGB behind later)
        # then each point predicted has a coordinate and I can map this point-based
    
    return classified_PFTs_path

