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

import rpy2.robjects as ro

from pathlib import Path
from initialize.lad import preprocess_lad

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
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

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
        # prepare_spatial_data_for_ic = ro.r('prepare_spatial_data_for_ic') 
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
                                                    extracted_features_csv[0],
                                                    ic_type, 
                                                    pcaInsteadOfWavelengths)
        log.info('Classified PFTs for FATES initializaiton stored in: '
                f'{classified_PFTs_path}')

        ### 3) Generate lad profiles on new plots 
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


def prepare_spatial_data_for_ic(site, year_inv, year_aop, ic_type, data_raw_aop_path, 
                                data_int_path, stacked_aop_path, data_final_path, ic_type_path,
                                use_case, n_plots, min_distance, plot_length, 
                                aggregate_from_1m_to_2m_res):
    # This function generates random plots across the remote sensing 
    # spatial extent and clips laz and rasters 
    generate_random_plots = ro.r('generate_random_plots') 
    clip_norm_laz_to_shp = ro.r('clip_norm_laz_to_shp') 
    extract_spectra_from_polygon = ro.r('extract_spectra_from_polygon') 
    
    ### Specify spatial extent of initial conditions
    if ic_type == "rs_inv_plots": 
        # use existing inventory plot extents
        plots_shp_path = os.path.join(data_int_path, site, str(year_inv), "inventory_plots", "plots.shp")
        laz_clipped_path = os.path.join(data_int_path, site, str(year_inv), "clipped_to_plots")
                
    elif ic_type == "rs_random_plots":
        
        plots_shp_path = os.path.join(ic_type_path, "random_plots.shp")
        laz_clipped_path = os.path.join(ic_type_path, "clipped_to_plots")

        if not os.path.exists(plots_shp_path):
            # Generate random plots across AOP spatial extent
            
            generate_random_plots(site, year_aop, year_inv, data_raw_aop_path, ic_type_path,
                                  n_plots, min_distance, plot_length)
            
        # Clip normalized point clouds to randomized plots
        if not os.path.exists(laz_clipped_path):
            print("2")
            laz_clipped_path = clip_norm_laz_to_shp(site, year_inv, data_int_path,
                                                    ic_type_path, plots_shp_path)
       
    extracted_features_csv = extract_spectra_from_polygon(site=site,
                                                           year=year_inv,
                                                           data_int_path=data_int_path,
                                                           data_final_path=data_final_path,
                                                           stacked_aop_path=stacked_aop_path,
                                                           shp_path=plots_shp_path,
                                                           use_case=use_case,
                                                           ic_type=ic_type,
                                                           ic_type_path=ic_type_path,
                                                           aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res)

    return (plots_shp_path, laz_clipped_path, extracted_features_csv)



# ais this function below is garbange and doesn't work and takes forever to run and lost me probably a full work day
# def extract_spectra_from_polygon(site, year, shp_path, data_int_path, data_final_path,
#                                 stacked_aop_path, use_case, ic_type, ic_type_path,
#                                 aggregate_from_1m_to_2m_res):
#     # Check if the shapefile is empty
#     shp_gdf = gpd.read_file(shp_path)
#     if shp_gdf.empty:
#         print(f"Shapefile is empty for {site} - {year}")
#         return None

#     # Define the output file path
#     training_data_dir = os.path.join(data_int_path, site, str(year), "training", 
#                                      f"{os.path.basename(shp_path).split('.')[0]}-extracted_features_inv.csv")
#     if os.path.exists(training_data_dir):
#         print(f"Output file already exists for {site} - {year}")
#         return training_data_dir
#     shapefile_description = os.path.basename(shp_path).split('.')[0]

#     # Specify destination for extracted features
#     if use_case == "train":
#         extracted_features_path = training_data_dir
#         extracted_features_filename = os.path.join(extracted_features_path,
#                                                   f"{shapefile_description}-extracted_features_inv.csv")
    
#     elif use_case == "predict":
#         extracted_features_path = ic_type_path
#         extracted_features_filename = os.path.join(extracted_features_path,
#                                                   f"{shapefile_description}-extracted_features.csv")
#     else:
#         print("need to specify use_case")
  
#     if not os.path.exists(extracted_features_filename): 
#         # Get the list of stacked AOP files
#         stacked_aop_files = [os.path.join(stacked_aop_path, f) for f in os.listdir(stacked_aop_path) if f.endswith('.tif')]

#         # Filter the stacked AOP files to only include those that intersect with the shapefile
#         shp_gdf['center_X'] = shp_gdf.centroid.x
#         shp_gdf['center_Y'] = shp_gdf.centroid.y
#         for file in stacked_aop_files: 
#             # read current tile of stacked AOP data
#             with rasterio.open(file) as src:
#                 stacked_aop_data = src.read()
#                 east_north_string = f"{int(src.transform[2])}_{int(src.transform[5])}"

#             # If csv file already exists, skip
#             if os.path.exists(os.path.join(extracted_features_path, f"extracted_features_{east_north_string}_{shapefile_description}.csv")):
#                 print("csv already exists... skipping to next coordinates")
#                 continue

#             # figure out which plots are within the current tile by comparing each
#             # X,Y coordinate to the extent of the current tile
#             shapes_in = shp_gdf[(shp_gdf['center_X'] >= src.bounds.left) &
#                                 (shp_gdf['center_X'] < src.bounds.right) &
#                                 (shp_gdf['center_Y'] >= src.bounds.bottom) &
#                                 (shp_gdf['center_Y'] < src.bounds.top)]
            
#             # if no polygons are within the current tile, skip to the next one
#             if shapes_in.empty:
#                 print("no shapes located within current tile... skipping to next shapefile")
#                 continue
#             else:
#                 if use_case == "train":
#                     print(f"Extracting {len(shapes_in)} trees in current tile, {east_north_string}")
#                 else:
#                     print(f"Extracting {len(shapes_in)} plots in current tile, {east_north_string}")
            
#             # Aggregate to 2m resolution from 1m
#             if aggregate_from_1m_to_2m_res:
#                 with rasterio.open(file) as src:
#                     stacked_aop_data = src.read()
#                     stacked_aop_data = np.mean(stacked_aop_data.reshape(-1, 2, 2), axis=(1, 2))

#             # clip the hyperspectral raster stack with the polygons within current tile.
#             shapes_in_gdf = gpd.GeoDataFrame(shapes_in, geometry='geometry')
#             with rasterio.open(file) as src:
#                 raster_data = src.read()
#                 crs = src.crs
#                 # extracted_spectra, _ = rasterio.mask.mask(src, shapes_in_gdf.geometry, crop=True)
                
#             # # Plot polygons of training data and RGB image of extracted pixels
#             # for group in shapes_in_gdf['groupID'].unique():
#             #     shapes_temp = shapes_in_gdf[shapes_in_gdf['groupID'] == group]
#             #     buff_temp = shapes_temp.buffer(10)
#             #     bbox_temp = buff_temp.total_bounds
#             #     with rasterio.open(file) as src:
#             #         cropped_rast, _ = rasterio.mask.mask(src, buff_temp.geometry, crop=True)
                
#             #     # save as png
#             #     plt.imsave(os.path.join(training_data_dir, f"group_{group}.png"), cropped_rast)
                            
#                 # Loop through each polygon in the shapefile
#                 pixels = []
#                 for index, row in shapes_in_gdf.iterrows():
#                     # Extract the pixels within the current polygon
#                     extracted_spectra, _ = rasterio.mask.mask(src, [row.geometry], crop=False)

#                     # Get the mask for the current polygon
#                     mask = extracted_spectra != src.nodata

#                     # Extract the pixel data for the current polygon
#                     pixel_values = raster_data[:, mask[0, :, :]]
#                     pixel_values = pixel_values.transpose((1, 0))   
                    
#                     # Create a DataFrame with the pixel data and the polygon name
#                     df = pd.DataFrame(pixel_values, columns=[f'band_{i+1}' for i in range(src.count)]) #ais bands stop at 426
#                     df['shapeID'] = row['PLOTID'] 

#                     # Add columns for the X and Y coordinates
#                     x_coords, y_coords = src.xy(mask[0, :, :].nonzero()[0],mask[0, :, :].nonzero()[1])
                                    
#                     # Convert the x and y coordinates to longitude and latitude
#                     east_coord, north_coords = rasterio.warp.transform(crs, 'EPSG:32611', x_coords, y_coords)
                    
#                     # Add columns for the longitude and latitude coordinates
#                     df['easting'] = east_coord
#                     df['northing'] = north_coords

#                     # Remove rows with NaN values
#                     df_no_nan = df.dropna(subset=['band_1'])
                    
#                     # Append the pixel DataFrame to the list
#                     pixels.append(df_no_nan)

#                 # Concatenate the list of DataFrames into a single DataFrame
#                 pixel_df = pd.concat(pixels, ignore_index=True)

#             # AIS later when I'm editing, witch from workign with a 2D dataframe to a 3D array. 
#             # ask cborg this: "in python, I load a multiband raster and a shapefile with 
#                 # multiple polygons. I want to mask the raster pixels for each polygon and add a layer 
#                 # labelling which polygon each unmasked pixel belongs to"

#             # # remove the geometry column to avoid issues when writing to csv later 
#             #     spectra_write <- merge(shapes_metadata,extracted_spectra,by="shapeID") %>% 
#             #         dplyr::select(shapeID, everything()) %>% 
#             #         dplyr::select(-geometry)
            
#             # write extracted spectra and other remote sensing data values to file
#             # ais these lines below here take a long time
#             pixel_df.to_csv(os.path.join(extracted_features_path, 
#                                                 f"extracted_features_{east_north_string}_{shapefile_description}.csv"), 
#                                                 index=False)

#     # combine all extracted features into a single .csv
#     paths_ls = [os.path.join(extracted_features_path, f) for f in os.listdir(extracted_features_path) if f.endswith('.csv')]
    
#     # refine the output csv selection
#     csvs = [f for f in paths_ls if f"000_{shapefile_description}.csv" in f]

#     # header_saved = False
#     # with open(extracted_features_filename,'wb') as fout:
#     #     for file in csvs:
#     #         with open(file, 'r') as fin:
#     #             # header = next(fin)
#     #             # if not header_saved:
#     #             #     fout.write(header)
#     #             #     header_saved = True
#     #             for line in fin:
#     #                 fout.write(line)

    
#     # combine all .csv data into a single data frame
#     spectra_all =  pd.concat(map(pd.read_csv, csvs), ignore_index=True) #10:59
#     # pd.concat([pd.read_csv(f) for f in csvs]) #10:14 ais this takes an incredibly long time
    
#     spectra_all['site'] = site
#     spectra_all['year'] = year
        
#     # write ALL the spectra to a single .csv file 
#     spectra_all.to_csv(extracted_features_filename, index=False)
        
#     # delete the individual csv files for each tile
#     for f in csvs:
#         os.remove(f)
    
#     return extracted_features_filename



import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

# def predict_plot_level_PFTs(site, year, data_int_path, data_final_path,
#                             stacked_aop_path, ic_type_path, rf_model_path, extracted_features_csv,
#                             ic_type, pcaInsteadOfWavelengths):
    
#     classified_PFTs_path = os.path.join(ic_type_path, "pfts_per_plot.csv")
    
#     if not os.path.exists(classified_PFTs_path):
#         print("Predicting plot-level PFTs")
        
#         # Load wavelengths and stacked AOP layer names
#         wavelengths = pd.read_csv(os.path.join(stacked_aop_path, "wavelengths.txt"))['wavelengths'].round(0).astype(int)
#         stacked_aop_layer_names = pd.read_csv(os.path.join(stacked_aop_path, "stacked_aop_layer_names.txt"))['stacked_aop_layer_names']
        
#         # Filter out unwanted wavelengths
#         filter_out_wavelengths = ro.r('filter_out_wavelengths') 
#         wavelength_lut = filter_out_wavelengths(wavelengths, stacked_aop_layer_names)
        
#         # Define feature names
#         featureNames = ["shapeID"] + wavelength_lut['xwavelength'].tolist() + [name for name in stacked_aop_layer_names if not name.isdigit()]
        
#         # Load extracted features CSV
#         extracted_features_df = pd.read_csv(extracted_features_csv)
        
#         # Select features of interest
#         features_df = extracted_features_df[featureNames]
        
#         # Remove rows with NAs
#         features_df_noNA = features_df.dropna()
        
#         # Perform PCA if specified
#         if pcaInsteadOfWavelengths:
#             features = apply_PCA(features_df_noNA, wavelength_lut, ic_type_path)
#         else:
#             features = features_df_noNA
        
#         # Load RF model
#         rf_model = pd.read_pickle(rf_model_path)
        
#         # Predict PFT ID
#         features['pred_PFT'] = rf_model.predict(features)
        
#         # Summarize as percentages per plot
#         features_pft_by_pct = features.groupby(['shapeID', 'pred_PFT']).size().reset_index(name='count')
#         features_pft_by_pct['pct'] = features_pft_by_pct['count'] / features_pft_by_pct.groupby('shapeID')['count'].transform('sum')
        
#         # Add "_central" to shapeID if ic_type is "rs_inv_plots"
#         if ic_type == "rs_inv_plots":
#             features_pft_by_pct['shapeID'] = features_pft_by_pct['shapeID'] + "_central"
        
#         # Write to CSV
#         features_pft_by_pct.to_csv(classified_PFTs_path, index=False)
    
#     return classified_PFTs_path


# # TEST CASE
# site='SJER',
# year_inv='2021',
# year_aop = '2021-06',
# use_case = "predict"
# ic_type = "rs_inv_plots" # field_inv_plots, rs_inv_plots, rs_random_plots
# aggregate_from_1m_to_2m_res = False
# px_thresh = 2
# ntree = 5000
# pcaInsteadOfWavelengths = True
# n_plots = 1000
# min_distance = 20
# plot_length = 20
# data_raw_aop_path = "/Users/AISpiers/dev/RS-PRISMATIC/inititalize/data/raw/aop",
# data_int_path="/Users/AISpiers/dev/RS-PRISMATIC/inititalize/data/intermediate", 
# data_final_path='/Users/AISpiers/dev/RS-PRISMATIC/inititalize/data/final',
# rf_model_path = "/Users/AISpiers/dev/RS-PRISMATIC/inititalize/data/intermediate/SOAP/2019/training/rf_tree_crowns_training/rf_model_tree_crowns-training.RData",
# stacked_aop_path = "/Users/AISpiers/dev/RS-PRISMATIC/inititalize/data/intermediate/SOAP/2019/stacked_aop",
# biomass_path='/Users/AISpiers/dev/RS-PRISMATIC/inititalize/data/intermediate/SOAP/2019/biomass',
                       
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
# # to do later ais (inititalize/hyperspectral.py)
# add argument name before input in every function call in main.py 
# save tifs into their own data product subfolders e.g. raw/aop/tif/chm
# combine lidar and hs aop downloads into one step? think on this
# add 06-plot_aop_imagery function for example tile - is this helpful -maybe combine with script 10 function?
    # even since we have the full tile saved in stacked aop?
# create function for script 10 - save plots in diagnostics folder?
# check through functions in 00 script