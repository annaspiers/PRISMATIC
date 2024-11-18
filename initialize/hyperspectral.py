import logging
import os
import sys
import whitebox
import zipfile
import ray

import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import rasterio
from rasterio.plot import show
import fiona
from sklearn.preprocessing import OneHotEncoder
# from keras.models import Sequential
# from keras.layers import Conv2D, MaxPooling2D, Flatten, Dense
import geopandas as gpd
from shapely.geometry import Point
import rpy2.robjects as ro
    from collections import Counter

# data preparation
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import ADASYN
from imblearn.over_sampling import SVMSMOTE
from imblearn.pipeline import Pipeline as imbpipeline
from sklearn.pipeline import Pipeline

# model selecation
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split    

# model evaluation
from sklearn.metrics import f1_score

# machine learning
from sklearn.ensemble import RandomForestClassifier

# Save files
import joblib
from joblib import dump, load

from utils.download_functions import download_aop_files
from utils.apply_brdf_corrections import move_red_yellow_subset, \
                                        convert_hdf_to_envi, \
                                        create_config, \
                                        implement_brdf_correction , \
                                        convert_envi_to_tif
from pathlib import Path

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = str(conda_env_path/'share'/'proj')

log = logging.getLogger(__name__)

def download_hyperspectral(site, year, data_raw_aop_path, hs_type):
    """Download raw, uncorrected hyperspectral data (raw and raster) for site-year
    Analog to 04-download_aop_imagery.R from https://github.com/earthlab/neon-veg

    Parameters
    ----------
    site : str
        Site name
    year : str
        Year of aop collection
    data_raw_aop_path : str
        Path to store the downloaded data
    hs_type : str
        "tile" (L3 1km x 1km tile) OR "flightline" (L1 flightline)

    Returns
    -------
    (str, str)
        Path to the result folder for hyperspectral uncorrected imagery (L3 tiles) ('*/hs_tile_h5' or '*/hs_flightline_h5')
        and the result raster folder ('*/tif')
    """
    path = Path(data_raw_aop_path)/site/year

    product_codes = ['DP3.30015.001', 'DP3.30010.001'] 
    file_type = 'tif' 
    for product_code in product_codes: 
        p = path/file_type
        download_aop_files(product_code,
                           site,
                           year,
                           str(p),
                           match_string=file_type,
                           check_size=False)

    product_code = 'DP3.30026.001'
    p = path/file_type
    download_aop_files(product_code,
                       site,
                       year,
                       str(p),
                       match_string='.zip',
                       check_size=False)
    file_type = 'tif'
    zip_files = [file for file in os.listdir(p) if file.endswith('.zip')] # get the list of files
    for zip_file in zip_files:  #for each zipfile
        with zipfile.ZipFile(Path(p/zip_file)) as item: # treat the file as a zip
            item.extractall(p)  # extract it in the working directory
            item.close()
    vi_error_files = [file for file in os.listdir(p) if file.endswith('_error.tif')]
    for vi_error_file in vi_error_files: 
        os.remove(Path(p/vi_error_file)) # remove error files

    if hs_type=="tile":
        product_code = 'DP3.30006.001'
        p = path/'hs_tile_h5'
    elif hs_type=="flightline": 
        product_code = 'DP1.30006.001'
        p = path/'hs_flightline_h5'
    else:
        print("must specify hs_type argument")
    file_type = 'h5'
    download_aop_files(product_code,
                       site,
                       year,
                       str(p),
                       match_string=file_type,
                       check_size=False)
    
    return str(p)

def generate_pft_reference(sites, data_raw_inv_path, data_int_path, trait_table_path):
    """ Generate a csv that translates each species into a PFT usable as training data
    """
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

    # Generate PFT reference
    generate_pft_reference = ro.r('generate_pft_reference')
    Path(data_int_path).mkdir(parents=True, exist_ok=True)
    pft_reference_path = generate_pft_reference(sites, data_raw_inv_path, data_int_path, trait_table_path)
    log.info('Generated PFT reference for: {sites}'
             f'{pft_reference_path}')
    #ais now where to feed in pft_reference csv path


def correct_flightlines(site, year_inv, year_aop, data_raw_aop_path, data_int_path):
    """Correct raw L1 NEON AOP flightlines.
        1) Apply BRDF/topographic corrections
            output format: *.envi 
        2)  Convert envi to tif
            output format: *.tif
        3) Crop flightlines to AOP tile bounds
            output format: *.tif
        4) Merge tiffs that overlap spatially
            output format: *.tif
    """
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

    flightline_h5_path = os.path.join(data_raw_aop_path,site,year_aop,"hs_flightline_h5")
    
    # 1) Apply BRDF/topographic corrections
    log.info(f'Applying BRDF/topo corrections for: {site} {year_inv}')
    # move_red_yellow_subset(site, flightline_h5_path) 
    # flightline_envi_path = convert_hdf_to_envi(flightline_h5_path, 
    #                                            os.path.join(data_int_path,site,year_inv,
    #                                                         'hs_flightline/'))
    flightline_envi_path = os.path.join(data_int_path,site,year_inv, 'hs_flightline/')
    # config_path = create_config(flightline_envi_path)
    # implement_brdf_correction(config_path) 

    # 2)  Convert corrected envi to tif
    log.info(f'Converting from envi to tif format for: {site} {year_inv}')
    convert_envi_to_tif(site, flightline_h5_path, flightline_envi_path) 

    # 3) Crop flightlines to AOP tile bounds
    log.info(f'Cropping flightline tifs to AOP tile bounds for: {site} {year_inv}')
    crop_flightlines_to_tiles = ro.r('crop_flightlines_to_tiles')
    # ???? = crop_flightlines_to_tiles(site, year_inv, data_raw_aop_path, flightline_envi_path) #ais

    # 4) Merge tiffs that overlap spatially
    log.info(f'Merging spatially overlapping tiffs for: {site} {year_inv}')
    

def create_tree_crown_polygons(site, year, data_raw_inv_path, data_int_path, biomass_path, 
                               pft_reference_path, px_thresh):
    """Create geospatial features (points, polygons with half the maximum crown diameter) 
        for every tree in the NEON woody vegetation data set. 
        Analog to 02-create_tree_features.R from https://github.com/earthlab/neon-veg 

        Generate polygons that intersect with independent pixels in the AOP data 
        Analog to 03-process_tree_features.R from https://github.com/earthlab/neon-veg 

        Extract features (remote sensing data) for each sample (pixel) within the 
        specified shapefile (containing polygons that correspond to trees at the NEON site)
        Analog to 07-extract_training_features.R from https://github.com/earthlab/neon-veg 
    """

    # Create features (points or polygons) for each tree 
    log.info(f'Creating tree crown training data features for: {site} {year}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))
    
    # Create tree polygons
    # create_tree_crown_polygons = ro.r('create_tree_crown_polygons') #swap out for manually identified crowns
    create_tree_crown_polygons = ro.r('sort_out_manual_delineations') #ais replace python function name with this later on?
    training_shp_path = create_tree_crown_polygons(site, year, data_raw_inv_path, data_int_path, 
                                                   biomass_path, pft_reference_path, px_thresh)
    log.info('Clipped tree crown polygons data saved at: '
             f'{training_shp_path}')
    
    return training_shp_path


def prep_aop_imagery(site, year, hs_type, hs_path, tif_path, data_int_path, use_tiles_w_veg):
    """ Prepare aop imagery for extracting descriptive features for classification by 
            (1) clean/mask out shadow and non-veg
            (2) stacking imagery 
     Analog to 05-prep_aop_imagery.R from https://github.com/earthlab/neon-veg 

     ais to do: Then plots one tile of the prepped imagery. Analog to 06-plot_aop_imagery.R

    Calls a function written in R

    Parameters
    ----------
    site : str
        Site name
    year : str
        Year of aop collection
    hs_path : str
        Path where either raw hyperspectral h5 files OR BRDF-corrected tifs are stored
    tif_path : str
        Path where the raw aop tif files are stored

    Returns
    -------
    (str)
        Path to the result folder for stacked HS data ready for classification
    """
    log.info(f'Prepping aop data for: {site} {year}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

    # Clean/ mask imagery

    # Stack imagery
    prep_aop_imagery = ro.r('prep_aop_imagery')
    stacked_aop_path = prep_aop_imagery(site, year, hs_type, hs_path, tif_path, data_int_path, 
                                        use_tiles_w_veg)
    log.info('Stacked AOP data and saved at: '
             f'{stacked_aop_path}')
    return stacked_aop_path


def extract_spectra_from_polygon(site, year, shp_path, data_int_path, data_final_path, stacked_aop_path, 
                         use_case, ic_type, aggregate_from_1m_to_2m_res):
    """Create geospatial features (points, polygons with half the maximum crown diameter) 
        for every tree in the NEON woody vegetation data set. 
        Analog to 02-create_tree_features.R from https://github.com/earthlab/neon-veg 

        Generate polygons that intersect with independent pixels in the AOP data 
        Analog to 03-process_tree_features.R from https://github.com/earthlab/neon-veg 

        Extract features (remote sensing data) for each sample (pixel) within the 
        specified shapefile (containing polygons that correspond to trees at the NEON site)
        Analog to 07-extract_training_features.R from https://github.com/earthlab/neon-veg 
    """

    # Create features (points or polygons) for each tree 
    log.info(f'Extracting spectra from tree crowns for: {site} {year}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))
    
    # Extract training data from AOP data with tree polygons
    extract_spectra_from_polygon = ro.r('extract_spectra_from_polygon')
    training_spectra_path = extract_spectra_from_polygon(site=site, 
                                                         year=year, 
                                                         shp_path=shp_path,
                                                         data_int_path=data_int_path, 
                                                         data_final_path=data_final_path, 
                                                         stacked_aop_path=stacked_aop_path, 
                                                         use_case=use_case, 
                                                         ic_type=ic_type,
                                                         aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res)
 
    # # Extract features (remote sensing data) for each sample (pixel) within the
    # # specified shapefile (containing polygons that correspond to trees at the NEON site)
    
    # # Check if the stacked AOP path is empty
    # if not os.listdir(stacked_aop_path):
    #     print(f"Cannot extract data because stacked_aop_path is empty for {site} - {year}")
    #     return None
    
    # # Define the training data directory
    # training_data_dir = os.path.join(data_int_path, site, year, "training")
    
    # # Get a description of the shapefile to use for naming outputs
    # shapefile_description = os.path.basename(shp_path).split('.')[0]
    
    # # Specify destination for extracted features
    # if use_case == "train":
    #     extracted_features_path = training_data_dir
    #     extracted_features_filename = os.path.join(extracted_features_path,
    #                                                f"{shapefile_description}-extracted_features_inv.csv")
    # elif use_case == "predict":
    #     extracted_features_path = ic_type_path
    #     extracted_features_filename = os.path.join(extracted_features_path,
    #                                                f"{shapefile_description}-extracted_features.csv")
    # else:
    #     print("Need to specify use_case")
    #     return None
    
    # # Check if the extracted features file already exists
    # if os.path.exists(extracted_features_filename):
    #     print(f"Extracted features file already exists for {site} - {year}")
    #     return extracted_features_filename
    
    # # Load the shapefile
    # shp_gdf = gpd.read_file(shp_path)
    # stacked_aop_list = os.listdir(stacked_aop_path)

    # # Create a list to store the pixel values and corresponding shape IDs
    # pixel_values_list = []
    # shape_ids_list = []

    # # Loop through each shape in the shapefile
    # for i, row in shp_gdf.iterrows():
    #     # Get the geometry of the current shape
    #     geometry = row.geometry
        
    #     # Loop through each raster in the folder
    #     for file in stacked_aop_list:
    #         if file.endswith('.tif'):
    #             # Load the current raster
    #             raster_path = os.path.join(stacked_aop_path, file)
    #             with rasterio.open(raster_path) as src:
    #                 # Check if the current shape overlaps with the current raster
    #                 if src.bounds.intersects(geometry.bounds):
    #                     # Extract the pixel values for the current shape from the current raster
    #                     pixel_values, _ = rasterio.mask.mask(src, [geometry], crop=True, all_touched=True)
                        
    #                     # Get the shape ID
    #                     shape_id = row['shape_id']  # Replace 'shape_id' with the actual column name in your shapefile
                        
    #                     # Append the pixel values and shape ID to the lists
    #                     pixel_values_list.append(pixel_values)
    #                     shape_ids_list.append(np.full(pixel_values.shape[1:], shape_id, dtype=np.int64))

    # # Convert the lists to numpy arrays
    # pixel_values_array = np.concatenate(pixel_values_list, axis=1)
    # shape_ids_array = np.concatenate(shape_ids_list, axis=0)

    # # Save the pixel values and shape IDs to a numpy file
    # output_path = os.path.join(training_data_dir,'training_data.npy')
    # np.savez(output_path, pixel_values=pixel_values_array, shape_ids=shape_ids_array)
    
    # # Filter to tiles containing veg to speed up the next for-loop
    # if use_case == "train" or ic_type == "rs_inv_plots":
    #     tiles_w_veg = pd.read_csv(os.path.join(training_data_dir, "tiles_w_veg.txt"), header=None)
    #     tiles_w_veg = tiles_w_veg[0].tolist()
    #     stacked_aop_list = [file for file in stacked_aop_list if any(tile in file for tile in tiles_w_veg)]
    
        
    # # Aggregate to 2m resolution from 1m
    # if aggregate_from_1m_to_2m_res:
    #     # This step is not implemented in the original R code
    #     pass
    
    # # Clip the hyperspectral raster stack with the polygons within current tile.
    # # The returned objects are data frames, each row corresponds to a pixel in the
    # # hyperspectral imagery. The ID number refers to which tree that the 
    # # the pixel belongs to. A large polygon will lead to many extracted pixels
    # # (many rows in the output data frame), whereas tree stem points will
    # # lead to a single extracted pixel per tree. 
    

    # # Plot polygons of training data and RGB image of extracted pixels
    # for group in shapes_in_gdf['groupID'].unique():
    #     shapes_temp = shapes_in_gdf[shapes_in_gdf['groupID'] == group]
    #     buff_temp = shapes_temp.buffer(10)
    #     bbox_temp = buff_temp.bounds
    #     with rasterio.open(stacked_aop_filename) as src:
    #         window = rasterio.windows.from_bounds(bbox_temp.left, bbox_temp.bottom, bbox_temp.right, bbox_temp.top, src.transform)
    #         data = src.read(1, window=window)
    #         plt.imshow(data, cmap='gray')
    #         plt.show()
    
    # # Merge the extracted spectra and other data values with the tree info 
    # shapes_metadata = shapes_in_gdf.drop(columns=['geometry'])
    
    # # Combine the additional data with each spectrum for writing to file.
    # # Remove the geometry column to avoid issues when writing to csv later 
    # spectra_write = pd.concat([shapes_metadata, pd.DataFrame(extracted_spectra)], axis=1)
    
    # # Write extracted spectra and other remote sensing data values to file 
    # spectra_write.to_csv(os.path.join(extracted_features_path, f"extracted_features_{east_north_string}_{shapefile_description}.csv"), index=False)
    
    # # Combine all extracted features into a single .csv
    # paths_ls = [os.path.join(extracted_features_path, file) for file in os.listdir(extracted_features_path) if file.endswith('.csv')]
    
    # # Refine the output csv selection 
    # csvs = [file for file in paths_ls if shapefile_description in file]
    
    # # Combine all .csv data into a single data frame 
    # spectra_all = pd.concat([pd.read_csv(file) for file in csvs], ignore_index=True)
    
    # # Add site and year columns
    # spectra_all['site'] = site
    # spectra_all['year'] = year
    
    # # Write ALL the spectra to a single .csv file 
    # spectra_all.to_csv(extracted_features_filename, index=False)
    
    # # Delete the individual csv files for each tile 
    # for file in csvs:
    #     os.remove(file)


    log.info('Spectral features for training data saved at: '
             f'{training_spectra_path}')
    
    return training_spectra_path





# def custom_train_test_split(grouped,  train_size=0.8, random_state=42): #X, y, 
#     # Get the class counts
#     # class_counts = y.value_counts()
#     class_counts = pd.concat(grouped)['pft'].value_counts() #grouped.size()

#     # Get the number of cedar training samples
#     cedar_train_count = int(class_counts['cedar'] * train_size)
    
#     # Split the data for each class
#     train_data = pd.DataFrame()
#     test_data = pd.DataFrame()
    
#     for class_label, group in grouped:
#         #shape_ids = group['shapeID'].unique()

#         if class_label == 'cedar':
#             class_train_data, class_test_data = train_test_split(group, train_size=train_size, random_state=random_state)
#         elif class_counts[class_label] > class_counts['cedar']:
#             class_train_data, class_test_data = train_test_split(group, train_size=cedar_train_count, random_state=random_state)
#         else:
#             class_train_data, class_test_data = train_test_split(group, train_size=train_size, random_state=random_state)
        
#         train_data = pd.concat([train_data, class_train_data])
#         test_data = pd.concat([test_data, class_test_data])
    
#     return train_data, test_data



def train_pft_classifier(sites, data_int_path, pcaInsteadOfWavelengths, ntree, 
                         randomMinSamples, independentValidationSet, balance_training_to_min_PFT):    
    """Train a Random Forest (RF) model to classify tree PFT using in-situ tree
        measurements for PFT labels and remote sensing data as descriptive features
        Analog to 08-classify_species.R from https://github.com/earthlab/neon-veg 

        Train random forest model and assess PFT classification accuracy.
        Model outputs will be written to a folder within the training directory 
        starting with "rf_" followed by a description of each shapefile 
        containing points or polygons per tree. 

        Evaluate classifier performance: Out-of-Bag accuracy, independent
        validation set accuracy, and Cohen's kappa. Generate confusion matrix
        to show the classification accuracy for each PFT. 
        Analog to 09-assess_accuracy.R from https://github.com/earthlab/neon-veg 

        ntree 
          RF tuning parameter, number of trees to grow. default value 500
    
        randomMinSamples 
          To reduce bias, this boolean variable indicates whether the same number 
          of samples are randomly selected per PFT class. Otherwise,
          all samples per class are used for training the classifier. 
    
        independentValidationSet=T
          if TRUE, keep separate set for validation
          ais this doesnt work if for F - troubleshoot this later

    """
    log.info(f'Training model on sites: {sites}')

    features_df = pd.DataFrame()
    for site in sites:
        site_year_paths = [os.path.join(data_int_path, site, year) for year in os.listdir(os.path.join(data_int_path, site))]
        for site_year_path in site_year_paths:

            if not os.path.exists(os.path.join(site_year_path, "stacked_aop/wavelengths.txt")):
                continue

            # Load data
            wavelengths = pd.read_csv(os.path.join(site_year_path, "stacked_aop/wavelengths.txt"))
            # Stacked AOP layer names
            stacked_aop_layer_names = pd.read_csv(os.path.join(site_year_path, "stacked_aop/stacked_aop_layer_names.txt"))
            # Labelled, half-diam crown polygons to be used for training/validation
            shapefile_description = os.path.basename(os.path.join(site_year_path, "training/tree_crowns_training.shp")).split('.')[0]
            # Csv file containing extracted features
            extracted_features_filename = os.path.join(site_year_path, "training/tree_crowns_training-extracted_features_inv.csv")
            if not os.path.exists(extracted_features_filename):
                continue

            # Directory to save the classification outputs 
            rf_output_dir = os.path.join(data_int_path, "rf_dir")
            if not os.path.exists(rf_output_dir):
                os.makedirs(rf_output_dir)
                os.makedirs(os.path.join(rf_output_dir,"validation"))
            
            # # Filter out unwanted wavelengths
            # wavelength_lut <- filter_out_wavelengths(wavelengths=wavelengths, 
            #                                          layer_names=stacked_aop_layer_names)
            
            # Prep extracted features csv for RF model
            extracted_features_X_df = pd.read_csv(extracted_features_filename)
            # Remove X from wavelength column names
            features_temp = extracted_features_X_df.rename(columns=lambda x: x[1:] if x.startswith('X') else x)                  
            # filter the data to contain only the features of interest 
            # features_temp = extracted_features_df[featureNames] 
            # features_temp.columns = featureNames[len(wavelengths):] + [f"X{i+1}" for i in range(len(wavelengths))] 
            
            # make sure that column names are the same so we can rbind them across site-year
            numeric_columns = [col for col in features_temp.columns if col.isdigit()]
            features_temp.columns = features_temp.columns.map(lambda x: f'X{numeric_columns.index(x)+1}' if x in numeric_columns else x)

            features_df = pd.concat([features_df, features_temp])
          

    nPCs = 9
    # remove the individual spectral reflectance bands from the training data
    features_noWavelengths = features_df.drop([col for col in features_df.columns if col.startswith("X")], axis=1)
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features_df[[col for col in features_df.columns if col.startswith("X")]])
    pca = PCA(n_components=nPCs)
    features_pca = pca.fit_transform(features_scaled) #features_scaled[[col for col in features_df_scaled.columns if col.startswith("X")]])
    features = pd.concat([features_noWavelengths.reset_index(drop=True), pd.DataFrame(features_pca, columns=[f"PC{i+1}" for i in range(nPCs)])],axis=1)
    features.drop(columns=['indvdID','center_X','center_Y','groupID','pixelNumber','eastingIDs','northingIDs'], inplace=True)
    
    # visualize PCA    
    plt.figure(figsize=(8, 6))
    for value in features['pft'].unique():
        plt.scatter(features.loc[features['pft'] == value, 'PC1'], features.loc[features['pft'] == value, 'PC2'], label=value)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PC1 vs PC2')
    plt.legend()
    plt.savefig(os.path.join(rf_output_dir, 'pc1_vs_pc2.png')) # plot PC1 vs PC2
    
    plt.figure(figsize=(8, 6))
    for value in features['pft'].unique():
        plt.scatter(features.loc[features['pft'] == value, 'PC2'], features.loc[features['pft'] == value, 'PC3'], label=value)
    plt.xlabel('PC2')
    plt.ylabel('PC3')
    plt.title('PC2 vs PC3')
    plt.legend()
    plt.savefig(os.path.join(rf_output_dir, 'pc2_vs_pc3.png')) # plot PC2 vs PC3

    # group the data by 'pft' and 'shapeID'
    grouped_features = features.groupby(['pft', 'shapeID'])

    # create a list of groups
    groups = [group for _, group in grouped_features]

    # split the groups into training and testing sets
    train_groups, test_groups = train_test_split(groups, test_size=0.2, random_state=42)

    train_df = pd.concat(train_groups)
    test_df = pd.concat(test_groups)
    
    # # never finished this
    # # Balance training data so that no group has more training data than the smallest target PFT
    # # in this case, the smallest target PFT is cedar
    # if balance_training_to_min_PFT:

    #     # how many pixels in training data for smallest target PFT?
    #     min_pft = 'cedar'
    #     min_pft_pixel_count = len(train_df[train_df['pft'] == min_pft])

    #     for pft_class in train_df['pft'].unique():
    #         class_df = train_df[train_df['pft'] == pft_class]

    #         # while number of pixels in this class are much larger than cedar 
    #         # move all pixels from one shape into testing
    #         while len(class_df) > (min_pft_pixel_count*1.4): 
    #             #ais I just chose 1.4 arbitrarily as an acceptable size

    #             shapeID_to_move = class_df['shapeID'].unique()[0]  # move the first shape_ID
    #             rows_to_move = train_df[train_df['shapeID'] == shapeID_to_move]
    #             test_df = pd.concat([test_df, rows_to_move])
    #             train_df = train_df.drop(rows_to_move.index)
                
                
    # create the training and testing DataFrames
    train_df = train_df[train_df['pft'] != "other_herb"]
    train_df_final = train_df.drop(["shapeID",'site','year'],axis=1)
    test_df = test_df[test_df['pft'] != "other_herb"]
    test_df_final = test_df.drop(["shapeID",'site','year'],axis=1)
    # check the class distribution in the training and testing sets
    print("Training set class distribution:")
    print(train_df['pft'].value_counts(normalize=True))
    print("Testing set class distribution:")
    print(test_df['pft'].value_counts(normalize=True))

    # create a bar plot to visualize the count of data points in training and testing datasets for each prediction class
    plt.figure(figsize=(10, 6))
    train_counts = train_df['pft'].value_counts()
    test_counts = test_df['pft'].value_counts()
    plt.bar(train_counts.index, train_counts.values, label='Training')
    plt.bar(test_counts.index, test_counts.values, label='Testing')
    plt.xlabel('Prediction Class')
    plt.ylabel('Count')
    plt.title('Count of Data Points in Training and Testing Datasets')
    plt.legend()
    plt.savefig(os.path.join(rf_output_dir, 'train_test_count.png'))
    
    # # create a bar plot to visualize the percentage of data points in training and testing datasets for each prediction class
    # # ais doesn't visualize the %s I want
    # plt.figure(figsize=(10, 6))
    # train_percentages = (train_df['pft'].value_counts() / len(train_df)) * 100
    # test_percentages = (test_df['pft'].value_counts() / len(test_df)) * 100
    # plt.bar(train_percentages.index, train_percentages.values, label='Training')
    # plt.bar(test_percentages.index, test_percentages.values, label='Testing')
    # plt.xlabel('Prediction Class')
    # plt.ylabel('Percentage')
    # plt.title('Percentage of Data Points in Training and Testing Datasets')
    # plt.legend()
    # plt.savefig(os.path.join(rf_output_dir, 'train_test_percentage.png'))

    if randomMinSamples:
        minSamples = train_df["pft"].value_counts().min()
        train_df = train_df.groupby("pft").apply(lambda x: x.sample(minSamples)).reset_index(drop=True)

    rf_model_path = os.path.join(rf_output_dir, f"rf_model_{shapefile_description}.joblib")
    if not os.path.exists(rf_model_path):
        rf_model = RandomForestClassifier(n_estimators=ntree, random_state=42, class_weight="balanced")
        rf_model.fit(train_df_final.drop("pft", axis=1), train_df_final["pft"])
        dump(rf_model, rf_model_path)
    else:
        rf_model = load(rf_model_path)

    y_pred_val = rf_model.predict(test_df_final.drop("pft", axis=1))
    accuracy_val = accuracy_score(test_df_final["pft"], y_pred_val)
    print("Validation Accuracy:", accuracy_val)
    print("Validation Classification Report:")
    print(classification_report(test_df_final["pft"], y_pred_val))
    print("Validation Confusion Matrix:")
    print(confusion_matrix(test_df_final["pft"], y_pred_val))
    
    # # Visualize classifications
    # # in the testing dataframe, group to identify which plants are near each other
    # # create a new column called "groupID"
    # vis_predict_df = pd.concat([test_df.drop(['index'],axis=1).reset_index(drop=True),pd.DataFrame(y_pred_val,columns=['pft_pred'])],axis=1,join='outer')
    # vis_predict_df['groupID'] = np.nan
    # group_id = 0
    # # loop through each row in the DataFrame
    # for i, row in vis_predict_df.iterrows():
    #     # check if the row has already been assigned a group ID
    #     if np.isnan(row['groupID']):
    #         group_id += 1
    #         vis_predict_df.loc[i, 'groupID'] = group_id

    #         # find all rows within 30m of the current row
    #         nearby_rows = vis_predict_df[(vis_predict_df['eastingIDs'] >= row['eastingIDs'] - 30) &
    #                         (vis_predict_df['eastingIDs'] <= row['eastingIDs'] + 30) &
    #                         (vis_predict_df['northingIDs'] >= row['northingIDs'] - 30) &
    #                         (vis_predict_df['northingIDs'] <= row['northingIDs'] + 30)]
    #         # assign the nearby rows to the same group ID
    #         vis_predict_df.loc[nearby_rows.index, 'groupID'] = group_id

    # # Create a colormap with a discrete color ramp
    # cmap = cm.get_cmap('tab10', len(vis_predict_df['pft'].unique()))
    # color_dict = {value: cmap(i) for i, value in enumerate(vis_predict_df['pft'].unique())}    
    # vis_predict_df['color'] = vis_predict_df['pft'].map(color_dict) # Add a new column with the corresponding colors

    # # Plot predictions
    # for i in vis_predict_df['groupID'].unique():
    #     vis_predict_temp = vis_predict_df[vis_predict_df['groupID']==i]
    #     fig, ax = plt.subplots(figsize=(10, 10))

    #     x = int(vis_predict_temp['eastingIDs'].unique())
    #     y = int(vis_predict_temp['northingIDs'].unique())
        
    #     xmin = x+1000
    #     xmax = x
    #     ymin = y
    #     ymax = y-1000
    #     for i, row in vis_predict_temp.iterrows():
    #         if (x+serial_index_to_coordinates(row['pixelNumber'],1000)[0]) < xmin:
    #             xmin = x+serial_index_to_coordinates(row['pixelNumber'],1000)[0]
    #         if (x+serial_index_to_coordinates(row['pixelNumber'],1000)[0]) > xmax:
    #             xmax = x+serial_index_to_coordinates(row['pixelNumber'],1000)[0]
    #         if (y-(1000-serial_index_to_coordinates(row['pixelNumber'],1000)[1])) < ymin:
    #             ymin = y-(1000-serial_index_to_coordinates(row['pixelNumber'],1000)[1])
    #         if (y-(1000-serial_index_to_coordinates(row['pixelNumber'],1000)[1])) > ymax:
    #             ymax = y-(1000-serial_index_to_coordinates(row['pixelNumber'],1000)[1])

    #     # search for the image that includes the specified easting and northing coordinate
    #     imagery_folder = os.path.join(data_int_path,row['site'],str(row['year']),"stacked_aop")
    #     for file in os.listdir(imagery_folder):
    #         if file.endswith('.tif'):
    #             imagery_path = os.path.join(imagery_folder, file)
    #             imagery = rasterio.open(imagery_path)
    #             if imagery.bounds.left == x and y == imagery.bounds.top:

    #                 # plot the remote sensing pixels used for validation with a 10 pixel buffer
    #                 window = rasterio.windows.Window(xmin-10, ymin-10, xmax-xmin+20,ymax-ymin+20)
    #                 data = imagery.read(window=window)
    #                 # multiband image
    #                 # specify which three bands to plot
    #                 bands = [9, 10, 11]  # e.g., plot bands 4, 3, and 2 as RGB                   
    #                 data = imagery.read(bands)  # read the specified bands                    
    #                 data = np.transpose(data, (1, 2, 0)) # transpose the data to the correct shape for plotting

    #                 # plot the data with RGB colors
    #                 ax.imshow(data, cmap='viridis', extent=(xmin-10, xmax+10, ymin-10, ymax+10))

    #                 # loop through each row in the DataFrame
    #                 for i, row in vis_predict_temp.iterrows():
    #                     x_temp = serial_index_to_coordinates(row['pixelNumber'],1000)[0]
    #                     y_temp = serial_index_to_coordinates(row['pixelNumber'],1000)[1]
                        
    #                     # plot the classified predictions on top with a small amount of transparency
    #                     ax.scatter(x, y, c=row['color'], alpha=0.5, label=row['pft']) 

    #                 break
         
    #     ax.legend() # add a legend
    #     # save the plot
    #     plt.savefig(os.path.join(rf_output_dir, "validation", row['site']+"_"+str(row['year'])+'_'+str(round(imagery.bounds.left))+'_'+str(round(imagery.bounds.top))+'.png'))
    
    log.info('Trained PFT classifier saved in this folder: '
             f'{rf_model_path}')
    
    return rf_model_path    


def serial_index_to_coordinates(serial_index, raster_width=1000):
    """
    Translates a serial index of a raster into coordinates of the raster.

    Parameters:
    serial_index (int): The serial index of the pixel in the raster.
    raster_width (int): The width of the raster.

    Returns:
    tuple: A tuple containing the x and y coordinates of the pixel in the raster.
    """
    x = (serial_index // raster_width)-1
    y = 1000 - (serial_index % raster_width) 
    return x, y





# def train_pft_classifier(sites, data_int_path,
#                          pcaInsteadOfWavelengths, ntree, randomMinSamples, independentValidationSet):

# # used originaly when I was training the RF in R

#     """Train a Random Forest (RF) model to classify tree PFT using in-situ tree
#         measurements for PFT labels and remote sensing data as descriptive features
#         Analog to 08-classify_species.R from https://github.com/earthlab/neon-veg 

#         Evaluate classifier performance: Out-of-Bag accuracy, independent
#         validation set accuracy, and Cohen's kappa. Generate confusion matrix
#         to show the classification accuracy for each PFT. 
#         Analog to 09-assess_accuracy.R from https://github.com/earthlab/neon-veg 
#     """
    
#     log.info(f'Training model on sites: {sites}')
#     r_source = ro.r['source']
#     r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))
    
#     # Train random forest model 
#     train_pft_classifier = ro.r('train_pft_classifier')
#     rf_model_path = train_pft_classifier(sites, data_int_path, pcaInsteadOfWavelengths, ntree, 
#                                          randomMinSamples, independentValidationSet)
#     log.info('Trained PFT classifier saved in this folder: '
#              f'{rf_model_path}')
    
#     return rf_model_path





# def train_pft_classifier2(sites, data_int_path, pcaInsteadOfWavelengths, ntree, 
#                          randomMinSamples, independentValidationSet):  

# # Started adapting from Fricker et al

#     # Load shapefiles and remote sensing imagery
#     shapefiles = []
#     for file in os.listdir('shapefiles'):
#         shapefiles.append(fiona.open(os.path.join('shapefiles', file), 'r'))

#     imagery = rasterio.open('imagery.tif', 'r')

#     # Extract pixels and create dataset
#     dataset = []
#     for shapefile in shapefiles:
#         for feature in shapefile:
#             pixels = []
#             for x, y in feature['geometry']['coordinates'][0]:
#                 pixels.append(imagery.read(1, window=rasterio.windows.Window(x, y, 1, 1)))
#             dataset.append((np.array(pixels), feature['properties']['pft']))

#     # Split dataset into training and testing data
#     train_data, test_data = train_test_split(dataset, test_size=0.2, random_state=42) 

#     # Normalize pixel values
#     scaler = StandardScaler()
#     for i, (pixels, label) in enumerate(train_data):
#         train_data[i] = (scaler.fit_transform(pixels), label)

#     for i, (pixels, label) in enumerate(test_data):
#         test_data[i] = (scaler.transform(pixels), label)

#     # One-hot encode labels
#     encoder = OneHotEncoder()
#     labels = [label for _, label in train_data]
#     encoder.fit(labels)
#     for i, (pixels, label) in enumerate(train_data):
#         train_data[i] = (pixels, encoder.transform([label]))

#     for i, (pixels, label) in enumerate(test_data):
#         test_data[i] = (pixels, encoder.transform([label]))

#     # TRAIN THE CNN
#     # Define CNN architecture
#     model = Sequential()
#     model.add(Conv2D(32, (3, 3), activation='relu', input_shape=(1, 1, 3)))
#     model.add(MaxPooling2D((2, 2)))
#     model.add(Conv2D(64, (3, 3), activation='relu'))
#     model.add(MaxPooling2D((2, 2)))
#     model.add(Conv2D(128, (3, 3), activation='relu'))
#     model.add(MaxPooling2D((2, 2)))
#     model.add(Flatten())
#     model.add(Dense(128, activation='relu'))
#     model.add(Dense(len(encoder.categories_[0]), activation='softmax'))

#     # Compile model
#     model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

#     # Train model
#     model.fit([pixels for pixels, _ in train_data], [label for _, label in train_data], epochs=10, batch_size=32)

#     # VISUALIZE

#     # Make predictions on testing data
#     predictions = model.predict([pixels for pixels, _ in test_data])

#     # Visualize predictions
#     for i, (pixels, label) in enumerate(test_data):
#         prediction = predictions[i]
#         pft = encoder.inverse_transform(prediction)
#         plt.imshow(imagery.read(1, window=rasterio.windows.Window(pixels[0, 0], pixels[0, 1], 1, 1)), cmap='gray')
#         plt.text(pixels[0, 0], pixels[0, 1], pft, fontsize=10, color='red')
#         plt.savefig(f'prediction_{i}.png')
#         plt.close()