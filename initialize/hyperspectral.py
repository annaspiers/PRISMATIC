import logging
import os
import sys
import whitebox
import zipfile
import ray

import rpy2.robjects as ro

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

def download_hyperspectral(site, year, data_raw_aop_path, hs_type="tile"):
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

    flightline_h5_path = data_raw_aop_path + "/" + site + "/" + year_aop + "/hs_flightline_h5/"

    # 1) Apply BRDF/topographic corrections
    log.info(f'Applying BRDF/topo corrections for: {site} {year_inv}')
    move_red_yellow_subset(site, flightline_h5_path) 
    create_config(flightline_h5_path)
    convert_hdf_to_envi(flightline_h5_path) 
    implement_brdf_correction(flightline_h5_path)

    # result_ids = convert_hdf_to_envi.remote(flightline_h5_path) 
    # sum = 0 
    # while len(result_ids): 
    #     done_id, result_ids = ray.wait(result_ids) 
    #     sum = implement_brdf_correction(sum, flightline_h5_path) 
    # print("result = ", sum) 

    # # 2)  Convert corrected envi to tif
    # log.info(f'Converting from envi to tif format for: {site} {year_inv}')
    # convert_envi_to_tif(flightline_h5_path)

    # 3) Crop flightlines to AOP tile bounds
    log.info(f'Cropping flightlines to AOP tile bounds for: {site} {year_inv}')
    
    # 4) Merge tiffs that overlap spatially
    log.info(f'Merging spatially overlapping tiffs for: {site} {year_inv}')
    

def create_tree_crown_polygons(site, year, data_raw_inv_path, data_int_path, biomass_path, px_thresh):
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
    create_tree_crown_polygons = ro.r('create_tree_crown_polygons')
    training_shp_path = create_tree_crown_polygons(site, year, data_raw_inv_path, data_int_path, biomass_path, px_thresh)
    log.info('Clipped tree crown polygons data saved at: '
             f'{training_shp_path}')
    
    return training_shp_path


def prep_aop_imagery(site, year, hs_type, hs_path, tif_path, data_int_path, use_tiles_w_veg):
    """Stack the AOP remote sensing imagery and prepare it for
     extracting descriptive features for classification. 
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
    prep_aop_imagery = ro.r('prep_aop_imagery')
    stacked_aop_path = prep_aop_imagery(site, year, hs_type, hs_path, tif_path, data_int_path,use_tiles_w_veg)
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
                                                         use_case="train", 
                                                         ic_type=ic_type,
                                                         aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res)
    log.info('Spectral features for training data saved at: '
             f'{training_spectra_path}')
    
    return training_spectra_path



def train_pft_classifier(site, year, stacked_aop_path, training_shp_path, training_spectra_path, data_int_path,
                         pcaInsteadOfWavelengths, ntree, randomMinSamples, independentValidationSet):
    """Train a Random Forest (RF) model to classify tree PFT using in-situ tree
        measurements for PFT labels and remote sensing data as descriptive features
        Analog to 08-classify_species.R from https://github.com/earthlab/neon-veg 

        Evaluate classifier performance: Out-of-Bag accuracy, independent
        validation set accuracy, and Cohen's kappa. Generate confusion matrix
        to show the classification accuracy for each PFT. 
        Analog to 09-assess_accuracy.R from https://github.com/earthlab/neon-veg 
    """
    
    log.info(f'Creating tree crown training data features for: {site} {year}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

    # Train random forest model 
    train_pft_classifier = ro.r('train_pft_classifier')
    rf_model_path = train_pft_classifier(site, year, stacked_aop_path, training_shp_path, training_spectra_path, data_int_path,
                                         pcaInsteadOfWavelengths, ntree, randomMinSamples, independentValidationSet)
    log.info('Trained PFT classifier saved in this folder: '
             f'{rf_model_path}')
    
    return rf_model_path

