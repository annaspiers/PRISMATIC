import logging
import os
import sys
import whitebox
import zipfile

import rpy2.robjects as ro

from utils.download_functions import download_aop_files
from pathlib import Path

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = str(conda_env_path/'share'/'proj')

log = logging.getLogger(__name__)
# wht = whitebox.WhiteboxTools()
# wht.set_verbose_mode(False)

# ais work on today -------------------------------------
# by end of today: successfully run create_training_data
# then: write function for 08-classify_pft, 09-assess_accuracy, 10-create_figures (and 01_main as needed)
# then: 11-predict_pft (and 01_main as needed)

# to do later
# save tifs into their own data product subfolders e.g. raw/aop/tif/chm
# combine lidar and hs aop downloads into one step? think on this
# add 06-plot_aop_imagery function for example tile - is this helpful 
    # even since we have the full tile saved in stacked aop?

def download_hs_L3_tiles(site, year, data_raw_aop_path):
    """Download raw, uncorrected hyperspectral data (raw and raster) for site-year
    Analog to 04-download_aop_imagery.R from https://github.com/earthlab/neon-veg

    Parameters
    ----------
    site : str
        Site name
    year : str
        Year of aop collection
    data_out_path : str
        Path to store the downloaded data

    Returns
    -------
    (str, str)
        Path to the result folder for hyperspectral uncorrected imagery (L3 tiles) ('*/hs_L3_path')
        and the result raster folder ('*/tif')
    """
    path = Path(data_raw_aop_path)/site/year
    
    product_code = 'DP3.30006.001'
    file_type = 'h5'
    p = path/'hs_L3_h5'
    download_aop_files(product_code,
                       site,
                       year,
                       str(p),
                       match_string=file_type,
                       check_size=False)
    
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
            os.remove(Path(p/zip_file)) #delete the zip_file too
    vi_error_files = [file for file in os.listdir(p) if file.endswith('_error.tif')]
    for vi_error_file in vi_error_files: 
        os.remove(Path(p/vi_error_file)) # remove error files
    
    return str(path/'hs_L3_tiles'), str(path/'tif')


def prep_aop_imagery(site, year, hs_L3_path, tif_path, data_out_path):
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
    hs_L3_path : str
        Path where raw hyperspectral h5 files are stored
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
    stacked_aop_path = prep_aop_imagery(site, year, hs_L3_path, tif_path, data_out_path)
    log.info('Downloaded inventory data saved at: '
             f'{stacked_aop_path}')
    return stacked_aop_path



def create_training_data(site, year, data_out_path, stacked_aop_path, use_case):
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
    training_shp_path = create_tree_crown_polygons(site, year, biomass_path, data_out_path)
    log.info('Clipped tree crown polygons data saved at: '
             f'{training_shp_path}')
    
    # Extract training data from AOP data with tree polygons
    extract_spectra_from_polygon = ro.r('extract_spectra_from_polygon')
    training_spectra_path = extract_spectra_from_polygon(site, year, data_out_path, stacked_aop_path, 
                                                         training_shp_path, use_case)
    log.info('Spectral features for training data saved at: '
             f'{training_spectra_path}')
    return training_spectra_path

# def train_pft_classifier(site, year, data_out_path):
#     """Create geospatial features (points, polygons with half the maximum crown diameter) 
#         for every tree in the NEON woody vegetation data set. 
#         Analog to 02-create_tree_features.R from https://github.com/earthlab/neon-veg 

#         Generate polygons that intersect with independent pixels in the AOP data 
#         Analog to 03-process_tree_features.R from https://github.com/earthlab/neon-veg 

#         Extract features (remote sensing data) for each sample (pixel) within the 
#         specified shapefile (containing polygons that correspond to trees at the NEON site)
#         Analog to 07-extract_training_features.R from https://github.com/earthlab/neon-veg 

#     """

#     # Create features (points or polygons) for each tree 
#     log.info(f'Creating tree crown training data features for: {site} {year}')
#     r_source = ro.r['source']
#     r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

#     # Create tree polygons
#     create_tree_crown_polygons = ro.r('create_tree_crown_polygons')
#     training_shp_path = create_tree_crown_polygons(site, year, biomass_path, data_out_path)
#     log.info('Clipped tree crown polygons data saved at: '
#              f'{training_shp_path}')
    
#     return ???

# create functions train_PFT_classifier and predict_PFT 
    # train_PFT_classifier()
    # Analog to 08, 09, 10 from https://github.com/earthlab/neon-veg
    # 
    # predict_PFT()
    # Analog to 11-predict_pft.R from neon-veg-SOAPpfts