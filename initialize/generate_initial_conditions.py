import logging
import os
import sys
import glob

import pandas as pd
import geopandas as gpd
import rasterio
from rasterio import mask
from rasterio.plot import show
from rasterio.features import geometry_mask
import numpy as np
from shapely.geometry import Point
from shapely.geometry import box
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

    log.info(f'Generating initial conditions for: {site} {year_inv} {ic_type}')
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
        classified_plot_PFTs_path = predict_pixel_PFTs(site=site, 
                                                       year_inv=year_inv, 
                                                       year_aop=year_aop,
                                                       data_raw_aop_path=data_raw_aop_path,
                                                       data_int_path=data_int_path,
                                                       stacked_aop_path=stacked_aop_path, 
                                                        ic_type=ic_type, 
                                                        plots_shp_path=plots_shp_path,
                                                        ic_type_path=ic_type_path, 
                                                        rf_model_path=rf_model_path,
                                                        pcaInsteadOfWavelengths=pcaInsteadOfWavelengths)
        # classified_plot_PFTs_path = predict_plot_level_PFTs(stacked_aop_path=stacked_aop_path, 
        #                                             ic_type=ic_type, 
        #                                             ic_type_path=ic_type_path, 
        #                                             rf_model_path=rf_model_path,
        #                                             extracted_features_filename=extracted_features_csv[0],
        #                                             pcaInsteadOfWavelengths=pcaInsteadOfWavelengths)
        log.info('Classified PFTs for FATES initializaiton stored in: '
                f'{classified_plot_PFTs_path}')

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
                                        classified_plot_PFTs_path=classified_plot_PFTs_path,
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
                    if filename.endswith("90pc_flux_ml.shp"):
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
       # if len(os.listdir(plots_laz_path)) != len(gpd.read_file(plots_shp_path)): 
       #     clip_norm_laz_to_shp(site, year_inv, data_int_path, ic_type_path, plots_shp_path)
        #ais sometimes thereare going to be some patches dropped - what to do?
       
    # ais tried changing this function to python, but it ran incredibly slowly and didn't work 
    # put too much time into trying to make it work (see scrap below)
    # extracted_features_csv = extract_spectra_from_polygon(site=site,
    #                                                        year=year_inv,
    #                                                        data_int_path=data_int_path,
    #                                                        data_final_path=data_final_path,
    #                                                        stacked_aop_path=stacked_aop_path,
    #                                                        shp_path=plots_shp_path,
    #                                                        use_case=use_case,
    #                                                        aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res,                                                           
    #                                                        ic_type=ic_type)
    extracted_features_csv = None

    return (plots_shp_path, plots_laz_path, extracted_features_csv)



def predict_plot_level_PFTs(data_int_path, stacked_aop_path, ic_type_path, rf_model_path, 
                            extracted_features_filename,
                            ic_type, pcaInsteadOfWavelengths):

    classified_plot_PFTs_path = os.path.join(ic_type_path,"pfts_per_plot.csv")

    if not os.path.exists(classified_plot_PFTs_path):

        wavelengths = pd.read_csv(os.path.join(stacked_aop_path, 
                                               "wavelengths.txt"))['wavelengths'].tolist()
        stacked_aop_layer_names = pd.read_csv(os.path.join(stacked_aop_path, 
                                                           "stacked_aop_layer_names.txt"))['stacked_aop_layer_names'].tolist()
        #ais remove the following rows earlier in the workflow so I don't need to do it here
        stacked_aop_layer_names = [x for x in stacked_aop_layer_names]
            
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
            
            # Load scaler and PCA used during training
            loaded_scaler = joblib.load(os.path.join(data_int_path, "rf_dir",'scaler.joblib'))
            loaded_pca = joblib.load(os.path.join(data_int_path, "rf_dir",'pca.joblib'))

            # PCA
            features_scaled = loaded_scaler.fit_transform(features_df[[col for col in features_df.columns if col.startswith("X")]])
            # nPCs = sum(1 for item in rf_model.feature_names_in_ if item.startswith('PC'))
            # pca = PCA(n_components=nPCs) 
            # #ais ^ I may have 9 PCs used in training the RF model, but to get to 99% variance for only rs_inv plots, 
            # # one PC may be sufficient, so I'm forcing the same number of PCs as model
            features_pca = loaded_pca.fit_transform(features_scaled) 
            features = pd.concat([features_noWavelengths.reset_index(drop=True), pd.DataFrame(features_pca, columns=[f"PC{i+1}" for i in range(nPCs)])],axis=1)
        else:
            features = features_df
            
        # Predict PFT ID
        features['pred_PFT'] = rf_model.predict(features.drop('shapeID',axis=1)) 
        
        # Summarize as percentages per plot
        features_pft_by_pct = features.groupby(['shapeID', 'pred_PFT']).size().reset_index(name='count')
        features_pft_by_pct['pct'] = features_pft_by_pct['count'] / features_pft_by_pct.groupby('shapeID')['count'].transform('sum')
        
        # Write to CSV
        features_pft_by_pct.to_csv(classified_plot_PFTs_path, index=False)

        # Notes w Nicola
        # save geometry of each pixel in addition to the coordinates of the plot (for plotting RGB behind later)
        # then each point predicted has a coordinate and I can map this point-based
    
    return classified_plot_PFTs_path


def read_raster(input_image):
    #author: Nicola Falco
    
    from osgeo import gdal
    
    # Tell GDAL to throw Python exceptions, and register all drivers
    gdal.UseExceptions()
    gdal.AllRegister()
    
    # Read in image
    img_gdal_obj = gdal.Open(input_image, gdal.GA_ReadOnly)
    
    # initialize img to zeros
    img = img_gdal_obj.ReadAsArray()
    
    return img, img_gdal_obj;



def image_vectorization_3d(img):
    #author: Nicola Falco
    
    new_shape = (img.shape[0], img.shape[1] * img.shape[2])
    img_vec = img[:, :, :].reshape(new_shape)
    img_vec = img_vec.T
    
    return img_vec



def array2gtiff(image_reference, map_img, format_img, savefile):
    #author: Nicola Falco

    from osgeo import gdal
    
    datatype = gdal.GDT_Float32
    
    img = gdal.Open(image_reference)
    cols = img.RasterXSize
    rows = img.RasterYSize

    driver = gdal.GetDriverByName(format_img)

    outDataRaster = driver.Create(savefile, cols, rows, 1, datatype)
    # sets same geotransform as input
    outDataRaster.SetGeoTransform(img.GetGeoTransform())
    # sets same projection as input
    outDataRaster.SetProjection(img.GetProjection())

    outDataRaster.GetRasterBand(1).WriteArray(map_img)
    outDataRaster.GetRasterBand(1).SetNoDataValue(-9999)

    outDataRaster.FlushCache()  # remove from memory
    del outDataRaster  # delete the data (not the actual geotiff)



def predict_pixel_PFTs(site, year_inv, year_aop, data_raw_aop_path, data_int_path, stacked_aop_path, plots_shp_path, ic_type_path, 
                            rf_model_path, ic_type, pcaInsteadOfWavelengths):

    classified_tiles_dir = os.path.join(data_int_path, site, year_inv, "rf_classified_tiles")
    if not os.path.exists(classified_tiles_dir):
        os.makedirs(classified_tiles_dir)

    # Load relevant data
    ic_type_plots_gpd = gpd.read_file(plots_shp_path)  
    stacked_aop_list = glob.glob(os.path.join(stacked_aop_path,"*.tif"))
    stacked_aop_layer_names = pd.read_csv(os.path.join(stacked_aop_path, 
                                                        "stacked_aop_layer_names.txt"))['stacked_aop_layer_names'].tolist()
    wavelengths = pd.read_csv(os.path.join(stacked_aop_path,  "wavelengths.txt"))['wavelengths'].tolist()
    rf_model = joblib.load(rf_model_path) 

    # Filter out unwanted wavelengths    
    wavelength_lut = filter_out_wavelengths(wavelengths=wavelengths, layer_names=stacked_aop_layer_names)

    # Create a mapping from unique strings to unique numeric values
    predicted_classes = ["grass_rock_bare", "pine", "cedar", "fir", "shrub", "oak"]
    string_to_number = {string: i for i, string in enumerate(predicted_classes)} 
    # Create a colormap
    cmap = plt.get_cmap('tab10', len(string_to_number))  # Using a colormap with enough colors
    cmap.set_under('lightgray')  # Color for NA values

    flipped_dict = {}
    for index, value in np.ndenumerate(arr):
        if value not in flipped_dict:
            flipped_dict[value] = []
        flipped_dict[value].append(index)
    return flipped_dict

    # Filter to stacked_aop that overlaps with ic_type shapes
    overlapping_coordinates = []
    for raster_path in stacked_aop_list:
        with rasterio.open(raster_path) as src:
           # Get the bounding box of the raster
            bounds = src.bounds  # returns (minx, miny, maxx, maxy)
            raster_box = box(bounds.left, bounds.bottom, bounds.right, bounds.top)

            # Check for overlap with polygons
            if ic_type_plots_gpd.geometry.intersects(raster_box).any():
                overlapping_coordinates.append(f'{int(bounds.left)}_{int(bounds.bottom)}')

    # Classify pixels in each overlapping tile
    for stacked_aop_filename in stacked_aop_list:

        # Check if the tif file overlaps with the polygons
        if any(coord in stacked_aop_filename for coord in overlapping_coordinates):
            
            # Load NEON tile
            tile_img, tile_gdal_obj = read_raster(stacked_aop_filename)
            [bands,rows,cols] = tile_img.shape

            #can I just use coord rather than all this below
            # Extract relevant info from tile for naming
            xmin = tile_gdal_obj.GetGeoTransform()[0] 
            ymin = tile_gdal_obj.GetGeoTransform()[3] - 1000 
            east_north_string = f"{int(xmin)}_{int(ymin)}" #ais get right indices

            # save classificaiton of tif file in intermediate data in new classified_pixels/ folder
            east_north_classified_tif_path = os.path.join(classified_tiles_dir,
                                                f"classified_{east_north_string}.tif")
                          
            # If CSV file already exists, skip
            if not os.path.exists(east_north_classified_tif_path):            
                # Image vectorization
                tile_vec = image_vectorization_3d(tile_img)

                # Convert the NumPy array to a Pandas DataFrame and assign column names
                tile_df = pd.DataFrame(tile_vec, columns=stacked_aop_layer_names)

                # features to use in the RF models
                rf_features_noWavelengths = [x for x in stacked_aop_layer_names if not x.isdigit()]
                rf_features_wavelengths = round(wavelength_lut['wavelength']).astype(int).astype(str).tolist() #wavelength_lut['xwavelength'].tolist()
                rf_features_all =  rf_features_noWavelengths + rf_features_wavelengths 
                
                # Filter image to selected bands
                tile_df_features = tile_df[rf_features_all] 
                
                # Remove rows with NaN values
                tile_df_features_noNA = tile_df_features.dropna()

                ### Predict PFT ID for FATES patches 
                print(f"Predicting pixel-level PFTs across tile: {east_north_string}")
                    
                if pcaInsteadOfWavelengths:

                    # Load scaler and PCA used during training
                    loaded_scaler = joblib.load(os.path.join(data_int_path, "rf_dir",'scaler.joblib'))
                    loaded_pca = joblib.load(os.path.join(data_int_path, "rf_dir",'pca.joblib'))
                    # loaded_pca.n_components = 6  

                    # PCA
                    features_scaled = loaded_scaler.fit_transform(tile_df_features_noNA[rf_features_wavelengths])
                    # nPCs = sum(1 for item in rf_model.feature_names_in_ if item.startswith('PC'))
                    # pca = PCA(n_components=nPCs) 
                    # #ais ^ I may have 9 PCs used in training the RF model, but to get to 99% variance for only rs_inv plots, 
                    # # one PC may be sufficient, so I'm forcing the same number of PCs as model
                    features_pca = loaded_pca.fit_transform(features_scaled) 
                    features = pd.concat([tile_df_features_noNA[rf_features_noWavelengths].reset_index(drop=True), 
                                    pd.DataFrame(features_pca, columns=[f"PC{i+1}" for i in range(features_pca.shape[1])])],axis=1)
                else:
                    features = tile_df_features_noNA
                
                # Predict PFT ID
                features['pred_PFT'] = rf_model.predict(features) 

                # Convert the array of strings to an array of numbers using the mapping
                pred_img_num = np.vectorize(string_to_number.get)(features['pred_PFT'])

                # Replace NaN values in the original DataFrame with predictions
                # Create a Series for predicted classes with the same index as the original DataFrame
                pred_PFT_withNAs = pd.Series(index=tile_df_features.index, dtype='object')
                pred_PFT_withNAs[tile_df_features_noNA.index] = pred_img_num  # Fill in predictions
                
                # reshape the vector to image shape
                pred_img = pred_PFT_withNAs.to_numpy().reshape(rows, cols)
                
                # im = plt.imshow(pred_img, cmap="jet")
                # plt.colorbar(im)
                # plt.show()
                
                #Save to geotiff
                format_img = "GTiff"    
                array2gtiff(stacked_aop_filename, pred_img, format_img, east_north_classified_tif_path)
    
    # Clip classified TIFFs to the ic_type plots and save them
    # Also plot each one side-by-side with RGB
        #ais I could integrate this for-loop below into the for-loop above, but keeping
        # them separate allows me to clip new plots
    for classified_tile_path in glob.glob(os.path.join(classified_tiles_dir, '*.tif')):

         # Check if the tif file overlaps with the polygons
        if any(coord in classified_tile_path for coord in overlapping_coordinates):
            
            with rasterio.open(classified_tile_path) as src:
                print(f"Clipping {classified_tile_path} to overlapping polygons...")

                # Filter to only polygons that overlap this tif
                classified_tile_bounds = src.bounds
                raster_box = box(classified_tile_bounds.left, classified_tile_bounds.bottom, classified_tile_bounds.right, classified_tile_bounds.top)
                overlapping_polygons = ic_type_plots_gpd[ic_type_plots_gpd.geometry.within(raster_box)]
                
                # Clip the raster to the extent of the polygons
                for _, polygon in overlapping_polygons.iterrows():

                    # Create a filename for the clipped raster
                    polygon_name = f"{polygon['PLOTID']}_{int(polygon['e_min'])}_{int(polygon['n_min'])}"  
                    clipped_plot_tif_path = os.path.join(ic_type_path, "clipped_to_plots", f"{polygon_name}_classified_clipped.tif")
                    
                    if not os.path.exists(clipped_plot_tif_path):   
                    
                        # Clip to the current polygon
                        out_image, out_transform = mask.mask(src, [polygon['geometry']], crop=True)
                        
                        with rasterio.open(clipped_plot_tif_path, 'w', driver='GTiff', height=out_image.shape[1],
                                        width=out_image.shape[2], count=1, dtype='int32',
                                        crs=src.crs, transform=out_transform) as dst:
                            dst.write(out_image.astype(rasterio.uint8))

                        # Clip RGB image to the polygon and plot predicted classes + RGB
                        for rgb_image_path in glob.glob(os.path.join(data_raw_aop_path,site,year_aop,'tif','*0_image.tif')):
                            rgb_filename = os.path.basename(rgb_image_path)
                            rgb_coords = rgb_filename.split('_')[-3:-1]  

                            # Check for overlap based on bounding box
                            if (int(rgb_coords[0]) == classified_tile_bounds[0] and int(rgb_coords[1]) == classified_tile_bounds[1]):

                                # Load corresponding RGB image
                                with rasterio.open(rgb_image_path) as src_rgb:
                                    out_image_rgb, out_transform_rgb = rasterio.mask.mask(src_rgb, [polygon['geometry']], crop=True)

                                # Plot and save image of clipped classificaiton and RGB rasters
                                fig, ax = plt.subplots(1, 2, figsize=(12, 6))

                                # Plot the classified raster
                                im_classified = ax[0].imshow(out_image.squeeze(), cmap=cmap, vmin=0, vmax=len(string_to_number)-1)  
                                # ax[0].imshow(out_image.squeeze(), cmap=cmap, vmin=0, vmax=len(string_to_number)-1)  
                                ax[0].set_title('Clipped Classified Raster')
                                ax[0].axis('off')

                                # Plot the RGB image
                                ax[1].imshow(out_image_rgb.transpose(1, 2, 0))  
                                ax[1].set_title('Clipped RGB Image')
                                ax[1].axis('off')

                                # Create a colorbar for the classified raster below the plots
                                cbar_classified = fig.colorbar(im_classified, ax=ax[0], orientation='horizontal', pad=0.1, aspect=50)
                                cbar_classified.set_ticks(list(string_to_number.keys()))
                                cbar_classified.set_ticklabels(list(string_to_number.values()))  # Set the labels to fruit names
                                cbar_classified.set_label('Predicted PFTs')

                                # Save the comparison plot
                                comparison_plot_path = os.path.join(ic_type_path, "clipped_to_plots", f"{polygon_name}_comparison.png")
                                plt.tight_layout()  # Adjust layout to prevent clipping of titles and labels
                                plt.savefig(comparison_plot_path, bbox_inches='tight')
                                plt.close(fig)  # Close the figure to free up memory

    return os.path.join(ic_type_path, "clipped_to_plots") 
