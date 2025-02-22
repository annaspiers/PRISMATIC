library(lidR)
library(terra)
library(neonUtilities)
options(stringsAsFactors = FALSE)


move_downloaded_files <- function(dir_out, dp_id, dp_name, file_pattern, delete_orig = FALSE, unzip = FALSE){
    # from Victoria Scholl: https://github.com/earthlab/neon-veg/blob/894bf3a9a22f1e0c3c9c6bbdfc4c35d30d88fdfa/00-supporting_functions.R#L372
  # moves AOP files downloaded using the neonUtilities::byTileAOP
  # function into a folder with an intuitive name for each 
  # remote sensing data type. 
  # 
  # Args
  #   dir_out - character string path to the main directory where
  #             files were downloaded.
  #             example: "data/data_raw
  #
  #   dp_id - character string with the data product ID 
  #          example: "DP3.30015.001" for the Canopy Height Model
  #
  #   dp_name - character string with a short name describing the data 
  #            this is the folder name where the files will be moved. 
  #            this folder will be within the dir_out folder. 
  #            example: "chm" for Canopy Height Model 
  # 
  #   file_pattern - character string used to identify downloaded files 
  #                  that will be moved by this function. 
  #                  example: "*CHM.tif$" for Canopy Height Model tiles
  #                           that include CHM and end with a .tif extension
  # 
  #   delete_orig - optional logical parameter to delete the original 
  #                 downloaded files. Default value of FALSE 
  # 
  #   unzip - optional logical parameter to unzip the downloaded files
  #           after they have been moved. Default value of FALSE. 
  # 
  # Example function call: 
  # move_downloaded_files(dir_out = "data/data_raw/NIWO_2017"
  #                       ,dp_id = "DP3.30026.001"
  #                       ,dp_name = "veg_indices"
  #                       ,file_pattern = "*VegIndices.zip$"
  #                       ,delete_orig = TRUE
  #                       ,unzip = TRUE)
  
  
  # list all filenames 
  download_path <- file.path(dir_out, dp_id)
  list_files <- list.files(path = download_path
                           ,pattern = file_pattern
                           ,recursive = TRUE
                           ,full.names = TRUE)
  
  # move files into a new folder with an intuitive name
  move_dir <- file.path(dir_out, dp_name)
  if (!dir.exists(move_dir)){dir.create(move_dir)}
  files_from <- list_files
  files_to <- paste(move_dir
                    ,sapply(stringr::str_split(list_files
                                               ,.Platform$file.sep)
                            ,tail, 1)
                    ,sep = .Platform$file.sep)
  # copy the downloaded files into the simplified directory
  file.copy(from = files_from
            ,to = files_to
            ,overwrite = TRUE)
  
  # If the copy was successful, delete original files in the nested directories
  if(delete_orig){
    if(all(file.exists(files_to))){
      unlink(download_path
             ,recursive = TRUE)
    }
  }
  
  # Unzip compressed files 
  if(unzip){
    for(zipFile in files_to){
      utils::unzip(zipFile, exdir = move_dir)
    }
  }
  
  print(paste("Downloaded files moved from", download_path, "to",move_dir))
}



download_by_tile <- function(dpID, site, year, eastings, northings, savepath,
                                buffer=0, check.size=F) {

    neonUtilities::byTileAOP(dpID=dpID, site=site, 
                 year=year, check.size=F,
                 buffer = 0,
                 easting=eastings,
                 northing=northings, 
                 savepath=savepath,
                 include.provisional = FALSE)
}

calc_leaf_area_density <- function(laz_path, dz=0.5, z0=1) {
    # Calculate LAD from voxelization
    laz_data = lidR::readLAS(laz_path)
    lad_profile = lidR::LAD(laz_data@data$Z,
                      dz=dz,
                      z0=z0)
    return(lad_profile)
}

normalize_lidR <- function(las_raw_path, normalized_output_path) {
    message(normalized_output_path)
    ctg <- readLAScatalog(las_raw_path)
    opt_output_files(ctg) <- normalized_output_path
    opt_laz_compression(ctg) <- TRUE
    output <- normalize_height(ctg, tin())
}


clip_lidar_to_polygon_lidR <- function(las_norm_path, shp_path, clipped_output_path) {

    ctg <- readLAScatalog(las_norm_path)
    opt_output_files(ctg) <- paste0(clipped_output_path,"/{plotID}")
    opt_laz_compression(ctg) <- TRUE

    plots_shp <- sf::st_read(shp_path)
    
    # o <- lidR::clip_roi(ctg, plots_shp)

    # Plots
    o <- clip_roi(ctg, plots_shp)
}


clip_raster_to_polygon_lidR <- function(tif_file,shp_path,output) {

    plots_shp <- sf::st_read(shp_path)
    r <- terra::rast(tif_file)
    x <- terra::crop (r,plots_shp)
    y <- terra::mask(x,plots_shp)
    write_sf(y,output)
    
}
