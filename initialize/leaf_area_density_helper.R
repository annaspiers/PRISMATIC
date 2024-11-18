library(lidR)
library(terra)
options(stringsAsFactors = FALSE)

calc_leaf_area_density <- function(laz_path, dz=0.5, z0=1) {
    # Calculate LAD from voxelization
    laz_data = readLAS(laz_path)
    lad_profile = LAD(laz_data@data$Z,
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
    for (idx in 1:nrow(plots_shp)) {
        o <- clip_roi(ctg, plots_shp[idx,])
    }
}


clip_raster_to_polygon_lidR <- function(tif_file,shp_path,output) {

    plots_shp <- sf::st_read(shp_path)
    r <- terra::rast(tif_file)
    x <- terra::crop (r,plots_shp)
    y <- terra::mask(x,plots_shp)
    write_sf(y,output)
    
}
