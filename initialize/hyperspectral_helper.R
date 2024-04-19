# This script is where packages are installed and parameters are defined for
# the other scripts within this workflow.

install.packages("ggbiplot", repos = c("https://cloud.r-project.org"))
library(ggbiplot)
library(stringr)
library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(caret)
library(tidyselect) #all_of
library(parallel) #mclapply



st_erase <- function(x, y) sf::st_difference(x, st_union(st_combine(y))) 
# ^ from https://r-spatial.github.io/sf/reference/geos_binary_ops.html
# keeps only columns from original sf object


list_tiles_with_veg <- function(veg_df, out_dir){
  # written by Scholl et al from https://github.com/earthlab/neon-veg
    
    # generate a list of 1km x 1km tiles containing field data
  # write list to a text file in the output directory provided
  
  print("Generating list of tiles containing veg stems...")
  
  # get easting, northing coordinates of all plants
  e <- NA
  n <- NA
  for (i in 1:dim(veg_df)[1]){
    easting <- floor(as.numeric(veg_df$easting[i])/1000)*1000
    northing <- floor(as.numeric(veg_df$northing[i])/1000)*1000
    e[i] <- easting
    n[i] <- northing
  }
  
  # find unique rows repesenting tiles
  easting_northings <- data.frame(as.character(e),as.character(n))
  colnames(easting_northings) <- c('e','n')
  tiles <- unique(easting_northings[,c('e','n')])
  
  # order by ascending tile coordinates 
  tiles <- tiles %>%
    arrange(e)
    
  tiles$e <- gsub("e+05","00000", tiles$e, fixed=T)
  tiles$n <- gsub("e+05","00000", tiles$n, fixed=T)
  
  # write to text file 
  tile_names <- paste(tiles$e, tiles$n, sep="_")
  tiles_file <- file(file.path(out_dir,"tiles_with_veg.txt"))
  writeLines(tile_names, tiles_file)
  close(tiles_file)
  
  return(tiles)
}



clip_overlap <- function(df, thresh) {
  # modified from Scholl et al from https://github.com/earthlab/neon-veg

  # Each polygon will be compared with the rest to determine which ones are
  # overlapping. For each overlapping pair, the taller polygon clips the
  # shorter one. If the clipped polygon is smaller than the specified area
  # threshold then it is deleted from subsequent analysis.
  #
  # Args
  #   df
  #     (Simple Features object) containing polygons to compare / clip.
  #
  #   thresh
  #     (integer) Area [m^2] threshold to remove small polygons.
  #

  message("\nChecking for overlap & clipping shorter polygons...")

    # how should be polygons be reordered?
    # for now, reorder then based on height.
    polys_ordered <- df[order(df$height,
                              decreasing = TRUE),] %>%
      dplyr::select(-area_m2) %>%
      # create unique ID
      tibble::rownames_to_column() %>%
      mutate(individualID_unique = paste(individualID,rowname,sep="_")) %>%
      relocate(individualID_unique, .before = geometry)

    # create a copy of the polygons to update with clipped/deleted entries
    polys_filtered <- polys_ordered

    # create an empty vector of polygon pairs that have been compared
    compared_pairs <- c()

    for (i in 1:nrow(polys_ordered)){

      individualID_unique <- polys_ordered$individualID_unique[i]

      # if this polygon was removed from the polys_filtered
      # data frame in a previous iteration, skip it
      if(sum(polys_filtered$individualID_unique==individualID_unique) == 0){
        next
      }

      # extract current vs. all other polygon from data frame
      current_poly <- polys_ordered[i,]
      other_polys <- polys_filtered[polys_filtered$individualID_unique!=current_poly$individualID_unique,]
      # other_polys_insubplot <- other_polys[other_polys$plotID.x == current_poly$plotID.x &
      #                                     other_polys$subplotID == current_poly$subplotID,] #ais added this


      # check for overlap between current polygon and all polygons
      # overlap <- raster::intersect(current_poly, other_polys) #crashes here - why?
      overlap_sf <- other_polys[sf::st_intersects(current_poly, other_polys)[[1]], ]
      n_overlap <- nrow(overlap_sf)

      if (n_overlap>0){
        for (o in 1:n_overlap){

          # if current polygon ID is not in filtered set
          if(sum(polys_filtered$individualID_unique==individualID_unique) == 0){
            break
          }

          # get height of current and test polygons that overlap
          current_poly.height <- current_poly$height

          test_poly <- overlap_sf[o,]
          test_poly.height <- test_poly$height

          # combine the ID's of the current and test polygons
          # to keep track of which pairs have been compared
          id.pair <- paste(current_poly$individualID_unique,
                           test_poly$individualID_unique,
                           sep = " ")

          # if polygon pair was already compared, skip to the next pair
          if (id.pair %in% compared_pairs) {

            next

          } else {

            # add to the list of polygon pairs that have been compared
            compared_pairs <- append(compared_pairs,id.pair)

            # add opposite combination of polygons
            id.pair.opp <- paste(test_poly$individualID_unique,
                                 current_poly$individualID_unique,
                                 sep = " ")
            compared_pairs <- append(compared_pairs,id.pair.opp)

          }

          # if test polygon is not in filtered set, skip to the next overlapping polygon
          # (because it has already been removed)
          if (sum(polys_filtered$individualID_unique==test_poly$individualID_unique) == 0){
            next
          }

          # if the area of one of the polygons is equivalent to zero, delete it and
          # skip to the next overlapping polygon.
          if (as.numeric(st_area(current_poly)) < 0.01){

            # delete the current polygon from the polygons_filtered data frame.
            polys_filtered <- polys_filtered[polys_filtered$individualID_unique!=current_poly$individualID_unique,]
            next

          } else if (as.numeric(st_area(test_poly)) < 0.01){
            # delete the test polygon from the polygons_filtered data frame if its
            # area is essentially 0 (using the criterion < 0.01 here). This step
            # was added for instances where after clipping, a small edge fragment remains.
            polys_filtered <- polys_filtered[polys_filtered$individualID_unique!=test_poly$individualID_unique,]
            next
          }

          # compare the heights of the polygons
          if (current_poly.height >= test_poly.height){
            
            #ais add case for when the two are the same height, choose larger stem
            # is this necessary? individualID_unique=="NEON.PLA.D17.SOAP.05666_87"
            
            # if current polygon is taller, clip the test polygon.
            clipped <- st_erase(test_poly, current_poly)

            # if the clipped region is NULL, (tree is obscured completely), 
            # delete from filtered df
            if (nrow(clipped) == 0){
              # remove test from filtered
              polys_filtered <- polys_filtered[polys_filtered$individualID_unique!=test_poly$individualID_unique,]
              
            } else{
              # get the index within the filtered polygons data frame
              # where the test polygon belongs
              j <- which(polys_filtered$individualID_unique == test_poly$individualID_unique)
              # check the area of the clipped test polygon. If it is greater than
              # or equal to the area threshold, replace it as the polygon
              # geometry for the entry matching the test individualID_unique in the
              # polygons_filtered data frame.
              if(as.numeric(st_area(clipped))  >= thresh){
                # replace the original polygon with the clipped polygon
                polys_filtered[j,] <- clipped
              } else {
                # otherwise, delete the test polygon from the polygons_filtered data frame.
                polys_filtered <- polys_filtered[polys_filtered$individualID_unique!=test_poly$individualID_unique,]
              }
            }

          } else {

            # otherwise, the test polygon is taller: clip the current polygon.
            clipped <-  st_erase(current_poly,test_poly)

            if (nrow(clipped) == 0){
              # remove current from filetered from filtered
              polys_filtered <- polys_filtered[polys_filtered$individualID_unique!=current_poly$individualID_unique,]
              
            } else {
              # get the index within the filtered polygons data frame
              # where the current polygon belongs.
              j <- which(polys_filtered$individualID_unique == current_poly$individualID_unique)
              
              if(as.numeric(st_area(clipped))  >= thresh){
                polys_filtered[j,] <- clipped
              } else {
                polys_filtered <- polys_filtered[polys_filtered$individualID_unique!=test_poly$individualID_unique,]
              }
            }
          }
        }
      }
    }

    # write final polygons to file after checking for overlap
    #writeOGR(polys_filtered, getwd(),
    #         paste(shp_filename),
    #         driver="ESRI Shapefile", overwrite_layer = TRUE)

    return(polys_filtered)

}

extract_ind_validation_set <- function(df, percentTrain=0.8) {
  # randomly select <percentTrain> of data from the clipped half diameter 
  # polygons for training, use the remaining samples for validation.
  # keep track of which pixelNumber, easting, northing
  # that these pixels correspond to. 
  
  # percentTrain: randomly select this amount of data for training, 
  # use the rest for validation
  
  #ais create alternatives to random split
  
  # remove any spectra with a height of 0 and remove any factors
  df_val <- df %>% 
    dplyr::filter(chm>0)#ais remove this step once I add in bare ground
  
  # randomly sample rows from this data set 
  set.seed(104)
  
  # the number of sampled rows is calculated based on 
  # percentTrain and the number of rows in the validation set. 
  # percentTrain may have a value like 0.80 (80% data used to train)
  train <- sample(nrow(df_val), 
                  percentTrain*nrow(df_val), 
                  replace = FALSE)
  trainSet <- df_val[train,] 
  valSet <- df_val[-train,]
    
  return (list(trainSet, valSet))
}

get_model_performance_metrics <- function(confusion_matrix) {
  # confusion_matrix is a (true, pred) square matrix
  true_positives <- diag(confusion_matrix)
  false_positives <- colSums(confusion_matrix) - diag(confusion_matrix)
  false_negatives <- rowSums(confusion_matrix) - diag(confusion_matrix)
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  f1 <- 2 * (precision * recall) / (precision + recall)
  tibble(idx = colnames(accuracy$confusion),
         precision = precision, 
         recall = recall, 
         f1 = f1)
}


plot_variable_importance <- function(rf_model, rf_output_dir) {
  # create a figure with two plots, MDA and MDG importance
  
  varImportance <- data.frame(randomForest::importance(rf_model)) 
  varImportance$feature <- rownames(varImportance)
  
  # rename the features for easier interpretation in the variable importance plot
  varImportance$feature_orig <- varImportance$feature
  varImportance$feature[varImportance$feature == "rgb_meanR"] = "RGB red mean"
  varImportance$feature[varImportance$feature == "rgb_meanG"] = "RGB green mean"
  varImportance$feature[varImportance$feature == "rgb_meanB"] = "RGB blue mean"
  varImportance$feature[varImportance$feature == "rgb_sdR"] = "RGB red sd"
  varImportance$feature[varImportance$feature == "rgb_sdG"] = "RGB green sd"
  varImportance$feature[varImportance$feature == "rgb_sdB"] = "RGB blue sd"
  #ais where did scholl et al create a combined variable rgb_mean_sd_R?
  
  varImportance <- varImportance %>% 
    dplyr::select(feature, MeanDecreaseAccuracy, MeanDecreaseGini, everything())
  varImportanceMDA <- varImportance %>% dplyr::arrange(desc(MeanDecreaseAccuracy))
  varImportanceMDG <- varImportance %>% dplyr::arrange(desc(MeanDecreaseGini))
  
  # MDA importance bar plot 
  mda_plot <- ggplot(data = varImportanceMDA, aes(x = reorder(feature, MeanDecreaseAccuracy), 
                                                  y = MeanDecreaseAccuracy
                                                  #,fill = MeanDecreaseAccuracy
  )) + 
    geom_bar(stat = 'identity', color = "black", linewidth = 0.1, width = 0.5, show.legend = FALSE) + 
    labs(x = "AOP-derived feature\n", 
         y = "Mean Decrease in Accuracy") +
    # x-axis label gets cut off otherwise after coord_flip
    ylim(0, max(varImportanceMDA$MeanDecreaseAccuracy) + 10) + 
    coord_flip() + 
    theme_bw() + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16)) +
  ggtitle("AOP-derived feature importance")
  
  # MDG importance bar plot 
  mdg_plot <- ggplot(data = varImportanceMDG, 
                     aes(x = reorder(feature, MeanDecreaseGini), 
                         y = MeanDecreaseGini)) + 
    geom_bar(stat = 'identity', color = "black", linewidth = 0.1, width = 0.5, show.legend = FALSE) + 
    labs(x = "",  # no y-axis label since it's the same as the MDA plot on the left
         y = "Mean Decrease in Gini") + 
    # x-axis label gets cut off otherwise after coord_flip
    ylim(0, max(varImportanceMDA$MeanDecreaseGini) + 10) + 
    coord_flip() + 
    theme_bw() + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16)) 
  
  # generate the plot to view in RStudio
  gridExtra::grid.arrange(mda_plot, 
                          mdg_plot, 
                          nrow = 1
                          # ,top = "Variable Importance" # don't need main title for figure in manuscript 
  )
  
  # write the plot to an image file
  var_imp_plot <- gridExtra::arrangeGrob(mda_plot, 
                              mdg_plot, 
                              nrow = 1) 
  ggsave(filename = file.path(rf_output_dir, "variable_importance.png"), 
         plot = var_imp_plot, width = 10, units = "in", dpi = 500)
}



filter_out_wavelengths <- function(wavelengths, layer_names){
    # define the "bad bands" wavelength ranges in nanometers, where atmospheric 
    # absorption creates unreliable reflectance values. 
    bad_band_window_1 <- c(1340, 1445)
    bad_band_window_2 <- c(1790, 1955)
    #ais ^hardcoded here probably isn't the best
    
    # remove the bad bands from the list of wavelengths 
    remove_bands <- wavelengths[(wavelengths > bad_band_window_1[1] & 
                                     wavelengths < bad_band_window_1[2]) | 
                                    (wavelengths > bad_band_window_2[1] & 
                                         wavelengths < bad_band_window_2[2])]
    
    # Make sure printed wavelengths and stacked AOP wavelengths match
    if (!all(paste0("X", as.character(base::round(wavelengths))) == 
             layer_names[1:length(wavelengths)])) {
        message("wavelengths do not match between wavelength.txt and the stacked imagery")
    }
    
    # create a LUT that matches actual wavelength values with the column names,
    # X followed by the base::rounded wavelength values. 
    # Remove the rows that are within the bad band ranges. 
    wavelength_lut <- data.frame(wavelength = wavelengths,
                                 xwavelength = paste0("X", as.character(base::round(wavelengths))),
                                 stringsAsFactors = FALSE) %>% 
        filter(!wavelength %in% remove_bands)
    
    return(wavelength_lut)
}



prep_features_for_RF <- function(extracted_features_filename, featureNames) {
    # read the spectra values extracted from the data cube 
    df <- read.csv(extracted_features_filename) %>%
        #dplyr::select(-rowname) %>%
        mutate(dfID = paste(eastingIDs, northingIDs, 
                             pixelNumber, sep = "_") ) %>%
        # remove gbase::round pixels ais
        filter(chm>0) %>% 
        # also reset the factor levels (in case there are dropped pft levels)
        droplevels()
    
    features_df <- df %>% 
        # filter the data to contain only the features of interest 
        dplyr::select(shapeID, all_of(featureNames))
    
    return(features_df)
}



apply_PCA <- function(df, wavelengths, output_dir, nPCs=3) {

  #ais find a way to automate determining nPCs
    # elbow plot
    # factoextra::fviz_eig(hs_pca)
    # features_final <- cbind(features, hs_pca$x[,1:nPCs]) %>% # add first n PCs to features df
    #     # create unique ID for each pixel in the input imagery,
    #     mutate(dfIDs = paste(plotID, pixelNumber, eastingIDs, northingIDs, 
    #                          sep = "_")) %>%
    #     dplyr::select(c(plotID, subplotID, chm, slope, aspect_cat, ARVI, EVI, NDVI,
    #                     PRI, SAVI, rgb_meanR, rgb_meanG, rgb_meanB, rgb_sdR, rgb_sdG,
    #                     rgb_sdB, rgb_mean_sd_R, rgb_mean_sd_G, rgb_mean_sd_B,
    #                     PC1, PC2, PC3, dfIDs, pixelNumber)) 
    # ais ^ find a way to automate identifying nPCs rather than hardcoding it

    # remove the individual spectral reflectance bands from the training data
    features_noWavelengths <- df %>% dplyr::select(-c(wavelengths$xwavelength))
    
    # PCA: calculate Principal Components 
    hs <- df %>% dplyr::select(c(wavelengths$xwavelength)) %>% as.matrix()
    hs_pca <- stats::prcomp(hs, center = TRUE, scale. = TRUE)
    # summary(hs_pca) ais print pc summary info somewhere?
    # add first n PCs to features data frame
    features <- cbind(features_noWavelengths, hs_pca$x[,1:nPCs]) 
    
    # visualize where each sample falls on a plot with PC2 vs PC1 
    ggbiplot::ggbiplot(hs_pca,
                       choices = 1:2, # which PCs to plot
                       obs.scale = 1, var.scale = 1, # scale observations & variables
                       var.axes=FALSE, # remove arrows
                       groups = df$pft, # color the points by PFT
                       ellipse = TRUE, # draw ellipse abase::round each group
                       circle = TRUE ) + # draw circle abase::round center of data set
        ggtitle("PCA biplot, PC1 and PC2") +
        scale_color_brewer(palette="Spectral") +
        theme_bw()
    # save to file
    ggplot2::ggsave(file.path(output_dir,"pcaPlot.pdf"))
    
    return(features)
}

stack_hyperspectral <- function(h5, out_dir) {
    # written by Scholl et al from https://github.com/earthlab/neon-veg

    # this function is called in prep_aop_imagery

    # This function creates a stacked raster object for the specified HDF5 
    # filename. 
    #
    # Args: 
    # h5
    #   character string filename of HDF5 file 
    # out_dir
    #   directory for output files where the wavelengths will be written to
    #   a text file for further analysis
    #
    # Returns:
    # s
    #   Stacked raster (collection of raster layers with the same spatial extent
    #   and resolution) containing nrows x ncols x nbands, based on the 
    #   resolution and number of bands in the HDF5 file. This stacked raster
    #   can then be clipped using Spatial vector layers. 
    #

    # list the contents of HDF5 file
    h5_struct <- rhdf5::h5ls(h5, all=T)
    
    # construct the string using "/Reflectance/Metadata/Coordinate_System",
    # without explicitly using a site code 
    crs_tag <- h5_struct$group[grepl("/Reflectance/Metadata/Coordinate_System", 
                                     h5_struct$group)][1] 
    
    # read coordinate reference system data
    crs_info <- rhdf5::h5read(h5, crs_tag)
    
    # convert "UTM" to lowercase "utm" for proper usage later
    crs_info$Proj4 <- sf::st_crs(chartr("UTM", "utm", crs_info$Proj4))
    
    # get attributes for the Reflectance dataset.
    # construct the string using "/Reflectance/Reflectance_Data"" 
    refl_tag <- paste0(h5_struct$group[grepl("/Reflectance", 
                                             h5_struct$group)][1],
                       "/Reflectance_Data")
    
    # read the reflectance metadata
    refl_info <- rhdf5::h5readAttributes(h5,refl_tag)
    
    # get the dimensions of the reflectance data
    n_rows <- refl_info$Dimensions[1]
    n_cols <- refl_info$Dimensions[2]
    n_bands <- refl_info$Dimensions[3]
    
    # print dimensions 
    # print(paste0("# Rows: ", as.character(n_rows)))
    # print(paste0("# Columns: ", as.character(n_cols)))
    # print(paste0("# Bands: ", as.character(n_bands)))
    
    # read the wavelengths of the hyperspectral image bands
    wavelength_tag <- paste0(h5_struct$group[grepl("/Reflectance/Metadata/Spectral_Data", 
                                                   h5_struct$group)][1],
                             "/Wavelength")
    wavelengths <- rhdf5::h5read(h5,
                                 wavelength_tag)
    
    # define spatial extent: extract resolution and origin coordinates
    map_info <- unlist(strsplit(crs_info$Map_Info, 
                                split = ", "))
    res_x <- as.numeric(map_info[6])
    res_y <- as.numeric(map_info[7])
    x_min <- as.numeric(map_info[4])
    y_max <- as.numeric(map_info[5])
    
    # calculate the maximum X and minimum Y values 
    x_max <- (x_min + (n_cols * res_x))
    y_min <- (y_max - (n_rows * res_y))
    tile_extent <- terra::ext(x_min, x_max, y_min, y_max) 
    
    # read reflectance data for all bands
    refl <- rhdf5::h5read(h5, refl_tag,
                          index = list(1:n_bands, 1:n_cols, 1:n_rows))
    
    # view and apply scale factor to convert integer values to reflectance [0,1]
    # and data ignore value
    scale_factor <- refl_info$Scale_Factor
    data_ignore <- refl_info$Data_Ignore_Value
    refl[refl == data_ignore] <- NA 
    refl_scaled <- refl / scale_factor
    
    # create georeferenced raster using band 1 
    # convert first band to matrix
    # transpose the image pixels for proper orientation to match
    # the other layers. create a raster for this band and assign
    # the CRS.
    # print("Transposing reflectance data for proper orientation")
    r1 <- terra::t(terra::rast(refl_scaled[1,,], crs = crs_info$Proj4))
    terra::ext(r1) <- tile_extent
    
    # start the raster stack with first band 
    s <- r1
    
    # loop through bands and create a giant rasterstack with 426 (n_bands) bands
    for(b in 2:n_bands){
        
        # create raster with current band
        
        r <- terra::t(terra::rast(refl_scaled[b,,], crs = crs_info$Proj4)) 
        terra::ext(r) <- tile_extent
        
        # add additional band to the stack with the addLayer function
        s <- c(s,r) 
        
    }
    
    # adjust the names for each layer in raster stack to correspond to wavelength
    names(s) <- base::round(wavelengths)
    
    # write wavelengths to a text file 
    # write the exact wavelengths to file for future use 
    write.table(data.frame(wavelengths = wavelengths),
                file.path(out_dir,"wavelengths.txt"),
                row.names=FALSE)
    
    # return the stacked hyperspectral data to clip with vector files 
    return(s)
}



prep_aop_imagery <- function(site, year, hs_L3_path, tif_path, data_int_path) {
    # modified from Scholl et al from https://github.com/earthlab/neon-veg
    
    # output directory for stacked AOP data
    stacked_aop_data_dir <- file.path({data_int_path}, {site}, {year}, "stacked_aop")
    if (!dir.exists(stacked_aop_data_dir)) {
      dir.create(stacked_aop_data_dir) 
    } 

    # ais change prep_aop_imagery to stack only 1km x 1km tiles
    # ais add plotting aop imagery from script 6
    
    # list the files in each data directory; filter results based on file type
    h5_ls <- list.files(path = file.path({hs_L3_path}),
                        pattern = "reflectance.h5", recursive=T, full.names=T)
    chm_ls <- list.files(path = file.path({tif_path}),
                         pattern = "_CHM.tif", recursive = T, full.names = T)
    slope_ls <- list.files(path = file.path({tif_path}),
                           pattern = "_slope.tif", recursive=T, full.names=T)
    aspect_ls <- list.files(path = file.path({tif_path}),
                            pattern = "_aspect.tif", recursive=T, full.names=T)
    rgb_ls <- list.files(path = file.path({tif_path}),
                         pattern = "_image.tif", recursive=T, full.names=T)
    savi_ls <- list.files(path = file.path({tif_path}),
                          pattern = "_SAVI.tif", recursive=T, full.names=T)
    pri_ls <- list.files(path = file.path({tif_path}),
                         pattern = "_PRI.tif", recursive=T, full.names=T)
    ndvi_ls <- list.files(path = file.path({tif_path}),
                          pattern = "_NDVI.tif", recursive=T, full.names=T)
    evi_ls <- list.files(path = file.path({tif_path}),
                         pattern = "_EVI.tif", recursive=T, full.names=T)
    arvi_ls <- list.files(path = file.path({tif_path}),
                          pattern = "_ARVI.tif", recursive=T, full.names=T)
    
    # create data cubes with AOP-derived features for each tile
    for (h5_path in h5_ls) {
        
        # each current hyperspectral tile must be read and stacked into a
        # georeferenced rasterstack (so it can be clipped with fpoint/polygon
        # shapefiles). The process of creating a rasterstack takes a while for
        # each tile, so after creating each rasterstack once, each object gets
        # written to a file.
        
        # Build up the rasterstack filename by parsing out the easting/northing
        # coordinates from the current h5 filename.
        
        # parse the UTM easting and northing values from the current h5 filename
        easting <- stringr::str_split(tail(stringr::str_split(h5_path, "/")[[1]],n=1),"_")[[1]][5]
        northing <- stringr::str_split(tail(stringr::str_split(h5_path, "/")[[1]],n=1),"_")[[1]][6]
        # combine them with an underscore; use this to find corresponding tiles
        # of various remote sensing data
        east_north_string <- paste0(easting,"_",northing)
        
        message(paste0("Stacking ",east_north_string,", AOP tile ",match(h5_path,h5_ls)," out of ",length(h5_ls)))

        # generate a filename for the stacked AOP data
        stacked_aop_data_filename = file.path(stacked_aop_data_dir,
                                              paste0("stacked_aop_data_",
                                                     east_north_string, ".tif"))
        if (file.exists(stacked_aop_data_filename)) {
            
            # if it exists, read that instead of re-generating the same rasterstack.
            message("stacked_aop_data already created for current tile.")
            
            # restore / read the rasterstack from file
            next
            
        } else {
            
            # if it doesn't exist, create the features from the aop data to file 
            # hyperspectral and lidar features ------------------------------------
            
            # read the corresponding remote sensing data layers for current tile
            if (length(grep(east_north_string, chm_ls, value=TRUE))==0) {
                message("CHM does not exist for current tile.")
            } else {
                chm <- terra::rast(grep(east_north_string, chm_ls, value=TRUE)) #problem line
                h5 <- stack_hyperspectral(h5_path, stacked_aop_data_dir)
                
                if (terra::ext(h5) != terra::ext(chm)) {
                    message("different extents of h5 and chm")
                    next
                }
                
                slope <- terra::rast(grep(east_north_string, slope_ls, value=TRUE))
                aspect <- terra::rast(grep(east_north_string, aspect_ls, value=TRUE))
                
                thresh_mat <- matrix(c(0,45,0, 45,135,1, 135,225,2, 225,315,3, 315,360,0),
                                     ncol=3, byrow=TRUE)
                aspect_cat <- terra::classify(x=aspect, rcl=thresh_mat)
                
                # Vegetation indices
                savi <- terra::rast(grep(east_north_string, savi_ls, value=TRUE))
                pri <- terra::rast(grep(east_north_string, pri_ls, value=TRUE))
                ndvi <- terra::rast(grep(east_north_string, ndvi_ls, value=TRUE))
                evi <- terra::rast(grep(east_north_string, evi_ls, value=TRUE))
                arvi <- terra::rast(grep(east_north_string, arvi_ls, value=TRUE))
                
                # set the raster name for each layer to be simply the name of the data 
                # (i.e. "aspect") as opposed to the full filename 
                # (i.e. ""NEON_D13_NIWO_DP3_452000_4431000_aspect")
                
                names(chm) <- "chm"
                names(slope) <- "slope"
                names(aspect_cat) <- "aspect_cat"
                names(savi) <- "savi"
                names(pri) <- "pri"
                names(ndvi) <- "ndvi"
                names(evi) <- "evi"
                names(arvi) <- "arvi"
                
                # RGB features
                # The RGB data tile has 10,000 x 10,000 pixels, 10cm res
                # All other layers tiles have 1,000 x 1,000 pixels, 1 meter res
                # Aggregate red, green, blue intensities within each coarser grid cell using
                # statistics such as mean and standard deviation.
                
                # rgb data has 3 bands. read each one individually
                # rgb data has 3 bands. read each one individually
                rgb_rast <- terra::rast(grep(east_north_string, rgb_ls, value = T))
                rgb_red <- rgb_rast[[1]]
                rgb_green <- rgb_rast[[2]]
                rgb_blue <- rgb_rast[[3]]
                
                # The "fact" parameter of the raster::aggregate function is the number of cells
                # in each direction (horizontal and vertically) to aggregate across.
                # Since the RGB data has a spatial resolution that is 1/10th of the 
                # other data layers (10cm compared to 1m), fact should be 10 to produce 
                # an output raster with 1000 x 1000 pixels. 
                # mean intensity per 1m x 1m grid cell 
                rgb_meanR <- terra::aggregate(rgb_red, fact = 10, fun = mean)
                rgb_meanG <- terra::aggregate(rgb_green, fact = 10, fun = mean)
                rgb_meanB <- terra::aggregate(rgb_blue, fact = 10, fun = mean)
                
                # standard deviation of intensity per 1m x 1m grid cell 
                rgb_sdR <- terra::aggregate(rgb_red, fact = 10, fun = sd)
                rgb_sdG <- terra::aggregate(rgb_green, fact = 10, fun = sd)
                rgb_sdB <- terra::aggregate(rgb_blue, fact = 10, fun = sd)
                
                # set the names of each layer to reflect the metric it contains
                names(rgb_meanR) <- "rgb_meanR" # mean intensity
                names(rgb_meanG) <- "rgb_meanG"
                names(rgb_meanB) <- "rgb_meanB"
                names(rgb_sdR) <- "rgb_sdR"     # standard deviation of intensity
                names(rgb_sdG) <- "rgb_sdG"
                names(rgb_sdB) <- "rgb_sdB"
                
                # stack up all the RGB features
                rgb_features <- c(rgb_meanR, rgb_meanG, rgb_meanB,
                                  rgb_sdR, rgb_sdG, rgb_sdB)
                
                # Create the pixel number grid as a layer to add to the data cube.
                # this one keeps track of individual pixel ID's to avoid duplicate
                # spectra being extracted. Basically, assign an integer ID to each
                # pixel in the 1000x1000 raster. This raster needs to have the same
                # dimensions, extent, and crs as the other layers so they can be stacked
                # together. create a vector of IDs from 1 to the number of pixels in one band
                pixelID <- 1:(dim(chm)[1]*dim(chm)[2]) #1:(nrow(s) * ncol(s))
                # add tile east, north coordinates 
                #pixelID <- paste(pixelID, east_north_string, sep="_")
                # reshape this 1D vector into a 2D matrix 
                dim(pixelID) <- c(dim(chm)[1],dim(chm)[2]) #c(nrow(s),ncol(s))
                # create a raster layer of pixel numbers 
                pixelNumbers <- terra::rast(pixelID, crs = terra::crs(chm)) 
                terra::ext(pixelNumbers) <- terra::ext(chm) 
                names(pixelNumbers) <- "pixelNumber"
                
                # Create similar layers to keep track of the tile where the pixel is located
                eastingID <- rep(as.numeric(easting), times = (dim(chm)[1] * dim(chm)[2]))
                northingID <- rep(as.numeric(northing), times = (dim(chm)[1] * dim(chm)[2]))
                # reshape to be two-dimensional
                dim(eastingID) <- c(dim(chm)[1],dim(chm)[2])
                dim(northingID) <- c(dim(chm)[1],dim(chm)[2])
                # create rasters to contain the easting and northing values
                eastingIDs <- terra::rast(eastingID, crs = terra::crs(chm))
                northingIDs <- terra::rast(northingID, crs = terra::crs(chm))
                # assign extent and CRS to match the other layers in the stack
                terra::ext(eastingIDs) <- terra::ext(chm)
                terra::ext(northingIDs) <- terra::ext(chm)
                names(eastingIDs) <- "eastingIDs"
                names(northingIDs) <- "northingIDs"
                
                # now, all the remote sensing data have been read in for the current
                # tile. add each one to the hyperspectral data stack along with the
                # layer to keep track of pixel number within the tile.
                stacked_aop_data <- c(h5, chm, slope, aspect_cat, savi,
                                      pri, ndvi, evi, arvi, 
                                      rgb_features, pixelNumbers, 
                                      eastingIDs, northingIDs)
                terra::crs(stacked_aop_data) <- terra::crs(chm)
                
                # save the stacked AOP data to file for easy clipping later
                # terra raster allows layer names are saved with the .tif
                terra::writeRaster(stacked_aop_data, file = stacked_aop_data_filename)
            }
        }
    }
    
    # Save stacked layer names separately as metadata
    write.table(data.frame(stacked_aop_layer_names = names(stacked_aop_data)),
                file.path(stacked_aop_data_dir,"stacked_aop_layer_names.txt"),
                row.names=FALSE)
    
    return("{stacked_aop_data_dir}")
}



create_tree_crown_polygons <- function(site, year, biomass_path, data_out_path, px_thresh) {
  # modified from Scholl et al from https://github.com/earthlab/neon-veg
  
  # Create geospatial features (points, polygons with half the maximum crown diameter) 
  # for every tree in the NEON woody vegetation data set that intersect with independent 
  # pixels in the AOP data 
  # Analog to 02-create_tree_features.R and 03-process_tree_features.R
  # from https://github.com/earthlab/neon-veg
  
  # px_thresh is the area threshold (# of of 1m pixels) to filter out small tree crowns
  
  # output directory for training data
  training_data_dir <- file.path({data_out_path}, {site}, {year}, "training")
  if (!dir.exists(training_data_dir)) {
    dir.create(training_data_dir)
  }
  
  live_trees_path <- file.path({biomass_path},"pp_veg_structure_IND_IBA_IAGB_live.csv")
  
  # Load cleaned inventory data for site and year
  # ais this is live trees only. should figure out how to incorporate dead trees
  veg_ind <- read.csv(live_trees_path) %>%
    dplyr::rename(easting = adjEasting, 
                  northing = adjNorthing,
                  # Rename columns so they're unique when truncated when saving shp 
                  lat = adjDecimalLatitude,
                  long = adjDecimalLongitude,
                  elev = adjElevation,
                  elevUncert = adjElevationUncertainty,
                  indStemDens = individualStemNumberDensity,
                  indBA = individualBasalArea,
                  dendroH = dendrometerHeight,
                  dendroGap = dendrometerGap,
                  dendroCond = dendrometerCondition,
                  bStemDiam = basalStemDiameter,
                  bStemMeasHgt = basalStemDiameterMsrmntHeight,
                  sciNameFull = scientificName,
                  sp_gen = scientific,
                  wood_dens1 = wood.dens,
                  wood_dens2 = wood_dens) #ais what is the source of columns in this df? NEON or our own changes?
  
  # Remove all rows with missing an easting or northing coordinate
  veg_has_coords <- veg_ind[complete.cases(veg_ind[, c("northing", "easting")]),]
  
  # Keep only the entries with crown diameter and tree height
  # since these measurements are needed to create and compare polygons.
  # ais future work: include ninetyCrownDiameter?
  veg_has_coords_size <- veg_has_coords[complete.cases(veg_has_coords$height) & 
                                          complete.cases(veg_has_coords$maxCrownDiameter),] %>%
    dplyr::filter(maxCrownDiameter < 50) %>%
    dplyr::filter(individualID != "NEON.PLA.D17.SOAP.04659") %>% #no tree here
    dplyr::filter(individualID != "NEON.PLA.D17.SOAP.05847") %>% #no tree here
    dplyr::filter(individualID != "NEON.PLA.D17.SOAP.04256")  #no tree here
    #"NEON.PLA.D17.SOAP.04972" #about 1/3 is ground
  
  # Assign a PFT to each individual
  # ais need to make this an exhaustive match for all sites, not just SOAP
  # ais use neon trait table from Marcos?
  veg_training <- veg_has_coords_size %>%
    left_join(pft_reference %>% dplyr::select(-growthForm) , by=join_by(siteID, taxonID)) %>%
    dplyr::mutate(pft =  ifelse(!is.na(pft),pft,ifelse(growthForm == "small shrub" | growthForm == "single shrub" ,"montane_shrub",
                                  ifelse(taxonID=="PIPO" | taxonID=="PILA" | taxonID=="PINUS", "pine", 
                                          ifelse(taxonID=="CADE27", "cedar", 
                                                ifelse(taxonID=="ABCO" | taxonID=="ABMA" , "fir", 
                                                        ifelse(taxonID=="QUCH2" | taxonID=="QUKE"| taxonID=="QUERC"|
                                                               taxonID=="QUWI2", "oak", 
                                                              ifelse(taxonID=="ARVIM" | taxonID=="CECU" | taxonID=="CEMOG" | 
                                                                          taxonID=="CEIN3" | taxonID=="RIRO" | taxonID=="FRCA6" | 
                                                                          taxonID=="RHIL" | taxonID=="AECA" | taxonID=="SANI4" | 
                                                                          taxonID=="CEANO", "montane_shrub", "other")))))))) 
  
  write.csv(veg_training, file = file.path(training_data_dir, "vst_training.csv"))
  
  coord_ref <- sf::st_crs(32611)
  #  ais ^ will need to change this from being hardcoded to finding
  # the UTM zone that any data are associated
  
  # Convert the data frame to an SF object 
  veg_training_pts_sf <- sf::st_as_sf(veg_training,
                                      coords = c("easting", "northing"),
                                      crs = coord_ref)
  # sf::st_write(obj = veg_training_pts_sf,
  #              dsn = file.path(training_data_dir, "veg_points_all.shp"),
  #              delete_dsn = TRUE)
  
  # Generate a list of AOP tiles that overlap with the inventory data
  tiles <- list_tiles_with_veg(veg_df = veg_training,
                               out_dir = training_data_dir)
  
  # # plot all mapped stem locations
  # ggplot2::ggplot() +
  #   ggplot2::geom_sf(data = mapped_stems_sf) +
  #   ggplot2::ggtitle("Map of Stem Locations")
  
  # Write shapefile with CIRCULAR POLYGONS for all mapped stems with 
  # height & crown diameter. Size: 1/2 Maximum crown diameter
  polys_half_diam <- sf::st_buffer(veg_training_pts_sf,
                                   dist = base::round((veg_training_pts_sf$maxCrownDiameter/4),
                                                digits = 1))
  # sf::st_write(obj = polys_half_diam,
  #     dsn = file.path(training_data_dir, "veg_polygons_half_diam.shp"),
  #     delete_dsn = TRUE)
  
  # Apply area threshold to remove small polygons (smaller than 4 1m pixels)---------------------------
  # 2. An area threshold is then applied to remove any small trees area values 
  # less than the area of four hyperspectral pixels. This threshold 
  # was selected with the coarser resolution of the hyperspectral and LiDAR data 
  # products in mind. By preserving the trees with larger crowns, it is believed 
  # that purer spectra will be extracted for them for the training data sets, 
  # as opposed to extracting mixed pixels which signal from smaller plants and 
  # backgbase::round materials or neighboring vegetation. 
  
  # Define the area threshold in units of [m^2]
  # Area of RGB pixel, 10cm by 10cm
  px_area_rgb <- 0.1 * 0.1 #[m^2]
  # Gridded LiDAR products and HS pixels are 1m x 1m
  px_area_hs <- 100 * px_area_rgb #[m^2]
  # Multiply area of 1 HS pixel by the number of pixels; 
  # as defined in the main script by the px_thresh variable
  thresh <- px_area_hs * px_thresh #[m^2]
  
  # Remove half-diam polygons with area < n hyperspectral pixels 
  polygons_thresh_half_diam <- polys_half_diam %>% 
    dplyr::mutate(area_m2 = sf::st_area(polys_half_diam)) %>% 
    dplyr::filter(as.numeric(area_m2) > thresh)
  
  # Clip shorter polygons with taller polygons -----------------------------
  
  polygons_clipped_half_diam <- clip_overlap(polygons_thresh_half_diam, thresh)
  
  # Check and fix/delete invalid geometries.
  # Remove invalid geometries if the reason
  # for invalidity is "Too few points in geometry component[453729.741 4433259.265]"
  # since there are too few points to create a valid polygon. 
  # Use the reason = TRUE paramater in sf::st_isvalid to see the reason. 
  polygons_clipped_half_diam_is_valid <- sf::st_is_valid(x = polygons_clipped_half_diam)
  
  # polygons_clipped_valid <- lwgeom::st_make_valid(polygons_clipped)
  polygons_clipped_half_diam_valid <- polygons_clipped_half_diam %>% 
    dplyr::filter(polygons_clipped_half_diam_is_valid)
  
  
  # Write shapefile with clipped tree crown polygons half the max crown diameter 
  training_shp_path <- file.path(training_data_dir, "tree_crowns_training.shp")
  sf::st_write(obj = polygons_clipped_half_diam_valid,
               dsn = training_shp_path,
               delete_dsn = TRUE)
  
  return(training_shp_path)
}



extract_spectra_from_polygon <- function(site, year, data_int_path, data_final_path,
                                        stacked_aop_path, shp_path, use_case, ic_type=NA, 
                                         ic_type_path=NA, aggregate_from_1m_to_2m_res) {
    # modified from Scholl et al from https://github.com/earthlab/neon-veg
    
    # Extract features (remote sensing data) for each sample (pixel) within the
    # specified shapefile (containing polygons that correspond to trees at the NEON site)
    # Analog to 07-extract_training_features.R from https://github.com/earthlab/neon-veg
    
    # ais make the countdown more useful (e.g., x/y plots remaining to extract)
    # ais when use_case="predict", make the temp csv files save to data_final_path instead of training_path 
    
    # use_case is "train" or "predict"

    assertthat::assert_that(length(list.files(stacked_aop_path))>0,
      msg=paste0("Cannot extract data because stacked_aop_path is empty for ",site," - ",year) )
    
    training_data_dir <- file.path({data_int_path}, {site}, {year}, "training")
    
    # clip data cube - extract features
    # get a description of the shapefile to use for naming outputs
    shapefile_description <- tools::file_path_sans_ext(basename(shp_path))
    
    # Specify destination for extracted features
    if (use_case == "train") {
        extracted_features_path <- file.path(training_data_dir)
        extracted_features_filename <- file.path(extracted_features_path,
                                                 paste0(shapefile_description, 
                                                        "-extracted_features_inv.csv"))
        
    } else if (use_case=="predict") {
        extracted_features_path <- file.path(ic_type_path)
        extracted_features_filename <- file.path(extracted_features_path,
                                                     paste0(shapefile_description, 
                                                            "-extracted_features.csv"))
    } else {
        message("need to specify use_case")
    }
    
    # Only run the whole function if the file does not exist
    if (!file.exists(extracted_features_filename)) {
        
        # Extract features (remote sensing data) for each sample (plot) within the 
        # specified shapefile (random 20m x 20m plot, FATES patch)
        
        # loop through stacked AOP data and extract features
        stacked_aop_list <- list.files(stacked_aop_path, 
                                       full.names = TRUE, pattern=".tif")
        
        shp_sf <- sf::st_read(shp_path)
        shp_coords <- shp_sf %>% # add columns for the center location of each plot
            sf::st_centroid() %>%  # get the centroids first for polygon geometries
            sf::st_coordinates() %>%
            as.data.frame()
        
        # add new columns for the plot center coordinates
        shp_sf$center_X <- shp_coords$X
        shp_sf$center_Y <- shp_coords$Y
        
        # create column to track shape ID 
        # if training, this is the tree crown boundary
        # if predicting, this is the plot boundary)
        if (use_case=="train") {
            shp_sf$shapeID <- paste0("tree_crown_", rownames(shp_sf))    
        } else if (use_case=="predict" ){
            if ("PLOTID" %in% colnames(shp_sf)) {
                shp_sf <- shp_sf %>%
                    dplyr::rename(shapeID=PLOTID)
            } else {
                shp_sf <- shp_sf %>%
                    dplyr::rename(shapeID=plotID)
            }
        }
        
        # Filter to tiles containing veg to speed up the next for-loop
        if (use_case=="train" | ic_type=="rs_inv_plots") {
            tiles_with_veg <- read.table(file.path(training_data_dir, "tiles_with_veg.txt")) %>%
                pull(V1)
              
            stacked_aop_list <- stacked_aop_list[stringr::str_detect(stacked_aop_list, tiles_with_veg %>% paste(collapse = "|"))]
        }
        
        # loop through AOP tiles 
        for (stacked_aop_filename in stacked_aop_list) { #12 index
            
            # read current tile of stacked AOP data 
            stacked_aop_data <- terra::rast(stacked_aop_filename) #old raster::stack
            
            # construct the easting northing string for naming outputs
            east_north_string <- paste0(stacked_aop_data$eastingIDs[1], "_",
                                        stacked_aop_data$northingIDs[1])
            
            # If csv file already exists, skip
            if (file.exists(file.path(training_data_dir,paste0("extracted_features_",
                                                               east_north_string, "_",shapefile_description,".csv")))) {
                next
            }
            
            # # Save spatial extent of tile
            # ext <- st_as_sf(as(extent(stacked_aop_data),"SpatialPolygons")) 
            # st_crs(ext) <- st_crs(stacked_aop_data) 
            # write_sf(ext,file.path(training_data_dir,"stacked_aop_extent",
            #                       paste0(east_north_string,".shp")))
            
            # figure out which plots are within the current tile by comparing each
            # X,Y coordinate to the extent of the current tile 
             shapes_in <- shp_sf %>% 
                dplyr::filter(center_X >= terra::ext(stacked_aop_data)[1] & #old raster::extent
                                  center_X < terra::ext(stacked_aop_data)[2] & 
                                  center_Y >= terra::ext(stacked_aop_data)[3] & 
                                  center_Y  < terra::ext(stacked_aop_data)[4]) 
            
            # if no polygons are within the current tile, skip to the next one
            if (nrow(shapes_in)==0){
                message("no shapes located within current tile... skipping to next shapefile")
                next
            } else {
                
                if (use_case == "train") {
                    message(paste0("Extracting ",nrow(shapes_in), 
                                   " trees in current tile, ",east_north_string))
                } else {
                    message(paste0("Extracting ",nrow(shapes_in), 
                                   " plots in current tile, ",east_north_string)) 
                } 
            }
            
            # Aggregate to 2m resolution from 1m
            if (aggregate_from_1m_to_2m_res == T) {
                stacked_aop_data <- raster::aggregate(stacked_aop_data,2) 
                # takes 5min for a 1kmx1km tile
            }
            
            # clip the hyperspectral raster stack with the polygons within current tile.
            # the returned objects are data frames, each row corresponds to a pixel in the
            # hyperspectral imagery. The ID number refers to which tree that the 
            # the pixel belongs to. A large polygon will lead to many extracted pixels
            # (many rows in the output data frame), whereas tree stem points will
            # lead to a single extracted pixel per tree. 
            
            shapes_in_spv <- terra::vect(shapes_in)
            extracted_spectra <- terra::extract(stacked_aop_data, shapes_in_spv, ID=TRUE) %>%
                left_join(data.frame(shapeID = shapes_in_spv$shapeID, 
                                     ID = 1:nrow(shapes_in_spv))) %>%
                dplyr::select(-ID)
            
            # VS-NOTE TO DO: ais
            # adjust this extract step to only get pixels WITHIN each tree polygon,
            # also try calculating the percentage that each pixel is within a polygon
            # and keep only pixels with > 50% overlap 
            
            # merge the extracted spectra and other data values with the tree info 
            shapes_metadata <- tibble(shapes_in)
            
            # combine the additional data with each spectrum for writing to file.
            # remove the geometry column to avoid issues when writing to csv later 
            spectra_write <- merge(shapes_metadata,
                                   extracted_spectra,
                                   by="shapeID") %>% 
                dplyr::select(shapeID, everything()) %>% 
                dplyr::select(-geometry)
            #ais find where subplotID is assigned
            
            # write extracted spectra and other remote sensing data values to file 
            write.csv(spectra_write, 
                      file = file.path(extracted_features_path,
                                       paste0("extracted_features_",
                                              east_north_string, "_",
                                              shapefile_description,
                                              ".csv")),
                      row.names = FALSE) 
        }  
        
        # combine all extracted features into a single .csv
        paths_ls <- list.files(extracted_features_path, full.names = TRUE)
        
        # refine the output csv selection 
        csvs <- paths_ls[grepl(paste0("*000_", shapefile_description, ".csv"), paths_ls)]
        
        # combine all .csv data into a single data frame 
        for (c in 1:length(csvs)){
          
            csv <- read.csv(csvs[c])
            
            if(c==1){
                spectra_all <- csv
            } else {
                spectra_all <- rbind(spectra_all, csv) #AIS check if extracted features CSV has already been made
            }
        }
        
        # write ALL the spectra to a single .csv file 
        write.csv(spectra_all, file=extracted_features_filename,
                  row.names = FALSE)
        
        # delete the individual csv files for each tile 
        file.remove(csvs)
    }
    
    return(extracted_features_filename)
}



train_pft_classifier <- function(site, year, stacked_aop_path, training_shp_path, training_spectra_path, 
                              data_out_path, pcaInsteadOfWavelengths, ntree, randomMinSamples, 
                              independentValidationSet) {
  # Train random forest model and assess PFT classification accuracy.
  #   Model outputs will be written to a folder within the training directory 
  #   starting with "rf_" followed by a description of each shapefile 
  #   containing points or polygons per tree. 

  #   ntree 
  #       RF tuning parameter, number of trees to grow. default value 500

  #   randomMinSamples 
  #       To reduce bias, this boolean variable indicates whether the same number 
  #       of samples are randomly selected per PFT class. Otherwise,
  #       all samples per class are used for training the classifier. 

  #   independentValidationSet=T
  #       if TRUE, keep separate set for validation
  #       ais this doesnt work if for F - troubleshoot this later
   
  
  # Load data
  # Hyperspectral wavelengths
  wavelengths <- read.csv(file.path({stacked_aop_path},"wavelengths.txt")) %>%
    pull(wavelengths)
  # Stacked AOP layer names
  stacked_aop_layer_names <- read.csv(file.path({stacked_aop_path},"stacked_aop_layer_names.txt")) %>%
    pull(stacked_aop_layer_names)
  # Labelled, half-diam crown polygons to be used for training/validation
  shapefile_description <- tools::file_path_sans_ext(basename(training_shp_path))
  # Csv file containing extracted features
  extracted_features_filename <- file.path(training_spectra_path)
  # Directory to save the classification outputs 
  rf_output_dir <- file.path({data_out_path}, {site}, {year}, "training", paste0("rf_",shapefile_description) )
  if (!dir.exists(file.path(rf_output_dir))) {
    dir.create(file.path(rf_output_dir))
  }
  
  # number of PCAs to keep
  nPCs <- 3  # ais automate this using the elbow method
  
  # Filter out unwanted wavelengths
  wavelength_lut <- filter_out_wavelengths(wavelengths=wavelengths, layer_names=stacked_aop_layer_names)
  
  # features to use in the RF models
  featureNames <- c(wavelength_lut$xwavelength,
                    stacked_aop_layer_names[!grepl("^X", stacked_aop_layer_names)],
                    "pft", "dfID")
  
  # Prep extracted features csv for RF model
  features_df <- prep_features_for_RF(extracted_features_filename, featureNames) 
  
    # alternative code to above for integrating all AOP/inventory data into one training csv
#     # features to use in the RF models
#     featureNames <- c(wavelength_lut$xwavelength,
#                       stacked_aop_layer_names[427:443],
#                       #stacked_aop_layer_names[!grepl("^X", stacked_aop_layer_names)],
#                       "pft", "dfID")
    
#     # Prep extracted features csv for RF model
#     features_df_carissa<- prep_features_for_RF(extracted_features_filename, featureNames) 
#     features_df_soap_2020<- prep_features_for_RF(extracted_features_filename, featureNames) 
#     features_df_soap_2019<- prep_features_for_RF(extracted_features_filename, featureNames) 
#     features_df_sjer_2019 <- prep_features_for_RF(extracted_features_filename, featureNames) 
    
#     colnames(features_df_sjer_2019) <- c("shpname",paste0("X",1:372),colnames(features_df_soap_2020)[374:392])
#     colnames(features_df_carissa)   <- c("shpname",paste0("X",1:372),colnames(features_df_soap_2020)[374:392])
#     colnames(features_df_soap_2020) <- c("shpname",paste0("X",1:372),colnames(features_df_soap_2020)[374:392])
#     colnames(features_df_soap_2019) <- c("shpname",paste0("X",1:372),colnames(features_df_soap_2020)[374:392])
    
#     features_df <- rbind(features_df_sjer_2019, features_df_carissa,
#                          features_df_soap_2020, features_df_soap_2019)
    
#     features_df_balanced <- features_df %>%
#         group_by(pft) %>%
#         dplyr::sample_n(size=min(table(features_df$pft))) %>% ungroup()
#     features_df<-features_df_balanced
    
  # ais see Scholl et al script 08 where they compare sampling bias of using half-diam
  # crowns instead of max diam crowns (search neonvegIDsForBothShapefiles variable)

  # perform PCA
  # remove the individual band reflectances.
  # this only needs to be done once, so check if the validationSet
  # already has a column named "PC1". 
  # testing whether PCA yields better accuracy than individual wavelength reflectance data
  if(pcaInsteadOfWavelengths == T){
        # remove the individual spectral reflectance bands from the training data
        features_noWavelengths <- features_df %>% dplyr::select(-starts_with("X"))
        
        # PCA: calculate Principal Components 
        hs <- data.matrix(features_df %>% dplyr::select(starts_with("X")))
        hs_pca <- stats::prcomp(hs, center = TRUE, scale. = TRUE)
        summ <- summary(hs_pca)
        summ$importance[2,]
        # summary(hs_pca) ais print pc summary info somewhere?
        # add first n PCs to features data frame
        features <- cbind(features_noWavelengths, hs_pca$x[,1:nPCs]) 
        
        # visualize where each sample falls on a plot with PC2 vs PC1 
        ggbiplot::ggbiplot(hs_pca,
                           choices = 1:2, # which PCs to plot
                           obs.scale = 1, var.scale = 1, # scale observations & variables
                           var.axes=FALSE, # remove arrows
                           groups = features_df$pft, # color the points by PFT
                           ellipse = TRUE, # draw ellipse abase::round each group
                           circle = TRUE ) + # draw circle abase::round center of data set
            ggtitle("PCA biplot, PC1 and PC2") +
            scale_color_brewer(palette="Spectral") +
            theme_bw()
        ggplot2::ggsave(file.path(data_int_path,"pcaPlot1vs2.pdf"))
        
        ggbiplot::ggbiplot(hs_pca,
                           choices = 2:3, # which PCs to plot
                           obs.scale = 1, var.scale = 1, # scale observations & variables
                           var.axes=FALSE, # remove arrows
                           groups = features_df$pft, # color the points by PFT
                           ellipse = TRUE, # draw ellipse abase::round each group
                           circle = TRUE ) + # draw circle abase::round center of data set
            ggtitle("PCA biplot, PC3 and PC3") +
            scale_color_brewer(palette="Spectral") +
            theme_bw()
        ggplot2::ggsave(file.path(data_int_path,"pcaPlot2vs3.pdf"))
        
        return(features)
        # features <- apply_PCA(df=features_df, wavelengths=wavelength_lut, output_dir=rf_output_dir)
    } else {
    features <- features_df
  }

  # Split labelled data into training and validation sets
  if(independentValidationSet){
        
    training_val_set <- extract_ind_validation_set(features)
    features_train <- training_val_set[[1]]
    features_val <- training_val_set[[2]]
    
    # Visualize the training and validation data
    train_val_breakdown <- features_val %>% dplyr::count(pft) %>% dplyr::rename(valSet=n) %>%
      left_join(features_train %>% dplyr::count(pft) %>% dplyr::rename(trainSet=n)) %>%
      left_join(features %>% dplyr::count(pft) %>% dplyr::rename(total=n)) 
    
    # Histograms
    train_val_byCount <- rbind(features_train,features_val) %>%
      mutate(train_val = ifelse(dfID %in% features_train$dfID,"train","val")) %>%
      ggplot(aes(x=pft, fill=train_val)) +
      ggtitle("Plot by count") +
      geom_bar() # by count
    train_val_byPerc <- rbind(features_train,features_val) %>%
      mutate(train_val = ifelse(dfID %in% features_train$dfID,"train","val")) %>%
      ggplot(aes(x=pft, fill=train_val)) +
      ggtitle("Plot by %") +
      geom_bar(position="fill") # by percent
    together <- cowplot::plot_grid(train_val_byCount,train_val_byPerc)
    ggsave(together, 
           file=file.path(rf_output_dir,"train_val_split.pdf"),
           width=9, height=6, units="in")
    
    # remove the pixelNumber, easting, and northing columns since they
    # are not input features to the train the classifier 
    features_train <- features_train %>% 
      dplyr::select(-c(pixelNumber, eastingIDs, northingIDs, dfID))
    
  } else { #if we are not using a validation set and training the model on the whole dataset
    features_train <- features
  }
  
  # if(randomMinSamples){
  #   # reduce number of samples per PFT to avoid classifier bias
  #   
  #   # count the minimum number of samples for a single class
  #   minSamples <- min(featureSummary$total)  
  #   print(paste0("Randomly selecting ",
  #                as.character(minSamples),
  #                " samples per PFT class to avoid classifier bias"))
  #   
  #   # isolate the samples per PFT
  #   taxon1 <- features[features$pft==pft_list[1],]
  #   taxon2 <- features[features$pft==pft_list[2],]
  #   taxon3 <- features[features$pft==pft_list[3],]
  #   #taxon4 <- features[features$pft==pft_list[4],]
  #   
  #   # keep random minSamples of each PFT; merge
  #   taxon1 <- taxon1[sample(nrow(taxon1), minSamples), ]
  #   taxon2 <- taxon2[sample(nrow(taxon2), minSamples), ]
  #   taxon3 <- taxon3[sample(nrow(taxon3), minSamples), ]
  #   #taxon4 <- taxon4[sample(nrow(taxon4), minSamples), ]
  #   
  #   features <- rbind(taxon1, taxon2, taxon3)#, taxon4)
  #   
  # } else{
  #   #print("Using all samples per class")
  # }
  
  # TRAIN RF CLASSIFIER using training set  ---------------------------------
  set.seed(104)
  # drop any rows with NA ais investigate this
  features_train_noNA <- na.omit(features_train) 
  rf_model <- randomForest::randomForest(as.factor(features_train_noNA$pft) ~ .,
                                         data=features_train_noNA, 
                                         importance=TRUE, 
                                         ntree=ntree) 
  
  # save RF model to file 
  rf_model_path = file.path(rf_output_dir, paste0("rf_model_",shapefile_description,".RData"))
  save(rf_model, file = rf_model_path)
  
  # y = predicted data; (horizontal axis)
  # x = observed data (true class labels) (vertical axis)
  
  accuracies <- rfUtilities::accuracy(x = rf_model$y,
                                    y = rf_model$predicted)
  
  # record each accuracy metric in the table for a final comparison.
  # base::round each value to the nearest decimal place 
  OA_OOB <- base::round(accuracies$PCC, 1) # Overall Accuracy
  K <- base::round(accuracies$kappa, 3) #Cohen's Kappa 
  
  # INDEPENDENT VALIDATION  -------------------------------------------------
  
  # predict PFT ID for validation set 
  if(independentValidationSet){
    features_val_noNA <- na.omit(features_val) 
    # ^ re-writing here after assigning in script 1. Doing to get PC columns for validaiton
    predValidation <- predict(rf_model, features_val_noNA, type = "class")
    confusionTable <- table(predValidation, features_val_noNA$pft)
    val_OA <- sum(predValidation == features_val_noNA$pft) / 
      length(features_val_noNA$pft)
  }
  
  # Other performance metrics
  # adapted from https://github.com/annaspiers/NEON-NIWO-misclass/blob/master/validation.R
  #validation data only
  # metrics_per_class <- get_model_performance_metrics(confusionTable) 
  # val_perf_metr <- rbind(metrics_per_class,
  #       metrics_per_class %>%
  #         summarize(idx = "macro_values",
  #                   precision=mean(precision, na.rm=T),
  #                   recall=mean(recall, na.rm=T), 
  #                   f1=mean(f1, na.rm=T))) %>%
  #   mutate(across(c(precision, recall, f1), \(x) base::round(x, 2) ))
  
  model_stats <- caret::confusionMatrix(data = rf_model$predicted, 
                                        reference = rf_model$y, 
                                        mode = "prec_recall")
  model_stats_byclass <- data.frame(t(model_stats$byClass)) %>% 
    tibble::rownames_to_column() %>% 
    dplyr::rename(metric=rowname, cedar=Class..cedar, 
                      other=Class..other, 
                      montane_shrub=Class..montane_shrub, 
                      oak=Class..oak, pine=Class..pine) %>%
    mutate(across(c(cedar, oak, pine, montane_shrub, other), \(x) base::round(x, 2) ))
    
  # write all relevant information to the textfile: -------------------------
  
  # open a text file to record the output results
    rf_output_file <- file(description = file.path(rf_output_dir,"rf_model_summaries.txt"), 
              open="w")

  # shapefile name
  write(shapefile_description, rf_output_file, append=TRUE)
  
  # RF model summary, OOB error rate 
  write("\nRaw RF model output with confusion matrix for training set: ", 
        rf_output_file, append=TRUE) 
  capture.output(rf_model, file = rf_output_file, append=TRUE)
  #alt to confusion matrix for training set is accuracies$confusion
  
  write("\n\nOverall Accuracy (OOB):", rf_output_file, append=TRUE) 
  write(base::round(accuracies$PCC/100,3), rf_output_file, append=TRUE)
  
  write("\nCohen's Kappa:", rf_output_file, append=TRUE) 
  capture.output(base::round(accuracies$kappa, 3), file = rf_output_file, append=TRUE)
  
  write("\nUser's Accuracy:", rf_output_file, append=TRUE) 
  capture.output(accuracies$users.accuracy, file = rf_output_file, append=TRUE)
  
  write("\nProducer's Accuracy:", rf_output_file, append=TRUE) 
  capture.output(accuracies$producers.accuracy, file = rf_output_file, append=TRUE)
  
  val_OA <- sum(predValidation == features_val_noNA$pft) / 
    length(features_val_noNA$pft)
  write("\n\nOverall Accuracy (IndVal):", rf_output_file, append=TRUE)
  write(base::round(val_OA,3), rf_output_file, append=TRUE)
  
  # write the accuracy summary data frame to file 
  write("\nConfusion matrix for validation set:", rf_output_file, append=TRUE) 
  capture.output(confusionTable, file = rf_output_file, append=TRUE)
  
  # recall, precision, F1
  
    write("\nOther performance metrics for validation data: ", rf_output_file, append=TRUE) #newline
    write(colnames(model_stats_byclass), rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[1,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[2,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[3,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[4,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[5,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[6,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[7,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[8,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[9,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[10,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    write(as.matrix(model_stats_byclass)[11,], rf_output_file, ncolumns=6, append=TRUE, sep="\t")
    
  # training/val data breakdown
  featureSummary <- data.frame(train_val_breakdown)
  colnames(featureSummary) <- c("pft","validation","training", "total")
  write("\n\nLabelled data split:", rf_output_file, append=TRUE) #newline
  capture.output(featureSummary, file = rf_output_file, append=TRUE)
  
  # features used to describe each sample (pixel)
  write("\ndescriptive features used to train this model: ", rf_output_file, append=TRUE) #newline
  write(colnames(features_train), rf_output_file, ncolumns=5, append=TRUE)
    
  # close the text file
  close(rf_output_file)

  # Save variable importance plots
  plot_variable_importance(rf_model, rf_output_dir)
  
  return(rf_model_path)
}
