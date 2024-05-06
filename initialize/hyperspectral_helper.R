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
library(parallel) #mclapply https://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")
library(rhdf5)


st_erase <- function(x, y) sf::st_difference(x, st_union(st_combine(y))) 
# ^ from https://r-spatial.github.io/sf/reference/geos_binary_ops.html
# keeps only columns from original sf object


list_tiles_w_veg <- function(veg_df, out_dir){
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
  tiles_file <- file(file.path(out_dir,"tiles_w_veg.txt"))
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
  # randomly select <percentTrain> of individuals (select all pixels from individaul) 
  # from outlined polygons for training, use the remaining samples for validation.
  # keep track of which pixelNumber, easting, northing
  # that these pixels correspond to. 
  
  # percentTrain: randomly select this amount of data for training, 
  # use the rest for validation
  
  #ais create alternatives to random split
  
  # remove any spectra with a height of 0 and remove any factors
  df_val <- df 
  
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
    if (!all(base::round(wavelengths) == 
             as.numeric(layer_names[1:length(wavelengths)]))) {
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
    #print("test stop")
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
    r1 <- terra::t(terra::rast(refl_scaled[1,,])) #, crs = crs_info$Proj4
    terra::crs(r1) <- "epsg:32611"
    terra::ext(r1) <- tile_extent
    
    # start the raster stack with first band 
    s <- r1
    
    # loop through bands and create a giant rasterstack with 426 (n_bands) bands
    for(b in 2:n_bands){
        
        # create raster with current band
        
        r <- terra::t(terra::rast(refl_scaled[1,,])) #, crs = crs_info$Proj4
        terra::crs(r) <- "epsg:32611"
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


generate_shadow_mask <- function() {
  # Automatic and Accurate Shadow Detection Using Near-Infrared Information
  # source: http://ivmp.unsw.edu.au/~dominicr/Publications/ruefenacht_pami_2014.pdf
  
  # ais I may need to use raw radiance rather than reflectance?
  ais_normalize <- function(x): #function to scale numbers to between 0 and 1
      x_norm = (x-np.min(x))/(np.max(x)-np.min(x))
      return(x_norm)          
  def nonlinear_map(x, alpha=14, beta=0.5, gamma=2.2):
      # Dark maps
      f = 1/(1 + np.exp(-alpha*(1 - pow(x, 1/gamma) - beta)))
      return(f)
  def gen_ratio_image(x, tau = 100):
    # VIS:NIR ratios
    f = np.zeros((t_k.shape[0],t_k.shape[1]))
    for r in range(t_k.shape[0]):
        for c in range(t_k.shape[1]):
            f[r,c] = (1/tau)*min(x[r,c,:].max(), tau)
    return(f)


  # Normalize 
  p_rnorm = ais_normalize(refl_clip_rxr[11,:,:])      #wavelength 438, band index 11
  p_gnorm = ais_normalize(refl_clip_rxr[33,:,:])      #wavelength 548, band index 33
  p_bnorm = ais_normalize(refl_clip_rxr[55,:,:])      #wavelength 659, band index 55
  p_nirnorm = ais_normalize(refl_clip_rxr[103,:,:])   #wavelength 899., band index 103          
  l_vis = (p_rnorm + p_gnorm + p_bnorm)/3 # brightness image for VIS bands

  # Dark maps
  d_vis = nonlinear_map(l_vis.as_numpy())
  d_nir = nonlinear_map(p_nirnorm.as_numpy())
  d = d_vis*d_nir

  # plot fig 4

  # (a) Visible image
  ep.plot_rgb(refl_clip_arr,rgb=(55, 33, 11), title="Clipped RGB Image", figsize=(2,1)) 

  # (b) Visible shadow candidate map DVIS 
  d_vis.plot.imshow(figsize=(1,1)) #ais how to add title?

  # (c) NIR shadow candidate map DNIR 
  d_nir.plot.imshow(figsize=(1,1)) #ais how to add title?

  # (d) Shadow candidate map D 
  d.plot.imshow(figsize=(1,1)) #ais how to add title?
  plt.show()
  # Even though the presence of a black object in the scene confounds DVIS , D is quite accurate 
  # thanks to DNIR . Images are tone mapped for better visibility.

  # VIS:NIR ratios
  t_k = np.array([p_rnorm/p_nirnorm, p_gnorm/p_nirnorm, p_bnorm/p_nirnorm]).transpose((1,2,0))
  t_arr = gen_ratio_image(t_k)
  t_xr = xr.DataArray(t_arr, coords=d.coords)

  # plot fig 5
  # (a) Visible image. 
  ep.plot_rgb(refl_clip_arr,rgb=(55, 33, 11), title="Clipped RGB Image", figsize=(4,4))
  plt.show()

  # (b) Near-infrared image. 
  refl_clip_rxr[103,:,:].plot.imshow(figsize=(4,3)) #ais how to add title?
  plt.show()

  # (c) Ratio image T that already outlines the shadows (images tone mapped for better visibility).
  t_xr.plot.imshow(figsize=(4,3)) #ais how to add title?
  plt.show()

  # Binary shadow mask
  u = (1 - d) * (1 - t_xr)
  # plot histogram of u
  nu = 1.6
  nbins = round(nu * (math.log2(u.shape[0]*u.shape[1]) + 1))
  plt.hist(u.values.flatten(), bins=nbins)
  plt.show()
  # The threshold θ is set at the location of the first valley in the histogram. 
  # the first valley is the smallest valued bin of the histogram where the two neighboring bins to the
  # left and the two to the right have larger, increasing values. If there
  # is no valley according to our definition, we gradually increase the
  # number of bins until such a valley is found.
  theta = 0.9 #manually specify til I like it
  u_bin = np.where(u > theta, 1, -9999)
  u_bin_xr = xr.DataArray(u_bin, coords=d.coords)

  # plot rgb image
  ep.plot_rgb(refl_clip_arr,rgb=(55, 33, 11), title="Clipped RGB Image", figsize=(4,4)) #refl_arr_mask
  plt.show()

  # plot binary shadow mask
  u_bin_xr.plot.imshow(figsize=(4,3)) #ais how to add title?
  plt.show()
  
}

generate_ndvi_mask <- function(tile, threshold=0.3) {
  vis = refl_clip_rxr[57]
  nir = refl_clip_rxr[89]
  ndvi = np.divide((nir-vis),(nir+vis))
}

prep_aop_imagery <- function(site, year, hs_type, hs_path, tif_path, data_int_path,
                            use_tiles_w_veg=FALSE) {
    # modified from Scholl et al from https://github.com/earthlab/neon-veg
    
    # output directory for stacked AOP data
    stacked_aop_data_dir <- file.path({data_int_path}, {site}, {year}, "stacked_aop")
    if (!dir.exists(stacked_aop_data_dir)) {
      dir.create(stacked_aop_data_dir) 
    } 

    # ais change prep_aop_imagery to stack only 1km x 1km tiles
    # ais add plotting aop imagery from script 6
    
    # list the files in each data directory; filter results based on file type
    if (hs_type=="tile") {
      hs_ls <- list.files(path = file.path({hs_path}),
                        pattern = "reflectance.h5", recursive=T, full.names=T)
    } else {
      hs_ls <- list.files(path = file.path({hs_path}),
                        pattern = "????corrected_tile.tif", recursive=T, full.names=T)
    }    
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
    
    # If user wants to stack only tiles overlapping inventory plots,
    # then filter list to only those
    if (use_tiles_w_veg == TRUE) {
        tiles_w_veg <- read.table(file.path(data_int_path,site,year,
                                            "training/tiles_w_veg.txt")) %>%
                pull(V1)
        hs_ls <- hs_ls[stringr::str_detect(hs_ls, tiles_w_veg %>% paste(collapse = "|"))]
    }
    
    # Stack imagery - create data cubes with AOP-derived features for each tile
    for (hs_path in hs_ls) {        
        # each current hyperspectral tile must be read and stacked into a
        # georeferenced rasterstack (so it can be clipped with fpoint/polygon
        # shapefiles). The process of creating a rasterstack takes a while for
        # each tile, so after creating each rasterstack once, each object gets
        # written to a file.

        # Build up the rasterstack filename by parsing out the easting/northing
        # coordinates from the current h5 filename.

        # parse the UTM easting and northing values from the current h5 filename
        easting <- stringr::str_split(tail(stringr::str_split(hs_path, "/")[[1]],n=1),"_")[[1]][5]
        northing <- stringr::str_split(tail(stringr::str_split(hs_path, "/")[[1]],n=1),"_")[[1]][6]
        east_north_string <- paste0(easting,"_",northing)

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
                if (hs_type=="tile") {
                  hs <- stack_hyperspectral(hs_path, stacked_aop_data_dir)
                } else {
                  hs <- terra::rast(grep(east_north_string, hs_ls, value=TRUE))
                }    
                
                if (terra::ext(hs) != terra::ext(chm)) {
                    message("different extents of hs and chm")
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

                ### Apply masks
                # Create chm mask to mask out non-veg
                # Create ndvi mask to mask out non-veg
                # Create shadow mask
                shadow_mask <- generate_shadow_mask(tile)
                
                # Apply masks
                refl_arr_mask = copy.copy(refl_clip_arr_t) # copy reflectance array
                for band in range(426) : # apply mask to copy of reflectance
                    # chm > 3m
                    refl_arr_mask[:,:,band] = np.where(chm_clip_rxr>3, refl_arr_mask[:,:,band], -9999) 
                    # haze_cloud_water masked out
                    refl_arr_mask[:,:,band] = np.where(haze_clip_rxr==5, refl_arr_mask[:,:,band], -9999) 
                    # weather quality == green
                    refl_arr_mask[:,:,band] = np.where(weather_clip_rxr[1,:,:]>0, refl_arr_mask[:,:,band], -9999) #second index of wq
                    # ndvi > 0.3
                    refl_arr_mask[:,:,band] = np.where(ndvi>0.3, refl_arr_mask[:,:,band], -9999) 
                    # Shadow binary mask == 1
                    refl_arr_mask[:,:,band] = np.where(u_bin_xr==1, refl_arr_mask[:,:,band], -9999) 
                  
                refl_arr_mask_nan = np.where(refl_arr_mask==-9999,np.nan,refl_arr_mask)

                fig, ax = plt.subplots(figsize=(5,5))
                ep.plot_rgb(np.transpose(refl_arr_mask_nan, (2,0,1)), rgb=(58, 34, 19), ax=ax, title="Masked RGB Image") 
                plt.show() 
                # visualize tile after each mask application

                # remove water vapor bands
                valid_band_range = [i for j in (range(0,191), range(212, 281), range(315,415)) for i in j] 
                refl_arr_mask_noH20vap = refl_arr_mask_nan[:, :, valid_band_range] 
                refl_arr_mask_noH20vap.shape
                # Filter out nan's, https://www.w3schools.com/python/numpy/numpy_array_filter.asp
                refl_l_num = []
                for band in range(refl_arr_mask_noH20vap.shape[2]):
                    temp_band = refl_arr_mask_noH20vap[:,:,band]
                    temp_band_noNaN = temp_band[~np.isnan(temp_band)]
                    refl_l_num.append(temp_band_noNaN)
                refl_arr_num = np.array(np.transpose(refl_l_num))
                refl_arr_num.shape # ? x 360

              # ais stack then PCA

                message(paste0("Stacking ",site," ",year," ",east_north_string,", AOP tile ",match(hs_path,hs_ls)," out of ",length(hs_ls)))


                # now, all the remote sensing data have been read in for the current
                # tile. add each one to the hyperspectral data stack along with the
                # layer to keep track of pixel number within the tile.
                stacked_aop_data <- c(hs, chm, slope, aspect_cat, savi,
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



create_tree_crown_polygons <- function(site, year, data_raw_inv_path, data_int_path, biomass_path, px_thresh) {
  # modified from Scholl et al from https://github.com/earthlab/neon-veg
  
  # Create geospatial features (points, polygons with half the maximum crown diameter) 
  # for every tree in the NEON woody vegetation data set that intersect with independent 
  # pixels in the AOP data 
  # Analog to 02-create_tree_features.R and 03-process_tree_features.R
  # from https://github.com/earthlab/neon-veg
  
  # px_thresh is the area threshold (# of of 1m pixels) to filter out small tree crowns
  
  # output directory for training data
  training_data_dir <- file.path({data_int_path}, {site}, {year}, "training")
  if (!dir.exists(training_data_dir)) {
    dir.create(training_data_dir)
  }
  
  live_trees_path <- file.path({biomass_path},"pp_veg_structure_IND_IBA_IAGB_live.csv")
    
  pft_reference <- read.csv(pft_reference_path))
  
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
  tiles <- list_tiles_w_veg(veg_df = veg_training,
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



extract_spectra_from_polygon <- function(site, year, shp_path, data_int_path, data_final_path,
                                        stacked_aop_path, use_case, ic_type=NA, 
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
        extracted_features_path <- file.path(ic_type_path) #ais need to change the folder where this is saved
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
        if (use_case=="train" | ic_type=="rs_inv_plots" ) {
            tiles_w_veg <- read.table(file.path(training_data_dir, "tiles_w_veg.txt")) %>%
                pull(V1)
              
            stacked_aop_list <- stacked_aop_list[stringr::str_detect(stacked_aop_list, tiles_w_veg %>% paste(collapse = "|"))]
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



train_pft_classifier <- function(sites, 
                                  stacked_aop_path, training_shp_path, training_spectra_path, 
                              data_int_path, pcaInsteadOfWavelengths, ntree, randomMinSamples, 
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

    features_df <- data.frame()
    for (site in sites) {
      site_year_paths <- list.dirs(file.path(data_int_path,site), recursive=F)
      for (site_year_path in site_year_paths) {
            
        # Load data
        wavelengths <- read.csv(file.path(site_year_path,"stacked_aop/wavelengths.txt")) %>%
            pull(wavelengths)
        # Stacked AOP layer names
        stacked_aop_layer_names <- read.csv(file.path(site_year_path,"stacked_aop/stacked_aop_layer_names.txt")) %>%
            pull(stacked_aop_layer_names)
        # Labelled, half-diam crown polygons to be used for training/validation
        shapefile_description <- tools::file_path_sans_ext(basename(file.path(site_year_path,"training/tree_crowns_training.shp")))
        # Csv file containing extracted features
        extracted_features_filename <- file.path(site_year_path,"training/tree_crowns_training-extracted_features_inv.csv")
        # Directory to save the classification outputs 
        rf_output_dir <- file.path(data_int_path, "rf_dir")
        if (!dir.exists(file.path(rf_output_dir))) {
            dir.create(file.path(rf_output_dir))
        }
        
        # Filter out unwanted wavelengths
        wavelength_lut <- filter_out_wavelengths(wavelengths=wavelengths, 
                                                layer_names=stacked_aop_layer_names)
        
        # features to use in the RF models
        featureNames <- c(wavelength_lut$xwavelength,
                          stacked_aop_layer_names[(length(wavelengths)+1):length(stacked_aop_layer_names)],
                          "pft", "dfID")
        
        # Prep extracted features csv for RF model
        features_temp <- prep_features_for_RF(extracted_features_filename, featureNames) 
        colnames(features_temp) <- c("shapeID",paste0("X",1:nrow(wavelength_lut)),
                                    featureNames[-(1:nrow(wavelength_lut))])
        #^ just make sure that column names are the same so we can rbind them across site-year
        # ais may want to add columns for site and year 
        
        features_df <- rbind(features_df, features_temp)
        }
    }
  
    # number of PCAs to keep
    nPCs <- 5  # ais automate this using the elbow method
    
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



generate_pft_reference <- function(sites, data_raw_inv_path, trait_table_path) {
  # Assign a PFT to each species

  # load data
  all_veg_structure <- data.frame()
  for (site in sites) {
    site_veg_structure <- readr::read_csv(file.path(data_raw_inv_path,site,"veg_structure.csv")) %>%
      mutate(growthForm_neon = ifelse(grepl("tree", growthForm), "tree",
                                    ifelse(grepl("sapling", growthForm), "tree",
                                            ifelse(grepl("shrub", growthForm), "shrub", growthForm)))) %>%
      dplyr::select(scientificName, taxonID, growthForm_neon, site=siteID)

    all_veg_structure <- rbind(all_veg_structure, site_veg_structure)
  }
  # ais how to group in Carissa's data here?

  trait_table <- readr::read_csv(trait_table_path) %>%
      dplyr::select(scientific, leafPhen=leaf.phen, growthForm_traittable=growth.form)

  # carissa <- readr::read_csv("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/carissa/neon_2023_field_work_final_ais_grouped_species.csv") %>%
  #     mutate(scientificName=ifelse(Species=="OTHER",`If other species: describe here`,Species)) %>%
  #     filter(!is.na(scientificName)) %>%
  #     #rename(taxonID=) %>%
  #     rename(siteID = `NEON Site`) %>%
  #     #rename(growthForm=) %>%
  #     dplyr::select(scientificName,  siteID) %>% #taxonID,growthForm
  #     distinct() %>%
  #     # Use NEON taxonID table taht I downloaded to raw/inventory
  #     # https://data.neonscience.org/taxonomic-lists?taxonTypeCode=PLANT
  #     mutate(taxonID = case_when(
  #         scientificName == "Broad-leaved lupine (Lupinus latifolia var. columbianus)" | 
  #             scientificName == "Lupinus sp." ~ "LUPINSPP",
  #         scientificName == "Abies concolor (Gord. & Glend.) Lindl. ex Hildebr. var. lowiana (Gord. & Glend.) Lemmon (white fir)" |
  #             scientificName == "Abies magnifica A. Murray bis (CA red fir)" ~ "ABIES",
  #         scientificName == "White leaf manzanita (arctostaphulos viscida)" |
  #             scientificName == "Arctostaphylos patula Greene (greenleaf manzanita)" |
  #             scientificName == "White leaf manzanita (arctostaphulos viscida)"  |
  #             scientificName == "Pinemat manzanita (Arctostaphylos nevadensis)"  |
  #             scientificName == "Arctostaphylos nevadensis\n\n\n" |
  #             scientificName == "Arctostaphylos patula Greene (greenleaf manzanita)" ~ "ARCTO3SPP",
  #         scientificName == "Chrysolepis sempervirens (Kellogg) Hjelmqvist (Sierran Chinkapin)" |
  #             scientificName == "Bush chinquapin (Chrysolepis (Castonopsis) sempervirens)" ~ "CHSE11",
  #         scientificName == "Arroyo or Lemmon\xd5s willow" |
  #             scientificName == "Lemmon\xd5s willow (Salix Lemmonii - see notes)" |
  #             scientificName == "Pacific Willow or Arroyo Willow" ~ "SALIX",
  #         scientificName == "Vaccinium sp." ~ "VACCI",
  #         scientificName == "Sambucus nigra L. (elderberry)" ~ "SANI4",
  #         scientificName == "Artemesia californica" ~ "ARCA11",
  #         scientificName == "coffee berry" ~ "PSYCH",
  #         scientificName == "Quercus douglasii Hook. & Arn. (blue oak)" ~ "QUDO",
  #         scientificName == "Pinus sabiniana Douglas ex Douglas (CA foothill pine)" ~ "PISA2",
  #         scientificName == "Ceanothus cuneatus (Hook.) Nutt. (buckbrush)" |
  #             scientificName == "Ceanothus? thyrsiflorus - blue blossom" ~ "CECU",
  #         scientificName == "Quercus wislizeni A. DC. (interior live oak)" ~ "QUWI2",
  #         scientificName == "Ceanothus leucodermis Greene (chaparral whitethorn)" ~ "CELE2",
  #         scientificName == "Aesculus californica (Spach) Nutt. (CA buckeye)" ~ "AECA",
  #         scientificName == "Rhamnus ilicifolia Kellogg (hollyleaf redberry)" ~ "RHIL",
  #         scientificName == "90% thimbleberry, 10% some kind of pea" ~ "RUBUS",
  #         scientificName == "Scouler\xd5s willow \nSalix scouleriana" ~ "SASC",
  #         scientificName == "Black oak \n\nQuecus kolleggii" |
  #             scientificName == "Quercus kelloggii Newberry (CA black oak)" ~ "QUKE",
  #         scientificName == "Incense cedar" |
  #             scientificName == "Calocedrus decurrens (Torr.) Florin (CA incense cedar)" ~ "CADE27",
  #         scientificName == "Quaking aspen (populous tremuloides)" ~ "POTR5",
  #         scientificName == "Broad-leaved lotus (Lotus crassifolius)" ~ "LOCR",
  #         scientificName == "Western juniper (Juniperous occidentalis)" ~ "JUOC",
  #         scientificName == "Black Cottonwood (populus balsamifera)" ~ "POBA2",
  #         scientificName == "Sequoiadendron giganteum" ~ "SEGI2",
  #         scientificName == "Chamaebatia foliolosa Benth. (mountain misery)" ~ "CHFO",
  #         scientificName == "Mountain alder (Alnus incana spp. tenuifolia)" ~ "ALIN2",
  #         scientificName == "Wax currant (Ribes cereum)" ~ "RICE",
  #         scientificName == "Pinus lambertiana Douglas (sugar pine)" ~ "PILA",
  #         scientificName == "Prunus emarginata (Douglas ex Hook.) D. Dietr. (bitter cherry)" ~ "PREM",
  #         scientificName == "Pinus jeffreyi Balf. (Jeffrey pine)" ~ "PIJE",
  #         scientificName == "Pinus contorta Douglas ex Loudon (lodgepole pine)" ~ "PICO",
  #         scientificName == "Fern" ~ "ASPLE",
  #         .default="2PLANT")) %>%
  #     mutate(TEAK=ifelse(siteID=="TEAK","TEAK",NA),
  #            SJER=ifelse(siteID=="SJER","SJER",NA)) %>%
  #     dplyr::select(-siteID)
  
  neon_sites <- all_veg_structure %>%
      arrange(taxonID) %>%
      distinct() %>%
      filter(!is.na(scientificName)) %>%
      mutate(scientific = word(scientificName, 1,2, sep=" ")) %>%
      dplyr::select(scientific, taxonID, SOAP, SJER, TEAK, growthForm_neon) %>%
      mutate(growthForm_neon = ifelse(taxonID=="CONU4","tree",
                                    ifelse(taxonID=="ARNE","shrub",growthForm_neon))) %>%
      # Otherwise get rid of rows with NA in growthForm_neon (these are duplicates)
      filter(!is.na(growthForm_neon)) %>%
      left_join(trait_table) %>%
      # By default use trait table growth form
      mutate(growthForm = ifelse(is.na(growthForm_traittable),
                                growthForm_neon,growthForm_traittable)) %>%
      # Filter out remaining unlikely growth forms - AIS search on Calscape
      mutate(growthForm = case_when(
          taxonID == "ABIES" |
              taxonID == "ABLO" |
              taxonID == "CADE27" ~ "tree",
          taxonID == "CECO" |
              taxonID == "FRCA6" |
              taxonID == "LOIN4" |
              taxonID == "RHIL" |
              taxonID == "RIBES" |
              taxonID == "RIRO" |
              taxonID == "TODI" ~ "shrub",
          taxonID == "AECA" |
              taxonID == "SALIX" ~ "shrub/tree",
          taxonID == "2PLANT" |
              taxonID == "2PLANT-S" |
              taxonID == "LUAL4" ~ NA,
          .default = growthForm  )) %>%
      dplyr::select(scientific,taxonID, SOAP, SJER, TEAK, growthForm) %>%
      distinct()
      
  # Add in Carissa's dataset
  neon_carissa <- neon_sites %>%
      full_join(carissa) %>%
      # Use NEON species name if it exists, otherwise Carissa's
      mutate(scientific = ifelse(is.na(scientific),scientificName,scientific)) %>%
      dplyr::select( scientific,   taxonID, SOAP, SJER, TEAK, growthForm) %>%
      distinct() %>% 
      tidyr::unite("siteID", c(SOAP,  SJER,  TEAK), na.rm = TRUE) %>%
      # Manually remove remaining duplicates
      mutate(siteID=ifelse(taxonID=="AECA", "SOAP_SJER", siteID),
            siteID=ifelse(taxonID=="QUWIF", "SJER_TEAK", siteID),
            
            siteID=ifelse(taxonID=="ARCTO3SPP", "SJER_TEAK", siteID),
            scientific=ifelse(taxonID=="ARCTO3SPP","Arctostaphylos sp.",scientific),
            
            siteID=ifelse(taxonID=="CECU", "SOAP_SJER_TEAK", siteID),
            scientific=ifelse(taxonID=="CECU","Ceanothus sp.",scientific),
            
            siteID=ifelse(taxonID=="CHFO", "SOAP_TEAK", siteID),
            scientific=ifelse(taxonID=="CHFO","Chamaebatia foliolosa",scientific),
            
            siteID=ifelse(taxonID=="LUPINSPP", "SJER_TEAK", siteID),
            scientific=ifelse(taxonID=="LUPINSPP","Lupinus sp.",scientific),
            
            siteID=ifelse(taxonID=="QUKE", "SJER_TEAK", siteID),
            scientific=ifelse(taxonID=="QUKE","Quercus kelloggii",scientific),
            
            siteID=ifelse(taxonID=="QUWI2", "SOAP_TEAK", siteID),
            scientific=ifelse(taxonID=="QUWI2", "Quercus wislizeni", scientific),
            
            scientific=ifelse(taxonID=="SANI4","Sambucus nigra",scientific),
            siteID=ifelse(scientific=="Sambucus nigra", "SOAP_SJER", siteID),
            taxonID=ifelse(scientific=="Sambucus nigra", "SANI4", taxonID),
            
            scientific=ifelse(taxonID=="SALE","Salix sp.",scientific)) %>%
      filter(!(taxonID == "CADE27" & siteID=="SOAP")) %>%
      filter(!(taxonID == "2PLANT" & siteID=="SJER")) %>%
      filter(!(taxonID == "RIRO" & siteID=="SOAP")) %>%
      filter(!(taxonID == "TODI" & siteID=="SJER")) %>%
      filter(!(taxonID == "CECU" & is.na(growthForm))) %>%
      filter(!(taxonID == "CHFO" & is.na(growthForm))) %>%
      filter(!(taxonID == "QUKE" & is.na(growthForm))) %>%
      filter(!(taxonID == "QUWI2" & is.na(growthForm))) %>%
      filter(!(taxonID == "SANI4" & is.na(growthForm))) %>%
      distinct()
      
      
  pft_assignment_df <- neon_carissa %>%
      mutate(pft = case_when(
          taxonID == "ABIES" | taxonID == "ABCO" |
              taxonID == "ABLO" |
              taxonID == "ABMA" ~ "fir", #white or red fir
          taxonID == "CADE27" | taxonID == "SEGI2" ~ "cedar", #incense cedar
          taxonID == "PICO" | taxonID == "PIJE" |
              taxonID == "PILA" | taxonID == "PINACE" |
              taxonID == "PINACE" | taxonID == "PIPO" |
              taxonID == "PISA2" |
              taxonID == "PINUS" ~ "pine", #Ponderosa pine
          taxonID == "QUCH2" | taxonID == "QUKE" |
              taxonID == "QUDO" | taxonID == "QUERC" |
              taxonID == "QUWI2" ~ "oak", #target:canyon live oak (evergreen), also Black oak(deciduous), 
          taxonID == "ARNE" | taxonID == "ARPA" |
              taxonID == "ARPA6" | taxonID == "ARVIM" |
              taxonID == "CEANO" | taxonID == "CECA" |
              taxonID == "CECO" | taxonID == "CECU" |
              taxonID == "CEIN3" | taxonID == "CELE2" |
              taxonID == "CHFO" | taxonID == "CHSE11" |
              taxonID == "CRSE11" | taxonID == "DAWR2" |
              taxonID == "FRCA6" | taxonID == "FRCAC7" |
              taxonID == "LOIN4" | taxonID == "RHAMNA" |
              taxonID == "RHIL" | taxonID == "PREM" |
              taxonID == "RICE" | taxonID == "RIBES" |
              taxonID == "RIRO" | taxonID == "RIVI3" |
              taxonID == "SALIX" | taxonID == "SANI4" |
              taxonID == "CEMOG" | taxonID == "ARCTO3SPP" |
              taxonID == "SEFL3" | taxonID == "TODI" |
              taxonID == "SASC" | taxonID == "AECA" |
              taxonID == "CONU4" | taxonID == "LUPINSPP" ~ "montane_shrub", #Ceanothus sp.
          .default = "other"
      ))
      
  # Print any scientific names that were not addressed in logic tree and instead assigned as 'other'
  message("The following rows were classified as PFT 'other'")
  message(pft_assignment_df %>% filter(pft=="other"))
      
  readr::write_csv(pft_assignment_df, "/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/pft_reference.csv")	
  return(pft_reference_path)
}
























crop_flightlines_to_tiles <- function(site,year_aop,data_raw_aop_path,
                                      hs_flightline_corrected_path) {
  # The to-sensor zenith angle of each pixel should already have been added as a 
  # new layer in the corrected flightline tiff in convert_envi_to_tif() (python)

  # just load in a flightline and see what metadata there are
  # to sensor zenith?

  # Load corrected flightline tiffs
  tif_list_all <- list.files(path = file.path(hs_flightline_corrected_path), pattern=".tif")
  tif_list <- tif_list_all[!grepl("*\\.tif.aux.xml", tif_list_all)]

  # Get HS data into tiles
  aoi_shp <- st_read(file.path(data_raw_aop_path,site,year_aop,"*_merged_tiles.shp")) %>%
      tibble::rownames_to_column() #used as unique ID for each tile
  st_crs(aoi_shp) <- st_crs(32611)
  # plot(aoi_shp$geometry)
  # plot(aoi_shp[85,],axes=T)

  # Check or create directory to save flightline tifs cropped to tiles
  flightline_crop_to_tiles_raw_path <- file.path(location,"flightline_crop_to_tiles_raw")
  if( !dir.exists(flightline_crop_to_tiles_raw_path) ) {
      dir.create(flightline_crop_to_tiles_raw_path)
  }
}
#need to incorporate following lines into code
if (!file.exists(file.path("data/processed/SOAP/2019/shp_cut_to_tifs_1km_tiles",
                           "shp_cut_to_tifs_1km_tiles.shp"))) {
    # Create sf polygon of flightline extents
    for (tif_ind in 1:length(tif_list)) {
        
        # read in flightline
        t <- stack(file.path(location,"flightlines",tif_list[tif_ind]))
        NAvalue(t) <- -9999
        crs(t) <- 32611
        print(paste0("read in tif ",tif_ind," out of ", length(tif_list)))
        # plot(t[[1]],add=T)
        
        shp_crop_temp <- st_crop(aoi_shp, t)
        
        if (tif_ind==1) {
            shp_of_flightlines <- shp_crop_temp
        } else {
            shp_of_flightlines <- rbind(shp_of_flightlines,
                                        shp_crop_temp)
        }
    }
    
    # Save sf polygon locally
    shp_cut_to_tifs <- st_crop(aoi_shp, shp_of_flightlines) 
    plot(shp_cut_to_tifs$geometry)
    shp_cut_to_tifs_1km_tiles <- shp_cut_to_tifs %>%
        mutate(area = st_area(shp_cut_to_tifs)) %>% #units m^2
        # Filter to only tiles 1km x 1km
        filter(as.numeric(area) > 990*990) 
    plot(shp_cut_to_tifs_1km_tiles$geometry,col="red",add=T)
    
    # Save locally
    dir.create("data/processed/SOAP/2019/shp_cut_to_tifs_1km_tiles")
    st_write(shp_cut_to_tifs_1km_tiles, file.path("data/processed/SOAP/2019/shp_cut_to_tifs_1km_tiles",
                                                  "shp_cut_to_tifs_1km_tiles.shp"))
} else {
    shp_cut_to_tifs_1km_tiles <- st_read(file.path("data/processed/SOAP/2019/shp_cut_to_tifs_1km_tiles",
                                                   "shp_cut_to_tifs_1km_tiles.shp"))
}


# Crop tifs to tiles
crop_tif_into_tiles <- function(t_l) {
    
    print(paste0("processing ",t_l))
    
    # read in flightline
    t <- stack(file.path(location,"flightlines",t_l))
    crs(t) <- 32611
    NAvalue(t) <- -9999
    #t[[1]][t[[1]] < -9990] <- NA
    #system.time(t_m <- mask(t,t[[1]])) 
    # ais the above steps takes way too long - how to load in the raster stack with 
        
    # Crop 1km tile shp to extent of tif
    shp_crop_temp <- st_crop(shp_cut_to_tifs_1km_tiles, t[[1]])
    
    # plot(aoi_shp$geometry)
    # plot(shp_cut_to_tifs_1km_tiles$geometry,col="red",add=T)
    # plot(shp_crop_temp$geometry,col="blue",add=T)
    # plot(t[[1]],add=T)

    # for each tile overlapping with tif
    # crop tif to NEON tiles
    for (tile_ind in 1:nrow(shp_crop_temp)) { #ais I need to parallelize this
        
        tile_temp <- shp_crop_temp[tile_ind,] 
        st_crs(tile_temp) <- 32611
        # plot(tile_temp$geometry,col="lightblue",add=T)
        
        # Skip if the clipped tif already exists
        if (file.exists(file.path(flightline_crop_to_tiles_raw_path,
                                  paste0("tile_",tile_temp$rowname,
                                         "_x_flightline_",
                                         stringr::str_extract(
                                             t_l,"\\d{8}_\\d{6}"),".tif")))){
            cat(paste0(t_l, " already clipped to tile ",tile_temp$rowname," (",
                       tile_ind,"/",nrow(shp_crop_temp),")"),"\n", 
                file = file.path(flightline_crop_to_tiles_raw_path,"log.txt"), 
                append = TRUE)
            next
        } else {
        
            # if tile does not intersect with tif, then skip tile
            r_mask <- mask(t[[1]], tile_temp)
            
            # Skip if there is no overlap between raster and tile
            if (r_mask@data@min <= 0 && r_mask@data@max <= 0){
                cat(paste0("no overlap between tile ",
                           tile_temp$rowname," (",tile_ind,"/",
                           nrow(shp_crop_temp),") for ", t_l),"\n", 
                    file = file.path(flightline_crop_to_tiles_raw_path,"log.txt"), 
                    append = TRUE)
                next
            } else {
    
                # crop tif to tile
                t_c <- raster::crop(t, tile_temp)
                for (band in 1:dim(t_c)[3]) {
                    t_c[[band]][t_c[[band]] < -9990] <- NA
                    print(band)
                }
                t_c_b <- brick(t_c)
                
                # save cropped tif locally to a new folder with the tile coordinates so 
                # the next script will be able to call the cropped flihgtlines within the 
                # same tile together
                writeRaster(t_c_b, filename=file.path(flightline_crop_to_tiles_raw_path,
                                                    paste0("tile_",tile_temp$rowname,
                                                           "_x_flightline_",
                                                           stringr::str_extract(
                                                               t_l,"\\d{8}_\\d{6}"))),
                            format="GTiff", overwrite=T, 
                            options=c("INTERLEAVE=BAND","COMPRESS=LZW"),
                            NAflag=-9999)
                
                cat(paste0("cropped + saved tile ",tile_temp$rowname," (",
                           tile_ind,"/",nrow(shp_crop_temp),") for ", t_l),"\n", 
                    file = file.path(flightline_crop_to_tiles_raw_path,"log.txt"), 
                    append = TRUE)
            }
        }
    }
    
    # delete the full flightline locally after cropping is done
    unlink(file.path(location,"flightlines",t_l))
    cat(paste0(t_l, " deleted from disk"),"\n", 
        file = file.path(flightline_crop_to_tiles_raw_path,"log.txt"), 
        append = TRUE)
}
#crop_tif_into_tiles(tif_list)
my_cores <- detectCores()-1
system.time(
    test <- mclapply(tif_list, FUN=crop_tif_into_tiles, #mclapply doesn't like plot fn
         mc.preschedule = T, mc.silent = F, mc.cores = my_cores)
)
# Sanity check
# # do these cropped tiles
# 
# # full flightline tiff
# # cropped tif to tile (arbitrary)
# t_full <- stack(file.path("notebooks/data_neoncyverse/SOAP_2019",
#                      "flightlines",
#                      "NEON_D17_SOAP_DP1_20190613_164442_reflectance_BRDF_topo_corrected.tif"))
# crs(t_full) <- 32611
# NAvalue(t_full) <- -9998
# t_full[[1]][t_full[[1]] < -9990] <- NA
# system.time(t_full_m <- mask(t_full,t_full[[1]])) 
# 
# # cropped tif to tile (arbitrary)
# t_tile <- stack(file.path("notebooks/data_neoncyverse/SOAP_2019",
#                      "flightline_crop_to_tiles_raw",
#                      "tile_35_x_flightline_20190613_164442.tif"))
# crs(t_tile) <- 32611
# NAvalue(t_tile) <- -9999



# later, when merging to final tiff, 
    # Merge flightiness by tile by choosing pixels with lowest angle of collection
    # ais do this in python using the mosaic_tif_files.py script
#save final tiled hs data in "dev/neon-veg-SOAPpfts/data/data_raw/"

