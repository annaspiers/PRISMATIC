# library(stringr)
# library(raster)
# library(dplyr)
# library(sf)

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

extract_ind_validation_set <- function(df, percentTrain=0.9) {
  # randomly select <percentTrain> of data from the clipped half diameter 
  # polygons for training, use the remaining samples for validation.
  # keep track of which pixelNumber, easting, northing
  # that these pixels correspond to. 
  
  # percentTrain: randomly select this amount of data for training, 
  # use the rest for validation
  
  #ais create alternatives to random split
  
  # remove any spectra with a height of 0 and remove any factors
  df_val <- df %>% 
    dplyr::filter(chm>0) %>% #ais remove this step once I add in bare ground
    base::droplevels() 
  
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
