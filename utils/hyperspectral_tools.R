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
                            dplyr::select(-area_m2)
  
  # create a copy of the polygons to update with clipped/deleted entries
  polys_filtered <- polys_ordered
  
  # create an empty vector of polygon pairs that have been compared 
  compared_pairs <- c()
  
  for (i in 1:nrow(polys_ordered)){ 

    individualID <- polys_ordered$individualID[i]
    
    # if this polygon was removed from the polys_filtered
    # data frame in a previous iteration, skip it 
    if(sum(polys_filtered$individualID==individualID) == 0){
      next
    }
    
    # extract current vs. all other polygon from data frame
    current_poly <- polys_ordered[i,] 
    other_polys <- polys_filtered[polys_filtered$individualID!=current_poly$individualID,]
    # other_polys_insubplot <- other_polys[other_polys$plotID.x == current_poly$plotID.x & 
    #                                     other_polys$subplotID == current_poly$subplotID,] #ais added this
    
    
    # check for overlap between current polygon and all polygons
    # overlap <- raster::intersect(current_poly, other_polys) #crashes here - why?
    overlap_sf <- other_polys[sf::st_overlaps(current_poly, other_polys)[[1]], ] 
    n_overlap <- nrow(overlap_sf)
    
    if(n_overlap>0){ 
      for (o in 1:n_overlap){
        
        # if current polygon ID is not in filtered set
        if(sum(polys_filtered$individualID==individualID) == 0){
          break
        }
                
        # get height of current and test polygons that overlap
        current_poly.height <- current_poly$height
        
        test_poly <- overlap_sf[o,]        
        test_poly.height <- test_poly$height
        
        # combine the ID's of the current and test polygons
        # to keep track of which pairs have been compared 
        id.pair <- paste(current_poly$individualID,
                         test_poly$individualID,
                         sep = " ")
        
        # if polygon pair was already compared, skip to the next pair
        if (id.pair %in% compared_pairs) {
          
          next
          
        } else { 
          
          # add to the list of polygon pairs that have been compared 
          compared_pairs <- append(compared_pairs,id.pair)
          
                    # add opposite combination of polygons
          id.pair.opp <- paste(test_poly$individualID, 
                                       current_poly$individualID,
                                       sep = " ")
          compared_pairs <- append(compared_pairs,id.pair.opp)
          
        }
        
        # if test polygon is not in filtered set, skip to the next overlapping polygon
        # (because it has already been removed)
        if(sum(polys_filtered$individualID==test_poly$individualID) == 0){
          next
        }
                
        # if the area of one of the polygons is equivalent to zero, delete it and 
        # skip to the next overlapping polygon. 
        if(as.numeric(st_area(current_poly)) < 0.01){ 
          
          # delete the current polygon from the polygons_filtered data frame. 
          polys_filtered <- polys_filtered[polys_filtered$individualID!=current_poly$individualID,]
          
          next
          
        } else if(as.numeric(st_area(test_poly)) < 0.01){ 
          # delete the test polygon from the polygons_filtered data frame if its
          # area is essentially 0 (using the criterion < 0.01 here). This step 
          # was added for instances where after clipping, a small edge fragment remains.
          polys_filtered <- polys_filtered[polys_filtered$individualID!=test_poly$individualID,]
          
          next
          
        }
        
        # compare the heights of the polygons
        if(current_poly.height > test_poly.height){
          
          # if current polygon is taller, clip the test polygon.
          clipped <- st_erase(test_poly, current_poly)
                   
          # if the clipped region is NULL, skip to the next polygon comparison.
          if(is.null(clipped)){
            #----
            # VS-NOTE: seeing if this fixes the issue with NEON.PLA.D13.NIWO.01626
            # where an engulfed polygon remains after the workflow...
            # delete the test polygon from the polygons_filtered data frame
            # in this case because the test polygon is taller and totally 
            # encompassess the current polygon 
            #polys_filtered <- polys_filtered[polys_filtered$individualID!=current_poly$individualID,]
            #-----
            next

          } else{ 
              
            # get the index within the filtered polygons data frame 
            # where the test polygon belongs. 
            j <- which(polys_filtered$individualID == test_poly$individualID)
          }
          
        } else {
          
          # otherwise, the test polygon is taller: clip the current polygon.
          clipped <-  st_erase(current_poly,test_poly)
                    
          if(is.null(clipped)){
            next

          } else {          
          
          # get the index within the filtered polygons data frame 
          # where the current polygon belongs. 
          j <- which(polys_filtered$individualID == current_poly$individualID)
          
          }
        }
        
        # if there is no clipped area, skip to the next overlap polygon
        if(nrow(clipped) == 0){
          next
          
        } else{
          
          # check the area of the clipped test polygon. If it is greater than
          # or equal to the area threshold, replace it as the polygon 
          # geometry for the entry matching the test individualID in the 
          # polygons_filtered data frame. 
          if(as.numeric(st_area(clipped))  >= thresh){
            
            # replace the original polygon with the clipped polygon
            polys_filtered[j,] <- clipped 
            
            
          } else{
            # otherwise, delete the test polygon from the polygons_filtered data frame. 
            polys_filtered <- polys_filtered[polys_filtered$individualID!=clipped$individualID,]
            
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
