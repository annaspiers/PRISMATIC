# install.packages("bvls")
# install.packages("assertthat")
# install.packages("rjson")
library(lidR)
library(sf)
library(randomForest)
library(ggplot2)
library(bvls) #bvls()
library(tidyr) #pivot_longer()
library(assertthat)
library(dplyr)


prepare_spatial_data_for_ic <- function(site, year_inv, year_aop, ic_type, data_raw_aop_path, 
                                        data_int_path, stacked_aop_path, data_final_path, ic_type_path,
                                        use_case, n_plots, min_distance, plot_length, 
                                        aggregate_from_1m_to_2m_res) {
    # This function generates random plots across the remote sensing 
    # spatial extent and clips laz and rasters 

    ### Specify spatial extent of initial conditions
    if (ic_type=="rs_inv_plots") { 
        #use existing inventory plot extents
        plots_shp_path <- file.path(data_int_path,site,year_inv,"inventory_plots", "plots.shp")
        laz_clipped_path <- file.path(data_int_path,site,year_inv,"clipped_to_plots")
        
    } else if (ic_type=="rs_random_plots") {
        
        plots_shp_path <- file.path(ic_type_path,"random_plots.shp")
        laz_clipped_path <- file.path(ic_type_path,"clipped_to_plots")

        if (!file.exists(plots_shp_path)) {
            #Generate random plots across AOP spatial extent
            generate_random_plots(site, year_aop, year_inv, data_raw_aop_path, ic_type_path,
                                  n_plots, min_distance, plot_length)
        }
            
        # Clip normalized point clouds to randomized plots
        if (!dir.exists(laz_clipped_path)) {
            laz_clipped_path <- clip_norm_laz_to_shp(site, year_inv, data_int_path,
                                                 ic_type_path, plots_shp_path)
        }
    }
    
    extracted_features_csv <- extract_spectra_from_polygon(site=site,
                                                           year=year_inv,
                                                           data_int_path=data_int_path,
                                                           data_final_path=data_final_path,
                                                           stacked_aop_path=stacked_aop_path,
                                                           shp_path=plots_shp_path,
                                                           use_case=use_case,
                                                           ic_type=ic_type,
                                                           ic_type_path=ic_type_path,
                                                           aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res)

    return(c(plots_shp_path, laz_clipped_path, extracted_features_csv))
}



clip_norm_laz_to_shp <- function(site, year_inv, data_int_path, ic_type_path,
                                 plots_shp_path) {
    # Function to clip normalized lidar point cloud to plots being  
    # used to generate initial conditions
    
    laz_clipped_random_path <- file.path(ic_type_path, "clipped_to_plots")
    
    # Create directory if it doesn't exist                                     
    if (dir.exists(laz_clipped_random_path)) {
        
    } else {
        dir.create(laz_clipped_random_path) 
    }
    
    # Only clip laz if directory is empty
    if (length(list.files(laz_clipped_random_path)) == 0 ) {
        
        las_ctg <- readLAScatalog(file.path(data_int_path, site, year_inv,
                                            "normalized_lidar_tiles"))
        opt_output_files(las_ctg) <- file.path(laz_clipped_random_path,
                                               "plot_{PLOTID}_{XCENTER}_{YCENTER}")
        opt_laz_compression(las_ctg) <- T
        
        plots_sf <- read_sf(plots_shp_path)

        # Extract the plots. Returns a list of extracted point clouds in R
        random_plots_ctg <- lidR::clip_roi(las_ctg, plots_sf)
        #ais can I parallelize this? ^
        #ais address differences in tile count between normalized_lidar_tiles, aop_tiles, and this

    } else {
       # jump to bottom and return path
        message("normalized laz already clipped to plots")
    }
    
    return(laz_clipped_random_path)
}



generate_random_plots <- function(site, year_aop, year_inv, data_raw_aop_path, ic_type_path,
                                  n_plots, min_distance, plot_length) {
    # Function to generate random points within a flight boundary
    
    random_plots_path <- file.path(ic_type_path, "random_plots.shp")
    
    if (file.exists(random_plots_path)) {
        message(paste0("file with ",n_plots," random plots already generated"))
    } else {
        aop_tiles_path <- list.files(file.path(data_raw_aop_path, site, year_aop, "shape/"),
                                     pattern="merged_tiles.shp", full.names=T)
        aop_tiles <- sf::st_read(aop_tiles_path)
        
        # Create or load a polygon
        # Limited by hyperspectral extent with "green" cloud conditions
        flight_bdry_sf <- aop_tiles$geometry
        
        #extent_temp <- read_sf("/Users/AISpiers/Downloads/QGIS_LDRD_StudyArea/NEON_raster_CHM/NEON_SOAP_A01_2019_Boundary.shp")
        #^ais what is this? Also, will need to filter this based on what tiles
        #we filter to based on the 'green' cloud cover flag, etc.
        
        points <- st_sfc(crs=st_crs(flight_bdry_sf))
        plots <- st_sfc(crs=st_crs(flight_bdry_sf))
        n <- 0
        
        while (n < n_plots) {
            point <- st_sample(flight_bdry_sf, size=1)
            
            # Check that point sufficiently far from existing points and polygon edge
            if (n == 0 ||
                as.numeric(min(st_distance(point, points))) >= min_distance) {
                if (as.numeric(min(st_geometry(obj = flight_bdry_sf) %>%
                                   st_cast(to = 'LINESTRING') %>%
                                   st_distance(y = point))) >= min_distance) {
                    
                    # Convert to 20m x 20m
                    plot <- point %>%
                        st_buffer(dist = plot_length/2) %>%
                        st_bbox() %>%
                        st_as_sfc() %>%
                        st_as_sf()
                    
                    plots <- rbind(plots, plot) %>%
                        as.data.frame %>%
                        sf::st_as_sf()
                    st_crs(plots) <- st_crs(flight_bdry_sf)
                    
                    # points <- rbind(points, point) %>%
                    #     as.data.frame %>%
                    #     sf::st_as_sf()
                    
                    n <- n+1
                }
            }
        }
        
        # Visualize 
        plot(flight_bdry_sf)
        plot(plots, pch=".", add=T)
        
        plots_final <- plots %>%
            tibble::rowid_to_column() %>%
            rename(PLOTID=rowid) %>%
            mutate(subplotID = "1") 
        
        # save as shp
        st_write(plots_final, random_plots_path, delete_layer = TRUE)
    }
    
    return(random_plots_path)
}



predict_plot_level_PFTs <- function(site, year, data_int_path, data_final_path,
                         stacked_aop_path, ic_type_path, rf_model_path, extracted_features_filename,
                         ic_type, pcaInsteadOfWavelengths) {

    # ais Commented out the following lines since I want to predict the inventory PFTs
    # across all inventory plots using RF classifier, not just use the pixels in 
    # tree crown polygons
    # if (ic_type=="rs_inv_plots") {
    # 
    #     features <- read.csv(extracted_features_filename)
    # 
    #     features_noNA <- na.omit(features) 
    # 
    #     # Summarize as percentages per plot
    #     features_pft_by_pct <- features_noNA %>% 
    #         #ais here, I should assign pixels with chm<1 to a bare ground PFT
    #         # but before I can do that, I need to change the system of equations for 
    #         # assigning PFTs across cohorts
    #         dplyr::group_by(shapeID, PFT) %>%
    #         dplyr::summarize(count = n()) %>%
    #         dplyr::group_by(shapeID)%>%
    #         dplyr::mutate(pct = count / sum(count))
    # 
    #     classified_PFTs_path <- file.path(ic_type_path,"pfts_per_plot_inv.csv")
    # 
    # } else if (ic_type=="rs_random_plots") {
    #     ### Set parameters and load data

    classified_PFTs_path <- file.path(ic_type_path,"pfts_per_plot.csv")

    if (!file.exists(classified_PFTs_path)){
        message("Predicting plot-level PFTs")

    wavelengths <- read.csv(file.path({stacked_aop_path},"wavelengths.txt")) %>%
        pull(wavelengths)
    # Stacked AOP layer names
    stacked_aop_layer_names <- read.csv(file.path({stacked_aop_path},"stacked_aop_layer_names.txt")) %>%
        pull(stacked_aop_layer_names)
       
    # Filter out unwanted wavelengths
    wavelength_lut <- filter_out_wavelengths(wavelengths=wavelengths, layer_names=stacked_aop_layer_names)
    
    # features to use in the RF models
    featureNames <- c("shapeID", wavelength_lut$xwavelength,
                      stacked_aop_layer_names[!grepl("^X", stacked_aop_layer_names)],
                      "dfID")
    
    # extracted_features_raw <- read.csv(extracted_features_filename)
    
    # Prep extracted features csv for RF model
    features_df <- prep_features_for_RF(extracted_features_filename, featureNames) 
    
    features_df_noNA <- na.omit(features_df) #ais why are there rows with NAs?
    
    # perform PCA
    if(pcaInsteadOfWavelengths == T){
        features <- apply_PCA(df=features_df_noNA, wavelengths=wavelength_lut, 
                              output_dir=ic_type_path)
    } else {
        features <- features_df_noNA
    }
    
    
    ### Predict PFT ID for FATES patches
    
    # load the current RF model
    load(rf_model_path)
        
    # Classify PFTs
    features$pred_PFT <- predict(rf_model, features, type = "class")
    # ^ takes < 3min
    
    # confusionTable <- table(predValidation, features_noNA$pft)
    # val_OA <- sum(predValidation == features_noNA$pft) / 
    #     length(features_noNA$pft)
    
    # Summarize as percentages per plot
    features_pft_by_pct <- features %>% 
        #ais here, I should assign pixels with chm<1 to a bare ground PFT
        # but before I can do that, I need to change the system of equations for 
        # assigning PFTs across cohorts
        dplyr::count(shapeID, pred_PFT) %>%
        dplyr::rename(count=n) %>%
        dplyr::group_by(shapeID)%>%
        dplyr::mutate(pct = count / sum(count))

    if (ic_type=="rs_inv_plots" ) {
        features_pft_by_pct$shapeID <- paste0(features_pft_by_pct$shapeID, "_central")
    }
    
    
    write.csv(features_pft_by_pct, classified_PFTs_path,
              row.names = FALSE)
    } 

    return(classified_PFTs_path)
}



divide_patch_into_cohorts <- function(laz_plot_paths, ic_type) {
    
    # Initialize dataframe
    by_patch_df <- data.frame()
    #height_step <- 0.5 #m #for height step in lad_csv
    
    for (c in 1:length(laz_plot_paths)) {
        print(c)
    
        # Load leaf area density files
        lad_json <- rjson::fromJSON(file=paste0(laz_plot_paths[c],"_lad.json"))
        lad_csv <- read.csv(paste0(laz_plot_paths[c],"_lad.csv"))
    
        if ( nrow(lad_csv)==0 | length(lad_json$layer_height)==0 ) {
            # In this case, the patch is flat and has no cohorts
            patch_temp <- data.frame()
    
            #ais but what about the cases where json has data but csv is empty?
               # ^ in these cases I got the warning "lad could not be calculated"
               # in a previous step. what's the difference?
    
        } else {
            # lad units are m2/m3
            #lad_json$cohort_idx <- length(lad_json$layer_height):1
            if (is.null(lad_csv$X)) {
                lad_csv <- lad_csv %>%
                    #cut off rows beyond highest layer_height
                    filter(z <= max(lad_json$layer_height)) #ais add 0.5 here?
            } else {
                lad_csv <- lad_csv %>%
                    dplyr::select(-c("X")) %>%
                    #cut off rows beyond highest layer_height
                    filter(z <= max(lad_json$layer_height)) #ais add 0.5 here?
            }
            
            #ais figure out how to vectorize this rather than for-loop
            for (l in 1:nrow(lad_csv)) {
                # identify the closest layer_height for each row
                # choosing reference height as top height in cohort
                lad_csv$cohort_height[l] <- min(lad_json$layer_height[
                    lad_json$layer_height >= lad_csv$z[l]])
    
                # assign which cohort number
                #lad_csv$cohort_idx[l] <- lad_json$cohort_idx[which(
                #    lad_json$layer_height == lad_csv$cohort_height[l])]
    
                # calculate the difference in height from one row to the next
                # (for calculating LAI)
                lad_csv$diff_z[l] <- ifelse(l==1,lad_csv$z[l],
                                            lad_csv$z[l]-lad_csv$z[l-1])
            }
    
            cohort_idx_df <- lad_csv %>%
                filter(lad>0) %>%
                distinct(cohort_height) %>%
                arrange(desc(cohort_height)) %>%
                dplyr::mutate(cohort_idx = row_number())
            
            if (grepl("central", basename(laz_plot_paths[c]))) {
                patch_temp <- lad_csv %>%
                    filter(lad>0) %>%
                    # assign patch # and cohort index
                    mutate(patch = basename(laz_plot_paths[c])) %>%
                    left_join(cohort_idx_df)
            } else {
                if (ic_type=="rs_inv_plots") {
                    patch_temp <- lad_csv %>%
                        filter(lad>0) %>%
                        # assign patch # and cohort index
                        mutate(patch = sub(".*_(.*?)\\.*", "\\1",
                                           laz_plot_paths[c])) %>%
                        left_join(cohort_idx_df)
                } else if (ic_type=="rs_random_plots") {
                    patch_temp <- lad_csv %>%
                        filter(lad>0) %>%
                        # assign patch # and cohort index
                        mutate(patch = sub(".*/plot_(\\d+)_.*", "\\1", 
                                           laz_plot_paths[c])) %>%
                        left_join(cohort_idx_df)
                }
                
            }
        }
    
        # assign things back into rows of cohort_df
        by_patch_df <- rbind(by_patch_df, patch_temp)
    }
    
    return(by_patch_df)
}



assign_pft_across_cohorts <- function(by_patch_df, allom_params, 
                                      pfts_by_cohort_wide) {
    
    breakdown_to_pft_df <- by_patch_df %>%
        # lai = sum(layer thickness (diff_z) * LAD of layer)
        # lai (layer 1) = stem density (layer 1) * ILA(as a funciton of z)
        # ILA = Bleaf * SLA , SLA is specific to height and PFT - for now just pft
        dplyr::mutate(lai = lad*diff_z) %>%
        dplyr::group_by(patch,cohort_idx,cohort_height) %>%
        dplyr::reframe(lai_cohort = sum(lai)) %>%
    
        # assign Hmax layer
        dplyr::mutate(layer = ifelse(cohort_height <= allom_params %>%
                                  filter(pft=="oak_PFT") %>% pull(Hmax),
                              "a", "b")) %>%
    
        # assign layer LAI
        dplyr::group_by(patch, layer) %>%
        dplyr::mutate(lai_layer = sum(lai_cohort)) %>%
        dplyr::ungroup() %>%
    
        # assign total LAI for patch
        dplyr::group_by(patch) %>%
        dplyr::mutate( lai_T = sum(lai_cohort) ) %>%
        dplyr::ungroup() %>%
    
        # percent of each pft per patch
        dplyr::left_join(pfts_by_cohort_wide) %>%
        # filter out rows where by_patch_df had a patch but pfts_by_cohort_wide doesn't
        drop_na() %>%
        
        # known fractions
        #f_o_b = fraction of oaks in layer b
        #f_o_a = fraction of oaks in layer a
        mutate(f_o = ifelse(layer=="b", 0, (p_o_T*lai_T)/lai_layer))
    
    # Estimate f_c and f_p for each layer
    message("Estimating PFT fractions across cohorts")
    
    # initialize
    params_est <- data.frame()
    for (p in unique(breakdown_to_pft_df$patch)) {
        temp_df <- breakdown_to_pft_df %>%
            dplyr::filter(patch == p)
    
        lai_a <- temp_df %>% dplyr::filter(layer=="a") %>% dplyr::distinct(lai_layer) %>% dplyr::pull(lai_layer)
        lai_b <- ifelse(nrow(temp_df %>% dplyr::filter(layer=="b"))==0, 0,
                        temp_df %>% dplyr::filter(layer=="b") %>%
                            dplyr::distinct(lai_layer) %>% dplyr::pull(lai_layer))
        lai_T <- temp_df %>% dplyr::distinct(lai_T) %>% dplyr::pull(lai_T)
        p_c_T <- temp_df %>% dplyr::distinct(p_c_T) %>% dplyr::pull(p_c_T)
        p_p_T <- temp_df %>% dplyr::distinct(p_p_T) %>% dplyr::pull(p_p_T)
        f_o_a <- temp_df %>% dplyr::filter(layer=="a") %>% dplyr::distinct(f_o) %>% pull(f_o)
    
        # system of equations
        # 1 - f_o_a = f_c_a + f_p_a... the same as 1 = f_c_a + f_p_a + f_o_a
        #         1 = f_c_b + f_p_b
        # p_c_T*lai_T = f_c_a*lai_a + f_c_b*lai_b
        # p_p_T*lai_T = f_p_a*lai_a + f_p_b*lai_b
    
        # also an equation, but not included above
        # p_o_T*lai_T = f_o_a*lai_a + f_o_b*lai_b
        # but since there is no oak in layer b, f_o_b=0
        # p_o_T*lai_T = f_o_a*lai_a
    
        # Write in matrix form
        # column order: f_c_a, f_c_b, f_p_a, f_p_b
        X <- matrix(c(1,0,1,0,
                      0,1,0,1,
                      lai_a,lai_b,0,0,
                      0,0,lai_a,lai_b), 4,4, byrow=TRUE)
        b_l <- rep(0,ncol(X))
        b_u <- rep(1,ncol(X))
        
        if (f_o_a > 1) { #minority of cases
            #e.g., "SOAP_021_central"
            # f_o_a should be between 0 and 1, but sometimes it is greater than 1
            # when f_o_a > 1, just set all of layer a to oak,
            # and split up layer b between pine and cedar
            f_o_a <- 1
            # which then forces f_p_a = f_c_a = 0, and we need to recalibrate.
            # since oak can only be in layer a, but the proportion of lai in
            # layer a is less than p_o_T, then we set p_o_T = lai_a/lai_T,
            # and recalibrate p_c_T and p_p_T accordingly
            rescal <- (1 - lai_a/lai_T)/(p_p_T + p_c_T)
            if (rescal!="Inf") {
                p_c_T <- rescal*p_c_T
                p_p_T <- rescal*p_p_T
            }
            # so now the system of equations above should balance out
    
            y <- matrix(c(1 - f_o_a, 1, p_c_T*lai_T, p_p_T*lai_T), 4, 1)
    
            sol_temp <- bvls(X, y, b_l, b_u, key=1, istate=c(-1,-3,2,4,2))
            dev_df <- data.frame(f_c_a=sol_temp$x[1], f_c_b=sol_temp$x[2],
                                 f_p_a=sol_temp$x[3], f_p_b=sol_temp$x[4],
                                 b_l_1=b_l[1], b_l_2=b_l[2],
                                 b_l_3=b_l[3], b_l_4=b_l[4],
                                 dev=sol_temp$deviance)
            i <- 1
            while (max(b_l) < 1) {
                ind_fixed <- which(dev_df[i,c(2,4)] == min(dev_df[i,c(2,4)]))
                b_l[ind_fixed] <- b_l[ind_fixed] + 0.01
                sol_temp <- bvls(X, y, b_l, b_u, key=1, istate=c(-1,-3,2,4,2))
    
                # Add next row
                dev_df <- rbind(dev_df,
                                c(x_1=sol_temp$x[1], x_2=sol_temp$x[2],
                                  x_3=sol_temp$x[3], x_4=sol_temp$x[4],
                                  b_l_1=b_l[1], b_l_2=b_l[2],
                                  b_l_3=b_l[3], b_l_4=b_l[4],
                                  dev=sol_temp$deviance))
                i <- i+1
            }
        } else { #majority of cases
            #e.g., "SOAP_002_central"
    
            y <- matrix(c(1 - f_o_a, 1, p_c_T*lai_T, p_p_T*lai_T), 4, 1)
    
            sol_temp <- bvls(X, y, b_l, b_u)
            dev_df <- data.frame(f_c_a=sol_temp$x[1], f_c_b=sol_temp$x[2],
                                 f_p_a=sol_temp$x[3], f_p_b=sol_temp$x[4],
                                 b_l_1=b_l[1], b_l_2=b_l[2],
                                 b_l_3=b_l[3], b_l_4=b_l[4],
                                 dev=sol_temp$deviance)
            i <- 1
            while (max(b_l) < 1) {
                ind_fixed <- which(dev_df[i,1:4] == min(dev_df[i,1:4]))
                b_l[ind_fixed] <- b_l[ind_fixed] + 0.01
                sol_temp <- bvls(X, y, b_l, b_u) #solve(X, y)
    
                # Add next row
                dev_df <- rbind(dev_df,
                                c(x_1=sol_temp$x[1], x_2=sol_temp$x[2],
                                  x_3=sol_temp$x[3], x_4=sol_temp$x[4],
                                  b_l_1=b_l[1], b_l_2=b_l[2],
                                  b_l_3=b_l[3], b_l_4=b_l[4],
                                  dev=sol_temp$deviance))
                i <- i+1
            }
    
        }
    
        params_temp <- cbind(p, f_o_a, dev_df %>%
                                 dplyr::filter(dev==min(dev_df$dev)) %>%
                                 dplyr::slice_head(n=1)  )
    
        params_est <- rbind(params_est ,params_temp)
    }
    
    params_est_f_o <- params_est %>%
        dplyr::select(patch=p, f_o_a) %>%
        mutate(f_o_b = 0) %>%
        tidyr::pivot_longer(cols=c(f_o_a,f_o_b), names_prefix = "f_o_", names_to="layer",
                     values_to="f_o")
    params_est_f_c <- params_est %>%
        dplyr::select(patch=p, f_c_a, f_c_b) %>%
        tidyr::pivot_longer(cols=c(f_c_a,f_c_b), names_prefix = "f_c_", names_to="layer",
                     values_to="f_c")
    params_est_f_p <- params_est %>%
        dplyr::select(patch=p, f_p_a, f_p_b) %>%
        tidyr::pivot_longer(cols=c(f_p_a,f_p_b), names_prefix = "f_p_", names_to="layer",
                     values_to="f_p")
    
    # Sanity check: system of equations
    # 1 = f_o_a + f_c_a + f_p_a and 1 = f_c_b + f_p_b
        test1 <- breakdown_to_pft_df %>%
        dplyr::select(-c(f_o)) %>%
        dplyr::left_join(params_est_f_o) %>%
        dplyr::left_join(params_est_f_c) %>%
        dplyr::left_join(params_est_f_p) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(should_be_1 = f_o + f_c + f_p) %>%
        dplyr::select(patch,layer,should_be_1) %>% distinct() %>% data.frame()
    # right column should be 1's
    assertthat::assert_that(round(sum(test1$should_be_1))==nrow(test1),
                            msg="pft fractions do not sum to 1 for at least one patch")
    
    # p_c_T*lai_T = f_c_a*lai_a + f_c_b*lai_b
    test2 <- breakdown_to_pft_df %>%
        dplyr::left_join(params_est_f_c) %>%
        dplyr::group_by(patch) %>%
        dplyr::reframe(cedar_lai_total = p_c_T*lai_T) %>% distinct() %>%
        dplyr::left_join(breakdown_to_pft_df %>%
                             dplyr::left_join(params_est_f_c) %>%
                             dplyr::group_by(patch,layer) %>%
                             dplyr::reframe(cedar_by_layer = f_c*lai_layer) %>%
                             dplyr::distinct() %>%
                             dplyr::group_by(patch) %>%
                             dplyr::summarize(cedar_tot = sum(cedar_by_layer)) ) %>%
        dplyr::mutate(diff=cedar_tot - cedar_lai_total)
    # right column should be as close to 0 as possible
    message(paste0("largest difference between allometrically and statistically derived cedar ",max(abs(test2$diff))))
    # assertthat::assert_that(max(abs(test2$diff)) < 0.1, #ais arbitrarily picked this threshold
    #                         msg="difference between allometrically and statistically derived cedar is too large for at least one patch")
    
    
    # p_p_T*lai_T = f_p_a*lai_a + f_p_b*lai_b
    test3 <- breakdown_to_pft_df %>%
        dplyr::left_join(params_est_f_p) %>%
        dplyr::group_by(patch) %>%
        dplyr::reframe(pine_lai_total = p_p_T*lai_T) %>% distinct() %>%
        dplyr::left_join(breakdown_to_pft_df %>%
                             dplyr::left_join(params_est_f_p) %>%
                             dplyr::group_by(patch,layer) %>%
                             dplyr::reframe(pine_by_layer = f_p*lai_layer) %>%
                             dplyr::distinct() %>%
                             dplyr::group_by(patch) %>%
                             dplyr::summarize(pine_tot = sum(pine_by_layer)) ) %>%
        dplyr::mutate(diff=pine_tot - pine_lai_total)
    # right column should be as close to 0 as possible
    message(paste0("largest difference between allometrically and statistically derived pine ",max(abs(test3$diff))))
    # assertthat::assert_that(max(abs(test3$diff)) < 0.1, #ais arbitrarily picked this threshold
    #                 msg="difference between allometrically and statistically derived cedar is too large for at least one patch")
    
    # ais are not all on par (diff is not 0), but good enough

    patch_all_df <- breakdown_to_pft_df %>%
            dplyr::select(-c(lai_layer)) %>%
            dplyr::left_join(params_est_f_c) %>%
            dplyr::left_join(params_est_f_p) %>%
            tidyr::pivot_longer(cols=c(f_o,f_c,f_p), names_to="pft", values_to="f") %>%
            dplyr::mutate(lai = f*lai_cohort) %>%
            dplyr::filter(lai > 0) %>%
            dplyr::group_by(patch) %>%
            dplyr::mutate(cohort_idx_new = row_number()) %>% 
            dplyr::ungroup() %>%
            dplyr::select(c(patch,pft,cohort_height,cohort_idx = cohort_idx_new,lai)) %>%
            dplyr::mutate(pft = case_when(
                pft == "f_p" ~ "pine_PFT",
                pft == "f_c" ~ "cedar_PFT",
                pft == "f_o" ~ "oak_PFT") ) %>%
            dplyr::left_join(allom_params) %>%
            dplyr::group_by(patch,pft,cohort_height,cohort_idx) %>%
            dplyr::reframe(lai = lai,
                    dbh = d1 * pmin(cohort_height,Hmax)^d2,
                    agb = a1 * dbh^a2,
                    leaf_biom  = x1 * pmin(cohort_height,Hmax)^x2,
                    n_stemdens = lai / (leaf_biom  * SLA)) %>%
            dplyr::ungroup() %>% distinct() %>% arrange(patch,cohort_idx) %>%
            dplyr::rename(lai_cohort=lai)
        
        # Sanity check: LAI
        breakdown_to_pft_df %>% group_by(patch) %>% 
            dplyr::summarize(lai_tot_a = sum(lai_cohort)) %>%
            dplyr::left_join(patch_all_df %>% group_by(patch) %>% 
            dplyr::summarize(lai_tot_b = sum(lai_cohort)) ) %>%
            dplyr::mutate(diff = lai_tot_a - lai_tot_b)
        #ais almost spot on (diff==0) - good enough!
    
    return(patch_all_df)
}



generate_cohort_patch_files <- function(site, year, data_int_path, biomass_path,
                                        ic_type, ic_type_path, plots_shp_path=NULL, classified_PFTs_path=NULL,
                                        lad_laz_clipped_path=NULL) {
    # Generate cohort (.css) and patch (.pss) files for FATES initialization
    
    allom_path <- file.path(data_int_path,site,"AllometricParameters_SOAPeqns_byMarcos.csv")
    allom_params <- read.csv(allom_path) #ais this doesn't exist for the other sites yet
    # ais how to access this data reproducibly? I just manually put this file in that folder
    
    if (ic_type=="field_inv_plots") {
        message("entered if statement")
        # Load cleaned individual-level data
        veg <- read.csv(file.path(biomass_path,"pp_veg_structure_IND_IBA_IAGB_live.csv"))%>% 
            rename(SLA_sytoan = SLA) %>% #to differentiate from SLA from allometry file
           #ais isn't this path name a global variable?, assign with this instead
            # already filtered to live plants
            
            # Filter to single stem if plant is an oak
            # veg %>% dplyr::count(individualID,pft) %>% filter(n>1)
            # for soap 2019, all individualIDs with > 1 count are shrub or oak
            
            dplyr::group_by(individualID) %>% 
            dplyr::slice(which.max(used_diameter)) %>% #508 expected
            ungroup() %>%
            
            #Assign pft
            mutate(pft = ifelse(growthForm == "small shrub" | 
                                    growthForm == "single shrub" ,"shrub_PFT",
                                ifelse(taxonID=="PIPO" | taxonID=="PILA" | taxonID=="PINUS", "pine_PFT", 
                                       
                                       ifelse(taxonID=="CADE27", "cedar_PFT", 
                                              
                                              ifelse(taxonID=="ABCO" | taxonID=="ABMA" , "fir_PFT", 
                                                     
                                                     ifelse(taxonID=="QUCH2" | taxonID=="QUKE"| taxonID=="QUERC", "oak_PFT", 
                                                            
                                                            ifelse(taxonID=="ARVIM" | taxonID=="CECU" | taxonID=="CEMOG" | 
                                                                       taxonID=="CEIN3" | taxonID=="RIRO" | taxonID=="FRCA6" | 
                                                                       taxonID=="RHIL" | taxonID=="AECA" | taxonID=="SANI4" | 
                                                                       taxonID=="CEANO", "shrub_PFT", "other"))))))) %>%
            # Filter out pft 'other'
            filter(pft != "other") %>%
            # Assign unknown plants to PFT - ais do this - all of them were shrubs anyway
            #filter(scientificName != "Unknown plant") %>% #20 plants
            left_join(allom_params) %>%
            dplyr::mutate(time = year, 
                   #subplotID_agg = ifelse(!is.na(subplotID),subplotID, pointID),
                   patch = paste0(plotID,"_central"), #paste0(plotID,"_",subplotID_agg),
                   index = row_number(),
                   dbh = used_diameter,
                   #ais or dbh should be used_diameter
                   n = 1/sampling_area,
                   agb = a1 * dbh^a2,
                   leaf_biom = x1 * pmin(height,Hmax)^x2,
                   lai = n * SLA * leaf_biom)
        # lai (layer 1) = stem density (layer 1) * ILA(as a funciton of z)
        # ILA = Bleaf * SLA , SLA is specific to height and PFT - for now just pft
        
        # Visualize / sanity check
        # by patch
        veg %>% group_by(patch) %>% summarize(sum_lai = sum(lai)) %>% 
            left_join( veg %>% 
                           group_by(patch, pft) %>% summarize(sum_lai = sum(lai)) %>% 
                           slice_max(sum_lai) %>% 
                           dplyr::select(-c(sum_lai)) %>% 
                           rename(pft_dom=pft) ) %>%
            ggplot(aes(sum_lai,fill=pft_dom)) + geom_histogram() + 
            xlab("LAI by patch (Inventory)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_lai.png"))) 
        
        veg %>% group_by(patch) %>% summarize(sum_n = sum(n)) %>% 
            left_join( veg %>% 
                           group_by(patch, pft) %>% summarize(sum_lai = sum(lai)) %>% 
                           slice_max(sum_lai) %>% 
                           dplyr::select(-c(sum_lai)) %>% 
                           rename(pft_dom=pft) ) %>%
            ggplot(aes(sum_n,fill=pft_dom)) + geom_histogram() + 
            xlab("Stem density (1/m2) by patch (Inventory)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_stemdens.png"))) 
        
        veg %>% group_by(patch) %>% summarize(sum_lb = sum(leaf_biom*n)) %>% 
            left_join( veg %>% 
                           group_by(patch, pft) %>% summarize(sum_lai = sum(lai)) %>% 
                           slice_max(sum_lai) %>% 
                           dplyr::select(-c(sum_lai)) %>% 
                           rename(pft_dom=pft) ) %>%
            ggplot(aes(sum_lb,fill=pft_dom)) + geom_histogram() + 
            xlab("Leaf biomass (kg/m2) by patch (Inventory)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_leafbiom.png"))) 
        
        veg %>% group_by(patch) %>% summarize(sum_agb = sum(agb*n)) %>% 
            left_join( veg %>% 
                           group_by(patch, pft) %>% summarize(sum_lai = sum(lai)) %>% 
                           slice_max(sum_lai) %>% 
                           dplyr::select(-c(sum_lai)) %>% 
                           rename(pft_dom=pft) ) %>%
            ggplot(aes(sum_agb,fill=pft_dom)) + geom_histogram() + 
            xlab("Aboveground Biomass (kg/m2) by patch (Inventory)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_agb.png"))) 
        
        # Create final files
        cohort_df_final <- veg %>%
            mutate(pft = case_when(
                pft == "pine_PFT" ~ 1,
                pft == "cedar_PFT" ~ 2,
                pft == "fir_PFT" ~ 3,
                pft == "shrub_PFT" ~ 4,
                pft == "oak_PFT" ~ 5 )) %>%
            dplyr::select(time, patch, index, dbh, height, pft, n)  %>%
            mutate(bdead = 0,
                   balive = 0,
                   avgRG = 0)
        
        # Patch file
        patch_df <- veg %>%
            distinct(patch) %>% arrange(patch) %>%
            mutate(area = 400/(400*nrow(veg %>% distinct(patch)))) %>% 
            #ais how to automate this to get 400 without saying 400?
            mutate( 
                time = year,
                trk = 0,
                age = 0,
                water = 0,
                fsc = 0, 
                stsc = 0, 
                stsl = 0, 
                ssc = 0, 
                psc = 0,
                msn = 0, 
                fsn = 0) %>%  
            relocate(area, .after=age) %>%
            relocate(patch, .after=time)
    
    } else { # then ic_type == rs_inv_plots or rs_random_plots
        
        # load plots shp
        plots_sf <- st_read(plots_shp_path)
        
        # load plots classified as PFT by percentage
        plot_by_pft_majority_raw <- read.csv(classified_PFTs_path)
        # ais why do some patches have no PFT classification? only a few e.g. random plot 162
        
        plot_by_pft_majority <- plot_by_pft_majority_raw %>%
            dplyr::rename(patch = shapeID) %>% 
            dplyr::mutate(patch = as.character(patch)) %>%
            group_by(patch) %>%
            # choose PFT with greatest coverage of patch, break tie with count 
            slice_max(tibble(pct, count), n = 1, with_ties = FALSE) 
        # to search for a specific patch: which(stringr::str_detect(laz_plot_paths,"_862"))
        
        pfts_by_cohort_wide <- plot_by_pft_majority_raw %>%
            dplyr::rename(pft = pred_PFT,
                   patch = shapeID) %>% 
            dplyr::mutate(patch = as.character(patch)) %>%
            dplyr::select(-c(count)) %>% 
            tidyr::pivot_wider(names_from = pft, values_from = pct) %>%
            dplyr::rename(p_c_T=cedar_PFT,
                   p_p_T=pine_PFT,
                   p_o_T=oak_PFT)
        pfts_by_cohort_wide[is.na(pfts_by_cohort_wide)] <- 0 
        
        lad_laz_plot_paths <- list.files(file.path(lad_laz_clipped_path), 
                                     pattern="*.laz", full.names = T) %>%
                tools::file_path_sans_ext()
        
        by_patch_df <- divide_patch_into_cohorts(lad_laz_plot_paths, ic_type)
        
        # Do patches in pft df match that in patch df?
        assertthat::assert_that(sum(unique(by_patch_df$patch) %in% #ais why do these two vectors have slightly different subsets in the rs_random_plots case?
                                    unique(plot_by_pft_majority$patch))>1,
                                msg="Patch names in pft df and patch df don't match")
        
        patch_all_df <- assign_pft_across_cohorts(by_patch_df, 
                                                         allom_params,
                                                         pfts_by_cohort_wide) 
        
        # Sanity check plots by patch
        patch_all_df %>% dplyr::group_by(patch,pft) %>% 
            dplyr::reframe(sum_lai = sum(lai_cohort)) %>% 
            ggplot(aes(sum_lai,fill=pft)) + geom_histogram() + 
            xlab("LAI by patch (lidar)") # LAI should be less than 10 for the most part
        ggsave(file.path(ic_type_path,paste0(ic_type,"_lai.png"))) 
        patch_all_df %>% dplyr::group_by(patch,pft) %>% 
            dplyr::reframe(sum_n = sum(n_stemdens)) %>% 
            ggplot(aes(sum_n,fill=pft)) + geom_histogram() + 
            xlab("Stem density (1/m2) by patch (lidar)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_stemdens.png"))) 
        patch_all_df %>% dplyr::group_by(patch,pft) %>% 
            dplyr::reframe(sum_lb = sum(leaf_biom *n_stemdens)) %>% 
            ggplot(aes(sum_lb,fill=pft)) + geom_histogram() + 
            xlab("Leaf biomass (kg/m2) by patch (lidar)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_leafbiom.png"))) 
        patch_all_df %>% dplyr::group_by(patch,pft) %>% 
            dplyr::reframe(sum_agb = sum(agb*n_stemdens)) %>% 
            ggplot(aes(sum_agb,fill=pft)) + geom_histogram() + 
            xlab("Aboveground biomass (kg/m2) by patch (lidar)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_agb.png"))) 
        
        #cohort_df_final[which(is.na(cohort_df_final$pft)),]
        # Cohort file
        cohort_df_final <- patch_all_df %>%
            dplyr::select(-c(lai_cohort,agb,leaf_biom )) %>%
            dplyr::rename(index = cohort_idx,
                   height= cohort_height,
                   n = n_stemdens) %>%
            dplyr::mutate(pft = case_when(
                pft == "pine_PFT" ~ 1,
                pft == "cedar_PFT" ~ 2,
                pft == "fir_PFT" ~ 3,
                pft == "shrub_PFT" ~ 4,
                pft == "oak_PFT" ~ 5 ),
                time = year,
                bdead = 0,
                balive = 0) %>%
            dplyr::mutate(patch = as.numeric(as.factor(patch))) %>%
            dplyr::relocate(time,patch) %>%
            dplyr::relocate(dbh, .after=index) %>%
            dplyr::relocate(height, .after=dbh) %>%
            dplyr::relocate(n, .after=pft) %>%
            dplyr::mutate(avgRG = 0)
        
        # Patch file
        patch_df <- data.frame( #each row is one subplot, or patch
            time = year,
            patch = as.character(1:nrow(plots_sf)),
            trk = 0,
            age = 0,
            water = 0,
            fsc = 0, #ais find realistic value
            stsc = 0, #ais find realistic value
            stsl = 0, #ais find realistic value
            ssc = 0, #ais find realistic value
            psc = 0,
            msn = 0, #ais find realistic value
            fsn = 0) %>%  #ais find realistic value
            dplyr::mutate(area = 1/nrow(plots_sf)) %>%
            dplyr::relocate(area, .after=age)
        
    } 
    
    # Save final initial conditions --------------------------------------------
    
    # wiki: https://github.com/EDmodel/ED2/wiki/Initial-conditions#files-types-and-
    # formats-for-nlied_init_mode6
    # FATES initial conditions are generated the same way as ED2
    
    if (ic_type == "field_inv_plots") {
        cohort_path <- file.path(ic_type_path,"cohort_ic_field_inv.css")
        patch_path <- file.path(ic_type_path,"patch_ic_field_inv.pss")
        
    } else if (ic_type=="rs_inv_plots") { 
        cohort_path <- file.path(ic_type_path,"cohort_ic_rs_inv_plots.css")
        patch_path <- file.path(ic_type_path,"patch_ic_rs_inv_plots.pss")
        
    } else if (ic_type == "rs_random_plots") {
        cohort_path <- file.path(ic_type_path,"cohort_ic_rs_random_plots.css")
        patch_path <- file.path(ic_type_path,"patch_ic_rs_random_plots.pss")
        
    } else {
        print("need to specify dataset type")
    }
    
    write.csv(cohort_df_final, cohort_path, row.names = FALSE)
    write.csv(patch_df, patch_path, row.names = FALSE)
    
    return(c(cohort_path, patch_path))
}
