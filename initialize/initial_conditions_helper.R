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


#ais replace this funciton below with sy-toan's clip_lidar_to_polygon_lidR
clip_norm_laz_to_shp <- function(site, year_inv, data_int_path, ic_type_path,
                                 plots_shp_path) {
    # Function to clip normalized lidar point cloud to plots being  
    # used to generate initial conditions
    
    laz_clipped_random_path <- file.path(ic_type_path, "clipped_to_plots")
    
    # Create directory if it doesn't exist                                     
    if (!dir.exists(laz_clipped_random_path)) {
        dir.create(laz_clipped_random_path) 
    }
    
    # Only clip laz if directory is empty
    if (length(list.files(laz_clipped_random_path)) == 0 ) {

        las_ctg <- lidR::readLAScatalog(file.path(data_int_path, site, year_inv,
                                            "normalized_lidar_tiles"))
        opt_output_files(las_ctg) <- file.path(laz_clipped_random_path,
                                               "plot_{PLOTID}_{XCENTER}_{YCENTER}")
        opt_laz_compression(las_ctg) <- TRUE 
        plots_sf <- sf::st_read(plots_shp_path)

        # Extract the plots. Returns a list of extracted point clouds in R
        random_plots_ctg <- lidR::clip_roi(las_ctg, plots_sf)
        #ais can I parallelize this? ^
        #ais address differences in tile count between normalized_lidar_tiles, aop_tiles, and this

    } else {
        # jump to bottom and return path
        message("normalized laz already clipped to plots")
    }
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




divide_patch_into_cohorts <- function(lad_laz_plot_paths, ic_type) {
    
    # Initialize dataframe
    by_patch_df <- data.frame()
    #height_step <- 0.5 #m #for height step in lad_csv
    
    for (c in 1:length(lad_laz_plot_paths)) {
        print(c)
        
        # Load leaf area density files
        lad_json <- rjson::fromJSON(file=paste0(lad_laz_plot_paths[c],"_lad.json"))
        lad_csv <- read.csv(paste0(lad_laz_plot_paths[c],"_lad.csv")) 
        
        if ( nrow(lad_csv)==0 | length(lad_json$layer_height)==0 ) {
            # In this case, the patch is flat and has no cohorts
            patch_temp <- data.frame() #ais change this so that an empty patch is added to patch file, even though no cohorts
            print(paste0("patch ",c," has nrow(lad_csv)==0 | length(lad_json$layer_height)==0"))
            # ais look through some of these cases with Marcos
            # note that the index c is not the same as the patch number
            # in some cases the csv has values, but json is empty - why?
            
            #ais but what about the cases where json has data but csv is empty?
            # ^ in these cases I got the warning "lad could not be calculated"
            # in a previous step. what's the difference?
            
        } else {
            # lad units are m2/m3
            if (is.null(lad_csv$X)) { #ais what is this column X?
                lad_csv <- lad_csv %>%
                    #cut off rows beyond highest layer_height
                    filter(z <= max(lad_json$layer_height)) #ais add 0.5 here?
            } else {
                lad_csv <- lad_csv %>%
                    dplyr::select(-c("X")) %>%
                    #cut off rows beyond highest layer_height
                    filter(z <= max(lad_json$layer_height)) 
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
            
            # ais when would 'central' be part of filename now?
            if (grepl("central", basename(lad_laz_plot_paths[c]))) { 
                patch_temp <- lad_csv %>%
                    filter(lad>0) %>%
                    # assign patch and cohort index
                    mutate(patch = basename(lad_laz_plot_paths[c])) %>%
                    left_join(cohort_idx_df)
            } else {
                if (ic_type=="rs_inv_plots") {
                    patch_temp <- lad_csv %>%
                        filter(lad>0) %>%
                        # assign patch # and cohort index
                        mutate(patch = sub(".*/(.*)$", "\\1", #".*_(.*?)\\.*", "\\1",
                                           lad_laz_plot_paths[c])) %>%
                        left_join(cohort_idx_df)
                } else if (ic_type=="rs_random_plots") {
                    patch_temp <- lad_csv %>%
                        filter(lad>0) %>%
                        # assign patch # and cohort index
                        mutate(patch = sub(".*/plot_(\\d+)_.*", "\\1", 
                                           lad_laz_plot_paths[c])) %>%
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

    pfts_by_cohort_wide$patch = as.character(pfts_by_cohort_wide$patch)

    # Sort out number of layers based on number of distinct Hmax's
    hmax_values = sort(unique(allom_params$Hmax))
    by_patch_df$layer <- base::cut(by_patch_df$cohort_height,
                             breaks = c(0,hmax_values), 
                             labels = c(letters[1:length(hmax_values)]))

    breakdown_to_pft_df <- by_patch_df %>%
        # lai = sum(layer thickness (diff_z) * LAD of layer)
        # lai (layer 1) = stem density (layer 1) * ILA(as a function of z)
        # ILA = Bleaf * SLA , SLA is specific to height and PFT - for now just pft
        dplyr::mutate(lai = lad*diff_z) %>%
        dplyr::group_by(patch,cohort_idx,cohort_height,layer) %>%
        dplyr::reframe(lai_cohort = sum(lai)) %>%
        dplyr::ungroup() %>%

        # assign layer LAI
        dplyr::group_by(patch, layer) %>%
        dplyr::mutate(lai_layer = sum(lai_cohort)) %>%
        dplyr::ungroup() %>%

        # assign total LAI for patch
        dplyr::group_by(patch) %>%
        dplyr::mutate( lai_patch_T = sum(lai_cohort) ) %>%
        dplyr::ungroup() %>%
        
        # percent of each pft per patch
        dplyr::left_join(pfts_by_cohort_wide) %>%
        # filter out rows where by_patch_df had a patch but pfts_by_cohort_wide doesn't
        tidyr::drop_na() %>%
        
        # known fractions
        #f_o_c = fraction of oaks in layer c
        #f_o_b = fraction of oaks in layer b
        #f_o_a = fraction of oaks in layer a
        mutate(f_s = ifelse(layer=="a", (p_s_T*lai_patch_T)/lai_layer, 0))
    
    # Estimate fraction of each PFT (f_p,f_c,f_o,f_f) for each layer
    message("Estimating PFT fractions across cohorts")
    
    # initialize
    params_est <- data.frame()
    for (p in unique(breakdown_to_pft_df$patch)) {
        
        temp_df <- breakdown_to_pft_df %>%
            dplyr::filter(patch == p)
        
        # when p="13", we have no cohorts in layer a, but we do have p_s_T - how to deal with this?
        # If a proportion of pixels are shrub, but there is no a layer, skip this patch
        if (sum(temp_df$p_s_T)>0 & nrow(temp_df %>% filter(layer=="a"))==0) {
            next
        } #ais talk to marcos about this
        
        # Isolate the lai total for each layer
        lai_a <- ifelse(nrow(temp_df %>% dplyr::filter(layer=="a"))==0, 0,
                        temp_df %>% dplyr::filter(layer=="a") %>%
                            dplyr::distinct(lai_layer) %>% dplyr::pull(lai_layer))
        lai_b <- ifelse(nrow(temp_df %>% dplyr::filter(layer=="b"))==0, 0,
                            temp_df %>% dplyr::filter(layer=="b") %>%
                                dplyr::distinct(lai_layer) %>% dplyr::pull(lai_layer))
        lai_c <- ifelse(nrow(temp_df %>% dplyr::filter(layer=="c"))==0, 0,
                            temp_df %>% dplyr::filter(layer=="c") %>%
                                dplyr::distinct(lai_layer) %>% dplyr::pull(lai_layer))
        lai_patch_T <- temp_df %>% dplyr::distinct(lai_patch_T) %>% dplyr::pull(lai_patch_T)
        p_p_T <- ifelse("p_p_T" %in% colnames(temp_df),temp_df %>% dplyr::distinct(p_p_T) %>% dplyr::pull(p_p_T),0)
        p_c_T <- ifelse("p_c_T" %in% colnames(temp_df),temp_df %>% dplyr::distinct(p_c_T) %>% dplyr::pull(p_c_T),0)
        p_f_T <- ifelse("p_f_T" %in% colnames(temp_df),temp_df %>% dplyr::distinct(p_f_T) %>% dplyr::pull(p_f_T),0)
        p_o_T <- ifelse("p_o_T" %in% colnames(temp_df),temp_df %>% dplyr::distinct(p_o_T) %>% dplyr::pull(p_o_T),0)
        p_s_T <- ifelse("p_s_T" %in% colnames(temp_df),temp_df %>% dplyr::distinct(p_s_T) %>% dplyr::pull(p_s_T),0)
            
        f_s_a <- ifelse(lai_a>0, temp_df %>% dplyr::filter(layer=="a") %>% dplyr::distinct(f_s) %>% dplyr::pull(f_s), 0)
        
        
        # system of equations
        # 1 - f_s_a = f_c_a + f_p_a + f_f_a + f_o_a  (f_s_a is known)
        #         1 = f_c_b + f_p_b + f_f_b + f_o_b
        #         1 = f_c_c + f_p_c + f_f_c
        # p_c_T*lai_patch_T = f_c_a*lai_a + f_c_b*lai_b + f_c_c*lai_c
        # p_p_T*lai_patch_T = f_p_a*lai_a + f_p_b*lai_b + f_p_c*lai_c
        # p_f_T*lai_patch_T = f_f_a*lai_a + f_f_b*lai_b + f_f_c*lai_c
        # p_o_T*lai_patch_T = f_o_a*lai_a + f_o_b*lai_b 
        
        # also an equation, but not included above
        # p_s_T*lai_patch_T = f_s_a*lai_a + f_s_b*lai_b
        # but since there is no shrub beyond layer a, f_o_b=f_o_c=0
        # p_s_T*lai_patch_T = f_s_a*lai_a
        
        # Write in matrix form
        # ais need to scale this out to include more PFTs too
        # column order: f_c_a,f_c_b,f_c_c,f_p_a,f_p_b,f_p_c,f_f_a,f_f_b,f_f_c,f_o_a,f_o_b
        X <- matrix(c(1,0,0,1,0,0,1,0,0,1,0, #a
                      0,1,0,0,1,0,0,1,0,0,1, #b
                      0,0,1,0,0,1,0,0,1,0,0, #c
                      lai_a,lai_b,lai_c,0,0,0,0,0,0,0,0,  #cedar
                      0,0,0,lai_a,lai_b,lai_c,0,0,0,0,0,  #pine
                      0,0,0,0,0,0,lai_a,lai_b,lai_c,0,0,  #fir
                      0,0,0,0,0,0,0,0,0,lai_a,lai_b),     #oak
                    7,11, byrow=TRUE)
        b_l <- rep(0,ncol(X))
        b_u <- rep(1,ncol(X))
        
        
        # if (f_s_a == 0) { #minority of cases
        #     #e.g., "patch 13"
        #     # sometimes there are no cohorts in layer a, so when this happens and
        #     # if p_s_T is not 0, we need to set f_s_a==0
            
        if (f_s_a > 1) { #minority of cases
            #e.g., "SOAP_021_central"
            # f_s_a should be between 0 and 1, but sometimes it is greater than 1
            # when f_s_a > 1, just set all of layer a to shrub,
            # and split up layer b between the rest of the PFTs 
            f_s_a <- 1
            # which then forces f_p_a = f_c_a = f_f_a = f_o_a = 0, and we need 
            # to recalibrate. Since shrub can only be in layer a, but the 
            # proportion of lai in layer a is less than p_s_T, then we set 
            # p_s_T = lai_a/lai_patch_T,and recalibrate the other PFTs accordingly
            rescal <- (1 - lai_a/lai_patch_T)/(p_p_T + p_c_T + p_f_T + p_o_T)
            if (rescal!="Inf") {
                if ("p_p_T" %in% colnames(temp_df)) { p_p_T <- rescal*p_p_T } 
                if ("p_c_T" %in% colnames(temp_df)) { p_c_T <- rescal*p_c_T } 
                if ("p_f_T" %in% colnames(temp_df)) { p_f_T <- rescal*p_f_T } 
                if ("p_o_T" %in% colnames(temp_df)) { p_o_T <- rescal*p_o_T } 
            }
            # so now the system of equations above should balance out
            
            # see system of equations above - the left hand side
            y <- matrix(c(1 - f_s_a, 1, 1, p_c_T*lai_patch_T, p_p_T*lai_patch_T,
                          p_f_T*lai_patch_T, p_o_T*lai_patch_T), 7, 1)
            
            # If key > 0, istate is as follows: the last contains the total
            # number of components at their bounds (the bound variables). 
            
            # The absolute values of the first nbound <- tail(istate,1) entries of
            # istate are the indices of these bound components of x. The sign of
            # istate[1:nbound] indicates whether x(abs(istate[1:nbound])) is at
            # its upper or lower bound. istate[1:nbound] is positive if the
            # component is at its upper bound, negative if the component is at
            # its lower bound. istate[(nbound+1):ncol(A)] contain the indices of
            # the components of x that are active (i.e. are expected to lie
            # strictly within their bounds). When key > 0, the routine initially
            # sets the active components to the averages of their upper and
            # lower bounds.
            
            sol_temp <- bvls(X, y, b_l, b_u) #, key=1, istate=c(-1,-3,2,4,2))
            dev_df <- data.frame(f_c_a=sol_temp$x[1], f_c_b=sol_temp$x[2], f_c_c=sol_temp$x[3], 
                                 f_p_a=sol_temp$x[4], f_p_b=sol_temp$x[5], f_p_c=sol_temp$x[6],
                                 f_f_a=sol_temp$x[7], f_f_b=sol_temp$x[8], f_f_c=sol_temp$x[9], 
                                 f_o_a=sol_temp$x[10], f_o_b=sol_temp$x[11], 
                                 b_l_1=b_l[1], b_l_2=b_l[2], b_l_3=b_l[3], 
                                 b_l_4=b_l[4], b_l_5=b_l[5], b_l_6=b_l[6],
                                 b_l_7=b_l[7], b_l_8=b_l[8], b_l_9=b_l[9], 
                                 b_l_10=b_l[10], b_l_11=b_l[11], 
                                 dev=sol_temp$deviance)
            i <- 1
            while (max(b_l) < 1) {
                # fix f_s_a=1
                ind_fixed <- which(dev_df[i,c(2,4)] == min(dev_df[i,c(2,4)]))
                b_l[ind_fixed] <- b_l[ind_fixed] + 0.01
                sol_temp <- bvls(X, y, b_l, b_u)#, key=1, istate=c(-1,-3,2,4,2))
                
                # Add next row
                dev_df <- rbind(dev_df,
                                c(x_1=sol_temp$x[1], x_2=sol_temp$x[2], x_3=sol_temp$x[3],
                                  x_4=sol_temp$x[4], x_5=sol_temp$x[5], x_6=sol_temp$x[6],
                                  x_7=sol_temp$x[7], x_8=sol_temp$x[8], x_9=sol_temp$x[9],
                                  x_10=sol_temp$x[10], x_11=sol_temp$x[11], 
                                  b_l_1=b_l[1], b_l_2=b_l[2], b_l_3=b_l[3], 
                                  b_l_4=b_l[4], b_l_5=b_l[5], b_l_6=b_l[6],
                                  b_l_7=b_l[7], b_l_8=b_l[8], b_l_9=b_l[9], 
                                  b_l_10=b_l[10], b_l_11=b_l[11], 
                                  dev=sol_temp$deviance))
                i <- i+1
            }
        } else { #majority of cases
            #e.g., "SOAP_002_central"
            
            y <- matrix(c(1 - f_s_a, 1, 1, p_c_T*lai_patch_T, p_p_T*lai_patch_T,
                          p_f_T*lai_patch_T, p_o_T*lai_patch_T), 7, 1)
            
            sol_temp <- bvls(X, y, b_l, b_u) 
            dev_df <- data.frame(f_c_a=sol_temp$x[1], f_c_b=sol_temp$x[2], f_c_c=sol_temp$x[3], 
                                 f_p_a=sol_temp$x[4], f_p_b=sol_temp$x[5], f_p_c=sol_temp$x[6],
                                 f_f_a=sol_temp$x[7], f_f_b=sol_temp$x[8], f_f_c=sol_temp$x[9], 
                                 f_o_a=sol_temp$x[10], f_o_b=sol_temp$x[11], 
                                 b_l_1=b_l[1], b_l_2=b_l[2], b_l_3=b_l[3], 
                                 b_l_4=b_l[4], b_l_5=b_l[5], b_l_6=b_l[6],
                                 b_l_7=b_l[7], b_l_8=b_l[8], b_l_9=b_l[9], 
                                 b_l_10=b_l[10], b_l_11=b_l[11], 
                                 dev=sol_temp$deviance)
            i <- 1
            while (max(b_l) < 1) {
                ind_fixed <- which(dev_df[i,1:4] == min(dev_df[i,1:4]))
                b_l[ind_fixed] <- b_l[ind_fixed] + 0.01
                sol_temp <- bvls(X, y, b_l, b_u) #solve(X, y)
                
                # Add next row
                dev_df <- rbind(dev_df,
                                c(x_1=sol_temp$x[1], x_2=sol_temp$x[2], x_3=sol_temp$x[3],
                                  x_4=sol_temp$x[4], x_5=sol_temp$x[5], x_6=sol_temp$x[6],
                                  x_7=sol_temp$x[7], x_8=sol_temp$x[8], x_9=sol_temp$x[9],
                                  x_10=sol_temp$x[10], x_11=sol_temp$x[11], 
                                  b_l_1=b_l[1], b_l_2=b_l[2], b_l_3=b_l[3], 
                                  b_l_4=b_l[4], b_l_5=b_l[5], b_l_6=b_l[6],
                                  b_l_7=b_l[7], b_l_8=b_l[8], b_l_9=b_l[9], 
                                  b_l_10=b_l[10], b_l_11=b_l[11], 
                                  dev=sol_temp$deviance))
                i <- i+1
            }
            
        }
        
        params_temp <- cbind(p, f_s_a, dev_df %>%
                                 dplyr::filter(dev==min(dev_df$dev)) %>%
                                 dplyr::slice_head(n=1)  )
        
        params_est <- rbind(params_est ,params_temp)
    }
    
    params_est_f_s <- params_est %>%
        dplyr::select(patch=p, f_s_a) %>%
        mutate(f_s_b = 0) %>%
        tidyr::pivot_longer(cols=c(f_s_a,f_s_b), names_prefix = "f_s_", names_to="layer",
                            values_to="f_s")
    params_est_f_o <- params_est %>%
        dplyr::select(patch=p, f_o_a, f_o_b) %>%
        tidyr::pivot_longer(cols=c(f_o_a,f_o_b), names_prefix = "f_o_", names_to="layer",
                            values_to="f_o")
    params_est_f_c <- params_est %>%
        dplyr::select(patch=p, f_c_a, f_c_b, f_c_c) %>%
        tidyr::pivot_longer(cols=c(f_c_a,f_c_b,f_c_c), names_prefix = "f_c_", names_to="layer",
                            values_to="f_c")
    params_est_f_p <- params_est %>%
        dplyr::select(patch=p, f_p_a, f_p_b, f_p_c) %>%
        tidyr::pivot_longer(cols=c(f_p_a, f_p_b, f_p_c), names_prefix = "f_p_", names_to="layer",
                            values_to="f_p")
    params_est_f_f <- params_est %>%
        dplyr::select(patch=p, f_f_a, f_f_b, f_f_c) %>%
        tidyr::pivot_longer(cols=c(f_f_a, f_f_b, f_f_c), names_prefix = "f_f_", names_to="layer",
                            values_to="f_f")
    
    # Sanity check: system of equations
    # 1 = f_o_a + f_c_a + f_p_a and 1 = f_c_b + f_p_b
    test1 <- breakdown_to_pft_df %>%
        dplyr::select(-c(f_s)) %>%
        dplyr::left_join(params_est_f_o) %>%
        dplyr::left_join(params_est_f_c) %>%
        dplyr::left_join(params_est_f_p) %>%
        dplyr::left_join(params_est_f_f) %>%
        dplyr::left_join(params_est_f_s) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(should_be_1 = f_o + f_c + f_p + f_f + f_s) %>%
        dplyr::select(patch,layer,should_be_1) %>% distinct() %>% 
        tidyr::drop_na() %>% data.frame()
    # right column should be 1's
    assertthat::assert_that(round(floor(sum(test1$should_be_1)))==nrow(test1),
                            msg="pft fractions do not sum to 1 for at least one patch")
    
    # p_c_T*lai_patch_T = f_c_a*lai_a + f_c_b*lai_b
    test2 <- breakdown_to_pft_df %>%
        dplyr::left_join(params_est_f_c) %>%
        dplyr::group_by(patch) %>%
        dplyr::reframe(cedar_lai_patch_Total = p_c_T*lai_patch_T) %>% distinct() %>%
        dplyr::left_join(breakdown_to_pft_df %>%
                             dplyr::left_join(params_est_f_c) %>%
                             dplyr::group_by(patch,layer) %>%
                             dplyr::reframe(cedar_by_layer = f_c*lai_layer) %>%
                             dplyr::distinct() %>%
                             dplyr::group_by(patch) %>%
                             dplyr::summarize(cedar_tot = sum(cedar_by_layer)) ) %>%
        tidyr::drop_na() %>% 
        dplyr::mutate(diff=cedar_tot - cedar_lai_patch_Total)
    # right column should be as close to 0 as possible
    message(paste0("largest LAI difference between allometrically and statistically derived cedar ",max(abs(test2$diff))))
    # assertthat::assert_that(max(abs(test2$diff)) < 0.1, #ais arbitrarily picked this threshold
    #                         msg="difference between allometrically and statistically derived cedar is too large for at least one patch")
    
    
    # p_p_T*lai_patch_T = f_p_a*lai_a + f_p_b*lai_b
    test3 <- breakdown_to_pft_df %>%
        dplyr::left_join(params_est_f_p) %>%
        dplyr::group_by(patch) %>%
        dplyr::reframe(pine_lai_patch_Total = p_p_T*lai_patch_T) %>% distinct() %>%
        dplyr::left_join(breakdown_to_pft_df %>%
                             dplyr::left_join(params_est_f_p) %>%
                             dplyr::group_by(patch,layer) %>%
                             dplyr::reframe(pine_by_layer = f_p*lai_layer) %>%
                             dplyr::distinct() %>%
                             dplyr::group_by(patch) %>%
                             dplyr::summarize(pine_tot = sum(pine_by_layer)) ) %>%
        tidyr::drop_na() %>% 
        dplyr::mutate(diff=pine_tot - pine_lai_patch_Total)
    # right column should be as close to 0 as possible
    message(paste0("largest LAI difference between allometrically and statistically derived pine ",max(abs(test3$diff))))
    # assertthat::assert_that(max(abs(test3$diff)) < 0.1, #ais arbitrarily picked this threshold
    #                 msg="difference between allometrically and statistically derived cedar is too large for at least one patch")
    
    # ais are not all on par (diff is not 0), but good enough
    
    patch_all_df <- breakdown_to_pft_df %>%
        dplyr::select(-c(lai_layer)) %>%
        dplyr::left_join(params_est_f_c) %>%
        dplyr::left_join(params_est_f_p) %>%
        dplyr::left_join(params_est_f_f) %>%
        dplyr::left_join(params_est_f_o) %>%
        tidyr::pivot_longer(cols=c(f_o,f_c,f_p,f_s,f_f), names_to="pft", values_to="f") %>%
        dplyr::mutate(lai = f*lai_cohort) %>%
        dplyr::filter(lai > 0) %>%
        dplyr::group_by(patch) %>%
        dplyr::mutate(cohort_idx_new = row_number()) %>% 
        dplyr::ungroup() %>%
        dplyr::select(c(patch,pft,cohort_height,cohort_idx = cohort_idx_new,lai)) %>%
        dplyr::mutate(pft = case_when(
            pft == "f_p" ~ "pine",
            pft == "f_f" ~ "fir",
            pft == "f_c" ~ "cedar",
            pft == "f_o" ~ "oak",
            pft == "f_s" ~ "shrub") ) %>%
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
        dplyr::summarize(lai_patch_Tot_a = sum(lai_cohort)) %>%
        dplyr::left_join(patch_all_df %>% group_by(patch) %>% 
                             dplyr::summarize(lai_patch_Tot_b = sum(lai_cohort)) ) %>%
        tidyr::drop_na() %>%
        dplyr::mutate(diff = lai_patch_Tot_a - lai_patch_Tot_b)
    #ais almost spot on (diff==0) - good enough!
    
    return(patch_all_df)
}

# old assign_pft_across_cohorts <- function(by_patch_df, allom_params, 
#                                       pfts_by_cohort_wide) {
#     print(colnames(by_patch_df))
#     breakdown_to_pft_df <- by_patch_df %>%
#         # lai = sum(layer thickness (diff_z) * LAD of layer)
#         # lai (layer 1) = stem density (layer 1) * ILA(as a funciton of z)
#         # ILA = Bleaf * SLA , SLA is specific to height and PFT - for now just pft
#         dplyr::mutate(lai = lad*diff_z) %>%
#         dplyr::group_by(patch,cohort_idx,cohort_height) %>%
#         dplyr::reframe(lai_cohort = sum(lai)) %>%
        
#         # assign Hmax layer
#         dplyr::mutate(layer = ifelse(cohort_height <= allom_params %>%
#                                          filter(pft=="oak") %>% pull(Hmax),
#                                      "a", "b")) %>%
        
#         # assign layer LAI
#         dplyr::group_by(patch, layer) %>%
#         dplyr::mutate(lai_layer = sum(lai_cohort)) %>%
#         dplyr::ungroup() %>%
        
#         # assign total LAI for patch
#         dplyr::group_by(patch) %>%
#         dplyr::mutate( lai_T = sum(lai_cohort) ) %>%
#         dplyr::ungroup() %>%
        
#         # percent of each pft per patch
#         dplyr::left_join(pfts_by_cohort_wide) %>%
#         # filter out rows where by_patch_df had a patch but pfts_by_cohort_wide doesn't
#         drop_na() %>%
        
#         # known fractions
#         #f_o_b = fraction of oaks in layer b
#         #f_o_a = fraction of oaks in layer a
#         mutate(f_o = ifelse(layer=="b", 0, (p_o_T*lai_T)/lai_layer))
    
#     # Estimate f_c and f_p for each layer
#     message("Estimating PFT fractions across cohorts")
    
#     # initialize
#     params_est <- data.frame()
#     for (p in unique(breakdown_to_pft_df$patch)) {
#         temp_df <- breakdown_to_pft_df %>%
#             dplyr::filter(patch == p)
        
#         lai_a <- temp_df %>% dplyr::filter(layer=="a") %>% dplyr::distinct(lai_layer) %>% dplyr::pull(lai_layer)
#         lai_b <- ifelse(nrow(temp_df %>% dplyr::filter(layer=="b"))==0, 0,
#                         temp_df %>% dplyr::filter(layer=="b") %>%
#                             dplyr::distinct(lai_layer) %>% dplyr::pull(lai_layer))
#         lai_T <- temp_df %>% dplyr::distinct(lai_T) %>% dplyr::pull(lai_T)
#         p_c_T <- temp_df %>% dplyr::distinct(p_c_T) %>% dplyr::pull(p_c_T)
#         p_f_T <- temp_df %>% dplyr::distinct(p_f_T) %>% dplyr::pull(p_f_T)
#         p_p_T <- temp_df %>% dplyr::distinct(p_p_T) %>% dplyr::pull(p_p_T)
#         f_o_a <- temp_df %>% dplyr::filter(layer=="a") %>% dplyr::distinct(f_o) %>% pull(f_o)
        
#         # system of equations
#         # 1 - f_o_a = f_c_a + f_p_a... the same as 1 = f_c_a + f_p_a + f_o_a
#         #         1 = f_c_b + f_p_b
#         # p_c_T*lai_T = f_c_a*lai_a + f_c_b*lai_b
#         # p_p_T*lai_T = f_p_a*lai_a + f_p_b*lai_b
        
#         # also an equation, but not included above
#         # p_o_T*lai_T = f_o_a*lai_a + f_o_b*lai_b
#         # but since there is no oak in layer b, f_o_b=0
#         # p_o_T*lai_T = f_o_a*lai_a
        
#         # Write in matrix form
#         # column order: f_c_a, f_c_b, f_p_a, f_p_b
#         X <- matrix(c(1,0,1,0,
#                       0,1,0,1,
#                       lai_a,lai_b,0,0,
#                       0,0,lai_a,lai_b), 4,4, byrow=TRUE)
#         b_l <- rep(0,ncol(X))
#         b_u <- rep(1,ncol(X))
        
#         # ais need to scale this out to include more PFTs too
#         # # column order: f_c_a, f_c_b, f_f_a, f_f_b,f_p_a, f_p_b
#         # X <- matrix(c(1,0,1,0,1,0,
#         #               0,1,0,1,0,1,
#         #               lai_a,lai_b,0,0,0,0,
#         #               0,0,lai_a,lai_b,0,0,
#         #               0,0,0,0,0,0,#
#         #               0,0,0,0,lai_a,lai_b), 6,6, byrow=TRUE)
#         # b_l <- rep(0,ncol(X))
#         # b_u <- rep(1,ncol(X))
        
#         if (f_o_a > 1) { #minority of cases
#             #e.g., "SOAP_021_central"
#             # f_o_a should be between 0 and 1, but sometimes it is greater than 1
#             # when f_o_a > 1, just set all of layer a to oak,
#             # and split up layer b between pine and cedar
#             f_o_a <- 1
#             # which then forces f_p_a = f_c_a = 0, and we need to recalibrate.
#             # since oak can only be in layer a, but the proportion of lai in
#             # layer a is less than p_o_T, then we set p_o_T = lai_a/lai_T,
#             # and recalibrate p_c_T and p_p_T accordingly
#             rescal <- (1 - lai_a/lai_T)/(p_p_T + p_c_T)
#             if (rescal!="Inf") {
#                 p_c_T <- rescal*p_c_T
#                 p_f_T <- rescal*p_f_T
#                 p_p_T <- rescal*p_p_T
#             }
#             # so now the system of equations above should balance out
            
#             y <- matrix(c(1 - f_o_a, 1, p_c_T*lai_T, p_p_T*lai_T), 4, 1)
            
#             sol_temp <- bvls(X, y, b_l, b_u, key=1, istate=c(-1,-3,2,4,2))
#             dev_df <- data.frame(f_c_a=sol_temp$x[1], f_c_b=sol_temp$x[2],
#                                  f_p_a=sol_temp$x[3], f_p_b=sol_temp$x[4],
#                                  b_l_1=b_l[1], b_l_2=b_l[2],
#                                  b_l_3=b_l[3], b_l_4=b_l[4],
#                                  dev=sol_temp$deviance)
#             i <- 1
#             while (max(b_l) < 1) {
#                 ind_fixed <- which(dev_df[i,c(2,4)] == min(dev_df[i,c(2,4)]))
#                 b_l[ind_fixed] <- b_l[ind_fixed] + 0.01
#                 sol_temp <- bvls(X, y, b_l, b_u, key=1, istate=c(-1,-3,2,4,2))
                
#                 # Add next row
#                 dev_df <- rbind(dev_df,
#                                 c(x_1=sol_temp$x[1], x_2=sol_temp$x[2],
#                                   x_3=sol_temp$x[3], x_4=sol_temp$x[4],
#                                   b_l_1=b_l[1], b_l_2=b_l[2],
#                                   b_l_3=b_l[3], b_l_4=b_l[4],
#                                   dev=sol_temp$deviance))
#                 i <- i+1
#             }
#         } else { #majority of cases
#             #e.g., "SOAP_002_central"
            
#             y <- matrix(c(1 - f_o_a, 1, p_c_T*lai_T, p_p_T*lai_T), 4, 1)
            
#             sol_temp <- bvls(X, y, b_l, b_u)
#             dev_df <- data.frame(f_c_a=sol_temp$x[1], f_c_b=sol_temp$x[2],
#                                  f_p_a=sol_temp$x[3], f_p_b=sol_temp$x[4],
#                                  b_l_1=b_l[1], b_l_2=b_l[2],
#                                  b_l_3=b_l[3], b_l_4=b_l[4],
#                                  dev=sol_temp$deviance)
#             i <- 1
#             while (max(b_l) < 1) {
#                 ind_fixed <- which(dev_df[i,1:4] == min(dev_df[i,1:4]))
#                 b_l[ind_fixed] <- b_l[ind_fixed] + 0.01
#                 sol_temp <- bvls(X, y, b_l, b_u) #solve(X, y)
                
#                 # Add next row
#                 dev_df <- rbind(dev_df,
#                                 c(x_1=sol_temp$x[1], x_2=sol_temp$x[2],
#                                   x_3=sol_temp$x[3], x_4=sol_temp$x[4],
#                                   b_l_1=b_l[1], b_l_2=b_l[2],
#                                   b_l_3=b_l[3], b_l_4=b_l[4],
#                                   dev=sol_temp$deviance))
#                 i <- i+1
#             }
            
#         }
        
#         params_temp <- cbind(p, f_o_a, dev_df %>%
#                                  dplyr::filter(dev==min(dev_df$dev)) %>%
#                                  dplyr::slice_head(n=1)  )
        
#         params_est <- rbind(params_est ,params_temp)
#     }
    
#     params_est_f_o <- params_est %>%
#         dplyr::select(patch=p, f_o_a) %>%
#         mutate(f_o_b = 0) %>%
#         tidyr::pivot_longer(cols=c(f_o_a,f_o_b), names_prefix = "f_o_", names_to="layer",
#                             values_to="f_o")
#     params_est_f_c <- params_est %>%
#         dplyr::select(patch=p, f_c_a, f_c_b) %>%
#         tidyr::pivot_longer(cols=c(f_c_a,f_c_b), names_prefix = "f_c_", names_to="layer",
#                             values_to="f_c")
#     params_est_f_p <- params_est %>%
#         dplyr::select(patch=p, f_p_a, f_p_b) %>%
#         tidyr::pivot_longer(cols=c(f_p_a,f_p_b), names_prefix = "f_p_", names_to="layer",
#                             values_to="f_p")
    
#     # Sanity check: system of equations
#     # 1 = f_o_a + f_c_a + f_p_a and 1 = f_c_b + f_p_b
#     test1 <- breakdown_to_pft_df %>%
#         dplyr::select(-c(f_o)) %>%
#         dplyr::left_join(params_est_f_o) %>%
#         dplyr::left_join(params_est_f_c) %>%
#         dplyr::left_join(params_est_f_p) %>%
#         dplyr::ungroup() %>%
#         dplyr::mutate(should_be_1 = f_o + f_c + f_p) %>%
#         dplyr::select(patch,layer,should_be_1) %>% distinct() %>% data.frame()
#     # right column should be 1's
#     assertthat::assert_that(round(sum(test1$should_be_1))==nrow(test1),
#                             msg="pft fractions do not sum to 1 for at least one patch")
    
#     # p_c_T*lai_T = f_c_a*lai_a + f_c_b*lai_b
#     test2 <- breakdown_to_pft_df %>%
#         dplyr::left_join(params_est_f_c) %>%
#         dplyr::group_by(patch) %>%
#         dplyr::reframe(cedar_lai_total = p_c_T*lai_T) %>% distinct() %>%
#         dplyr::left_join(breakdown_to_pft_df %>%
#                              dplyr::left_join(params_est_f_c) %>%
#                              dplyr::group_by(patch,layer) %>%
#                              dplyr::reframe(cedar_by_layer = f_c*lai_layer) %>%
#                              dplyr::distinct() %>%
#                              dplyr::group_by(patch) %>%
#                              dplyr::summarize(cedar_tot = sum(cedar_by_layer)) ) %>%
#         dplyr::mutate(diff=cedar_tot - cedar_lai_total)
#     # right column should be as close to 0 as possible
#     message(paste0("largest difference between allometrically and statistically derived cedar ",max(abs(test2$diff))))
#     # assertthat::assert_that(max(abs(test2$diff)) < 0.1, #ais arbitrarily picked this threshold
#     #                         msg="difference between allometrically and statistically derived cedar is too large for at least one patch")
    
    
#     # p_p_T*lai_T = f_p_a*lai_a + f_p_b*lai_b
#     test3 <- breakdown_to_pft_df %>%
#         dplyr::left_join(params_est_f_p) %>%
#         dplyr::group_by(patch) %>%
#         dplyr::reframe(pine_lai_total = p_p_T*lai_T) %>% distinct() %>%
#         dplyr::left_join(breakdown_to_pft_df %>%
#                              dplyr::left_join(params_est_f_p) %>%
#                              dplyr::group_by(patch,layer) %>%
#                              dplyr::reframe(pine_by_layer = f_p*lai_layer) %>%
#                              dplyr::distinct() %>%
#                              dplyr::group_by(patch) %>%
#                              dplyr::summarize(pine_tot = sum(pine_by_layer)) ) %>%
#         dplyr::mutate(diff=pine_tot - pine_lai_total)
#     # right column should be as close to 0 as possible
#     message(paste0("largest difference between allometrically and statistically derived pine ",max(abs(test3$diff))))
#     # assertthat::assert_that(max(abs(test3$diff)) < 0.1, #ais arbitrarily picked this threshold
#     #                 msg="difference between allometrically and statistically derived cedar is too large for at least one patch")
    
#     # ais are not all on par (diff is not 0), but good enough
    
#     patch_all_df <- breakdown_to_pft_df %>%
#         dplyr::select(-c(lai_layer)) %>%
#         dplyr::left_join(params_est_f_c) %>%
#         dplyr::left_join(params_est_f_p) %>%
#         tidyr::pivot_longer(cols=c(f_o,f_c,f_p), names_to="pft", values_to="f") %>%
#         dplyr::mutate(lai = f*lai_cohort) %>%
#         dplyr::filter(lai > 0) %>%
#         dplyr::group_by(patch) %>%
#         dplyr::mutate(cohort_idx_new = row_number()) %>% 
#         dplyr::ungroup() %>%
#         dplyr::select(c(patch,pft,cohort_height,cohort_idx = cohort_idx_new,lai)) %>%
#         dplyr::mutate(pft = case_when(
#             pft == "f_p" ~ "pine",
#             pft == "f_c" ~ "cedar",
#             pft == "f_o" ~ "oak") ) %>%
#         dplyr::left_join(allom_params) %>%
#         dplyr::group_by(patch,pft,cohort_height,cohort_idx) %>%
#         dplyr::reframe(lai = lai,
#                        dbh = d1 * pmin(cohort_height,Hmax)^d2,
#                        agb = a1 * dbh^a2,
#                        leaf_biom  = x1 * pmin(cohort_height,Hmax)^x2,
#                        n_stemdens = lai / (leaf_biom  * SLA)) %>%
#         dplyr::ungroup() %>% distinct() %>% arrange(patch,cohort_idx) %>%
#         dplyr::rename(lai_cohort=lai)
    
#     # Sanity check: LAI
#     breakdown_to_pft_df %>% group_by(patch) %>% 
#         dplyr::summarize(lai_tot_a = sum(lai_cohort)) %>%
#         dplyr::left_join(patch_all_df %>% group_by(patch) %>% 
#                              dplyr::summarize(lai_tot_b = sum(lai_cohort)) ) %>%
#         dplyr::mutate(diff = lai_tot_a - lai_tot_b)
#     #ais almost spot on (diff==0) - good enough!
    
#     return(patch_all_df)
# }


generate_pcss <- function(site, year, data_int_path, biomass_path,
                                        ic_type, ic_type_path, plots_shp_path=NULL, classified_PFTs_path=NULL,
                                        lad_laz_clipped_path=NULL) {
    # Generate cohort (.css) and patch (.pss) files for FATES initialization
    
    allom_path <- file.path(data_int_path,"AllometricParameters_SOAPeqns_byMarcos_ais.csv") 
    #ais what to do about this for the otehr sites? move the file to SJER and TEAK manually
    allom_params <- read.csv(allom_path) #ais this doesn't exist for the other sites yet
    # ais how to access this data reproducibly? I just manually put this file in that folder
    
    if (ic_type=="field_inv_plots") {
        # Load cleaned individual-level data
        veg <- read.csv(biomass_path)%>% 
            rename(SLA_sytoan = SLA) %>% #to differentiate from SLA from allometry file
            
            # If an individual is multistem, select only the largest stem (usually oak)
            dplyr::group_by(individualID) %>% 
            dplyr::slice(which.max(used_diameter)) %>% 
            ungroup() %>%
            
            #Assign pft
            # ais how does this work for non-NEON inventory data?
            mutate(pft = ifelse(growthForm == "small shrub" | 
                                    growthForm == "single shrub" ,"shrub",
                                ifelse(taxonID=="PIPO" | taxonID=="PILA" | taxonID=="PINUS", "pine", 
                                       
                                       ifelse(taxonID=="CADE27", "cedar", 
                                              
                                              ifelse(taxonID=="ABCO" | taxonID=="ABMA" , "fir", 
                                                     
                                                     ifelse(taxonID=="QUCH2" | taxonID=="QUKE"| taxonID=="QUERC", "oak", 
                                                            
                                                            ifelse(taxonID=="ARVIM" | taxonID=="CECU" | taxonID=="CEMOG" | 
                                                                       taxonID=="CEIN3" | taxonID=="RIRO" | taxonID=="FRCA6" | 
                                                                       taxonID=="RHIL" | taxonID=="AECA" | taxonID=="SANI4" | 
                                                                       taxonID=="CEANO", "shrub", "other"))))))) %>%
            # Filter out pft 'other'
            filter(pft != "other") %>%
            # Assign unknown plants to PFT - ais do this - all of them were shrubs anyway
            #filter(scientificName != "Unknown plant") %>% #20 plants
            left_join(allom_params) %>%
            dplyr::mutate(time = year, 
                          patch = ifelse(is.na(subplotID),paste0(plotID,"_central"), paste0(plotID,"_",subplotID)),
                          index = row_number(),
                          dbh = used_diameter,
                          n = 1/sampling_area,
                          agb = a1 * dbh^a2,
                          leaf_biom = x1 * pmin(height,Hmax)^x2, #ais can we calculate leaf_biomass from the neon biomass column?
                          lai = n * SLA * leaf_biom)
            # lai (layer 1) = stem density (layer 1) * ILA(as a funciton of z)
            # ILA = Bleaf * SLA , SLA is specific to height and PFT - for now just pft

        # Visualize / sanity check
        # ais check that i'm generating these plots the way i intend
        # by patch
        veg %>% group_by(patch) %>% summarize(sum_lai = sum(lai,na.rm=T)) %>% 
            left_join( veg %>% 
                        group_by(patch, pft) %>% summarize(sum_lai = sum(lai,na.rm=T)) %>% 
                        slice_max(sum_lai) %>% 
                        dplyr::select(-c(sum_lai)) %>% 
                        rename(pft_dom=pft) ) %>%
            ggplot(aes(sum_lai,fill=pft_dom)) + geom_histogram() + 
            xlab("LAI by patch (Inventory)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_lai.png"))) 

        veg %>% group_by(patch) %>% summarize(sum_n = sum(n,na.rm=T)) %>% 
            left_join( veg %>% 
                        group_by(patch, pft) %>% summarize(sum_lai = sum(lai,na.rm=T)) %>% 
                        slice_max(sum_lai) %>% 
                        dplyr::select(-c(sum_lai)) %>% 
                        rename(pft_dom=pft) ) %>%
            ggplot(aes(sum_n,fill=pft_dom)) + geom_histogram() + 
            xlab("Stem density (1/m2) by patch (Inventory)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_stemdens.png"))) 

        veg %>% group_by(patch) %>% summarize(sum_lb = sum(leaf_biom*n,na.rm=T)) %>% 
            left_join( veg %>% 
                        group_by(patch, pft) %>% summarize(sum_lai = sum(lai,na.rm=T)) %>% 
                        slice_max(sum_lai) %>% 
                        dplyr::select(-c(sum_lai)) %>% 
                        rename(pft_dom=pft) ) %>%
            ggplot(aes(sum_lb,fill=pft_dom)) + geom_histogram() + 
            xlab("Leaf biomass (kg/m2) by patch (Inventory)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_leafbiom.png"))) 

        veg %>% group_by(patch) %>% summarize(sum_agb = sum(agb*n,na.rm=T)) %>% 
            left_join( veg %>% 
                        group_by(patch, pft) %>% summarize(sum_lai = sum(lai,na.rm=T)) %>% 
                        slice_max(sum_lai) %>% 
                        dplyr::select(-c(sum_lai)) %>% 
                        rename(pft_dom=pft) ) %>%
            ggplot(aes(sum_agb,fill=pft_dom)) + geom_histogram() + 
            xlab("Aboveground Biomass (kg/m2) by patch (Inventory)")
        ggsave(file.path(ic_type_path,paste0(ic_type,"_agb.png"))) 

        # Create final files
        cohort_df_initial <- veg %>%
            mutate(pft = case_when(
                pft == "pine" ~ 1,
                pft == "cedar" ~ 2,
                pft == "fir" ~ 3,
                pft == "shrub" ~ 4,
                pft == "oak" ~ 5 )) %>%
            dplyr::select(time, patch, index, dbh, height, pft, n)  %>%
            mutate(bdead = 0,
                balive = 0,
                avgRG = 0)

        # Patch file
        patch_df_initial <- veg %>%
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
        plots_sf <- sf::st_read(plots_shp_path)
        
        # load plots classified as PFT by percentage
        classified_PFTs <- read.csv(classified_PFTs_path)
        # ais why do some patches have no PFT classification? only a few e.g. random plot 15, 162
        
        pfts_by_cohort_wide_raw <- classified_PFTs %>%
            dplyr::rename(pft = pred_PFT,
                          patch = shapeID) %>% 
            dplyr::mutate(patch = patch) %>%
            dplyr::select(-c(count)) %>% 
            tidyr::pivot_wider(names_from = pft, values_from = pct)
        pfts_by_cohort_wide_raw[is.na(pfts_by_cohort_wide_raw)] <- 0 
        # ^ ais what does this line do? 
        
        #ais change the code below so taht I save somewhere the names of all the
        # Rename columns
        if ("pine" %in% colnames(pfts_by_cohort_wide_raw)) {
            pfts_by_cohort_wide_raw <- pfts_by_cohort_wide_raw %>% dplyr::rename(p_p_T=pine)
        } 
        if ("cedar" %in% colnames(pfts_by_cohort_wide_raw)) {
            pfts_by_cohort_wide_raw <- pfts_by_cohort_wide_raw %>% dplyr::rename(p_c_T=cedar)
        } 
        if ("fir" %in% colnames(pfts_by_cohort_wide_raw)) {
            pfts_by_cohort_wide_raw <- pfts_by_cohort_wide_raw %>% dplyr::rename(p_f_T=fir)
        } 
        if ("oak" %in% colnames(pfts_by_cohort_wide_raw)) {
            pfts_by_cohort_wide_raw <- pfts_by_cohort_wide_raw %>% dplyr::rename(p_o_T=oak)
        } 
        if ("shrub" %in% colnames(pfts_by_cohort_wide_raw)) {
            pfts_by_cohort_wide_raw <- pfts_by_cohort_wide_raw %>% dplyr::rename(p_s_T=shrub)
        }
                
        pfts_by_cohort_wide <- pfts_by_cohort_wide_raw %>%
            # Keep only these columns: patch + focal PFTs
            dplyr::select(patch, starts_with("p_")) %>%
            # sum across the p_ columns
            dplyr::mutate(pft_sum = rowSums(pfts_by_cohort_wide %>% dplyr::select(-patch))) %>%
            # Get rid of any patches that have no target PFT classification
            dplyr::filter(pft_sum>0) #ais i want to make sure that bareground/rock is represented 
        # as taking up space in patches
        
        # After removing non-PFT columns, now scale PFT percentage columns 
        # so that they sum to 1 
        for (c in colnames(pfts_by_cohort_wide)) {
            if (grepl("p_", c)) {
                # divide each p_ column value by sum 
                pfts_by_cohort_wide[,c] = pfts_by_cohort_wide[,c]/pfts_by_cohort_wide$pft_sum 
                #ais doublecheck math - number/pft_sum = ?/100
            }
        }
        # Drop pft_sum column, no longer of use
        pfts_by_cohort_wide <- pfts_by_cohort_wide %>% dplyr::select(-pft_sum)
        
        lad_laz_plot_paths <- list.files(file.path(lad_laz_clipped_path), 
                                         pattern="*.laz", full.names = T) %>%
            tools::file_path_sans_ext()
        
        by_patch_df <- divide_patch_into_cohorts(lad_laz_plot_paths, ic_type)
        
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
        
        #cohort_df_initial[which(is.na(cohort_df_initial$pft)),]
        # Cohort file
        cohort_df_initial <- patch_all_df %>%
            dplyr::select(-c(lai_cohort,agb,leaf_biom )) %>%
            dplyr::rename(index = cohort_idx,
                          height= cohort_height,
                          n = n_stemdens) %>%
            dplyr::mutate(pft = case_when(
                pft == "pine" ~ 1,
                pft == "cedar" ~ 2,
                pft == "fir" ~ 3,
                pft == "shrub" ~ 4,
                pft == "oak" ~ 5 ),
                time = year,
                bdead = 0,
                balive = 0) %>%
            dplyr::relocate(time,patch) %>%
            dplyr::relocate(dbh, .after=index) %>%
            dplyr::relocate(height, .after=dbh) %>%
            dplyr::relocate(n, .after=pft) %>%
            dplyr::mutate(avgRG = 0)
        
        # Patch file
        patch_df_initial <- data.frame( #each row is one subplot, or patch
            time = year,
            patch = sort(unique(cohort_df_initial$patch)), 
            # ^ 1:nrow(plots_sf), not every patch from plots_sf is characterized by cohorts
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
            dplyr::mutate(area = 1/length(unique(cohort_df_initial$patch))) %>% # nrow(plots_sf))
            dplyr::relocate(area, .after=age)
        
    } 
    
    # Validate css and pss
    pcss_list <- validate_pss_css(year, cohort_df_initial, patch_df_initial, ic_type, allom_params)
    cohort_df_final <- pcss_list[[1]]
    patch_df_final <- pcss_list[[2]]
    
    # Save final initial conditions --------------------------------------------
    
    # wiki: https://github.com/EDmodel/ED2/wiki/Initial-conditions#files-types-and-
    # formats-for-nlied_init_mode6
    # FATES initial conditions are generated the same way as ED2
    
    if (ic_type == "field_inv_plots") {
        cohort_path <- file.path(ic_type_path,"ic_field_inv.css")
        patch_path <- file.path(ic_type_path,"ic_field_inv.pss")
        
    } else if (ic_type=="rs_inv_plots") { 
        cohort_path <- file.path(ic_type_path,"ic_rs_inv_plots.css")
        patch_path <- file.path(ic_type_path,"ic_rs_inv_plots.pss")
        
    } else if (ic_type == "rs_random_plots") {
        cohort_path <- file.path(ic_type_path,"ic_rs_random_plots.css")
        patch_path <- file.path(ic_type_path,"ic_rs_random_plots.pss")
        
    } else {
        print("need to specify dataset type")
    }
        
    readr::write_delim(x=cohort_df_final,file=cohort_path,delim=" ",append=FALSE,quote="none")
    readr::write_delim(x=patch_df_final,file=patch_path,delim=" ",append=FALSE,quote="none")
    
    return(c(cohort_path, patch_path))
}




validate_pss_css <- function(year, i_css, i_pss, ic_type, allom_params){
        
    initial_list = 
        tibble::tribble( ~desc               , ~in_suffix         , ~out_suffix       , ~year,              ~colour  , ~f_net_def
                         , "Inventory (Field)" , "ic_inv"           , "inventory_neon"  , as.numeric(year), "#CC8829", 1.0
                         , "Inventory (RS)"    , "ic_rs_over_inv"   , "rsinventory_neon", as.numeric(year), "#47B2AA", 1.0
                         , "Entire site (RS)"  , "ic_rs_over_random", "prismatic_neon"  , as.numeric(year), "#AC5CE5", 1.0 # 2.018285
        )
    pft_lookup = allom_params  %>%
        dplyr::mutate(fates_id = case_when(pft=="pine" ~ 1, pft=="cedar" ~ 2,
                                           pft=="fir" ~ 3, pft=="shrub" ~ 4,
                                           pft=="oak" ~ 5, .default = NA),
                      colour = case_when(pft=="pine" ~ "#66CCAE", pft=="cedar" ~ "#1B9E77",
                                         pft=="fir" ~ "#095941", pft=="shrub" ~ "#D95F02",
                                         pft=="oak" ~ "#7570B3", .default = "other"),
                      desc = case_when(pft=="pine" ~ "Pine", pft=="cedar" ~ "Cedar",
                                       pft=="fir" ~ "Fir", pft=="shrub" ~ "Shrub",
                                       pft=="oak" ~ "Oak", .default = "other"))
    
    #   Loop through output types, and fix files.
    i = ifelse(ic_type == "field_inv_plots", 1,
               ifelse(ic_type == "rs_inv_plots", 2,
                      ifelse(ic_type == "rs_random_plots", 3, 0 )))
    
    #   Load information.
    io_desc   = initial_list$desc      [i]
    io_year   = initial_list$year      [i]
    f_net_def = initial_list$f_net_def [i]
    message(" + Process data for ",tolower(io_desc)," initialisation.")    
    
    #   Validate files
    message("   - Validate files.")
    is_pcss_match = all(i_css$patch %in% i_pss$patch)
    if (! is_pcss_match){
        invalid_patch = i_css$patch[! (i_css$patch %in% i_pss$patch)]
        
        message("---~---"                                                                 )
        message(" FATAL ERROR!"                                                           )
        message("---~---"                                                                 )
        message("  The following patches exist in cohort files but not in the patch file:")
        for (p in seq_along(invalid_patch)){
            message(" + ",invalid_patch[p],".")
        }#end for (p in seq_along(invalid_patch))
        message("---~---"                                                                 )
        stop(" All patches in cohort file must appear in the patch file.")
    }#end if (! is_pcss_match)
    
    
    #   Validate and standardise PFTs.
    message("   - Validate PFTs.")
    if (is.character(i_css$pft)){
        is_pft_fine = all(i_css$pft %in% pft_lookup$as_name)
        
        #   Report invalid PFTs.
        if (! is_pft_fine){
            invalid_pft = sort(unique(i_css$pft[! (i_css$pft %in% pft_lookup$as_name)]))
            
            message("---~---"                                                           )
            message(" FATAL ERROR!"                                                     )
            message("---~---"                                                           )
            message("  The following PFTs exist in cohort files but are not recognised:")
            for (p in seq_along(invalid_pft)){
                message(" + ",invalid_pft[p],".")
            }#end for (p in seq_along(invalid_pft))
            message("---~---"                                                           )
            stop(" All PFTs in cohort file must match a value in \"pft_lookup\"."    )
        }
        
        #   Substitute names with the FATES id
        i_css$pft = pft_lookup$fates_id[match(i_css$pft,pft_lookup$as_name)]
        
    }else{
        is_pft_fine = all(i_css$pft %in% pft_lookup$fates_id)
        if (! is_pft_fine){
            invalid_pft = sort(unique(i_css$pft[! (i_css$pft %in% pft_lookup$fates_id)]))
            
            message("---~---"                                                           )
            message(" FATAL ERROR!"                                                     )
            message("---~---"                                                           )
            message("  The following PFTs exist in cohort files but are not recognised:")
            for (p in seq_along(invalid_pft)){
                message(" + ",invalid_pft[p],".")
            }#end for (p in seq_along(invalid_pft))
            message("---~---"                                                           )
            stop(" All PFTs in cohort file must match a value in \"pft_lookup\"."    )
        }#end if (! is_pft_fine)
    }#end if (is.character(icss$pft))
    
    
    #   Validate and standardise PFTs.
    message("   - Validate cohorts.")
    if (any(i_css$n == 0.)){
        #---~---
        #   Find invalid.
        #---~---
        sel_invalid = which(i_css$n == 0.)
        i_invalid   = i_css[sel_invalid,,drop=""]
        #---~---
        
        
        #---~---
        #   Report invalid PFTs.
        #---~---
        message ("---~---"                                                            )
        message (" FATAL ERROR!"                                                      )
        message ("---~---"                                                            )
        message ("   There are cohorts with zero stem density, and this is not valid.")
        message (" The first few lines are printed here. Check variable \"i_invalid\"")
        message (" for the full list."                                                )
        message ("---~---"                                                            )
        print(i_invalid)
        message ("---~---"                                                            )
        stop(" Cohort cannot have zero density."                                   )
    }
    
    #   Make sure that the patch area adds to 1
    message("   - Standardise patch area.")
    i_pss$area = i_pss$area / sum(i_pss$area)
    
    #   Standardise names...
    message("   - Standardise names.")
    if (! "nplant" %in% names(i_css)) i_css = i_css %>% rename( nplant = n     )
    if ("index"    %in% names(i_css)) i_css = i_css %>% rename( cohort = index )
    
    #   For patches and cohorts, we use either the original name or the hexadecimal code.
    # if (is.numeric(i_pss$patch )) i_pss$patch  = sprintf("0x%3.3X",i_pss$patch )
    # if (is.numeric(i_css$patch )) i_css$patch  = sprintf("0x%3.3X",i_css$patch )
    if (is.numeric(i_pss$patch )) i_pss$patch  = as.character(i_pss$patch)
    if (is.numeric(i_css$patch )) i_css$patch  = as.character(i_css$patch)
    if (is.numeric(i_css$cohort)) i_css$cohort = sprintf("0x%3.3X",i_css$cohort)
    
    #   Apply correction factor if needed.
    i_css$nplant = i_css$nplant * f_net_def
    
    
    #   Discard patches with tiny populations.
    no_css = i_css %>% filter(i_css$nplant <  1.e-8)
    i_css  = i_css %>% filter(i_css$nplant >= 1.e-8)
    
    #   Sort patches and cohorts
    i_pss = i_pss %>% arrange(patch)
    i_css = i_css %>% arrange(patch,-dbh)
    
    #   Create output patch and cohort structures.
    message("   - Create output structures.")
    o_pss = tibble::tibble( time  = sprintf("%4.4i"  , io_year    )
                            , patch = sprintf("%s"     , i_pss$patch)
                            , trk   = sprintf("%5i"    , 2L         )
                            , age   = sprintf("%6.1f"  , 0.         )
                            , area  = sprintf("%16.14f", i_pss$area )
                            , water = sprintf("%5i    ", 0L         )
                            , fsc   = sprintf("%10.5f" , 0.         )
                            , stsc  = sprintf("%10.5f" , 0.         )
                            , stsl  = sprintf("%10.5f" , 0.         )
                            , ssc   = sprintf("%10.5f" , 0.         )
                            , psc   = sprintf("%10.5f" , 0.         )
                            , msn   = sprintf("%10.5f" , 0.         )
                            , fsn   = sprintf("%10.5f" , 0.         )
    )#end tibble::tibble
    
    #   Create output cohort structure.
    o_css = tibble::tibble( time   = sprintf("%4.4i"  , io_year     )
                            , patch  = sprintf("%s"     , i_css$patch )
                            , cohort = sprintf("%s"     , i_css$cohort)
                            , dbh    = sprintf("%9.3f"  , i_css$dbh   )
                            , height = sprintf("%9.3f"  , i_css$height)
                            , pft    = sprintf("%5i"    , i_css$pft   )
                            , nplant = sprintf("%16.10f", i_css$nplant)
                            , bdead  = sprintf("%9.3f"  , 0.          )
                            , balive = sprintf("%9.3f"  , 0.          )
                            , avgRg  = sprintf("%9.3f"  , 0.          )
    )#end tibble:tibble    
    
    #   Write the output files
    message("   - Write formatted files.")
    return(list(o_css, o_pss))
}