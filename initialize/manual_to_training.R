# Clean tree crown polygons that were manually delineated and labelled in imagery with inventory data
# Save to reference shapes to respective site/year folder in raw/


# SJER 
    # 2017 NEON inventory

    # 2021 NEON inventory
    # 2023 carissa/anna fieldwork
    # 2024 carissa/anna fieldwork

# SOAP 
    # 2019 NEON inventory

    # 2021 NEON inventory


# TEAK 
    # 2021 NEON inventory
    # 2023 carissa/anna fieldwork
    # 2024 carissa/anna fieldwork

library(dplyr)
library(sf)
library(ggplot2)
library(forcats)

source("initialize/inventory_helper.R")

# Save NEON data tot heir respective folders -------------------------------------------------------------------

# SJER 2017 NEON inventory

# SOAP 2019 NEON inventory

# SOAP 2021 NEON inventory


# Combine 2021-2024 inventory data into 2021 - SJER 

# SJER 2021 NEON inventory
SJER_2021_neon_manual   <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/NEON - SJER - 2021/my_shapes.shp")
st_crs(SJER_2021_neon_manual)[[1]]=="WGS 84 / UTM zone 11N"
SJER_2021_neon_ref      <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/NEON - SJER - 2021/tree_crowns_training.shp")
st_crs(SJER_2021_neon_ref)[[1]]=="WGS 84 / UTM zone 11N"

# SJER 2023 carissa/anna fieldwork
SJER_2023_c_manual      <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/Carissa - SJER - 2023/my_shapes.shp") #manually exported from QGIS
st_crs(SJER_2023_c_manual)[[1]]=="WGS 84 / UTM zone 11N"
SJER_2023_c_ref         <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/Carissa - SJER - 2023/tgis_merged_species_32611.shp")
st_crs(SJER_2023_c_ref)[[1]]=="WGS 84 / UTM zone 11N"

# SJER 2024 carissa/anna fieldwork
SJER_2024_c_a_manual    <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/Carissa Anna - SJER - 2024/my_shapes.shp") #manually exported from QGIS
st_crs(SJER_2024_c_a_manual)[[1]]=="WGS 84 / UTM zone 11N"
SJER_2024_c_a_ref1      <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/Carissa Anna - SJER - 2024/SJER_trimble_carissa/trimble_merged_PLY_SJER.shp") %>%
    st_transform(32611)
st_crs(SJER_2024_c_a_ref1)[[1]]=="EPSG:32611"
SJER_2024_a_ref2        <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/Carissa Anna - SJER - 2024/Anna_V1_ais-point.shp")
st_crs(SJER_2024_a_ref2)[[1]]=="WGS 84 / UTM zone 11N"
SJER_2024_c_a_ref3      <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/Carissa Anna - SJER - 2024/SJER_Garmin_carissa/garmin_merged_PT_SJER.shp") %>%
    st_transform(32611)
st_crs(SJER_2024_c_a_ref3)[[1]]=="EPSG:32611"



# plot manual vs ref data
plot(SJER_2021_neon_manual$geometry,col="blue")
plot(SJER_2021_neon_ref$geometry,col="red", add=T)

plot(SJER_2023_c_manual$geometry,col="blue")
plot(SJER_2023_c_ref$geometry,col="red", add=T)

plot(SJER_2024_c_a_manual$geometry,col="blue")
plot(SJER_2024_c_a_ref1$geometry,col="red", add=T)
plot(SJER_2024_a_ref2$geometry,col="yellow", add=T)
plot(SJER_2024_c_a_ref3$geometry,col="green", add=T)

# Merge manual and reference data for each year 
# Here we get the right PFT for each shape in each year
SJER_2021_neon <- SJER_2021_neon_manual %>% select(-id) %>% 
    left_join(data.frame(SJER_2021_neon_ref) %>% select(indvdID,species=sp_gen))

SJER_2023_c <- SJER_2023_c_manual %>% 
    left_join(data.frame(SJER_2023_c_ref) %>% select(Name,species)) %>%
    rename(indvdID=Name)

SJER_2024_c_a <- SJER_2024_c_a_manual %>%  
    left_join(data.frame(SJER_2024_c_a_ref1 %>% 
                             dplyr::rename(site_type=`Site type`,
                                           species_other=`If other s`)) %>% 
                  select(Name,Species,species_other,site_type)) %>% 
    mutate(species=ifelse(!is.na(species_ai),species_ai,
                        ifelse(site_type=="Grassland",site_type,
                        ifelse(Species=="OTHER",species_other, 
                               ifelse(!is.na(Species),Species,site_type))))) %>%
    mutate(species=ifelse(is.na(species),"grass",species)) %>%
    tibble::rownames_to_column() %>%
    # assign arbitrary names to crowns delineated by me 
    # (not assiciated with carissa's dataset)
    mutate(Name = ifelse(is.na(Name), paste0(rowname,"_SJER24_ais"), Name)) %>%
    select(indvdID=Name, species)
    
# Merge shapefiles for all years since we're using only 2021 imagery for now

SJER_2021 <- SJER_2021_neon %>%
    rbind(SJER_2023_c) %>%
    rbind(SJER_2024_c_a) %>%
    mutate(pft = match_species_to_pft(species, indvdID)) %>%
    # Get rid of any rows with no PFT. No longer using the 'other_herb' PFT
    filter(indvdID!="A08") %>%
    select(c(indvdID,pft))

write_sf(SJER_2021,"/Users/AISpiers/Downloads/crown_delin_datasets/SJER_2021_merged.shp")

# Combine 2021-2024 inventory data into 2021 - TEAK --------------------------------------------------------------------

# TEAK 2021 NEON inventory
TEAK_2021_neon_manual   <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/NEON - TEAK - 2021/my_shapes.shp")
st_crs(TEAK_2021_neon_manual)[[1]]=="WGS 84 / UTM zone 11N"
TEAK_2021_neon_ref      <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/NEON - TEAK - 2021/tree_crowns_training.shp")
st_crs(TEAK_2021_neon_ref)[[1]]=="WGS 84 / UTM zone 11N"

# TEAK 2023 carissa/anna fieldwork
TEAK_2023_c_a_manual      <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/Carissa - TEAK - 2023/my_shapes.shp") #manually exported from QGIS
st_crs(TEAK_2023_c_a_manual)[[1]]=="WGS 84 / UTM zone 11N"
TEAK_2023_c_a_ref         <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/Carissa - TEAK - 2023/tgis_merged_species.shp") %>%
    st_transform(32611)
st_crs(TEAK_2023_c_a_ref)[[1]]=="EPSG:32611"

# TEAK 2024 carissa/anna fieldwork
TEAK_2024_c_a_manual    <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/Carissa - TEAK - 2024/my_shapes.shp") #manually exported from QGIS
st_crs(TEAK_2024_c_a_manual)[[1]]=="WGS 84 / UTM zone 11N"
TEAK_2024_c_a_ref       <- read_sf("/Users/AISpiers/Downloads/crown_delin_datasets/Carissa - TEAK - 2024/TEAK_tgis_carissa/tgis_merged_PLY_TEAK.shp") %>%
    st_transform(32611)
st_crs(TEAK_2024_c_a_ref)[[1]]=="EPSG:32611"


# plot manual vs ref data
plot(TEAK_2021_neon_manual$geometry,col="blue")
plot(TEAK_2021_neon_ref$geometry,col="red", add=T)

plot(TEAK_2023_c_a_manual$geometry,col="blue")
plot(TEAK_2023_c_a_ref$geometry,col="red", add=T)

plot(TEAK_2024_c_a_manual$geometry,col="blue")
plot(TEAK_2024_c_a_ref$geometry,col="red", add=T)


# Merge manual and reference data for each year 
# Here we get the right PFT for each shape in each year
TEAK_2021_neon <- TEAK_2021_neon_manual %>% select(-c(id,notes)) %>% #03011, 02218
    left_join(data.frame(TEAK_2021_neon_ref) %>% select(indvdID,species=sp_gen))

TEAK_2023_c_a <- TEAK_2023_c_a_manual %>% 
    left_join(data.frame(TEAK_2023_c_a_ref) %>% select(Name,species)) %>%
    mutate(species = case_when(
        grepl("C03",Name,ignore.case=T) ~ "oak", 
        grepl("C036",Name,ignore.case=T) ~ "willow", 
        grepl("A172",Name,ignore.case=T) ~ "pine", 
        grepl("A101",Name,ignore.case=T) ~ "manzanita",
        .default = species )) %>%
    filter(!is.na(species)) %>%
    #Remove where there is lotus in the species name - too few for training data
    filter(!grepl("lotus", species, ignore.case = TRUE)) %>%
    rename(indvdID=Name) 
    

TEAK_2024_c_a <- TEAK_2024_c_a_manual %>% 
    rename(species=species_ai) %>%
    # assign arbitrary names to crowns delineated by me 
    # (not assiciated with carissa's dataset)
    tibble::rownames_to_column() %>%
    mutate(indvdID = ifelse(is.na(id), paste0(rowname,"_TEAK24_ais"), id)) %>%
    filter(!is.na(species)) %>%
    select(indvdID, species)


# Merge shapefiles for all years since we're using only 2021 imagery for now

TEAK_2021 <- TEAK_2021_neon %>%
    rbind(TEAK_2023_c_a) %>%
    rbind(TEAK_2024_c_a) %>%
    mutate(pft = match_species_to_pft(species, indvdID)) %>%
    select(c(indvdID,pft))

write_sf(TEAK_2021,"/Users/AISpiers/Downloads/crown_delin_datasets/TEAK_2021_merged.shp")


















# # # Identify which plots are in burn scar -----------------------------------
# # 
# # # Read the shapefiles
# # plots_shp <- st_read("/Users/AISpiers/Downloads/rs_random_plots/random_plots.shp")
# # creek_fire_shp <- st_read("/Users/AISpiers/Downloads/marcos_qgis/CreekFire/2020_CreekFire_Perimeter_EPSG4326.shp") %>%
# #     st_transform(crs=32611)
# # 
# # plot(creek_fire_shp$geometry)
# # plot(plots_shp$geometry,fill="red",add=T)
# # # Check for overlaps
# # overlaps <- st_intersects(plots_shp, creek_fire_shp, sparse = FALSE)
# # 
# # # Extract the overlapping shapes
# # # Get indices of overlapping geometries in plots_shp
# # overlapping_indices <- which(rowSums(overlaps) > 0)
# # 
# # # Extract the overlapping shapes
# # overlapping_plots <- plots_shp[overlapping_indices, ]
# # 
# # plot(creek_fire_shp$geometry)
# # plot(overlapping_plots$geometry,fill="red",add=T)
# # 
# # # Optionally, write the overlapping shapes to a new shapefile
# # st_write(overlapping_plots, "/Users/AISpiers/Downloads/rs_random_plots/Inputs/plots_overlapping_fire.shp")





# # Identify which PFTs overlap with the tower footprints -------------------

# # Towers
# twr_sjer = sf::st_read("/Users/AISpiers/Downloads/marcos_qgis/AmerifluxSites/SJER_90pc_flux_ml.shp")
# twr_soap = sf::st_read("/Users/AISpiers/Downloads/marcos_qgis/AmerifluxSites/SOAP_90pc_flux_ml.shp")
# twr_teak = sf::st_read("/Users/AISpiers/Downloads/marcos_qgis/AmerifluxSites/TEAK_90pc_flux_ml.shp")

# # Veg data
# sjer_neon_veg = readr::read_csv("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/SJER/veg_structure.csv")
# sjer_neon_plots = readr::read_csv("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/SJER/plot_sampling_effort.csv") %>%
#     st_as_sf(coords=c("easting","northing")) %>%
#     select(plotID) %>% distinct() %>%
#     left_join(sjer_neon_veg %>% select(plotID, individualID, taxonID, scientificName) %>% distinct()) %>%
#     mutate(pft=match_species_to_pft(taxonID)) %>% mutate(site="SJER") %>%
#     left_join(SJER_2021 %>% as_tibble() %>% select(individualID=indvdID, pft))
# st_crs(sjer_neon_plots) = 32611
# sjer_intersect = st_intersection(sjer_neon_plots,twr_sjer)

# soap_neon_veg = readr::read_csv("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/SOAP/veg_structure.csv")
# soap_neon_plots = readr::read_csv("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/SOAP/plot_sampling_effort.csv")  %>%
#     st_as_sf(coords=c("easting","northing")) %>%
#     select(plotID) %>% distinct %>% mutate(site="SOAP") %>%
#     left_join(soap_neon_veg %>%  select(plotID, individualID, taxonID, scientificName) %>% distinct()) %>%
#     mutate(pft=match_species_to_pft(taxonID))
# st_crs(soap_neon_plots) = 32611
# soap_intersect = st_intersection(soap_neon_plots,twr_soap)

# teak_neon_veg = readr::read_csv("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/TEAK/veg_structure.csv")
# teak_neon_plots = readr::read_csv("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/TEAK/plot_sampling_effort.csv") %>%
#     st_as_sf(coords=c("easting","northing")) %>%
#     select(plotID) %>% distinct() %>% mutate(site="TEAK") %>%
#     left_join(teak_neon_veg %>% select(plotID, individualID, taxonID, scientificName) %>% distinct()) %>%
#     mutate(pft=match_species_to_pft(taxonID)) %>% 
#     left_join(TEAK_2021 %>% as_tibble() %>% select(individualID=indvdID, pft))
# st_crs(teak_neon_plots) = 32611
# teak_intersect = st_intersection(teak_neon_plots,twr_teak)


# # Identify overlap spatially
# plot(twr_sjer$geometry)
# plot(sjer_neon_plots$geometry, add=T)
# plot(sjer_intersect$geometry, col="red", add=T)

# plot(twr_soap$geometry)
# plot(soap_neon_plots$geometry, add=T)
# plot(soap_intersect$geometry, col="red", add=T)

# plot(twr_teak$geometry)
# plot(teak_21_24$geometry, add=T)
# plot(teak_intersect$geometry, col="red", add=T)

# # Plot species/PFTs
# rbind(sjer_intersect, soap_intersect, teak_intersect) %>%
#     ggplot( aes(x = forcats::fct_infreq(scientificName), fill=pft)) +  # Reorder strings by Count in descending order
#     geom_bar(stat = "count") +
#     labs(title = "Count of PFTs in tower footprint",
#          x="species", y = "Count") +
#     theme_minimal() + facet_wrap(vars(site)) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) 





# #2) read in inventory csv
# live_trees_path <- file.path("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/intermediate/SJER/2021/biomass/pp_veg_structure_IND_IBA_IAGB_live.csv")
# pft_reference <- read.csv("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/intermediate/pft_reference.csv")

# # Load cleaned inventory data for site and year
# # ais this is live trees only. should figure out how to incorporate dead trees
# veg_ind <- read.csv(live_trees_path) %>%
#     dplyr::rename(easting = adjEasting, 
#                   northing = adjNorthing,
#                   # Rename columns so they're unique when truncated when saving shp 
#                   lat = adjDecimalLatitude,
#                   long = adjDecimalLongitude,
#                   elev = adjElevation,
#                   elevUncert = adjElevationUncertainty,
#                   indStemDens = individualStemNumberDensity,
#                   indBA = individualBasalArea,
#                   dendroH = dendrometerHeight,
#                   dendroGap = dendrometerGap,
#                   dendroCond = dendrometerCondition,
#                   bStemDiam = basalStemDiameter,
#                   bStemMeasHgt = basalStemDiameterMsrmntHeight,
#                   sciNameFull = scientificName,
#                   sp_gen = scientific,
#                   wood_dens1 = wood.dens,
#                   wood_dens2 = wood_dens) #ais what is the source of columns in this df? NEON or our own changes?

# # Remove all rows with missing an easting or northing coordinate
# veg_has_coords <- veg_ind[complete.cases(veg_ind[, c("northing", "easting")]),]

# # Keep only the entries with crown diameter and tree height
# # since these measurements are needed to create and compare polygons.
# # ais future work: include ninetyCrownDiameter?
# veg_has_coords_size <- veg_has_coords[complete.cases(veg_has_coords$height) & 
#                                           complete.cases(veg_has_coords$maxCrownDiameter),] %>%
#     dplyr::filter(maxCrownDiameter < 50) %>%
#     dplyr::filter(individualID != "NEON.PLA.D17.SOAP.04659") %>% #no tree here
#     dplyr::filter(individualID != "NEON.PLA.D17.SOAP.05847") %>% #no tree here
#     dplyr::filter(individualID != "NEON.PLA.D17.SOAP.04256")  #no tree here
# #"NEON.PLA.D17.SOAP.04972" #about 1/3 is ground

# # Assign a PFT to each individual
# # ais need to make this an exhaustive match for all sites, not just SOAP
# # ais use neon trait table from Marcos?
# veg_training <- veg_has_coords_size %>%
#     left_join(pft_reference %>% dplyr::select(-growthForm) , by=join_by(siteID, taxonID)) %>%
#     dplyr::mutate(pft = case_when(
#         grepl("CADE27",taxonID,ignore.case=T) ~ "cedar",
#         grepl("ABCO",taxonID,ignore.case=T) |
#             grepl("ABMA",taxonID,ignore.case=T) ~ "fir",
#         grepl("PIPO",taxonID,ignore.case=T) |
#             grepl("PILA",taxonID,ignore.case=T) |
#             grepl("PISA2",taxonID,ignore.case=T) |
#             grepl("PINUS",taxonID,ignore.case=T) ~ "pine",
#         grepl("QUCH2",taxonID,ignore.case=T)  |
#             grepl("QUDO",taxonID,ignore.case=T) |
#             grepl("QUKE",taxonID,ignore.case=T) |
#             grepl("QUERC",taxonID,ignore.case=T) |
#             grepl("QUWI2",taxonID,ignore.case=T) ~ "oak",
#         grepl("SANI4",taxonID,ignore.case=T) |
#             grepl("AECA",taxonID,ignore.case=T) |
#             grepl("FRCA6",taxonID,ignore.case=T) |
#             grepl("FRCAC7",taxonID,ignore.case=T) |
#             grepl("RIRO",taxonID,ignore.case=T) |
#             grepl("RIRO",taxonID,ignore.case=T) |
#             grepl("LUAL4",taxonID,ignore.case=T) |
#             grepl("CEMOG",taxonID,ignore.case=T) |
#             grepl("CELE2",taxonID,ignore.case=T) |
#             grepl("CECU",taxonID,ignore.case=T) |
#             grepl("CEIN3",taxonID,ignore.case=T) |
#             grepl("ARVIM",taxonID,ignore.case=T) |
#             grepl("RHIL",taxonID,ignore.case=T) |
#             grepl("CEANO",taxonID,ignore.case=T) ~ "shrub",
#         .default="other" ))

# #3) link inventory csv and manual shapes
# SJER_2021_neon_manual %>% 
#     st_transform(32611) %>%
#     #  ais ^ will need to change this from being hardcoded to defining globally
#     left_join(veg_training %>% rename(indvdID=individualID)) %>%
#     dplyr::select(c(indvdID,pft,sciNameFull))

