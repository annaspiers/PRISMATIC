# install.packages("neonUtilities")
# install.packages("neonOS")
# install.packages("devtools")
# install.packages("glue")
library(neonUtilities)
library(neonOS)
#devtools::install_github("NEONScience/NEON-geolocation/geoNEON", dependencies=TRUE)
library(geoNEON)
library(glue)
options(stringsAsFactors = FALSE)

# set working directory
# adapt directory path for your system

download_veg_structure_data <- function(site, data_path) {
    wd <- glue("{data_path}/{site}")
    if (!dir.exists(wd)) {
        dir.create(wd)
    } else {
        print("dir exists")
    }
    setwd(wd)
    veglist <- neonUtilities::loadByProduct(dpID = "DP1.10098.001",
                            site = site,
                            package = "basic",
                            check.size = FALSE)
    vegmap <- geoNEON::getLocTOS(veglist$vst_mappingandtagging,
                            "vst_mappingandtagging")
    veg <- neonOS::joinTableNEON(veglist$vst_apparentindividual,
                         vegmap,
                         name1 = "vst_apparentindividual",
                         name2 = "vst_mappingandtagging")
    df <- data.frame(veg)
    df_perplotperyear <- data.frame(veglist$vst_perplotperyear)
    write.csv(df, glue("{wd}/veg_structure.csv"), row.names = FALSE)
    write.csv(df_perplotperyear, glue("{wd}/plot_sampling_effort.csv"),
              row.names = FALSE)
    return("wd")
}



# Function to match species to PFT
match_species_to_pft <- function(species = NULL, indvdID = NULL, taxonID = NULL, 
                                 growthForm = NULL) {
    if (!is.null(species)) species <- tolower(species)
    if (!is.null(indvdID)) indvdID <- tolower(indvdID)
    if (!is.null(taxonID)) taxonID <- tolower(taxonID)
    if (!is.null(growthForm)) growthForm <- tolower(growthForm)
    
    # Initialize pft as a character vector
    pft <- character(length(species))  # Preallocate a character vector for pft
    
    # Define a named list of patterns for each PFT
    pft_patterns <- list(
        oak = c("oak", "quercus", "quch2", "quke", "querc", "quwi2", "qudo",
                "alder", "aspen", "cottonwood", "populus", "eucalyptus"), # other broadleaf
        pine = c("pine", "pinus", "pico", "pije", "pinace",
                 "pisa2", "pipo", "pila", "pisa2"),
        fir = c("fir", "abies", "abco", "ablo", "abma"), 
        cedar = c("cedar", "cypress", "sequoia", 
                  "cade27", "segi2", "juniper"), 
        shrub = c("ceanothus", "wax", "cherry", "berry", "willow", "misery", 
                  "chinkapin", "chinquapin", "ribes", "alder", "manzan", 
                  "lupin", "buck", "leuco", "coffee", "rhamnus", 
                  "oleander", "lonicera", "kekiella", "olive", 
                  "aeca", "aece", "arcto3spp", "arne", "arpa", "arpa6", "arvim", 
                  "ceano", "ceca", "ceco", "cecu", "cein3", "cele2", "cemog", "chfo", 
                  "chse11", "conu4", "crse11", "dawr2", "frca6", "frcac7", "loin4", "lual4", 
                  "lupinspp", "prem", "rhil", "rhamna", "ribes", "rice", "riro", "rivi3", 
                  "rhila", "salix", "sani4", "sasc", "sefl3", "todi"), 
        grass_rock_bare = c("rock", "bareground", "grass")
    )    
    
    # Check for growthForm first
    for (i in seq_along(species)) {
        if (!is.null(growthForm) && growthForm[i] %in% c("small shrub", "single shrub")) {
            pft[i] <- "shrub"
            next  # Skip to the next iteration if growthForm matches
        }
        
        # Check for matches in the provided columns
        for (pft_name in names(pft_patterns)) {
            if ((any(grepl(paste(pft_patterns[[pft_name]], collapse = "|"), species[i], ignore.case = TRUE))) ||
                (any(grepl(paste(pft_patterns[[pft_name]], collapse = "|"), indvdID[i], ignore.case = TRUE))) ||
                (any(grepl(paste(pft_patterns[[pft_name]], collapse = "|"), taxonID[i], ignore.case = TRUE)))) {
                pft[i] <- pft_name  # Assign the matched PFT
                break  # Exit the loop once a match is found
            } else {
                pft[i] <- NA  # Default case if no matches are found
            }
        }
    }
    
    return(pft)  # Return the character vector of PFTs
}
