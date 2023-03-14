install.packages("neonUtilities")
install.packages("neonOS")
devtools::install_github("NEONScience/NEON-geolocation/geoNEON")
library(sp)
library(neonUtilities)
library(neonOS)
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
    veglist <- loadByProduct(dpID = "DP1.10098.001",
                            site = site,
                            package = "basic",
                            check.size = FALSE)
    vegmap <- getLocTOS(veglist$vst_mappingandtagging,
                            "vst_mappingandtagging")
    veg <- joinTableNEON(veglist$vst_apparentindividual,
                         vegmap,
                         name1 = "vst_apparentindividual",
                         name2 = "vst_mappingandtagging")
    df <- data.frame(veg)
    df_perplotperyear <- data.frame(veglist$vst_perplotperyear)
    write.csv(df, glue("{wd}/veg_structure.csv"), row.names = FALSE)
    write.csv(df_perplotperyear, glue("{wd}/plot_sampling_effort.csv"),
              row.names = FALSE)
    return(glue("raw data for {site} downloaded"))
}
