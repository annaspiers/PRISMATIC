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

download_veg_structure_data <- function(site) {
    wd <- getwd()
    wd <- glue("{wd}/data/{site}")
    setwd(wd)
    if (!dir.exists(wd)) {
        dir.create(wd)
    } else {
        print("dir exists")
    }
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
    write.csv(df, glue("{wd}/veg_structure.csv"), row.names = FALSE)
    return(glue("{wd}/veg_structure.csv"))
}
