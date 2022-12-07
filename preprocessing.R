library(sp)
library(neonUtilities)
library(neonOS)
library(geoNEON)
library(glue)
options(stringsAsFactors=F)

# set working directory
# adapt directory path for your system

download_veg_structure_data <- function(site) {
    wd <- glue("./data/{site}")
    setwd(wd)
    if (!dir.exists(wd)){
        dir.create(wd)
    }
    else{
        print("dir exists")
    }
    veglist <- loadByProduct(dpID="DP1.10098.001",
                            site="SOAP",
                            package="basic",
                            check.size = FALSE)
    vegmap <- getLocTOS(veglist$vst_mappingandtagging,
                            "vst_mappingandtagging")
    veg <- joinTableNEON(veglist$vst_apparentindividual,
                        vegmap,
                        name1="vst_apparentindividual",
                        name2="vst_mappingandtagging")
    df <- data.frame(veg)
    write.csv(df, wd, row.names=False)
    return(glue("raw data for {site} downloaded"))
}

