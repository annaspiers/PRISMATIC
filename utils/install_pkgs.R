if (!require("ggbiplot")){
    install.packages("ggbiplot")
    library(ggbiplot)
}

if (!require("neonOS")){
    install.packages("neonOS")
    library(neonOS)
}

if (!require("geoNEON")){
    install.packages("devtools")
    devtools::install_github("NEONScience/NEON-geolocation/geoNEON")
    library(geoNEON)
}