if (!require("ggbiplot")){
    install.packages("ggbiplot")
    library(ggbiplot)
}

if (!require("neonOS")){
    install.packages("neonOS", repos = c("https://cloud.r-project.org"))
    library(neonOS)
}

if (!require("geoNEON")){
    install.packages("devtools")
    devtools::install_github("NEONScience/NEON-geolocation/geoNEON")
    library(geoNEON)
}