if (!require(ggbiplot)){
    install.packages("ggbiplot")
}

if (!require(neonOS)){
    install.packages("neonOS")
}

if (!require(geoNEON)){
    install.packages("devtools")
    devtools::install_github("NEONScience/NEON-geolocation/geoNEON")
}