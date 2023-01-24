# preprocessing
Data acquisition and preprocessing for RS-PRISMATIC

# Installation
## Prerequisite
```
apt update
apt upgrade
apt install libcurl4-openssl-dev
apt install libbz2-dev liblzma-dev
apt install libfontconfig1-dev
apt install libharfbuzz-dev libfribidi-dev
apt install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```
## Install R packages
```
install.packages("neonUtilities")
install.packages("neonOS")
install.packages("sp")
install.packages("raster")
install.packages("devtools")
devtools::install_github("NEONScience/NEON-geolocation/geoNEON")
```
## Install conda environment
```
conda create -f environment.yml
```