<img src="docs/img/logo.png" align="right" width="25%"/>

[![Python Package using Conda](https://github.com/RS-PRISMATIC/preprocessing/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/RS-PRISMATIC/preprocessing/actions/workflows/python-package-conda.yml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/RS-PRISMATIC/preprocessing/HEAD)

# PRISMATIC-Initialization
Generating initial conditions for FATES 

# Install
```
conda env create -f environment.yml
```

# Troubleshoots
For MacOS M1/M2 users:
 - It does not work with R arch arm64, so you will need to reinstall R arch x86_64, follow this guide [here](https://github.com/rpy2/rpy2/issues/900#issuecomment-1499431341).
 - `Library not loaded: /opt/X11/lib/libX11.6.dylib`: run this command: `brew install xquartz --cask`
    - If you need to install `brew`, run this command: `$ /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"`

GDAL, OSR:
 - To import GDAL, use: `from osgeo import gdal, osr`

Cannot install leafR:
 - sudo aptitude install libgdal-dev

# Usage

`conf/paths/paths.yaml`: you may need to update the path.
- `data_raw_inv_path`: where to save inventory data.
- `data_raw_aop_path`: where to save raw NEON Airborne Observatory Platform (AOP) data.
- `data_int_path`: where to save processed data.
- `data_final_path`: where to save final-formatted data.

`conf/sites/sites.yaml`: you may need to add other sites if you want to process those.

Data procesing workflow:
[Workflow diagram of functions](https://drive.google.com/file/d/1Ttap0vm3rWWv8yI-nDyoWKjv-z10Mz5l/view?usp=sharing)
- For each site, download NEON data: 
    - `download_lidar`: discrete return lidar
    - `download_hyperspectral`: tiles or flightlines
    - `download_veg_structure_data`: woody veg data
    - `download_polygons`:  polygons for woody veg plots
    - `download_trait_table`: table of trait values for species
- For all sites at once, `generate_pft_reference`: generates table linking species
- For one site at a time, process intermediate data products:
    - `prep_veg_structure`: organize woody veg data into year and sampling effort
    - `prep_polygons`: identify subplot associated with woody veg data
    - `normalize_laz`: normalize laz tiles
    - `clip_lidar_by_plots`: clip the laz/tif files given plots generated from `prep_polygons`
    - `prep_lad`: generate leaf area density profile stratified by size class
    - `prep_biomass`: calculate biomass for all invidivuals in NEON data
    - if using hyperspectral flightlines rather than tiles
        - `correct_flightlines`: apply BRDF and topographic corrections and convert corrected flightlines into tiles
    - `prep_manual_training_data`: clean manually delineated tree crowns with plant functional type (PFT) labels
    - `prep_aop_imagery`: prepare NEON AOP-derived rasters for PFT classifier
    - `extract_spectra_from_polygon`: Extract features (remote sensing data) for each pixel within the specified shapefile (tree crown for training, plot for prediction)
- For all sites at once, `train_pft_classifier`: train PFT classifier
- For a single target site and year, `generate_initial_conditions`: generate FATES initial conditions (cohort and patch files) by combining forest structure and composition

We generate FATES intital conditions in three types: 
- `ic_type == field_inv_plots`: initialization from NEON forest inventory plots
- `ic_type == rs_inv_plots`: initialization from *remote sensing data over* NEON forest inventory plots
- `ic_type == rs_random_plots`: initialization from remote sensing data over *plots randomly generated across entire NEON site*

The final result is at `data_final_path/site/year/ic_type`.

```
# run all sites with default params
python main.py

# force prep_biomass to rerun for all sites/years
python main.py sites.global_run_params.force_rerun.prep_biomass=True

# run only SJER
python main.py sites.global_run_params.run=SJER

# run SJER and SOAP
python main.py sites.global_run_params.run='[SJER, SOAP]'

# run only SJER, force to rerun prep_biomass on SJER
python main.py sites.global_run_params.run=SJER sites.SJER.2019.force_rerun.prep_biomass=True
```
