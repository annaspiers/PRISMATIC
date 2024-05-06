<img src="docs/img/logo.png" align="right" width="25%"/>

[![Python Package using Conda](https://github.com/RS-PRISMATIC/preprocessing/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/RS-PRISMATIC/preprocessing/actions/workflows/python-package-conda.yml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/RS-PRISMATIC/preprocessing/HEAD)

# PRISMATIC-Preprocessing
Data acquisition and preprocessing for PRISMATIC

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
- `data_raw_inv_path`: where to save processed data.
- `data_raw_aop_path`: where to save raw NEON Airborne Observatory Platform (AOP) data.
- `data_int_path`: where to save processed data.
- `data_final_path`: where to save final-formatted data.

`conf/sites/sites.yaml`: you may need to add other sites if you want to process those.

List of processes for each site:
[Workflow diagram of functions](https://drive.google.com/file/d/1Ttap0vm3rWWv8yI-nDyoWKjv-z10Mz5l/view?usp=sharing)
- `download_lidar`: download lidar data from NEON
- `download_veg_structure_data`: download vegetation data from NEON
- `preprocess_veg_structure_data`: process vegetation data and sampling effort
- `download_polygons`: download polygons data from NEON
- `preprocess_polygons`: process polygons data
- `normalize_laz`: normalize laz files
- `clip_lidar_by_plots`: clip the laz/tif files given plots in processed vegetation structure and save to output
- `preprocess_biomass`: process biomass and save to output
- `preprocess_lad`: process Leaf Area Density and save to output
- `download_hyperspectral`: download imaging spectroscopy data from NEON
- `prep_aop_imagery`: prepare NEON AOP imagery for plant functional type (PFT) classifier
- `create_training_data`: generate training data for PFT classifier
- `train_pft_classifier`: train PFT classifier
- `generate_initial_conditions`: generate FATES initial conditions (cohort and patch files)

We generate FATES intital conditions in three types: 
- `ic_type == field_inv_plots`: initialization from NEON forest inventory plots
- `ic_type == rs_inv_plots`: initialization from *remote sensing data over* NEON forest inventory plots
- `ic_type == rs_random_plots`: initialization from remote sensing data over *plots randomly generated across entire NEON site*

The final result is at `data_final_path/site/year/ic_type`.

```
# run all sites with default params
python main.py

# force preprocess_biomass to rerun for all sites/years
python main.py sites.global_run_params.force_rerun.preprocess_biomass=True

# run only SJER
python main.py sites.global_run_params.run=SJER

# run SJER and SOAP
python main.py sites.global_run_params.run='[SJER, SOAP]'

# run only SJER, force to rerun preprocess_biomass on SJER
python main.py sites.global_run_params.run=SJER sites.SJER.2019.force_rerun.preprocess_biomass=True
```
