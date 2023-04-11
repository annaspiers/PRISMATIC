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
# Todo
 - [x] Data acquisition and data processing of inventory data
 - [x] Data acquisition and data processing of NEON plots
 - [x] Data acquisition and data processing of LiDAR data
 - [x] Encode allometric equations
 - [x] Calculate individual/plot-level biomass given stem diameter/basal stem diameter
 - [x] Seperate plots into 400m2 subplots
 - [x] Setup GitHub actions
 - [ ] Add unit tests and integration tests
 - [x] Add diagnostics
 - [x] Add logging
 - [x] Refactor configuration files

# Usage

`conf/paths/paths.yaml`: you may need to update the path.
- `data_path`: where to save processed data.
- `lidar_path`: where to download/load lidar files to.

`conf/sites/sites.yaml`: you may need to add other sites if you want to process those.

List of processes for each site:
- `download_lidar`: download lidar data from NEON
- `download_veg_structure_data`: download vegetation data from NEON
- `preprocess_veg_structure_data`: process vegetation data and sampling effort
- `download_polygons`: download polygons data from NEON
- `preprocess_polygons`: process polygons data
- `normalize_laz`: normalize laz files
- `clip_lidar_by_plots`: clip the laz/tif files given plots in processed vegetation structure and save to output
- `preprocess_biomass`: process biomass and save to output

The final result is at `data_path/site/year/output`.

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
