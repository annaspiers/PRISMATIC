<img src="docs/img/logo.png" align="right" width="25%"/>

[![Python Package using Conda](https://github.com/RS-PRISMATIC/preprocessing/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/RS-PRISMATIC/preprocessing/actions/workflows/python-package-conda.yml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/RS-PRISMATIC/preprocessing/HEAD)

# PRISMATIC-Preprocessing
Data acquisition and preprocessing for PRISMATIC

# Install
```
conda create -f environment.yml
```

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
