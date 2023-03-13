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
```
python main.py
```
