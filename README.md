<img src="docs/img/logo.png" align="right" width="25%"/>

[![Python Package using Conda](https://github.com/RS-PRISMATIC/preprocessing/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/RS-PRISMATIC/preprocessing/actions/workflows/python-package-conda.yml)

# PRISMATIC-Preprocessing
Data acquisition and preprocessing for RS-PRISMATIC

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
 - [ ] Seperate plots into 400m2 subplots
 - [x] Setup GitHub actions
 - [ ] Add unit tests and integration tests
 - [ ] Add logging
 - [ ] Refactor configuration files

# Usage
```
python run.py
```