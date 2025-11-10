# Grassland–Agriculture Modelling Code

This repository contains R scripts used to reproduce the statistical analyses for a study examining how surrounding land use and agricultural intensity influence grassland plant diversity at the continental scale.

## Repository structure
/scripts
├── 01_data_preparation.R
├── 02_model_fitting.R
├── 03_model_validation.R
├── 04_predictions_visualisation.R

Each script can be run independently but assumes that input data (vegetation plots, land-cover layers, and derived covariates) are available.

## Requirements

- **R version:** ≥ 4.3  
- **Core packages:**  
  `sf`, `spdep`, `INLA`, `ggplot2`, `data.table`, `raster`, `terra`, `mgcv`, `glmmTMB`  
  (additional dependencies are listed at the start of each script)

## Usage

Scripts are designed to be executed sequentially from 01 to 04.  
Paths to input and output directories should be set at the top of each script under the variable `base_dir`.

Example:
```r
base_dir <- "path/to/local/project/"
source("scripts/01_data_preparation.R")
