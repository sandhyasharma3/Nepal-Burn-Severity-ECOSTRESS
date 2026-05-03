# Nepal Burn Severity — ECOSTRESS & Spatial Predictors

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.1-blue.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19990957.svg)](https://doi.org/10.5281/zenodo.19990957)


## Overview

This repository contains the analysis code for studying **wildfire burn severity** across Nepal using NASA [ECOSTRESS](https://ecostress.jpl.nasa.gov/) thermal data (ESI, ET, WUE) combined with topographic and spatial predictors. Two Random Forest models are compared:

- **Annual model** — ECOSTRESS annual and spring-summer composites
- **AM/PM model** — ECOSTRESS instantaneous acquisitions split by ISS overpass time (AM vs PM)

Both models predict dNBR (Sentinel-2 derived) and include stratified analyses by land cover class, elevation, slope, and aspect across Nepal's 2021 fire season.

## Analysis Workflow

The analysis spans three platforms in sequence:

Step 1 — ArcGIS Pro  (preprocessing)
  └── Fire perimeter clipping, ECOSTRESS raster reprojection & alignment,
      predictor layer preparation, physiographic zone extraction

Step 2 — Google Earth Engine  (burn severity)
  └── Sentinel-2 cloud masking (QA60 + MSK_CLASSI + SCL),
      pre/post NBR composites, dNBR calculation,
      burn severity classification & export to Drive
      → see gee/Nepal_dNBR.js

Step 3 — R  (statistical modeling & figures)
  └── Random Forest modeling, stratified analyses, distribution plots
      → see scripts/


## Repository Structure

```
Nepal-Burn-Severity-ECOSTRESS/
├── gee/
│   ├── Nepal_dNBR.js              # Sentinel-2 dNBR + severity classes (GEE)
│   └── README_gee.md              # Instructions for running in GEE
├── scripts/
│   ├── PCNM_spatial_vectors.R     # Spatial eigenvectors (run on HPC before RF scripts)
│   ├── Annual_RF.R                # RF model: annual ECOSTRESS composites
│   ├── AmPm_RF.R                  # RF model: AM vs PM ECOSTRESS acquisitions
│   └── Severity_Distribution.R    # Halfeye distribution plots by severity class
├── LICENSE
└── README.md
```

> **Data not included.** Raw rasters, sampling point files, and PCNM outputs are not
> tracked due to file size. See [Data Access](#data-access) below.

## Data Access

All raster datasets except Sentinel-2 were downloaded via **NASA Earthdata AppEEARS**
([appeears.earthdatacloud.nasa.gov](https://appeears.earthdatacloud.nasa.gov/)).
An Earthdata login (free) is required. Sentinel-2 imagery was accessed directly
through Google Earth Engine.

| Dataset | AppEEARS Product Name | Variables used |

| ECOSTRESS | ECO_L3T_JET_002 | ESI, ET, WUE (annual, spring-summer, AM, PM) |
| SRTM DEM (30 m) | SRTMGL1_NC.003 | Elevation, Slope, Aspect |
| MODIS Land Cover | MCD12Q1.061 | Land cover classes |
| Sentinel-2 SR | via [Google Earth Engine](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR) (not AppEEARS) | B8A, B12 → dNBR |
| Nepal boundary | User GEE asset | ROI for GEE and spatial sampling |


## R Requirements

All R analyses were run in **R ≥ 4.1**. Install required packages with:

```r
install.packages(c(
  "sf", "terra", "dplyr", "tidyverse",
  "vegan", "adespatial",          # for PCNM spatial eigenvectors
  "randomForest",
  "ggplot2", "ggdist", "ggpointdensity",
  "scales", "patchwork", "cowplot"
))
```

### Running the R scripts

> **Important:** `PCNM_spatial_vectors.R` is computationally intensive (distance matrix
> on 40,000 points) and was run on an HPC cluster (SLURM). Update the `raster_dir` and
> `boundary_path` at the top of the script for your system before running.
> The two RF scripts and the figure script can be run locally once PCNM outputs exist.

Set the `rdir` path at the top of each script to point to your local raster directory.
Scripts should be run in this order:

```r
# Step 1 — Generate PCNM spatial eigenvectors (HPC recommended)
source("scripts/PCNM_spatial_vectors.R")

# Step 2 — Random Forest models (run either or both, independently)
source("scripts/Annual_RF.R")
source("scripts/AmPm_RF.R")

# Step 3 — Severity distribution figures
source("scripts/Severity_Distribution.R")
```

---

## GEE Instructions

See [`gee/README_gee.md`](gee/README_gee.md) for step-by-step instructions on running
the burn severity script in Google Earth Engine.

---

## Outputs

| Script | Key outputs |
|---|---|
| `Annual_RF.R` | `RF_final_metrics.csv`, `RF_final_varimp.csv`, `Figure3.tiff`, `Table1_LC_stratified.csv`, `Table2_Topo_stratified.csv`, `TableS1_LCxTopo.csv` |
| `AmPm_RF.R` | `RF_final_metrics.csv`, `RF_final_varimp.csv`, `Figure4.tiff`, `Table3_LC_stratified_ampm.csv`, `Table4_Topo_stratified_ampm.csv`, `TableS2_LCxTopo_ampm.csv` |
