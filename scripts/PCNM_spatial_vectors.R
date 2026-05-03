Pcnm spatial vectors · R
# =============================================================================
# Author : Sandhya Sharma
# Year   : 2026
#
# PCNM Spatial Eigenvectors — Nepal 2021 Fire Season
#
# Generates Principal Coordinates of Neighbour Matrices (PCNM) vectors
# for ~40,000 systematic grid points across Nepal, using a 150 km
# distance threshold. Output is used as spatial predictors in the
# Random Forest models (Annual_RF.R and AmPm_RF.R).
#
# NOTE: Computing a distance matrix on 40,000 points is memory-intensive.
#       This script was run on an HPC cluster (SLURM). Running locally
#       requires sufficient RAM (>= 32 GB recommended).
# =============================================================================
 
library(terra)
library(sf)
library(dplyr)
library(vegan)
library(ggplot2)
library(adespatial)
 
# -----------------------------------------------------------------------------
# USER-DEFINED PATHS — update before running
# -----------------------------------------------------------------------------
 
raster_dir    <- "path/to/your/raster/directory"   # <-- SET THIS
boundary_path <- "path/to/New_boundary.shp"         # <-- SET THIS
 
# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
 
n.points      <- 40000
n.points.name <- "40k"
radius        <- 150000    # 150 km in metres
radius.name   <- "150km"
 
# Output files (saved into raster_dir)
points_file    <- file.path(raster_dir, paste0("systematic_points_", n.points.name, ".rds"))
pcnm_file_full <- file.path(raster_dir, paste0("pcnm_full_", n.points.name, "_",
                                                radius.name, "_radius.rds"))
 
# -----------------------------------------------------------------------------
# 1. Load boundary
# -----------------------------------------------------------------------------
 
boundary      <- st_read(boundary_path)
boundary_vect <- vect(boundary)
 
# -----------------------------------------------------------------------------
# 2. Generate systematic grid points (~40k) clipped to boundary
# -----------------------------------------------------------------------------
 
area_m2 <- expanse(boundary_vect, unit = "m")
spacing  <- sqrt(area_m2 / n.points)
 
grid_sf <- st_make_grid(boundary, cellsize = spacing, what = "centers") %>%
  st_as_sf() %>%
  st_intersection(boundary)
 
grid_points <- st_cast(grid_sf, "POINT")
saveRDS(grid_points, file = points_file)
 
cat("Generated points:", nrow(grid_points), "\n")
 
# -----------------------------------------------------------------------------
# 3. Compute PCNM with 150 km threshold
# PCNMs are computed on all points before any NA removal so that
# eigenvectors capture the full spatial extent of the study area.
# -----------------------------------------------------------------------------
 
coords    <- st_coordinates(grid_points)
dist_mat  <- dist(coords)
 
cat("PCNM start time:", as.character(Sys.time()), "\n")
pcnm_obj <- pcnm(dist_mat, threshold = radius)
cat("PCNM finish time:", as.character(Sys.time()), "\n")
 
# -----------------------------------------------------------------------------
# 4. Save outputs
# -----------------------------------------------------------------------------
 
saveRDS(pcnm_obj, file = pcnm_file_full)
 
cat("Saved full PCNM object ->", pcnm_file_full, "\n")
cat("Positive eigenvalues:", sum(pcnm_obj$values > 0), "\n")
 
