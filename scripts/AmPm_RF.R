# =============================================================================
# Author : Sandhya Sharma
# Year   : 2026
#
# Burn Severity Prediction - AM/PM ECOSTRESS Model
# Random Forest Analysis: Nepal 2021 Fire Season
#
# Workflow:
#   1. Load data, extract rasters, build modelling table
#   2. Preliminary RF with all 20 PCNMs -> identify top 3 by importance
#   3. Final RF with top 3 PCNMs (80/20 train/test split)
#   4. Stratified RF by land cover class
#   5. Stratified RF by topographic class
#   6. Stratified RF by land cover x topographic class
#
# AM = ECOSTRESS instantaneous acquisitions before noon (ISS overpass time)
# PM = ECOSTRESS instantaneous acquisitions after noon (ISS overpass time)
# Both derived from ECO_L3T_JET_002_PTJPLSMinst product
# =============================================================================

library(sf)
library(terra)
library(dplyr)
library(randomForest)
library(ggplot2)
library(ggpointdensity)
library(scales)
library(patchwork)

# -----------------------------------------------------------------------------
# paths
# -----------------------------------------------------------------------------

rdir <- "D:/Users/sharm201/Desktop/RESEARCH2NDHALF/standardized_rasters/clean_rasters"

pts_file  <- file.path(rdir, "Am_pmDATA_100k_clean.rds")
pcnm_file <- file.path(rdir, "newpcnms", "pcnm_true_150km_radius_am_pm.rds")
burn_file <- file.path(rdir, "dNBRfiltered.tif")
out_dir   <- file.path(rdir, "ampm_RF_results")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

min_n_lc      <- 100
min_n_topo    <- 50
min_n_lc_topo <- 50

burn_threshold <- -0.1

lulc_codes <- c(
  "1"  = "Evergreen Needleleaf",
  "2"  = "Evergreen Broadleaf",
  "4"  = "Deciduous Broadleaf",
  "5"  = "Mixed Forest",
  "8"  = "Woody Savanna",
  "9"  = "Savanna",
  "10" = "Grassland"
)

aspect_levels <- c("North", "Northeast", "East", "Southeast",
                   "South", "Southwest", "West", "Northwest")
elev_levels   <- c("Low", "Medium", "High", "Very High")
slope_levels  <- c("Gentle", "Medium", "Steep", "Very Steep")

okabe_ito <- c(
  "Spatial"      = "#0072B2",
  "Hydrological" = "#009E73",
  "Topographic"  = "#D55E00"
)

# -----------------------------------------------------------------------------
# load sampling points
# -----------------------------------------------------------------------------

pts_raw <- readRDS(pts_file)

if (class(pts_raw)[1] == "sf") {
  coords <- st_coordinates(pts_raw)
  pts    <- as.data.frame(st_drop_geometry(pts_raw))
  pts$X  <- coords[, 1]
  pts$Y  <- coords[, 2]
  pts_vect <- vect(pts[, c("X", "Y")], geom = c("X", "Y"),
                   crs = st_crs(pts_raw)$wkt)
} else {
  pts      <- as.data.frame(pts_raw)
  pts_vect <- vect(pts[, c("X", "Y")], geom = c("X", "Y"),
                   crs = "EPSG:32644")
}

cat("Points loaded:", nrow(pts), "\n")

# -----------------------------------------------------------------------------
# load PCNM vectors
# calculated on all 40,000 points before NA removal
# -----------------------------------------------------------------------------

pcnm_raw <- readRDS(pcnm_file)

if (!is.null(pcnm_raw$vectors)) {
  pcnm_mat <- pcnm_raw$vectors
} else {
  pcnm_mat <- as.matrix(pcnm_raw)
}

pcnm20 <- as.data.frame(pcnm_mat[, 1:20])
colnames(pcnm20) <- paste0("PCNM", 1:20)

cat("PCNM rows:", nrow(pcnm20), "- Points rows:", nrow(pts), "\n")

# -----------------------------------------------------------------------------
# load rasters and extract values
# AM and PM variables from instantaneous ECOSTRESS product
# separated by ISS overpass time before this script
# -----------------------------------------------------------------------------

raster_files <- c(
  file.path(rdir, "Aspect_deg.tif"),
  file.path(rdir, "Elevation.tif"),
  file.path(rdir, "Slope.tif"),
  burn_file,
  file.path(rdir, "LC.tif"),
  file.path(rdir, "ESI_AM.tif"),
  file.path(rdir, "ESI_PM.tif"),
  file.path(rdir, "ET_AM.tif"),
  file.path(rdir, "ET_PM.tif"),
  file.path(rdir, "WUE_AM.tif"),
  file.path(rdir, "WUE_PM.tif")
)

rstack <- rast(raster_files)
names(rstack) <- c("Aspect_deg", "Elevation", "Slope", "dNBR", "LC",
                   "ESI_AM", "ESI_PM", "ET_AM", "ET_PM", "WUE_AM", "WUE_PM")

if (!same.crs(pts_vect, rstack)) {
  pts_vect <- project(pts_vect, crs(rstack))
}

vals <- extract(rstack, pts_vect, ID = FALSE)

cat("NA count per raster layer:\n")
print(colSums(is.na(vals)))

# -----------------------------------------------------------------------------
# build modelling table
# -----------------------------------------------------------------------------

dat <- cbind(pts[, c("X", "Y")], pcnm20, vals)

dat$Aspect_deg <- dat$Aspect_deg %% 360

dat$Elevation_cat <- cut(dat$Elevation,
                         breaks = c(-Inf, 1000, 2500, 5000, Inf),
                         labels = elev_levels, right = FALSE)

dat$Slope_cat <- cut(dat$Slope,
                     breaks = c(-Inf, 5, 15, 30, Inf),
                     labels = slope_levels, right = FALSE)

dat$Aspect_dir <- cut(dat$Aspect_deg,
                      breaks = c(-Inf, 22.5, 67.5, 112.5, 157.5,
                                 202.5, 247.5, 292.5, 337.5, Inf),
                      labels = c("North", "Northeast", "East", "Southeast",
                                 "South", "Southwest", "West", "Northwest", "North"),
                      right = FALSE)

dat$Aspect_dir[dat$Aspect_deg >= 337.5 | dat$Aspect_deg < 22.5] <- "North"
dat$Aspect_dir    <- factor(dat$Aspect_dir, levels = aspect_levels)
dat$Elevation_cat <- factor(dat$Elevation_cat, levels = elev_levels)
dat$Slope_cat     <- factor(dat$Slope_cat, levels = slope_levels)

dat <- dat[!is.na(dat$dNBR) & dat$dNBR >= burn_threshold, ]
names(dat)[names(dat) == "dNBR"] <- "newdNBR"

cat("Rows after dNBR filter:", nrow(dat), "\n")

dat$LC <- as.character(dat$LC)
dat <- dat[dat$LC %in% names(lulc_codes), ]
dat$LC_short <- factor(lulc_codes[dat$LC], levels = unname(lulc_codes))

hydro_vars <- c("ESI_AM", "ESI_PM", "ET_AM", "ET_PM", "WUE_AM", "WUE_PM")
topo_vars  <- c("Aspect_dir", "Elevation", "Slope")
pcnm_all   <- paste0("PCNM", 1:20)

keep_cols <- c("newdNBR", "LC_short", "Elevation_cat", "Slope_cat",
               "Aspect_dir", topo_vars, hydro_vars, pcnm_all)

dat <- dat[complete.cases(dat[, keep_cols]), ]

cat("Final modelling rows:", nrow(dat), "\n")

saveRDS(dat, file.path(out_dir, "modelling_table.rds"))
write.csv(dat, file.path(out_dir, "modelling_table.csv"), row.names = FALSE)

# =============================================================================
# STAGE 1 - preliminary RF with all 20 PCNMs
# used only to identify the top 3 most important PCNMs
# results from this stage are not reported in the paper
# =============================================================================

predictors_stage1 <- c(topo_vars, hydro_vars, pcnm_all)

df_stage1 <- dat[complete.cases(dat[, c("newdNBR", predictors_stage1)]), ]
df_stage1$Aspect_dir <- factor(df_stage1$Aspect_dir, levels = aspect_levels)

set.seed(123)
train_idx1 <- sample(nrow(df_stage1), floor(0.8 * nrow(df_stage1)))

set.seed(123)
rf_stage1 <- randomForest(
  x = df_stage1[train_idx1, predictors_stage1],
  y = df_stage1$newdNBR[train_idx1],
  ntree = 500,
  importance = TRUE
)

cat("Stage 1 OOB R2:", round(tail(rf_stage1$rsq, 1), 3), "\n")

imp_stage1 <- importance(rf_stage1, type = 1)[, "%IncMSE"]
pcnm_imp   <- sort(imp_stage1[pcnm_all], decreasing = TRUE)
top3_pcnms <- names(pcnm_imp)[1:3]

cat("Top 3 PCNMs:", paste(top3_pcnms, collapse = ", "), "\n")

write.csv(data.frame(PCNM = names(pcnm_imp), IncMSE = pcnm_imp),
          file.path(out_dir, "Stage1_PCNM_importance.csv"), row.names = FALSE)

# =============================================================================
# STAGE 2 - final RF with top 3 PCNMs, 80/20 train/test split
# =============================================================================

predictors_final <- c(topo_vars, hydro_vars, top3_pcnms)

df_model <- dat[complete.cases(dat[, c("newdNBR", predictors_final)]), ]
df_model$Aspect_dir <- factor(df_model$Aspect_dir, levels = aspect_levels)

set.seed(123)
train_idx <- sample(nrow(df_model), floor(0.8 * nrow(df_model)))
df_train  <- df_model[train_idx, ]
df_test   <- df_model[-train_idx, ]

cat("Train n:", nrow(df_train), "| Test n:", nrow(df_test), "\n")

set.seed(123)
rf_final <- randomForest(
  x = df_train[, predictors_final],
  y = df_train$newdNBR,
  ntree = 500,
  importance = TRUE
)

oob_r2 <- round(tail(rf_final$rsq, 1), 3)
cat("Final model OOB R2:", oob_r2, "\n")

test_pred <- predict(rf_final, newdata = df_test[, predictors_final])
test_r2   <- round(cor(df_test$newdNBR, test_pred)^2, 3)
test_rmse <- round(sqrt(mean((df_test$newdNBR - test_pred)^2)), 4)

cat("Test R2:", test_r2, "| RMSE:", test_rmse, "\n")

metrics <- data.frame(
  Metric = c("OOB_R2", "Test_R2", "Test_RMSE", "Train_n", "Test_n"),
  Value  = c(oob_r2, test_r2, test_rmse, nrow(df_train), nrow(df_test))
)
write.csv(metrics, file.path(out_dir, "RF_final_metrics.csv"), row.names = FALSE)
saveRDS(rf_final, file.path(out_dir, "RF_final_model.rds"))

# variable importance
imp_final <- importance(rf_final, type = 1)[, "%IncMSE"]

imp_df <- data.frame(
  Variable = names(imp_final),
  IncMSE = as.numeric(imp_final),
  stringsAsFactors = FALSE
)
imp_df <- imp_df[order(imp_df$IncMSE, decreasing = TRUE), ]

imp_df$Group <- ifelse(grepl("^PCNM", imp_df$Variable), "Spatial",
               ifelse(grepl("ESI|ET|WUE", imp_df$Variable), "Hydrological",
               "Topographic"))

imp_df$Label <- imp_df$Variable
imp_df$Label[imp_df$Variable == "Aspect_dir"] <- "Aspect"
imp_df$Label[imp_df$Variable == "Elevation"]  <- "Elevation"
imp_df$Label[imp_df$Variable == "Slope"]      <- "Slope"
imp_df$Label[imp_df$Variable == "ESI_AM"]     <- "ESI (AM)"
imp_df$Label[imp_df$Variable == "ESI_PM"]     <- "ESI (PM)"
imp_df$Label[imp_df$Variable == "ET_AM"]      <- "ET (AM)"
imp_df$Label[imp_df$Variable == "ET_PM"]      <- "ET (PM)"
imp_df$Label[imp_df$Variable == "WUE_AM"]     <- "WUE (AM)"
imp_df$Label[imp_df$Variable == "WUE_PM"]     <- "WUE (PM)"

write.csv(imp_df, file.path(out_dir, "RF_final_varimp.csv"), row.names = FALSE)

# figures
df_test_pred <- data.frame(
  Observed  = df_test$newdNBR,
  Predicted = as.numeric(test_pred)
)
lims <- range(c(df_test_pred$Observed, df_test_pred$Predicted), na.rm = TRUE)

p_imp <- ggplot(imp_df, aes(x = reorder(Label, IncMSE), y = IncMSE, fill = Group)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = okabe_ito) +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "% Increase in MSE", fill = "Predictor group") +
  theme(legend.position = "bottom")

p_pred <- ggplot(df_test_pred, aes(x = Observed, y = Predicted)) +
  geom_pointdensity(size = 0.8) +
  scale_color_viridis_c(option = "cividis", trans = "log10",
                        labels = label_scientific(), name = "Point density") +
  geom_abline(linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "steelblue", linewidth = 1) +
  coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
  theme_bw(base_size = 13) +
  labs(
    x = "Observed dNBR",
    y = "Predicted dNBR",
    subtitle = paste0("Training OOB R\u00b2 = ", oob_r2,
                      " | Test R\u00b2 = ", test_r2,
                      " | N = ", nrow(df_model))
  )

fig4 <- (p_imp + labs(tag = "A")) | (p_pred + labs(tag = "B"))

ggsave(file.path(out_dir, "Figure4.tiff"), fig4,
       device = "tiff", dpi = 600, width = 14, height = 7,
       units = "in", compression = "lzw")

ggsave(file.path(out_dir, "Figure4.png"), fig4,
       dpi = 300, width = 14, height = 7, units = "in")

# =============================================================================
# stratified RF by land cover class
# =============================================================================

cat("\n--- LC stratified RF ---\n")

predictors_strat <- c(hydro_vars, top3_pcnms)
df_lc <- dat[complete.cases(dat[, c("newdNBR", "LC_short", predictors_strat)]), ]

lc_results <- list()

for (lc in levels(df_lc$LC_short)) {

  d <- df_lc[df_lc$LC_short == lc, ]
  cat(lc, "| n =", nrow(d), "\n")

  if (nrow(d) < min_n_lc) next

  is_const <- sapply(d[, predictors_strat], function(x) length(unique(na.omit(x))) <= 1)
  preds    <- predictors_strat[!is_const]

  set.seed(123)
  rf  <- randomForest(x = d[, preds], y = d$newdNBR, ntree = 500, importance = TRUE)
  imp <- sort(importance(rf, type = 1)[, "%IncMSE"], decreasing = TRUE)

  lc_results[[lc]] <- data.frame(
    LC_short = lc,
    n        = nrow(d),
    OOB_R2   = round(tail(rf$rsq, 1), 3),
    Top1     = names(imp)[1],
    Top2     = names(imp)[2],
    Top3     = names(imp)[3],
    stringsAsFactors = FALSE
  )
}

lc_table <- do.call(rbind, lc_results)
lc_table  <- lc_table[order(lc_table$OOB_R2, decreasing = TRUE), ]
write.csv(lc_table, file.path(out_dir, "Table3_LC_stratified_ampm.csv"), row.names = FALSE)
print(lc_table)

# =============================================================================
# stratified RF by topographic class
# =============================================================================

cat("\n--- Topography stratified RF ---\n")

df_topo <- dat[complete.cases(dat[, c("newdNBR", "Elevation_cat", "Slope_cat",
                                      "Aspect_dir", predictors_strat)]), ]
topo_results <- list()

for (topo_var in c("Elevation_cat", "Slope_cat", "Aspect_dir")) {

  classes <- levels(droplevels(df_topo[[topo_var]]))

  for (cls in classes) {

    d <- df_topo[df_topo[[topo_var]] == cls, ]
    cat(topo_var, "=", cls, "| n =", nrow(d), "\n")

    if (nrow(d) < min_n_topo) next

    is_const <- sapply(d[, predictors_strat], function(x) length(unique(na.omit(x))) <= 1)
    preds    <- predictors_strat[!is_const]

    set.seed(123)
    rf  <- randomForest(x = d[, preds], y = d$newdNBR, ntree = 500, importance = TRUE)
    imp <- sort(importance(rf, type = 1)[, "%IncMSE"], decreasing = TRUE)

    key <- paste(topo_var, cls, sep = "_")
    topo_results[[key]] <- data.frame(
      Topography = topo_var,
      Class      = cls,
      n          = nrow(d),
      OOB_R2     = round(tail(rf$rsq, 1), 3),
      Top1       = names(imp)[1],
      Top2       = names(imp)[2],
      Top3       = names(imp)[3],
      stringsAsFactors = FALSE
    )
  }
}

topo_table <- do.call(rbind, topo_results)
topo_table  <- topo_table[order(topo_table$OOB_R2, decreasing = TRUE), ]
write.csv(topo_table, file.path(out_dir, "Table4_Topo_stratified_ampm.csv"), row.names = FALSE)
print(topo_table)

# =============================================================================
# stratified RF by land cover x topographic class (Table S2 - Supplement)
# =============================================================================

cat("\n--- LC x Topography stratified RF ---\n")

df_lc_topo <- dat[complete.cases(dat[, c("newdNBR", "LC_short", "Elevation_cat",
                                         "Slope_cat", "Aspect_dir",
                                         predictors_strat)]), ]
lc_topo_results <- list()

for (topo_var in c("Elevation_cat", "Slope_cat", "Aspect_dir")) {
  for (lc in levels(df_lc_topo$LC_short)) {
    for (cls in levels(droplevels(df_lc_topo[[topo_var]]))) {

      d <- df_lc_topo[df_lc_topo$LC_short == lc &
                      df_lc_topo[[topo_var]] == cls, ]

      if (nrow(d) < min_n_lc_topo) next

      is_const <- sapply(d[, predictors_strat], function(x) length(unique(na.omit(x))) <= 1)
      preds    <- predictors_strat[!is_const]
      if (length(preds) == 0) next

      set.seed(123)
      rf  <- randomForest(x = d[, preds], y = d$newdNBR, ntree = 500, importance = TRUE)
      imp <- sort(importance(rf, type = 1)[, "%IncMSE"], decreasing = TRUE)

      key <- paste(lc, topo_var, cls, sep = "_")
      lc_topo_results[[key]] <- data.frame(
        LC_short   = lc,
        Topography = topo_var,
        Class      = cls,
        n          = nrow(d),
        OOB_R2     = round(tail(rf$rsq, 1), 3),
        Top1       = names(imp)[1],
        Top2       = names(imp)[2],
        Top3       = names(imp)[3],
        stringsAsFactors = FALSE
      )
    }
  }
}

lc_topo_table <- do.call(rbind, lc_topo_results)
lc_topo_table  <- lc_topo_table[order(lc_topo_table$OOB_R2, decreasing = TRUE), ]
write.csv(lc_topo_table, file.path(out_dir, "TableS2_LCxTopo_ampm.csv"), row.names = FALSE)

cat("LC x Topo AM/PM done. Top rows:\n")
print(head(lc_topo_table, 10))

cat("\nAll outputs saved to:", out_dir, "\n")
