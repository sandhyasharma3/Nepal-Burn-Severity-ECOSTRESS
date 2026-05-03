# Author : Sandhya Sharma
# Year   : 2026
#
# Burn Severity Distribution — Halfeye Plot
# Nepal 2021 Fire Season
#
# Produces a half-eye (raincloud) plot showing the distribution of dNBR
# values across burn severity classes (Unburned, Low, Moderate, High).
# =============================================================================
 
library(tidyverse)
library(ggplot2)
library(ggdist)
library(terra)
library(cowplot)
 
# -----------------------------------------------------------------------------
# USER-DEFINED PATHS — update before running
# -----------------------------------------------------------------------------
 
raster_dir <- "path/to/your/raster/directory"   # <-- SET THIS
burn_file  <- file.path(raster_dir, "dNBRfiltered.tif")
out_file   <- file.path(raster_dir, "ridgeline_plot.png")
 
# -----------------------------------------------------------------------------
# 1. Load dNBR raster and convert to data frame
# -----------------------------------------------------------------------------
 
burn_rast <- rast(burn_file)
 
plot_df <- as.data.frame(burn_rast, xy = FALSE) %>%
  drop_na() %>%
  filter(dNBRfiltered > -0.1)
 
# -----------------------------------------------------------------------------
# 2. Classify burn severity
# Standard USGS dNBR thresholds
# -----------------------------------------------------------------------------
 
plot_df <- plot_df %>%
  mutate(severity_class = case_when(
    dNBRfiltered < 0.1  ~ "Unburned",
    dNBRfiltered < 0.27 ~ "Low",
    dNBRfiltered < 0.44 ~ "Moderate",
    TRUE                ~ "High"
  )) %>%
  mutate(severity_class = factor(severity_class,
                                 levels = c("Unburned", "Low", "Moderate", "High")))
 
# -----------------------------------------------------------------------------
# 3. Plot
# -----------------------------------------------------------------------------
 
p1 <- plot_df %>%
  ggplot(aes(y = severity_class, x = dNBRfiltered, fill = severity_class)) +
  stat_halfeye(alpha = 0.8) +
  scale_fill_manual(values = c(
    "Unburned" = "#fed98e",
    "Low"      = "#fe9929",
    "Moderate" = "#cc4c02",
    "High"     = "#993404"
  )) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x    = "dNBR",
    y    = NULL,
    fill = "Burn Severity"
  ) +
  theme(legend.position = "none")
 
print(p1)
 
# -----------------------------------------------------------------------------
# 4. Save
# -----------------------------------------------------------------------------
 
ggsave(
  out_file,
  plot   = p1,
  width  = 18,
  height = 12,
  units  = "cm",
  dpi    = 300
)
 
cat("Saved:", out_file, "\n")
