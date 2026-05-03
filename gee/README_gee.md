Running the GEE Script
Prerequisites

A free Google Earth Engine account
Access to the Nepal boundary asset (projects/ee-phdsandhya3/assets/Country_Boundary)
(If reproducing independently, substitute your own Nepal boundary shapefile uploaded as a GEE asset)

Steps

Go to code.earthengine.google.com
Create a new script and paste the contents of Nepal_dNBR.js
Update the Folder variable (line ~135) to the name of your Google Drive folder where outputs should be saved
Click Run — the map will display dNBR and severity class layers
In the Tasks panel (top right), click Run next to each export task to write the rasters to your Drive

Outputs
FileDescriptionNepal_dnbr.tifContinuous dNBR raster at 20 m resolutionNepal_BurnSeverity.tifClassified severity (0 = unburned, 1 = low, 2 = moderate, 3 = high)
Date Windows Used
PeriodStartEndPre-fire2020-01-012020-03-28Post-fire2021-06-012021-06-30
Cloud Masking
The script applies a robust three-layer mask:

QA60 — cloud and cirrus bits
MSK_CLASSI — opaque and cirrus classification (newer S2 collections)
SCL — cloud shadow (class 3) and snow/ice (class 11)
