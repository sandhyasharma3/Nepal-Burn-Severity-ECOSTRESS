// ============================================================
// Nepal dNBR (pre - post) — Sentinel-2 burn severity
// Author: Sandhya Sharma
// Year:   2026
//
// Cloud masking handles both:
//   - QA60 version (older S2 collections)
//   - MSK_CLASSI version (newer S2 SR collections)
// ============================================================
 
// ---------------------------
// 1) ROI
// ---------------------------
var roiFC = ee.FeatureCollection('projects/ee-phdsandhya3/assets/Country_Boundary');
var roi   = roiFC.geometry();
 
Map.setCenter(84.0, 28.3, 6);
Map.addLayer(roiFC, {}, 'Country_Boundary');
 
// ---------------------------
// 2) Date windows
// ---------------------------
var preStart  = '2020-01-01';
var preEnd    = '2020-03-28';
var postStart = '2021-06-01';
var postEnd   = '2021-06-30';
 
// ---------------------------
// 3) Cloud mask function
// Applies QA60, MSK_CLASSI, and SCL masks where available
// ---------------------------
function maskS2(img) {
  var bandNames = img.bandNames();
 
  // QA60 branch
  var qaMask = ee.Image(
    ee.Algorithms.If(
      bandNames.contains('QA60'),
      (function() {
        var qa = img.select('QA60');
        var cloudBitMask  = 1 << 10;
        var cirrusBitMask = 1 << 11;
        return qa.bitwiseAnd(cloudBitMask).eq(0)
                 .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
      })(),
      ee.Image(1)
    )
  );
 
  // MSK_CLASSI branch
  var mskMask = ee.Image(
    ee.Algorithms.If(
      bandNames.contains('MSK_CLASSI_OPAQUE'),
      img.select('MSK_CLASSI_OPAQUE').eq(0)
        .and(img.select('MSK_CLASSI_CIRRUS').eq(0)),
      ee.Image(1)
    )
  );
 
  // Combined cloud mask
  var cloudMask = qaMask.and(mskMask);
 
  // SCL cleanup: remove cloud shadow (3) and snow/ice (11)
  var sclMask = img.select('SCL')
    .neq(3)
    .and(img.select('SCL').neq(11));
 
  return img.updateMask(cloudMask).updateMask(sclMask);
}
 
// ---------------------------
// 4) Image collection
// ---------------------------
var s2 = ee.ImageCollection('COPERNICUS/S2_SR')
  .filterBounds(roi)
  .filterDate(preStart, postEnd)
  .map(maskS2)
  .select(['B8A', 'B12']);
 
var preCol  = s2.filterDate(preStart, preEnd);
var postCol = s2.filterDate(postStart, postEnd);
 
print('Pre-fire images:', preCol.size());
print('Post-fire images:', postCol.size());
 
// ---------------------------
// 5) Median composites
// ---------------------------
var preMed  = preCol.median().divide(10000);
var postMed = postCol.median().divide(10000);
 
// ---------------------------
// 6) NBR and dNBR
// NBR = (NIR - SWIR) / (NIR + SWIR)
// dNBR = pre-NBR - post-NBR
// ---------------------------
function nbr(img) {
  var nir  = img.select('B8A');
  var swir = img.select('B12');
  return nir.subtract(swir)
            .divide(nir.add(swir))
            .rename('NBR');
}
 
var dNBR = nbr(preMed)
  .subtract(nbr(postMed))
  .rename('dNBR');
 
// ---------------------------
// 7) Burn severity classification
// Based on USGS dNBR thresholds
// 0 = unburned, 1 = low, 2 = moderate, 3 = high
// ---------------------------
var severity = ee.Image(0)
  .where(dNBR.gte(0.10).and(dNBR.lt(0.27)), 1)
  .where(dNBR.gte(0.27).and(dNBR.lt(0.66)), 2)
  .where(dNBR.gte(0.66), 3)
  .rename('severity');
 
// ---------------------------
// 8) Visualization
// ---------------------------
Map.addLayer(
  dNBR.clip(roi),
  {min: -0.2, max: 0.8, palette: ['#1a9850', '#ffffbf', '#d73027']},
  'dNBR'
);
 
Map.addLayer(
  severity.clip(roi),
  {min: 0, max: 3, palette: ['#2c7bb6', '#abd9e9', '#fdae61', '#d7191c']},
  'Severity Class'
);
 
// ---------------------------
// 9) Export to Google Drive
// Update 'folder' to your Drive folder name
// ---------------------------
var exportFolder = 'YOUR_DRIVE_FOLDER_NAME';
 
Export.image.toDrive({
  image:       dNBR.clip(roi),
  description: 'Nepal_dnbr',
  folder:      exportFolder,
  region:      roi,
  scale:       20,
  maxPixels:   1e13
});
 
Export.image.toDrive({
  image:       severity.clip(roi).toInt16(),
  description: 'Nepal_BurnSeverity',
  folder:      exportFolder,
  region:      roi,
  scale:       20,
  maxPixels:   1e13
});
 
