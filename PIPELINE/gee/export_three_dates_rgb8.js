// =====================================================
// 10 km² ROI + 3 single-day images (BEFORE / ON / AFTER)
// Export: 8-bit RGB "PNG-ready" GeoTIFFs (convert to PNG later)
// FIX: visualize() outputs have bands [vis-red, vis-green, vis-blue],
//      so Map.addLayer must NOT use rgbVis on those.
// =====================================================

// -------------------- SETTINGS --------------------
var lat = 7.2699;
var lon = 80.5938;

var targetDate = ee.Date('2025-11-30');

var preLookbackDays = 30;
var postForwardDays = 10;

var cloudProbThresh = 40;       // cloud prob threshold for cloudFrac metric
var maxCloudFrac = 0.02;        // allow up to 2% cloudy pixels in ROI (tighten to 0.01 if you want)
var minValidFrac = 0.995;       // require 99.5% valid coverage (tighten to 0.999 for stricter)

var onFallbackNearestDays = 3;  // if no exact image on targetDate, nearest within ±N days

var cloudTieTolerance = 3;      // mean cloud prob tie tolerance

// (Optional) inset export region to reduce tile-edge artifacts.
// If you will crop later anyway, you can keep it small.
var exportBufferMeters = 30;    // try 30; if still black corners, use 80–200

// Export settings
var outFolder = 'SL_pairs';
var outScale = 10;
var outCRS = 'EPSG:32644';      // UTM zone 44N

// Visualization stretch for SR bands
var rgbVis = {bands: ['B4','B3','B2'], min: 0, max: 3000, gamma: 1.1};

// -------------------- 10 km² square ROI (true meters using UTM) --------------------
// Area = 10 km^2 => side = sqrt(10) km ≈ 3.162 km => half-side ≈ 1581 m
var squareAreaKm2 = 10;
var halfSideMeters = Math.sqrt(squareAreaKm2 * 1e6) / 2;

var utm = ee.Projection(outCRS);
var roiUtm = ee.Geometry.Point([lon, lat]).transform(utm, 1).buffer(halfSideMeters).bounds();
var roiLL = roiUtm.transform('EPSG:4326', 1);

var exportUtm = roiUtm.buffer(-exportBufferMeters).bounds();
var exportLL = exportUtm.transform('EPSG:4326', 1);

Map.centerObject(roiLL, 13);
Map.addLayer(roiLL, {color: 'red'}, 'ROI (10 km^2)', false);
Map.addLayer(exportLL, {color: 'yellow'}, 'Export region (inset)', false);

// -------------------- DATA: Sentinel-2 SR + s2cloudless --------------------
var s2sr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(roiLL);
var s2cp = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY').filterBounds(roiLL);

var joined = ee.ImageCollection(ee.Join.saveFirst('cloudprob').apply({
  primary: s2sr,
  secondary: s2cp,
  condition: ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'})
}));

function addCloudProb(img) {
  var cp = ee.Image(img.get('cloudprob'));
  cp = ee.Algorithms.If(
    cp,
    ee.Image(cp).select('probability'),
    ee.Image.constant(100).rename('probability')
  );
  return img.addBands(ee.Image(cp).rename('probability'))
            .copyProperties(img, ['system:time_start', 'CLOUDY_PIXEL_PERCENTAGE']);
}

// Metrics inside export region:
// roiCloud  = mean cloud probability (lower is better)
// cloudFrac = fraction of pixels with prob > cloudProbThresh (lower is better)
// validFrac = fraction of valid pixels from B2 mask (higher is better)
function addMetrics(img) {
  var prob = img.select('probability');

  var roiCloud = ee.Number(prob.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: exportUtm,
    scale: 20,
    bestEffort: true,
    maxPixels: 1e13
  }).get('probability'));
  roiCloud = ee.Number(ee.Algorithms.If(roiCloud, roiCloud, 100));

  var cloudFrac = ee.Number(prob.gt(cloudProbThresh).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: exportUtm,
    scale: 20,
    bestEffort: true,
    maxPixels: 1e13
  }).get('probability'));
  cloudFrac = ee.Number(ee.Algorithms.If(cloudFrac, cloudFrac, 1));

  var validFrac = ee.Number(img.select('B2').mask().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: exportUtm,
    scale: 10,
    bestEffort: true,
    maxPixels: 1e13
  }).get('B2'));
  validFrac = ee.Number(ee.Algorithms.If(validFrac, validFrac, 0));

  return img.set('roiCloud', roiCloud)
            .set('cloudFrac', cloudFrac)
            .set('validFrac', validFrac);
}

function qualityFilter(col) {
  return col
    .filter(ee.Filter.gte('validFrac', minValidFrac))
    .filter(ee.Filter.lte('cloudFrac', maxCloudFrac));
}

// Select from a date window: strict quality -> fallback valid-only -> fallback all
// Then choose min roiCloud, tie-break by time preference
function selectMinCloudWithTimePref(dateStart, dateEnd, preferLatest) {
  var colAll = joined.filterDate(dateStart, dateEnd)
    .map(addCloudProb)
    .map(addMetrics);

  print('Window:', dateStart.format('YYYY-MM-dd'), '→', dateEnd.format('YYYY-MM-dd'),
        'All images:', colAll.size());

  var colStrict = qualityFilter(colAll);
  var colValid  = colAll.filter(ee.Filter.gte('validFrac', minValidFrac));

  var col = ee.ImageCollection(ee.Algorithms.If(
    colStrict.size().gt(0), colStrict,
    ee.Algorithms.If(colValid.size().gt(0), colValid, colAll)
  ));

  var empty = ee.Image().updateMask(ee.Image(0));

  var minCloud = ee.Number(col.aggregate_min('roiCloud'));
  var best = col.filter(ee.Filter.lte('roiCloud', minCloud.add(cloudTieTolerance)));

  // preferLatest => descending time_start, else ascending
  var sorted = best.sort('system:time_start', !preferLatest);
  var chosen = ee.Image(sorted.first());

  return ee.Image(ee.Algorithms.If(colAll.size().gt(0), chosen, empty));
}

// ON: exact day if possible else nearest (closest day first, then cloud)
function selectOnDayExactOrNearest(tDate, nearestDays) {
  var exactAll = joined.filterDate(tDate, tDate.advance(1, 'day'))
    .map(addCloudProb).map(addMetrics);

  var exactStrict = qualityFilter(exactAll);
  var exactValid  = exactAll.filter(ee.Filter.gte('validFrac', minValidFrac));

  var exact = ee.ImageCollection(ee.Algorithms.If(
    exactStrict.size().gt(0), exactStrict,
    ee.Algorithms.If(exactValid.size().gt(0), exactValid, exactAll)
  ));

  var nearestAll = joined
    .filterDate(tDate.advance(-nearestDays, 'day'), tDate.advance(nearestDays + 1, 'day'))
    .map(addCloudProb).map(addMetrics)
    .map(function(img) {
      var d = ee.Number(ee.Date(img.get('system:time_start')).difference(tDate, 'day')).abs();
      var key = d.multiply(1000).add(ee.Number(img.get('roiCloud'))); // closeness dominates
      return img.set('sortKey', key).set('absDiffDays', d);
    });

  var nearestStrict = qualityFilter(nearestAll);
  var nearestValid  = nearestAll.filter(ee.Filter.gte('validFrac', minValidFrac));

  var nearest = ee.ImageCollection(ee.Algorithms.If(
    nearestStrict.size().gt(0), nearestStrict,
    ee.Algorithms.If(nearestValid.size().gt(0), nearestValid, nearestAll)
  ));

  var exactSize = exactAll.size();
  var exactChosen = ee.Image(exact.sort('roiCloud').first());
  var nearestChosen = ee.Image(nearest.sort('sortKey').first());

  return ee.Image(ee.Algorithms.If(exactSize.gt(0), exactChosen, nearestChosen))
    .set('usedExactOnDay', exactSize.gt(0));
}

// -------------------- WINDOWS --------------------
var preStart  = targetDate.advance(-preLookbackDays, 'day');
var preEnd    = targetDate; // end exclusive => strictly before

var postStart = targetDate.advance(1, 'day'); // strictly after
var postEnd   = targetDate.advance(postForwardDays + 1, 'day');

// -------------------- PICK IMAGES --------------------
var beforeImg = selectMinCloudWithTimePref(preStart, preEnd, true);     // latest among best cloud
var onImg     = selectOnDayExactOrNearest(targetDate, onFallbackNearestDays);
var afterImg  = selectMinCloudWithTimePref(postStart, postEnd, false);  // earliest among best cloud

// -------------------- 8-bit RGB outputs (PNG-ready) --------------------
// These images now have bands: [vis-red, vis-green, vis-blue]
var before8 = beforeImg.visualize(rgbVis).clip(exportUtm);
var on8     = onImg.visualize(rgbVis).clip(exportUtm);
var after8  = afterImg.visualize(rgbVis).clip(exportUtm);

// ---- DISPLAY ----
// IMPORTANT: do NOT pass rgbVis to these visualized images
Map.addLayer(before8, {min: 0, max: 255}, 'BEFORE (RGB8)');
Map.addLayer(on8,     {min: 0, max: 255}, 'ON (RGB8)');
Map.addLayer(after8,  {min: 0, max: 255}, 'AFTER (RGB8)');

// (Optional) also view raw SR using rgbVis (bands B4/B3/B2 exist here):
Map.addLayer(beforeImg.clip(exportUtm), rgbVis, 'BEFORE raw (SR)', false);
Map.addLayer(onImg.clip(exportUtm),     rgbVis, 'ON raw (SR)', false);
Map.addLayer(afterImg.clip(exportUtm),  rgbVis, 'AFTER raw (SR)', false);

// -------------------- PRINT METADATA --------------------
function printInfo(label, img) {
  print(label,
    'date:', ee.Date(img.get('system:time_start')).format('YYYY-MM-dd'),
    'roiCloud:', img.get('roiCloud'),
    'cloudFrac:', img.get('cloudFrac'),
    'validFrac:', img.get('validFrac'),
    'scene cloud%:', img.get('CLOUDY_PIXEL_PERCENTAGE')
  );
}
printInfo('BEFORE', beforeImg);
printInfo('ON', onImg);
printInfo('AFTER', afterImg);
print('ON usedExactOnDay?:', onImg.get('usedExactOnDay'));

// -------------------- PNG THUMB LINKS --------------------
print('BEFORE PNG (thumb):', before8.getThumbURL({region: exportLL, dimensions: 512, format: 'png'}));
print('ON PNG (thumb):',     on8.getThumbURL({region: exportLL, dimensions: 512, format: 'png'}));
print('AFTER PNG (thumb):',  after8.getThumbURL({region: exportLL, dimensions: 512, format: 'png'}));

// -------------------- EXPORTS --------------------
var dateDict = ee.Dictionary({
  before: ee.Date(beforeImg.get('system:time_start')).format('YYYYMMdd'),
  on:     ee.Date(onImg.get('system:time_start')).format('YYYYMMdd'),
  after:  ee.Date(afterImg.get('system:time_start')).format('YYYYMMdd')
});

dateDict.evaluate(function(d) {
  Export.image.toDrive({
    image: before8,
    description: 'BEFORE_RGB8_' + d.before,
    folder: outFolder,
    fileNamePrefix: 'BEFORE_RGB8_' + d.before + '_pre' + preLookbackDays + 'd',
    region: exportUtm,
    scale: outScale,
    crs: outCRS,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    skipEmptyTiles: true,
    formatOptions: {cloudOptimized: true}
  });

  Export.image.toDrive({
    image: on8,
    description: 'ON_RGB8_' + d.on,
    folder: outFolder,
    fileNamePrefix: 'ON_RGB8_' + d.on + '_target',
    region: exportUtm,
    scale: outScale,
    crs: outCRS,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    skipEmptyTiles: true,
    formatOptions: {cloudOptimized: true}
  });

  Export.image.toDrive({
    image: after8,
    description: 'AFTER_RGB8_' + d.after,
    folder: outFolder,
    fileNamePrefix: 'AFTER_RGB8_' + d.after + '_post' + postForwardDays + 'd',
    region: exportUtm,
    scale: outScale,
    crs: outCRS,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    skipEmptyTiles: true,
    formatOptions: {cloudOptimized: true}
  });
});

