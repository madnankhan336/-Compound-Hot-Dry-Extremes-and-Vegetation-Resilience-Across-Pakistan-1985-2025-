/*
Google Earth Engine script
Project: Compound Hot–Dry Extremes and Vegetation Resilience Across Pakistan (1985–2025)

Purpose
- Builds monthly Landsat NDVI composites for Pakistan
- Uses Landsat 5 TM, Landsat 7 ETM+, Landsat 8 OLI, and Landsat 9 OLI-2
- Applies QA masking
- Computes monthly P75 NDVI composites
- Aggregates to a common analysis grid
- Optionally fills short gaps (<= 2 months)
- Exports one 12-band GeoTIFF per year

Important note
In Google Earth Engine, your shapefile must first be uploaded as an Asset.
Then place the asset path in CUSTOM_BOUNDARY_ASSET below.
*/

var START_YEAR = 1985;
var END_YEAR = 2025;

// ==============================
// USER SETTINGS
// ==============================
var CRS = 'EPSG:4326';
var SCALE = 10000;              // analysis grid resolution in metres
var MAX_GAP = 2;                // only short-gap filling up to 2 months
var MASK_WATER = true;
var FILL_SHORT_GAPS = true;

var EXPORT_FOLDER = 'Pakistan_Landsat_Monthly_NDVI';
var EXPORT_PREFIX = 'pakistan_landsat_monthly_ndvi_';

// Set to true when you upload your official Pakistan boundary shapefile as a GEE asset
var USE_CUSTOM_BOUNDARY = true;

// Placeholder: replace with your own uploaded asset path
var CUSTOM_BOUNDARY_ASSET = 'users/your_username/your_pakistan_boundary_asset';

// ==============================
// STUDY AREA
// ==============================
var pakistanFC = USE_CUSTOM_BOUNDARY
  ? ee.FeatureCollection(CUSTOM_BOUNDARY_ASSET)
  : ee.FeatureCollection('FAO/GAUL/2015/level0')
      .filter(ee.Filter.eq('ADM0_NAME', 'Pakistan'));

var pakistan = pakistanFC.geometry();

Map.centerObject(pakistan, 6);
Map.addLayer(pakistan, {color: 'red'}, 'Pakistan boundary');

// ==============================
// HELPER FUNCTIONS
// ==============================
function scaleSR(image) {
  var optical = image.select('SR_B.*').multiply(0.0000275).add(-0.2);
  return image.addBands(optical, null, true);
}

function maskSR(image) {
  var qa = image.select('QA_PIXEL');

  var mask = qa.bitwiseAnd(1 << 0).eq(0)   // fill
    .and(qa.bitwiseAnd(1 << 1).eq(0))      // dilated cloud
    .and(qa.bitwiseAnd(1 << 2).eq(0))      // cirrus
    .and(qa.bitwiseAnd(1 << 3).eq(0))      // cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))      // cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0));     // snow/ice

  if (MASK_WATER) {
    mask = mask.and(qa.bitwiseAnd(1 << 7).eq(0));  // water
  }

  var radMask = image.select('QA_RADSAT').eq(0);

  return image.updateMask(mask).updateMask(radMask);
}

function prepL57(image) {
  image = maskSR(scaleSR(image));
  return image.normalizedDifference(['SR_B4', 'SR_B3'])
    .rename('NDVI')
    .copyProperties(image, ['system:time_start']);
}

function prepL89(image) {
  image = maskSR(scaleSR(image));
  return image.normalizedDifference(['SR_B5', 'SR_B4'])
    .rename('NDVI')
    .copyProperties(image, ['system:time_start']);
}

function toGrid(image) {
  return image
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 65535
    })
    .reproject({
      crs: CRS,
      scale: SCALE
    })
    .clip(pakistan);
}

function monthComposite(d, col) {
  d = ee.Date(d);

  var ndvi = col
    .filterDate(d, d.advance(1, 'month'))
    .select('NDVI')
    .reduce(ee.Reducer.percentile([75]))
    .rename('NDVI');

  return toGrid(ndvi)
    .set('system:time_start', d.millis())
    .set('year', d.get('year'))
    .set('month', d.get('month'))
    .set('date_str', d.format('YYYY-MM'));
}

// ==============================
// SHORT GAP FILLING
// ==============================
function fillAt(i, list, n) {
  i = ee.Number(i);
  var img = ee.Image(list.get(i));

  var edge = i.lt(MAX_GAP).or(i.gte(ee.Number(n).subtract(MAX_GAP)));

  return ee.Image(ee.Algorithms.If(edge, img, (function () {
    var p1 = ee.Image(list.get(i.subtract(1)));
    var p2 = ee.Image(list.get(i.subtract(2)));
    var n1 = ee.Image(list.get(i.add(1)));
    var n2 = ee.Image(list.get(i.add(2)));

    var prev = p1.select('NDVI').unmask(p2.select('NDVI'));
    var next = n1.select('NDVI').unmask(n2.select('NDVI'));

    var tPrev = ee.Image.constant(ee.Number(p1.get('system:time_start'))).toDouble()
      .updateMask(p1.select('NDVI').mask())
      .unmask(
        ee.Image.constant(ee.Number(p2.get('system:time_start'))).toDouble()
          .updateMask(p2.select('NDVI').mask())
      );

    var tNext = ee.Image.constant(ee.Number(n1.get('system:time_start'))).toDouble()
      .updateMask(n1.select('NDVI').mask())
      .unmask(
        ee.Image.constant(ee.Number(n2.get('system:time_start'))).toDouble()
          .updateMask(n2.select('NDVI').mask())
      );

    var tNow = ee.Image.constant(ee.Number(img.get('system:time_start'))).toDouble();

    var interp = prev.add(
      next.subtract(prev).multiply(
        tNow.subtract(tPrev).divide(tNext.subtract(tPrev))
      )
    );

    var fillMask = img.select('NDVI').mask().not()
      .and(prev.mask())
      .and(next.mask());

    return img.select('NDVI')
      .where(fillMask, interp)
      .rename('NDVI')
      .copyProperties(img, ['system:time_start', 'year', 'month', 'date_str']);
  })()));
}

function fillShortGaps(col) {
  var list = col.toList(col.size());
  var n = col.size();

  var out = ee.List.sequence(0, n.subtract(1)).map(function(i) {
    return fillAt(i, list, n);
  });

  return ee.ImageCollection.fromImages(out).sort('system:time_start');
}

// ==============================
// EXPORT HELPERS
// ==============================
function yearStack(col, year) {
  year = ee.Number(year);

  var ycol = col
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .sort('system:time_start');

  var stack = ycol.toBands().toFloat();

  var names = ee.List.sequence(1, 12).map(function(m) {
    return ee.String('NDVI_')
      .cat(year.format('%04d'))
      .cat('_')
      .cat(ee.Number(m).format('%02d'));
  });

  return stack.rename(names).clip(pakistan);
}

function exportYear(col, year) {
  var yearStr = ee.Number(year).format('%04d').getInfo();

  Export.image.toDrive({
    image: yearStack(col, year),
    description: EXPORT_PREFIX + yearStr,
    folder: EXPORT_FOLDER,
    fileNamePrefix: EXPORT_PREFIX + yearStr,
    region: pakistan,
    crs: CRS,
    scale: SCALE,
    maxPixels: 1e13
  });
}

// ==============================
// LOAD LANDSAT COLLECTIONS
// ==============================
var startDate = ee.Date.fromYMD(START_YEAR, 1, 1);
var endDate = ee.Date.fromYMD(END_YEAR + 1, 1, 1);

var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(pakistan)
  .filterDate('1985-01-01', '2013-06-01')
  .map(prepL57);

var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterBounds(pakistan)
  .filterDate('1999-01-01', endDate)
  .map(prepL57);

var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(pakistan)
  .filterDate('2013-01-01', endDate)
  .map(prepL89);

var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
  .filterBounds(pakistan)
  .filterDate('2021-01-01', endDate)
  .map(prepL89);

var landsatNDVI = l5.merge(l7).merge(l8).merge(l9).sort('system:time_start');

print('Merged Landsat NDVI collection:', landsatNDVI.limit(5));

// ==============================
// BUILD MONTHLY NDVI SERIES
// ==============================
var nMonths = endDate.difference(startDate, 'month');

var monthStarts = ee.List.sequence(0, nMonths.subtract(1)).map(function(m) {
  return startDate.advance(m, 'month');
});

var monthlyNDVI = ee.ImageCollection.fromImages(
  monthStarts.map(function(d) {
    return monthComposite(d, landsatNDVI);
  })
).sort('system:time_start');

print('Monthly NDVI collection (raw):', monthlyNDVI.limit(3));

var monthlyOut = FILL_SHORT_GAPS ? fillShortGaps(monthlyNDVI) : monthlyNDVI;

print('Monthly NDVI collection (final):', monthlyOut.limit(3));

// ==============================
// QUICK VISUAL CHECK
// ==============================
Map.addLayer(
  ee.Image(monthlyOut.first()),
  {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green']},
  'First month NDVI'
);

// ==============================
// EXPORT ONE 12-BAND STACK PER YEAR
// ==============================
ee.List.sequence(START_YEAR, END_YEAR).getInfo().forEach(function(year) {
  exportYear(monthlyOut, year);
});
