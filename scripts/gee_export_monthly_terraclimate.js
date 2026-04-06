/*
Purpose
- Exports monthly TerraClimate variables for Pakistan
- Keeps the same boundary placeholder and analysis grid logic as the Landsat script
- Exports one multiband GeoTIFF per year
- Variables exported:
    tmmx  = monthly maximum temperature
    pr    = precipitation
    vpd   = vapour pressure deficit
    soil  = soil moisture
    pet   = potential evapotranspiration
    aet   = actual evapotranspiration
    wb    = climatic water balance (pr - pet)

Important note
Upload your official Pakistan boundary shapefile as a GEE Asset first,
then replace CUSTOM_BOUNDARY_ASSET below with your own asset path.
*/

var START_YEAR = 1985;
var END_YEAR = 2025;

// ==============================
// USER SETTINGS
// ==============================
var CRS = 'EPSG:4326';
var SCALE = 10000;   // keep consistent with the Landsat export grid

var EXPORT_FOLDER = 'Pakistan_TerraClimate_Monthly';
var EXPORT_PREFIX = 'pakistan_terraclimate_monthly_';

// Set to true to use your uploaded official boundary asset
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
// VARIABLE SETTINGS
// ==============================
var OUT_BANDS = ['tmmx', 'pr', 'vpd', 'soil', 'pet', 'aet', 'wb'];

/*
TerraClimate scaling
- tmmx is commonly stored with scale factor 0.1 -> degrees C
- vpd is commonly stored with scale factor 0.01 -> kPa
- pr, soil, pet, aet are treated here in mm/month
- wb is computed as pr - pet in mm/month

If your preferred workflow uses native stored values instead of scaled units,
adjust the scaling below accordingly.
*/

// ==============================
// HELPER FUNCTIONS
// ==============================
function scaleTerraClimate(image) {
  var tmmx = image.select('tmmx').multiply(0.1).rename('tmmx');
  var pr   = image.select('pr').rename('pr');
  var vpd  = image.select('vpd').multiply(0.01).rename('vpd');
  var soil = image.select('soil').rename('soil');
  var pet  = image.select('pet').rename('pet');
  var aet  = image.select('aet').rename('aet');

  var wb = pr.subtract(pet).rename('wb');

  return ee.Image.cat([tmmx, pr, vpd, soil, pet, aet, wb])
    .copyProperties(image, ['system:time_start']);
}

function toGrid(image) {
  return image
    .resample('bilinear')
    .reproject({
      crs: CRS,
      scale: SCALE
    })
    .clip(pakistan);
}

function prepareMonthlyImage(image) {
  var d = ee.Date(image.get('system:time_start'));

  return toGrid(scaleTerraClimate(image))
    .set('system:time_start', d.millis())
    .set('year', d.get('year'))
    .set('month', d.get('month'))
    .set('date_str', d.format('YYYY-MM'));
}

function makeYearBandNames(year) {
  year = ee.Number(year);

  var months = ee.List.sequence(1, 12);

  var bandNames = months.map(function(m) {
    m = ee.Number(m);
    var mm = m.format('%02d');
    var yyyy = year.format('%04d');

    return ee.List([
      ee.String('tmmx_').cat(yyyy).cat('_').cat(mm),
      ee.String('pr_').cat(yyyy).cat('_').cat(mm),
      ee.String('vpd_').cat(yyyy).cat('_').cat(mm),
      ee.String('soil_').cat(yyyy).cat('_').cat(mm),
      ee.String('pet_').cat(yyyy).cat('_').cat(mm),
      ee.String('aet_').cat(yyyy).cat('_').cat(mm),
      ee.String('wb_').cat(yyyy).cat('_').cat(mm)
    ]);
  }).flatten();

  return bandNames;
}

function yearStack(col, year) {
  year = ee.Number(year);

  var ycol = col
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .sort('system:time_start')
    .select(OUT_BANDS);

  var stack = ycol.toBands().toFloat();
  var bandNames = makeYearBandNames(year);

  return stack.rename(bandNames).clip(pakistan);
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
// LOAD TERRACLIMATE
// ==============================
var startDate = ee.Date.fromYMD(START_YEAR, 1, 1);
var endDate = ee.Date.fromYMD(END_YEAR + 1, 1, 1);

var terraclimate = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
  .filterBounds(pakistan)
  .filterDate(startDate, endDate)
  .map(prepareMonthlyImage)
  .sort('system:time_start');

print('Prepared TerraClimate collection:', terraclimate.limit(3));

// ==============================
// QUICK VISUAL CHECK
// ==============================
var firstImage = ee.Image(terraclimate.first());

Map.addLayer(
  firstImage.select('tmmx'),
  {min: 10, max: 40, palette: ['blue', 'cyan', 'yellow', 'red']},
  'First month tmmx'
);

Map.addLayer(
  firstImage.select('pr'),
  {min: 0, max: 200, palette: ['white', 'lightblue', 'blue', 'darkblue']},
  'First month precipitation'
);

// ==============================
// OPTIONAL METADATA TABLE EXPORT
// ==============================
var monthlyMetadata = ee.FeatureCollection(terraclimate.map(function(img) {
  return ee.Feature(null, {
    date_str: img.get('date_str'),
    year: img.get('year'),
    month: img.get('month')
  });
}));

Export.table.toDrive({
  collection: monthlyMetadata,
  description: 'pakistan_terraclimate_monthly_metadata',
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'pakistan_terraclimate_monthly_metadata',
  fileFormat: 'CSV'
});

// ==============================
// EXPORT ONE YEAR AT A TIME
// Each year will contain 84 bands:
// 7 variables × 12 months
// ==============================
ee.List.sequence(START_YEAR, END_YEAR).getInfo().forEach(function(year) {
  exportYear(terraclimate, year);
});
