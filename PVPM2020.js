// State scale analysis of RICE vpm model
// Loading the data
// MOD9A1 data has been loaded before as modTerra
// Images for the study site arkansasRice2020

//Load the modis images
//Growing season of 2020
var MOD09A1Collection = ee.ImageCollection("MODIS/006/MOD09A1").filterDate("2020-01-01", "2020-12-31");

// Function to extract bitwise
function bitwiseExtract(value, fromBit, toBit) {
  if (toBit === undefined) toBit = fromBit;
  var maskSize = ee.Number(1).add(toBit).subtract(fromBit);
  var mask = ee.Number(1).leftShift(maskSize).subtract(1);
  return value.rightShift(fromBit).bitwiseAnd(mask);
}

// Masking function
// Function to filter the images having cloud, cloud shadow, aerosol;
var maskMOD09A1Clouds = function (image) {
  var qa = image.select('StateQA');
  var cloudState = bitwiseExtract(qa, 0, 1); 
  // var cloudShadowState = bitwiseExtract(qa, 2);
  // var cirrusState = bitwiseExtract(qa, 8, 9);
  var aerosolQuantity = bitwiseExtract(qa, 6, 7); 
  var mask = cloudState.eq(0) // Clear
    // .and(cloudShadowState.eq(0)) // No cloud shadow
    // .and(cirrusState.eq(0)) // No cirrus
    .and(aerosolQuantity.lte(1)); // No aerosol quantity
  var maskedImage = image.updateMask(mask);
  return maskedImage; 
}

// filter and cloud-mask image collection 
var maskmod = MOD09A1Collection.map(maskMOD09A1Clouds)

// scaling converting the modis band by the required scaling
var addscaleb1 = function(image) {
  var scaleb1= image.expression('float(b("sur_refl_b01")/10000)').rename('scaleb1');
  return image.addBands(scaleb1)
};

var modisb1 = maskmod.map(addscaleb1);


var addscaleb2 = function(image) {
  var scaleb2= image.expression('float(b("sur_refl_b02")/10000)').rename('scaleb2');
  return image.addBands(scaleb2)
};

var modisb2 = modisb1.map(addscaleb2);


var addscaleb3 = function(image) {
  var scaleb3= image.expression('float(b("sur_refl_b03")/10000)').rename('scaleb3');
  return image.addBands(scaleb3)
};
var modisb3 = modisb2.map(addscaleb3);

    
// NDVI Formula 
var addNDVI = function(image) {
  var ndvi = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).rename('NDVI');
  return image.addBands(ndvi);
};

var withNDVI = modisb3.map(addNDVI);


// EVI
var addEVI = function(image) {
  var evi= image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'RED': image.select("scaleb1"),
      'NIR': image.select("scaleb2"),
      'BLUE': image.select('scaleb3')}).rename('EVI');
  return image.addBands(evi)
};

var withevi = withNDVI.map(addEVI);


// LSWI
var addLSWI = function(image) {
  var lswi = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b06']).rename('LSWI');
  return image.addBands(lswi);
};

var withLSWI = withevi.map(addLSWI);

// 
//max LSWI of whole season
var maxLSWI = withLSWI.reduce(ee.Reducer.max());

// var vis_param = {bands: ['LSWImax'], gamma: 1.6};
// Map.addLayer(maxLSWI, vis_param);


// //LSWImax
var addLSWImax = function(image) {
  var lswimax = maxLSWI.select("LSWI_max").rename('LSWImax');
  return image.addBands(lswimax);
};

var final_col = withLSWI.map(addLSWImax);
// print(withLSWImax, "withLSWImax")
// Map.addLayer(withLSWImax, null, "withLSWImax")


// //LSWImax collection or mapping function only applies for one image and all image have same values.So we merged two collections. Since LSWImax was constant across different image collection
// ////////////////
// // Merging LSWImax with the other collection
// var mod1 = withLSWI
// var mod2 = withLSWImax.select('LSWImax')
// // Use an equals filter to define how the collections match.
// var filter = ee.Filter.equals({
//   leftField: 'system:index',
//   rightField: 'system:index'
// });

// // Create the join.
// var simpleJoin = ee.Join.simple();

// // Apply join
// var mod1join = ee.ImageCollection(simpleJoin.apply(mod1, mod2, filter))
// var mod2join = ee.ImageCollection(simpleJoin.apply(mod2, mod1, filter))

// print('Joined', mod1join, mod2join)

// var final_col = mod1join.map(function(img){

//   // Create a collection with 1 image
//   var temp = ee.ImageCollection(ee.List([img]));

//   // Apply join to collection 2
//   // Resulting collection will have 1 image with exact same date as img
//   var join = simpleJoin.apply(mod2join, temp, filter);

//   // Get resulting image
//   var i2 = ee.Image(join.first())

//   return img.addBands(i2)
// })

// print('final_col', final_col)
// Map.addLayer(final_col, null, "final_col")

// adding FAPAR

// //FAPAR function
var addFAPAR = function(image) {
  var evi = image.select("EVI")
  var fapar = (evi.subtract(0.1)).multiply(1.25).rename('FAPAR')
  return image.addBands(fapar);
};

var withFAPAR = final_col.map(addFAPAR);


// ////////////////////////////////////
// GAFILLING EVI
var fapar1 = ee.Image(withFAPAR.first())
Map.addLayer(fapar1, {bands: ['FAPAR'],
    }, "fapar1")

var mod_evi = withFAPAR.select(['EVI'])
// print(mod_evi)


//fill the gaps in the original MODIS NDVI time series by linear interpolation
var interl_m = require('users/Yang_Chen/GF-SG:Interpolation_v1');
var frame  = 8*4; 
var nodata = -9999; 
var mod_ndvi_interp = interl_m.linearInterp(mod_evi, frame, nodata);
// print("mod_ndvi_interp", mod_ndvi_interp)
var mod_ndvi_interp0 = mod_ndvi_interp.select(['MOD_NDVI_INTER']);
// print("mod_ndvi_interp0", mod_ndvi_interp0)


//Smooth the MODIS interpolation time series by the Savitzky–Golay filter
var sg_filter = require('users/Yang_Chen/GF-SG:SG_filter_v1');
//Reduce the residual noise in the synthesized NDVI time series by the weighted SG filter
var list_trend_sgCoeff = ee.List([-0.070588261,-0.011764720,0.038009040,0.078733027,0.11040724,0.13303168,0.14660634,
0.15113123,0.14660634,0.13303168,0.11040724,0.078733027,0.038009040,-0.011764720,-0.070588261]);   //suggested trend parameter:(7,7,0,2)
var list_sg_coeff = ee.List([0.034965038,-0.12820521,0.069930017,0.31468537,0.41724950,0.31468537,
0.069930017,-0.12820521,0.034965038]);   //suggested parameters of SG:(4,4,0,5)
var syn_series_sg = sg_filter.sg_filter_chen(mod_ndvi_interp0,list_trend_sgCoeff,list_sg_coeff);

// print("syn_series_sg", syn_series_sg)

// Map.addLayer(mod_evi, {bands: ['EVI'],
//     }, "mod_evi")
// Map.addLayer(syn_series_sg, {bands: ['MOD_NDVI_SG'],
//     }, "MOD_NDVI_SG")

function renameBandsETM(image) {
    var bands = ['MOD_NDVI_SG'];
    var new_bands = ['EVI_SG'];
    return image.select(bands).rename(new_bands);
}
var syn_series_sg_evirenamed = syn_series_sg
  .map(renameBandsETM)
// print(syn_series_sg_evirenamed)


// Gapfilled EVI is in another image collection
// join withfapar and gapfilledEVI 

var mod1 = withFAPAR
var mod2 = syn_series_sg_evirenamed.select('EVI_SG')
// Use an equals filter to define how the collections match.
var filter = ee.Filter.equals({
  leftField: 'system:index',
  rightField: 'system:index'
});

// Create the join.
var simpleJoin = ee.Join.simple();

// Applt join
var mod1join = ee.ImageCollection(simpleJoin.apply(mod1, mod2, filter))
var mod2join = ee.ImageCollection(simpleJoin.apply(mod2, mod1, filter))

// print('Joined', mod1join, mod2join)

var final_col = mod1join.map(function(img){

  // Create a collection with 1 image
  var temp = ee.ImageCollection(ee.List([img]));

  // Apply join to collection 2
  // Resulting collection will have 1 image with exact same date as img
  var join = simpleJoin.apply(mod2join, temp, filter);

  // Get resulting image
  var i2 = ee.Image(join.first())

  return img.addBands(i2)
})

// print('final_col', final_col)
///

// gapfill lswi data
var mod_lswi = final_col.select(['LSWI'])


//fill the gaps in the original MODIS NDVI time series by linear interpolation
var interl_m = require('users/Yang_Chen/GF-SG:Interpolation_v1');
var frame  = 8*7; 
var nodata = -9999; 
var mod_ndvi_interp = interl_m.linearInterp(mod_lswi, frame, nodata);
var mod_ndvi_interp0 = mod_ndvi_interp.select(['MOD_NDVI_INTER']);

//Smooth the MODIS interpolation time series by the Savitzky–Golay filter
var sg_filter = require('users/Yang_Chen/GF-SG:SG_filter_v1');
//Reduce the residual noise in the synthesized NDVI time series by the weighted SG filter
var list_trend_sgCoeff = ee.List([-0.070588261,-0.011764720,0.038009040,0.078733027,0.11040724,0.13303168,0.14660634,
0.15113123,0.14660634,0.13303168,0.11040724,0.078733027,0.038009040,-0.011764720,-0.070588261]);   //suggested trend parameter:(7,7,0,2)
var list_sg_coeff = ee.List([0.034965038,-0.12820521,0.069930017,0.31468537,0.41724950,0.31468537,
0.069930017,-0.12820521,0.034965038]);   //suggested parameters of SG:(4,4,0,5)
var syn_series_sg = sg_filter.sg_filter_chen(mod_ndvi_interp0,list_trend_sgCoeff,list_sg_coeff);

// print("syn_series_sg", syn_series_sg)

// Map.addLayer(mod_lswi, {bands: ['LSWI'],
//     }, "mod_lswi")
// Map.addLayer(syn_series_sg, {bands: ['MOD_NDVI_SG'],
//     }, "MOD_NDVI_SG_LSWI")

function renameBandsETMLSWI(image) {
    var bands = ['MOD_NDVI_SG'];
    var new_bands = ['LSWI_SG'];
    return image.select(bands).rename(new_bands);
}
var syn_series_sg_lswirenamed = syn_series_sg
  .map(renameBandsETMLSWI)
// print("syn_series_sg_lswirenamed", syn_series_sg_lswirenamed)

// adding LSWIsg to the collection
// join with lswimax (syn_series_sg_lswirenamed) and 

var mod1 = final_col
var mod2 = syn_series_sg_lswirenamed.select('LSWI_SG')
// Use an equals filter to define how the collections match.
var filter = ee.Filter.equals({
  leftField: 'system:index',
  rightField: 'system:index'
});

// Create the join.
var simpleJoin = ee.Join.simple();

// Applt join
var mod1join = ee.ImageCollection(simpleJoin.apply(mod1, mod2, filter))
var mod2join = ee.ImageCollection(simpleJoin.apply(mod2, mod1, filter))

// print('Joined', mod1join, mod2join)

var final_col = mod1join.map(function(img){

  // Create a collection with 1 image
  var temp = ee.ImageCollection(ee.List([img]));

  // Apply join to collection 2
  // Resulting collection will have 1 image with exact same date as img
  var join = simpleJoin.apply(mod2join, temp, filter);

  // Get resulting image
  var i2 = ee.Image(join.first())

  return img.addBands(i2)
})

// print('final_col', final_col)


///
// Calculate FAPAR SG
var addFAPARsg = function(image) {
  var evi = image.select("EVI_SG")
  var faparsg = (evi.subtract(0.1)).multiply(1.25).rename('FAPAR_sg')
  return image.addBands(faparsg);
};

var final_colwithFAPAR = final_col.map(addFAPARsg);
// print("final_colwithFAPAR", final_colwithFAPAR)


//Calculate WS

var addWs = function(image) {
  var lswi = image.select('LSWI_SG');
  var lswimax = image.select('LSWImax');
  var ws = lswi.add(1).divide(lswimax.add(1)).rename('Ws')
return image.addBands(ws);
};

var finalCol_withWs = final_colwithFAPAR.map(addWs);
// print(finalCol_withWs, "finalCol_withWs")

var finalCol_withWs_WS = finalCol_withWs.select(['LSWI_SG', 'LSWImax', 'Ws'])
// Map.addLayer(finalCol_withWs_WS, {bands: ['LSWI_SG', 'LSWImax', 'Ws']}, "finalCol_withWs_WS")


// map add layer
// Map.addLayer(arkansasRice2020) 
var palettes = require('users/gena/packages:palettes');

// Mosaic the visualization layers and display (or export).
var image = ee.Image(final_colwithFAPAR.first());
var imageRGB = image.visualize({bands: ['EVI_SG'], 
  max: 1,
  palette: palettes.misc.tol_rainbow[7]
});

var mosaic = ee.ImageCollection([imageRGB]).mosaic();
// Map.addLayer(mosaic, {}, 'mosaic');

// Display a clipped version of the mosaic.
// Map.addLayer(mosaic.clip(arkansasRice2020));

var addLUE = function(image) {
  var LUE = (ee.Image.constant(0.42)).rename('LUE')
  return image.addBands(LUE);
};

var withLUE = finalCol_withWs.map(addLUE);
// print("withLUE", withLUE)

// Map.addLayer(withLUE) 

var finalCol_withWs_WS = finalCol_withWs.select(['LSWI_SG', 'LSWImax', 'Ws'])


// NCEP data
var startDate = ee.Date('2020-01-01')
var endDate = ee.Date('2020-12-31')
var collection = ee.ImageCollection("NOAA/CFSV2/FOR6H")
var collection = collection.filterDate('2020-01-01', '2020-12-31')

// Renaming the bands
function renameBandsETMdswr(image) {
    var bands = ['Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average', 'Temperature_height_above_ground', 'Maximum_temperature_height_above_ground_6_Hour_Interval'];
    var new_bands = ['DSWR', 'Temp', 'Tmax'];
    return image.select(bands).rename(new_bands);
}
var collection = collection.map(renameBandsETMdswr)
// print(collection)



// 6 hourly to daily
var numberOfDays = endDate.difference(startDate, 'days')
var daily1 = ee.ImageCollection(
  ee.List.sequence(0, numberOfDays.subtract(1))
    .map(function (dayOffset) {
      var start = startDate.advance(dayOffset, 'days')
      var end = start.advance(1, 'days')
      return collection
        .filterDate(start, end)
        .mean()
        .set('system:time_start', start.millis());
    })
);

var eightday = daily1.filterDate('2020-01-01', "2020-12-31")
// print("eightday", eightday)


var startDate = ee.Date('2020-01-01')
var endDate = ee.Date('2020-12-31')
var dayOffsets = ee.List.sequence(
  0, 
  endDate.difference(startDate, 'days').subtract(1),
  8 // Single day every week
)

var weeklyMeans = ee.ImageCollection.fromImages(
  dayOffsets.map(function(dayOffset) {
    var start = startDate.advance(dayOffset, 'day')
    var end = start.advance(8, 'day')
    return eightday
      .filterDate(start, end)
      .mean()
      .set('system:time_start', start.millis());
  })  
);

// print('weeklyMeans', weeklyMeans)
var TempSWR = weeklyMeans.select(['Tmax', 'Temp', 'DSWR'])



///Combine two image collections
var mod1 = withLUE
var mod2 = TempSWR

var filter = ee.Filter.equals({
  leftField: 'system:time_start',
  rightField: 'system:time_start'
});

// Create the join.
var simpleJoin = ee.Join.inner();

// Inner join
var innerJoin = ee.ImageCollection(simpleJoin.apply(mod1, mod2, filter))

var joined = innerJoin.map(function(feature) {
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'));
})

// print('Joined', joined)

var joined_swr = joined.select(['DSWR'])
// Map.addLayer(joined_swr, {bands: ['DSWR']}, "joined_swr")


// Calculate PAR
// Calculate FAPAR SG
var addpar = function(image) {
  var dswr = image.select("DSWR")
  var dswr_par = dswr.multiply(0.9)
  var par = dswr_par.multiply(2.02).multiply(0.0864).rename('par')
  return image.addBands(par);
};

var withpar = joined.map(addpar);
// print("withpar", withpar)


var joined_par = withpar.select(['par'])
// Map.addLayer(joined_par, {bands: ['par']}, "joined_par")

// Temp Deg to celsius
var tmaxcelsius = function(image) {
  var tmax = image.select("Tmax")
  var tmaxcel = tmax.subtract(273.15).rename('tmaxcel')
  return image.addBands(tmaxcel);
};

var withcel = withpar.map(tmaxcelsius);

var tempcelsius = function(image) {
  var temp = image.select("Temp")
  var tempcel = temp.subtract(273.15).rename('tempcel')
  return image.addBands(tempcel);
};

var withcel = withcel.map(tempcelsius);
//print("withcel", withcel)

// mean temperature

var tmean = function(image) {
  var tmax = image.select("tmaxcel")
  var temp = image.select("tempcel")
  var tmean = tmax.add(temp).divide(2). rename("tmean")
  return image.addBands(tmean);
};

var withcel = withcel.map(tmean);
//print(withcel)

var withcel_temp = withcel.select(['tmaxcel', 'tempcel', 'tmean'])
//Map.addLayer(withcel_temp, {bands: ['tmaxcel', 'tempcel', 'tmean']}, "withcel_temp")


// Ts calculation

var Ts = function(image) {
  var tmean = image.select("tmean")
  var ts = (tmean.add(1).multiply(tmean.subtract(48))).divide((tmean.add(1).multiply(tmean.subtract(48))).subtract((tmean.subtract(30)).pow(2))).rename("ts")
  return image.addBands(ts);
};

var withcel = withcel.map(Ts);
//print(withcel)


var withcel_ts = withcel.select(['ts', 'tmean'])
//Map.addLayer(withcel_ts, {bands: ['ts', 'tmean']}, "withcel_ts")

var withcel_fapar = withcel.select(['FAPAR_sg'])
//Map.addLayer(withcel_fapar, {bands: ['FAPAR_sg']}, "withcel_fapar")

var withcel_LUE = withcel.select(['LUE'])
//Map.addLayer(withcel_LUE, {bands: ['LUE']}, "withcel_LUE")

// GPPVPM calculation

var VPM = function(image) {
  var par = image.select("par")
  var ts = image.select("ts")
  var Ws = image.select("Ws")
  var FAPAR_sg = image.select("FAPAR_sg")
  var LUE = image.select("LUE")
  var gpp = (par.multiply(ts).multiply(Ws).multiply(FAPAR_sg).multiply(LUE)).rename("gpp")
  return image.addBands(gpp);
};

var withcel = withcel.map(VPM);
//print("withcelVPM",withcel)


var withcel_gpp = withcel.select(['gpp'])
//Map.addLayer(withcel_gpp, {bands: ['gpp']}, "withcel_gpp")

// Image collection reduction
// Compute the median in each band, each pixel.
// Band names are B1_median, B2_median, etc.
var mean = withcel_gpp.reduce(ee.Reducer.mean());

// print(mean)




// PVPM calculation starts here
var bandSubset = withcel.select(['EVI_SG'], ['EVI_SG'])
// add doy
var addDate = function(image){
  var doy = image.date().getRelative('day', 'year');
  var doyBand = ee.Image.constant(doy).uint16().rename('doy')
  return image.addBands(doyBand);
};

var floatcollection = function(image){
  image.float();
  return image;
};

var bandSubset = bandSubset.map(addDate)
var bandSubset = bandSubset.map(floatcollection)
//print(bandSubset, "bandSubset")

// Define time range
// FULL TS
var startyear = 2020;
var endyear = 2020;

var startmonth = 1  
var endmonth = 12

var startday = 1
var endday = 31

var startdate = ee.Date.fromYMD(startyear,startmonth,startday);
var enddate = ee.Date.fromYMD(endyear,endmonth,endday);

// create list for years
var years = ee.List.sequence(startyear,endyear);
// create list for months
var months = ee.List.sequence(1,12);
// create list for pentads for each month
var pentads = ee.List.sequence(1,26,5,6);


var maxevi = bandSubset.reduce(ee.Reducer.max(bandSubset.first().bandNames().size()))
//print(maxevi, "maxevi")
//print(bandSubset, "bandSubset")

var max = ee.ImageCollection(
    years.map(function (y) {
      var start = ee.Date.fromYMD(y, startmonth,startday);
      var stop = ee.Date.fromYMD (y,endmonth,endday);
      var x = bandSubset.filterDate(start, stop)
      var w = x.qualityMosaic('EVI_SG').select('doy', 'EVI_SG')
    return w.set({'year': y})
}).flatten());
//print(max,'max')
//Map.addLayer(max, {min: 1, max: 365}, 'max')


var bandSubset = bandSubset.select(['EVI_SG', "doy"], ['EVI_SG', "doy"])


// Before
// Find the day of min before day of max for each year
var min_before = ee.ImageCollection(
    years.map(function (y) {
      var start = ee.Date.fromYMD(y, startmonth,startday);
      var stop = ee.Date.fromYMD (y,endmonth,endday);
      var x = bandSubset.filterDate(start, stop)
      var z = ee.Image(max
        .filterMetadata('year', 'equals', y).first());
      var w = x.map(
        function(img) {
        var date = img.select('doy')
        var k = img.updateMask(z.gt(date))
      return k
}).select('EVI_SG', 'doy').reduce(ee.Reducer.min(2)).rename('EVI_SG_min_before','doy')
  return w
}).flatten());
//print(min_before,'min_before')
//Map.addLayer(min_before, {min: 1, max: 365}, 'min_before')






// After
// Find the day of min after of max for each year
var min_after = ee.ImageCollection(
    years.map(function (y) {
      var start = ee.Date.fromYMD(y, startmonth,startday);
      var stop = ee.Date.fromYMD (y,endmonth,endday);
      var x = bandSubset.filterDate(start, stop)
      var z = ee.Image(max
        .filterMetadata('year', 'equals', y).first());
      var w = x.map(
        function(img) {
        var date = img.select('doy')
        var k = img.updateMask(z.lt(date))
      return k
}).select('EVI_SG', 'doy').reduce(ee.Reducer.min(2)).rename('EVI_SG_min_after','doy')
  return w
}).flatten());
//print(min_after,'min_after')
//Map.addLayer(min_after, {min: 1, max: 365}, 'min_after')



Map.addLayer(arkansasRice2020)


var bandSubsetonlyEVI = bandSubset.select(['EVI_SG'], ['EVI_SG'])
//Map.addLayer(bandSubsetonlyEVI, null, 'bandSubsetonlyEVI')

// renaming the max collection evi to max evi
var max = max.select(['EVI_SG', 'doy'],['EVI_SG_max', "doy_max"])
//print(max, "max")
//print(eviwithdop, "eviwithdop")

// converting the collection to image
var maximg = ee.Image(max.first());
//print("maximg", maximg)
//

//  add max minbefore and min after with evidop
var adddmaxminpoints = function(image) {
  var maxminpoints = maximg.select("EVI_SG_max", "doy_max");
  return image.addBands(maxminpoints);
};

var eviwithdopmaxmin = withcel.map(adddmaxminpoints);
//print("eviwithdopmaxmin", eviwithdopmaxmin)

//Map.addLayer(eviwithdopmaxmin, null, 'eviwithdopmaxmin')


// min before
// converting the collection to image
var min_beforeimg = ee.Image(min_before.first());
//print("min_beforeimg", min_beforeimg)



//  add max minbefore and min after with evidop
var adddminbeforepoints = function(image) {
  var minbeforepoints = min_beforeimg.select("EVI_SG_min_before");
  return image.addBands(minbeforepoints);
};

var eviwithdopmaxmin = eviwithdopmaxmin.map(adddminbeforepoints);
//print("eviwithdopmaxmin", eviwithdopmaxmin)

//Map.addLayer(eviwithdopmaxmin, null, 'eviwithdopmaxmin')

// min after
// converting the collection to image
var min_afterimg = ee.Image(min_after.first());
//print("min_afterimg", min_afterimg)



//  add max minbefore and min after with evidop
var adddmin_afterpoints = function(image) {
  var minafterpoints = min_afterimg.select("EVI_SG_min_after");
  return image.addBands(minafterpoints);
};

var eviwithdopmaxmin = eviwithdopmaxmin.map(adddmin_afterpoints);
//print("eviwithdopmaxmin", eviwithdopmaxmin)

//Map.addLayer(eviwithdopmaxmin, null, 'eviwithdopmaxmin')







var eviwithdopmaxminimg = ee.Image(eviwithdopmaxmin.first());

var addthresholdparams = function(image) {
  var evipoints = eviwithdopmaxminimg.select("EVI_SG_max", "EVI_SG_min_before", "EVI_SG_min_after", "doy_max");
  return image.addBands(evipoints);
};

var eviwithdop = bandSubset.map(addthresholdparams);
// print(eviwithdop, "eviwithdop")
// Map.addLayer(eviwithdop, null, "eviwithdop")



// calculate g1 and g2
var calculateg1 = function(image) {
  var g1band= image.expression(
    '(EVI_SG_max)-(EVI_SG_min_before)', {
      'EVI_SG_max': image.select("EVI_SG_max"),
      'EVI_SG_min_before': image.select("EVI_SG_min_before")}).rename('g1');
  return image.addBands(g1band)
};
var eviwithdop = eviwithdop.map(calculateg1);
// print(eviwithdop)

var calculateg2 = function(image) {
  var g2band= image.expression(
    '(EVI_SG_max)-(EVI_SG_min_after)', {
      'EVI_SG_max': image.select("EVI_SG_max"),
      'EVI_SG_min_after': image.select("EVI_SG_min_after")}).rename('g2');
  return image.addBands(g2band)
};
var eviwithdop = eviwithdop.map(calculateg2);
//print(eviwithdop)
//Map.addLayer(eviwithdop, null, "eviwithdop")


// Calculate the DOP
var calculateDOP = function(image) {
  var DOP= image.expression(
    '-95+ (EVI_SG_min_after  *139.33)+ (76.2488 * g1) + (0.7153 * doy_max)', {
      'EVI_SG_min_after': image.select("EVI_SG_min_after"),
      'g1': image.select("g1"),
      'doy_max': image.select("doy_max"),
      
    }).rename('DOP');
  return image.addBands(DOP)
};
var eviwithdop = eviwithdop.map(calculateDOP);

// Calculate the DOH
var calculateDOH = function(image) {
  var DOH= image.expression(
    '-17.6151 + (EVI_SG_min_after*110.4337) +(66.0053 *g1) + (1.0242 * doy_max)', {
      'EVI_SG_min_after': image.select("EVI_SG_min_after"),
      'g1': image.select("g1"),
      'doy_max': image.select("doy_max"),
      
    }).rename('DOH');
  return image.addBands(DOH)
};
var eviwithdop = eviwithdop.map(calculateDOH);
//Map.addLayer(eviwithdop, null, "eviwithdop")


// need to reduce the eviwithdop image collection to an image
// Reduce the collection.
var eviwithdopmean = eviwithdop.reduce(ee.Reducer.mean());
//Map.addLayer(eviwithdopmean, null, "eviwithdopmean")


//print(final_col, "finalcolbeforedopdoh")
// add dop and doh to withcel main collection
var adddopdoh = function(image) {
  var adddopdohbands = eviwithdopmean.select("DOP_mean", "DOH_mean");
  return image.addBands(adddopdohbands);
};

var final_col = eviwithdopmaxmin.map(adddopdoh);
//print(final_col, "finalcolafterdopdoh")
var final_col = final_col.map(addDate)

// Calculate DAP
var calculateDAP = function(image) {
  var DAP= image.expression(
    'DayOfYear - DOP_mean', {
      'DayOfYear': image.select("doy"),
      'DOP_mean': image.select("DOP_mean"),
      
    }).rename('DAP');
  return image.addBands(DAP)
};
var final_col = final_col.map(calculateDAP);
//Map.addLayer(final_col, null, "final_col")

//Modeled LUEMax

// subtract Dop from DOY
var addDAP = function(image) {

  var DAP = ee.Image().expression(
    "((i.doy - i.DOP_mean) < 0.0) ? 0.0"+
    ":i.doy - i.DOP_mean", { i: image}
  ).rename('DAP_1')
  
  return image.addBands(DAP)  
 }

var withDAPfromDOP = final_col.map(addDAP)
//print(withDAPfromDOP, "withDAPfromDOP")
//Map.addLayer( withDAPfromDOP)


// subtract Dop from DOY
var addDAP2 = function(image) {

  var DAP = ee.Image().expression(
    "((i.DOH_mean - i.doy) < 0.0) ? 0.0"+
    ":i.doy - i.DOP_mean", { i: image}
  ).rename('DAP_2')
  
  return image.addBands(DAP)  
 }
var withDAPfromDOP = withDAPfromDOP.map(addDAP2)
 
// subtract Dop from DOY
var addDAP3 = function(image) {

  var DAP3 = ee.Image().expression(
    "((i.DAP_2) < 0.0) ? 0.0"+
    ":i.DAP_2", { i: image}
  ).rename('DAP_3')
  
  return image.addBands(DAP3)  
 }

var withDAPfromDOP = withDAPfromDOP.map(addDAP3)
//print(withDAPfromDOP, "withDAPfromDOP")
//Map.addLayer( withDAPfromDOP)

// y_pred_modarrhenius<-(0.071*((-3537*exp((-1665*(x-69))/(x*8.14*69)))/(-3537-(-1665*(1-exp((-3537*(x-69))/(x*8.14*69)))))))



var pLUEmax = function(image) {
  var a = ee.Image().expression(
    'i.DAP_3*8.14*69', {i: image}
  )
  var b = ee.Image().expression(
    '-1665*(i.DAP_3-69)', {i: image}
  )
  var m = ee.Image().expression(
    '(b/a)', {i: image, a:a, b:b}
  )
  var c = ee.Image().expression(
    '-3537*(exp(m))', {i: image, m:m}
  )
  var d = ee.Image().expression(
    'i.DAP_3*8.14*69', {i: image}
  )
  var e = ee.Image().expression(
    '-3537*(i.DAP_3-69)', {i: image}
  )
  var f = ee.Image().expression(
    '1-(exp(e/d))', {i: image, d:d, e:e}
  )
  var g = ee.Image().expression(
    '-3537-(-1665*f)', {i: image, f:f}
  )
  var pLUEmax = ee.Image().expression(
    '0.071*(c/g)', {a: a, b: b, c: c, d:d, e:e, f:f, g:g, m:m}
  ).rename('pLUEmax')
  
  return image.addBands(pLUEmax)  
 }

var withpLUEmaxmodarr = withDAPfromDOP.map(pLUEmax)
//print(withpLUEmaxmodarr, 'withpLUEmaxmodarr')
//Map.addLayer( withpLUEmaxmodarr)



// LUEmax
var pLUEmax2 = function(image) {

  var LUEmax2 = ee.Image().expression(
    "((i.DAP_3) == 0.0) ? 0.0"+
    ":i.pLUEmax", { i: image}
  ).rename('pLUEmax2')
  
  return image.addBands(LUEmax2)  
 }

var withDAPfromDOP = withpLUEmaxmodarr.map(pLUEmax2)
//print(withDAPfromDOP, "withDAPfromDOP")
//Map.addLayer( withDAPfromDOP)



// GPPVPM calculation

var PVPMfunction = function(image) {
  var par = image.select("par")
  var ts = image.select("ts")
  var Ws = image.select("Ws")
  var FAPAR_sg = image.select("FAPAR_sg")
  var LUE = image.select("pLUEmax2")
  var gpppvpm = (par.multiply(ts).multiply(Ws).multiply(FAPAR_sg).multiply(LUE).multiply(12.011)).rename("gpppvpm")
  return image.addBands(gpppvpm);
};

var withcel = withDAPfromDOP.map(PVPMfunction);
//print(withcel)


var withcel_gpp = withcel.select(['gpppvpm'])
//Map.addLayer(withcel_gpp, {bands: ['gpppvpm']}, "withcel_gpp")



// Define the chart and print it to the console.
var chart =
    ui.Chart.image
        .series({
          imageCollection: withcel_gpp,
          region: arkansasRice2020,
          reducer: ee.Reducer.mean(),
          scale: 500,
          xProperty: 'system:time_start'
        })
        .setSeriesNames(['gpppvpm'])
        .setOptions({
          title: 'Average Vegetation Index Value by Date for Forest',
          hAxis: {title: 'Date', titleTextStyle: {italic: false, bold: true}},
          vAxis: {
            title: 'Vegetation index (x1e4)',
            titleTextStyle: {italic: false, bold: true}
          },
          lineWidth: 5,
          colors: ['e37d05'],
          curveType: 'function'
        });
print(chart);

var withcel_gpp = withcel.select(['gpppvpm'])


// Image collection reduction
// Compute the median in each band, each pixel.
// Band names are B1_median, B2_median, etc.
var mean = withcel_gpp.reduce(ee.Reducer.mean());
// Export the image, specifying scale and region.
Export.image.toDrive({
  image: mean.clipToCollection(arkansasRice2020),
  description: 'arkansasRice2020PVPM',
  scale: 500,
  region: geometry
});

// code if the minimum code does not work
// var z = ee.Image(max
//         .filterMetadata('year', 'equals', 2020).first());
// print(z, "z")

// var w = bandSubset.map(
//     function(img) {
//     var date = img.select('doy')
//     var k = z.updateMask(img.gt(date))
//     return k})
// print(w, "w")
// Map.addLayer(w, null, "w")

// print(bandSubset.first().date(), "Bandsubset first date")
// print(bandSubset.aggregate_max('system:index')) 





