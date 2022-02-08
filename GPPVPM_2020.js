// State scale analysis of RICE vpm model
// Loading the data
// MOD9A1 data has been loaded before as modTerra
// Images for the study site arkansasRice2020

//Load the modis images
//Growing season of 2020
var MOD09A1Collection = ee.ImageCollection("MODIS/006/MOD09A1").filterDate("2020-04-01", "2020-10-30");

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

var withLSWImax = withLSWI.map(addLSWImax);

// print("withLSWImax", withLSWImax)
// print("maxLSWI", maxLSWI)
// Map.addLayer( withLSWImax)


//LSWImax collection or mapping function only applies for one image and all image have same values.So we merged two collections. Since LSWImax was constant across different image collection
////////////////
// Merging LSWImax with the other collection
var mod1 = withLSWI
var mod2 = withLSWImax.select('LSWImax')
// Use an equals filter to define how the collections match.
var filter = ee.Filter.equals({
  leftField: 'system:index',
  rightField: 'system:index'
});

// Create the join.
var simpleJoin = ee.Join.simple();

// Apply join
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
var frame  = 8*4; 
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
Map.addLayer(arkansasRice2020) 
var palettes = require('users/gena/packages:palettes');

// Mosaic the visualization layers and display (or export).
var image = ee.Image(final_colwithFAPAR.first());
var imageRGB = image.visualize({bands: ['EVI_SG'], 
  max: 1,
  palette: palettes.misc.tol_rainbow[7]
});

var mosaic = ee.ImageCollection([imageRGB]).mosaic();
Map.addLayer(mosaic, {}, 'mosaic');

// Display a clipped version of the mosaic.
Map.addLayer(mosaic.clip(arkansasRice2020));

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

var eightday = daily1.filterDate('2020-05-08', "2020-12-31")
print("eightday", eightday)


var startDate = ee.Date('2020-05-08')
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

print('weeklyMeans', weeklyMeans)
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

print('Joined', joined)

var joined_swr = joined.select(['DSWR'])
Map.addLayer(joined_swr, {bands: ['DSWR']}, "joined_swr")


// Calculate PAR
// Calculate FAPAR SG
var addpar = function(image) {
  var dswr = image.select("DSWR")
  var par = dswr.multiply(2.02).multiply(0.0864).rename('par')
  return image.addBands(par);
};

var withpar = joined.map(addpar);
print("withpar", withpar)


var joined_par = withpar.select(['par'])
Map.addLayer(joined_par, {bands: ['par']}, "joined_par")

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
print("withcel", withcel)

// mean temperature

var tmean = function(image) {
  var tmax = image.select("tmaxcel")
  var temp = image.select("tempcel")
  var tmean = tmax.add(temp).divide(2). rename("tmean")
  return image.addBands(tmean);
};

var withcel = withcel.map(tmean);
print(withcel)

var withcel_temp = withcel.select(['tmaxcel', 'tempcel', 'tmean'])
Map.addLayer(withcel_temp, {bands: ['tmaxcel', 'tempcel', 'tmean']}, "withcel_temp")


// Ts calculation

var Ts = function(image) {
  var tmean = image.select("tmean")
  var ts = (tmean.add(1).multiply(tmean.subtract(48))).divide((tmean.add(1).multiply(tmean.subtract(48))).subtract((tmean.subtract(30)).pow(2))).rename("ts")
  return image.addBands(ts);
};

var withcel = withcel.map(Ts);
print(withcel)


var withcel_ts = withcel.select(['ts', 'tmean'])
Map.addLayer(withcel_ts, {bands: ['ts', 'tmean']}, "withcel_ts")

var withcel_fapar = withcel.select(['FAPAR_sg'])
Map.addLayer(withcel_fapar, {bands: ['FAPAR_sg']}, "withcel_fapar")

var withcel_LUE = withcel.select(['LUE'])
Map.addLayer(withcel_LUE, {bands: ['LUE']}, "withcel_LUE")

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
print(withcel, "withcel")


var withcel_gpp = withcel.select(['gpp'])
Map.addLayer(withcel_gpp, {bands: ['gpp']}, "withcel_gpp")

// Image collection reduction
// Compute the median in each band, each pixel.
// Band names are B1_median, B2_median, etc.
var mean = withcel_gpp.reduce(ee.Reducer.mean());

print(mean)

// The output is an Image.  Add it to the map.
var palettes = require('users/gena/packages:palettes');
var palette = palettes.misc.tol_rainbow[7];
var vis_param = {bands: ['gpp_mean'], min: 0, max: 25, palette: palette};
Map.addLayer(mean.clip(arkansasRice2020), vis_param);

// Reduce the region. The region parameter is the Feature geometry.
// var meanDictionary = mean.reduceRegion({
//   reducer: ee.Reducer.mean(),
//   geometry: arkansasRice2020.geometry(),
//   scale: 500,
//   maxPixels: 1e9
// });

// The result is a Dictionary.  Print it.
//print(meanDictionary);

// Export the image, specifying scale and region.
Export.image.toDrive({
  image: mean.clipToCollection(arkansasRice2020),
  description: 'arkansasRice2020',
  scale: 500,
  region: geometry
});

 // Plot gpp ---------------------------------------------------------------------------------------------
var gpp2020 = ui.Chart.image.seriesByRegion({
  imageCollection: withcel_gpp,
  regions: arkansasRice2020,
  reducer: ee.Reducer.mean(), //type of reduction. See ee.Reducers for other kinds of reductions
  scale: 50000, //spatial scale of MODIS product
  seriesProperty: 'NAME'  //property of roi to display in map
})
  .setOptions({
    title: 'MODIS NDSI',
    vAxis: {title: 'NDSI', maxValue: 15, minValue: -1},
    hAxis: {title: 'date', format: 'MM-yy', gridlines: {count: 12}},
  })

print(gpp2020)


// Define a Point object.

var point = ee.Geometry.Point(-122.082, 35.8105);


var band = 'EVI_SG';



var how = ee.Reducer.mean();
var where = ee.FeatureCollection([ee.Feature(arkansasRice2020)]);

print(ui.Chart.image.doySeriesByRegion(withcel, band, where));



  
// add dapmaxevi = 80

var DAPconstant = function(image) {
  var DAPconstant = (ee.Image.constant(80)).rename('DAPconstant')
  return image.addBands(DAPconstant);
};


var DAPconstant = withcel.map(DAPconstant);
print("DAPconstant", DAPconstant)
  
  
  //adding unmasked DOY
  //create date band (day of  year) for each image in collection
var addDate = function(image){
  var doy = image.date().getRelative('day', 'year');
  var doyBand = ee.Image.constant(doy).uint16().rename('doy')
  doyBand = doyBand.updateMask(image.select('EVI_SG').mask())
  
  return image.addBands(doyBand);
};
  
var DAPconstant = DAPconstant.map(addDate);
print(DAPconstant);  
  
  // band subsetting
  var DAPconstantSubset = DAPconstant.select(['EVI_SG','LSWImax', 'DAPconstant', 'doy'], ['EVI_SG','LSWImax', 'DAPconstant', 'DayOfYear'])
  print(DAPconstantSubset)
  
// Max evi pixel and the DAP = 80
var maxevi = DAPconstantSubset.reduce(ee.Reducer.max(DAPconstantSubset.first().bandNames().size())).rename(['EVImax', 'LSWImax', 'DAP', 'DOY'])
print(maxevi, "maxevi")
Map.addLayer(maxevi) 





// var addDOP = function(image) {
//   var DOP = image.subtract(['DOY', 'DAP']).rename('DOP');
//   return image.addBands(DOP);
// };


// var withDOP = addDOP(maxevi)
// print(withDOP)


var addDOP = function(image) {

  var DOP = ee.Image().expression(
    'i.DOY - i.DAP', {i: image}
  ).rename('DOP')
  
  return image.addBands(DOP)  
 }

var withDOPimage = addDOP(maxevi)
print(withDOPimage)

var withDOP = addDOP(maxevi).select('DOP')
Map.addLayer(withDOP)
 
// merge DOP with the collection

// create DOP collection
// //Add DOP to the
var addDOPcollection = function(image) {
  var adddopone = withDOP.select("DOP").rename('DOP');
  return image.addBands(adddopone);
};

var withonlyDOPcollection = DAPconstant.map(addDOPcollection);

print("withonlyDOPcollection", withonlyDOPcollection)
Map.addLayer( withonlyDOPcollection)
Map.addLayer( DAPconstant)

// subtract Dop from DOY
var addDAP = function(image) {

  var DAP = ee.Image().expression(
    "((i.doy - i.DOP) == 0.0) ? 1.0"+
    ":i.doy - i.DOP", { i: image}
  ).rename('DAP')
  
  return image.addBands(DAP)  
 }

var withDAPfromDOP = withonlyDOPcollection.map(addDAP)
print(withDAPfromDOP, "withDAPfromDOP")
Map.addLayer( withDAPfromDOP)




// y_pred_modarrhenius<-(0.0583*((-3709*exp((-1665*(x-68.05))/(x*8.14*68.05)))/(-3709-(-1665*(1-exp((-3709*(x-68.05))/(x*8.14*68.05)))))))



var pLUEmax = function(image) {
  var a = ee.Image().expression(
    'i.DAP*8.14*68.05', {i: image}
  )
  var b = ee.Image().expression(
    '-1665*(i.DAP-68.05)', {i: image}
  )
  var m = ee.Image().expression(
    '(b/a)', {i: image, a:a, b:b}
  )
  var c = ee.Image().expression(
    '-3709*(exp(m))', {i: image, m:m}
  )
  var d = ee.Image().expression(
    'i.DAP*8.14*68.05', {i: image}
  )
  var e = ee.Image().expression(
    '-3709*(i.DAP-68.05)', {i: image}
  )
  var f = ee.Image().expression(
    '1-(exp(e/d))', {i: image, d:d, e:e}
  )
  var g = ee.Image().expression(
    '-3709-(-1665*f)', {i: image, f:f}
  )
  var pLUEmax = ee.Image().expression(
    '0.055*(c/g)', {a: a, b: b, c: c, d:d, e:e, f:f, g:g, m:m}
  ).rename('pLUEmax')
  
  return image.addBands(pLUEmax)  
 }

var withpLUEmaxmodarr = withDAPfromDOP.map(pLUEmax)
print(withpLUEmaxmodarr, 'withpLUEmaxmodarr')
Map.addLayer( withpLUEmaxmodarr)


// GPP PVPM
var PVPM = function(image) {
  var par = image.select("par")
  var ts = image.select("ts")
  var Ws = image.select("Ws")
  var FAPAR_sg = image.select("FAPAR_sg")
  var pLUE = image.select("pLUEmax")
  var gpp_pvpm = (par.multiply(ts).multiply(Ws).multiply(FAPAR_sg).multiply(pLUE)).multiply(ee.Number(12)).rename("gpp_pvpm")
  return image.addBands(gpp_pvpm);
};

var PVPMgpp = withpLUEmaxmodarr.map(PVPM);
print(PVPMgpp, "PVPMgpp")


var withcel_gpp = PVPMgpp.select(['gpp_pvpm'])
Map.addLayer(withcel_gpp, {bands: ['gpp_pvpm']}, "withcel_gpp")

 


// // Negative DAP values
// function zeroDAPconversion(image){
//   var proxy = 0.0
//   image = image.unmask(proxy)
//   var currentHS = image.expression(
// "(b('DAP') == proxy) ? 1.0"+
// ":1*DAP",
// {"DAP" : image, "proxy": proxy});
//   return currentHS;
// }
 
 
// var withDAPfromDOP = withDAPfromDOP.map(zeroDAPconversion);
// print("withDAPfromDOP", withDAPfromDOP)
// Map.addLayer( withDAPfromDOP)