# Krige one descriptor per image, returning prediction and variance as GeoTIFFs.
# Polygon mask is applied after kriging, not used to invent missing observations.
krige_descriptor <- function(data, boundary, parameter, out_dir, image_id,
                             target_grid_n = 100L, criterion = 'RMSE_CV',
                             models = c('Sph','Exp','Gau','Mat'), utm_epsg = NULL) {
  for (p in c('sf', 'terra', 'gstat')) if (!requireNamespace(p, quietly = TRUE)) stop('Instale ', p)
  require_columns(data, c('longitude','latitude',parameter), 'dados espaciais')
  pts <- data[is.finite(data[[parameter]]) & is.finite(data$longitude) & is.finite(data$latitude), ]
  if (nrow(pts) < 10) stop('Menos de 10 partições válidas')
  pts_geo <- sf::st_as_sf(pts, coords = c('longitude','latitude'), crs = 4326)
  if (is.character(boundary)) boundary <- sf::st_read(boundary, quiet = TRUE)
  if (is.na(sf::st_crs(boundary))) stop('Limite espacial sem CRS: não atribuir WGS84 automaticamente')
  if (is.null(utm_epsg)) {
    lon <- mean(pts$longitude); lat <- mean(pts$latitude)
    zone <- max(1L, min(60L, floor((lon + 180) / 6) + 1L))
    utm_epsg <- if (lat >= 0) 32600 + zone else 32700 + zone
  }
  pts_utm <- sf::st_transform(pts_geo, utm_epsg)
  boundary_utm <- sf::st_transform(boundary, utm_epsg)
  fit <- fit_best_variogram(pts_utm, parameter, models, criterion)
  e <- terra::ext(terra::vect(boundary_utm))
  template <- terra::rast(ncols = target_grid_n, nrows = target_grid_n, ext = e,
                          crs = sf::st_crs(boundary_utm)$wkt)
  grid_xy <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
  grid <- sf::st_as_sf(data.frame(x = grid_xy[, 1], y = grid_xy[, 2]), coords = c('x','y'), crs = utm_epsg)
  f <- stats::as.formula(paste0('`', parameter, '` ~ 1'))
  prediction <- gstat::krige(f, as(pts_utm, 'Spatial'), as(grid, 'Spatial'), model = fit$model)
  pred <- template; variance <- template
  terra::values(pred) <- prediction$var1.pred
  terra::values(variance) <- prediction$var1.var
  pred <- terra::mask(pred, terra::vect(boundary_utm))
  variance <- terra::mask(variance, terra::vect(boundary_utm))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  prefix <- paste0(gsub('[^A-Za-z0-9_-]', '_', image_id), '_', gsub('[^A-Za-z0-9_-]', '_', parameter))
  pred_path <- file.path(out_dir, paste0(prefix, '_prediction_utm.tif'))
  var_path <- file.path(out_dir, paste0(prefix, '_variance_utm.tif'))
  terra::writeRaster(pred, pred_path, overwrite = TRUE)
  terra::writeRaster(variance, var_path, overwrite = TRUE)
  # WGS84 copies are optional convenience outputs; UTM versions are the computation's source.
  geo_pred <- terra::project(pred, 'EPSG:4326', method = 'bilinear')
  geo_var <- terra::project(variance, 'EPSG:4326', method = 'bilinear')
  terra::writeRaster(geo_pred, file.path(out_dir, paste0(prefix, '_prediction_wgs84.tif')), overwrite = TRUE)
  terra::writeRaster(geo_var, file.path(out_dir, paste0(prefix, '_variance_wgs84.tif')), overwrite = TRUE)
  write_semicolon(fit$diagnostics, file.path(out_dir, paste0(prefix, '_model_comparison.csv')))
  list(prediction = pred_path, variance = var_path, selected_model = fit$selected,
       diagnostics = fit$diagnostics, utm_epsg = utm_epsg)
}
