# Join by partition ID; do not assume row order or infer longitude from projected x/y.
join_descriptors_to_centroids <- function(centroids, descriptors, include_unsaved = FALSE,
                                          target_crs = 'EPSG:4326') {
  if (!requireNamespace('sf', quietly = TRUE)) stop('Instale sf')
  if (is.character(centroids) && length(centroids) == 1L) centroids <- read_semicolon(centroids)
  if (is.character(descriptors) && length(descriptors) == 1L) descriptors <- read_semicolon(descriptors)
  require_columns(centroids, c('partition', 'centroid_x', 'centroid_y', 'crs'), 'centroides')
  require_columns(descriptors, 'partition', 'descritores')
  centroids$partition <- partition_id(centroids$partition)
  descriptors$partition <- partition_id(descriptors$partition)
  if (anyDuplicated(centroids$partition) || anyDuplicated(descriptors$partition))
    stop('Identificadores duplicados; vinculação espacial seria ambígua')
  if ('salvo' %in% names(centroids) && !include_unsaved)
    centroids <- centroids[centroids$salvo %in% TRUE, , drop = FALSE]
  missing_location <- setdiff(descriptors$partition, centroids$partition)
  missing_result <- setdiff(centroids$partition, descriptors$partition)
  if (length(missing_location) || length(missing_result))
    stop('Partições não correspondem. Sem coordenadas: ', paste(missing_location, collapse = ','),
         '; sem descritores: ', paste(missing_result, collapse = ','))
  if (length(unique(centroids$crs)) != 1L || anyNA(centroids$crs))
    stop('CRS misto ou ausente nos centroides')
  centroids$centroid_x <- as_numeric_checked(centroids$centroid_x, 'centroid_x')
  centroids$centroid_y <- as_numeric_checked(centroids$centroid_y, 'centroid_y')
  if (any(!is.finite(centroids$centroid_x)) || any(!is.finite(centroids$centroid_y))) stop('Centroides inválidos')
  matching <- match(centroids$partition, descriptors$partition)
  merged <- cbind(centroids, descriptors[matching, setdiff(names(descriptors), 'partition'), drop = FALSE])
  sf_points <- sf::st_as_sf(merged, coords = c('centroid_x', 'centroid_y'), crs = unique(centroids$crs), remove = FALSE)
  lonlat <- sf::st_coordinates(sf::st_transform(sf_points, target_crs))
  merged$longitude <- lonlat[, 1]; merged$latitude <- lonlat[, 2]
  merged
}
