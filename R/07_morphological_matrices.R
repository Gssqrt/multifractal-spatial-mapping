# Combine co-registered kriging outputs into spatial descriptor matrix.
rasters_to_matrix <- function(paths, parameter_names, out_csv = NULL) {
  if (!requireNamespace('terra', quietly = TRUE)) stop('Instale terra')
  if (length(paths) != length(parameter_names) || !length(paths) || anyDuplicated(parameter_names))
    stop('paths e parameter_names inválidos')
  rs <- lapply(paths, terra::rast)
  ref <- rs[[1L]]
  for (i in seq_along(rs)[-1L]) {
    if (!terra::compareGeom(ref, rs[[i]], stopOnError = FALSE))
      stop('Rasters não alinhados; reprojete/reamostre explicitamente antes da matriz')
  }
  stack <- do.call(c, rs); names(stack) <- parameter_names
  result <- terra::as.data.frame(stack, xy = TRUE, na.rm = FALSE)
  names(result)[1:2] <- c('x', 'y')
  if (!is.null(out_csv)) write_semicolon(result, out_csv)
  result
}
