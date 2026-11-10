# Partition a georeferenced, one-band raster into numbered tiles.
# It does not assign/guess a CRS: the input MUST already be georeferenced.
partition_image <- function(infile, out_dir, nx = 17L, ny = 17L,
                            drop_border = TRUE, binarize = TRUE,
                            foreground_values = NULL, background_value = 0,
                            skip_empty = TRUE, overwrite = FALSE) {
  if (!requireNamespace('terra', quietly = TRUE)) stop('Instale terra')
  r <- terra::rast(infile)
  if (terra::nlyr(r) != 1L) stop('Imagem deve ter exatamente uma banda')
  if (!nzchar(terra::crs(r))) stop('Raster sem CRS; georreferencie antes de particionar')
  if (nx < 1 || ny < 1 || nx > terra::ncol(r) || ny > terra::nrow(r)) stop('nx/ny inválidos')
  nx <- as.integer(nx); ny <- as.integer(ny)
  base <- tools::file_path_sans_ext(basename(infile))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  h <- if (drop_border) terra::nrow(r) %/% ny else ceiling(terra::nrow(r) / ny)
  w <- if (drop_border) terra::ncol(r) %/% nx else ceiling(terra::ncol(r) / nx)
  if (drop_border) {
    e <- terra::ext(r)
    resolution <- terra::res(r)
    r <- terra::crop(r, terra::ext(e[1], e[1] + w * nx * resolution[1],
                                    e[4] - h * ny * resolution[2], e[4]), snap = 'near')
  }
  manifest <- vector('list', nx * ny)
  k <- 0L
  for (iy in seq_len(ny)) for (ix in seq_len(nx)) {
    k <- k + 1L
    r1 <- (iy - 1L) * h + 1L; r2 <- min(iy * h, terra::nrow(r))
    c1 <- (ix - 1L) * w + 1L; c2 <- min(ix * w, terra::ncol(r))
    # Exact pixel indices prevent accumulating floating-point coordinate errors.
    cells <- as.vector(outer(seq.int(r1, r2), seq.int(c1, c2),
                             function(rr, cc) (rr - 1L) * terra::ncol(r) + cc))
    xy <- terra::xyFromCell(r, cells)
    rx <- terra::res(r)[1]; ry <- terra::res(r)[2]
    xmin <- min(xy[, 1]) - rx / 2; xmax <- max(xy[, 1]) + rx / 2
    ymin <- min(xy[, 2]) - ry / 2; ymax <- max(xy[, 2]) + ry / 2
    tile <- terra::crop(r, terra::ext(xmin, xmax, ymin, ymax), snap = 'near')
    if (binarize) {
      if (is.null(foreground_values)) {
        tile <- terra::ifel(is.na(tile), NA, terra::ifel(tile != background_value, 255, 0))
      } else {
        # Explicit foreground values recommended for masks where 255 is nodata.
        tile <- terra::ifel(is.na(tile), NA,
                            terra::ifel(tile %in% foreground_values, 255, 0))
      }
    }
    v <- terra::values(tile, mat = FALSE)
    valid <- v[!is.na(v)]
    reason <- if (!length(valid)) 'all_NA' else if (binarize && !any(valid == 255)) 'all_background' else NA_character_
    save <- !(skip_empty && !is.na(reason))
    id <- sprintf('%0*d', max(2L, nchar(as.character(nx * ny))), k)
    filename <- file.path(out_dir, paste0(base, '_', id, '.tif'))
    if (save) terra::writeRaster(tile, filename, overwrite = overwrite,
                                datatype = if (binarize) 'INT1U' else terra::datatype(r),
                                gdal = c('PROFILE=BASELINE', 'COMPRESS=LZW', 'TILED=NO'))
    manifest[[k]] <- data.frame(image_id = base, partition = id,
      linha_idx = iy, coluna_idx = ix, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      centroid_x = (xmin + xmax) / 2, centroid_y = (ymin + ymax) / 2,
      arquivo = if (save) filename else NA_character_, salvo = save,
      motivo = reason, crs = terra::crs(r), stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, manifest)
  write_semicolon(df, file.path(out_dir, paste0(base, '_centroides.csv')))
  df
}
