# Synthetic demonstration; not actual satellite imagery or FracLac measurements.
root <- normalizePath(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])), '..'), mustWork = TRUE)
for (f in sort(list.files(file.path(root,'R'), pattern = '\\.R$', full.names = TRUE))) source(f)
for (p in c('sf','terra')) if (!requireNamespace(p, quietly = TRUE)) stop('Instale ',p)
set.seed(20261110)
dir.create(file.path(root,'examples','input'), recursive = TRUE, showWarnings = FALSE)
out <- file.path(root,'examples','input')
# 7 x 7 synthetic tiles, each with 5 q values (for illustrative geostatistical fitting).
ids <- sprintf('%02d', 1:49)
coords <- expand.grid(x = seq(-39.8, -39.2, length.out = 7), y = seq(-10.8, -10.2, length.out = 7))
manifest <- data.frame(image_id='synthetic_scene', partition=ids,
  centroid_x=coords$x, centroid_y=coords$y, crs='EPSG:4326', salvo=TRUE)
write_semicolon(manifest, file.path(out,'centroids.csv'))
write_semicolon(data.frame(partition=ids), file.path(out,'partition_order.csv'))
q <- seq(-2,2,by=1)
spec <- do.call(rbind, lapply(seq_along(ids), function(i) {
  alpha <- 1.1 + (0.065 + 0.001*i + 0.007*sin(i/4))*q + 0.02*i + 0.004*q^2
  data.frame(q=q, a=alpha, f.a.=1.75 - (q*0.10)^2 + (0.001*i + 0.015*cos(i/5))*q - (i %% 3)*0.04,
             d.q.=1.3 - q*0.06 + i*0.007)
}))
write_semicolon(spec, file.path(out,'cleaned_spectra.csv'))
poly <- sf::st_as_sfc(sf::st_bbox(c(xmin=-40, ymin=-11, xmax=-39, ymax=-10), crs=sf::st_crs(4326)))
sf::st_write(sf::st_sf(name='synthetic_boundary', geometry=poly),
             file.path(out,'boundary.gpkg'), delete_dsn=TRUE, quiet=TRUE)
cat('Arquivos sintéticos criados em examples/input; use rows_per_partition=5\n')
