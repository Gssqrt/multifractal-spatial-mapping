#!/usr/bin/env Rscript
# Usage: Rscript scripts/run_workflow.R cleaned_spectra.csv partition_order.csv centroids.csv boundary.gpkg output_dir [rows_per_partition]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) stop('Uso: Rscript scripts/run_workflow.R espectros.csv ordem.csv centroides.csv limite.gpkg saida/ [219]')
root <- normalizePath(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])), '..'), mustWork = TRUE)
for (f in sort(list.files(file.path(root, 'R'), pattern = '\\.R$', full.names = TRUE))) source(f)
sp <- read_semicolon(args[1]); order_df <- read_semicolon(args[2]); cent <- read_semicolon(args[3])
require_columns(order_df, 'partition', 'ordem')
rows <- if (length(args) >= 6L) as.integer(args[6]) else 219L
indexed <- index_spectra(sp, order_df$partition, rows_per_partition = rows)
write_semicolon(indexed, file.path(args[5], '01_indexed_spectra.csv'))
desc <- extract_descriptors(indexed)
write_semicolon(desc, file.path(args[5], '02_multifractal_descriptors.csv'))
spatial <- join_descriptors_to_centroids(cent, desc)
write_semicolon(spatial, file.path(args[5], '03_spatial_descriptors.csv'))
for (p in c('d.alpha', 'df.alpha', 'alpha.0')) {
  krige_descriptor(spatial, args[4], p, file.path(args[5], 'kriging'),
                   image_id = if ('image_id' %in% names(spatial)) unique(spatial$image_id)[1] else 'image')
}
cat('Concluído; resultados em: ', args[5], '\n')
