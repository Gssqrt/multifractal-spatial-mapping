# Attach partition IDs to cleaned FracLac numerical spectra.
# The 219 rows/partition reported for the study is a configurable data convention.
# This function deliberately DOES NOT discard textual lines: the numeric table
# must already be cleaned and ordered exactly as the processed tiles.
index_spectra <- function(spectrum, partition_order, rows_per_partition = 219L,
                          verify_q_sequence = TRUE) {
  if (is.character(spectrum) && length(spectrum) == 1L) spectrum <- read_semicolon(spectrum)
  if (is.character(partition_order) && length(partition_order) == 1L && file.exists(partition_order))
    partition_order <- read_semicolon(partition_order)$partition
  partition_order <- partition_id(partition_order)
  if (anyDuplicated(partition_order)) stop('partition_order contém identificadores duplicados')
  rows_per_partition <- as.integer(rows_per_partition)
  if (is.na(rows_per_partition) || rows_per_partition < 1L) stop('rows_per_partition inválido')
  if ('partition' %in% names(spectrum)) stop('Espectro já tem coluna partition; não sobrescrever')
  expected <- length(partition_order) * rows_per_partition
  if (nrow(spectrum) != expected) stop('Linhas do espectro: ', nrow(spectrum),
     '; esperado: ', expected, '. Verifique remoção das linhas textuais, ordem e contagem por partição.')
  if (verify_q_sequence && 'q' %in% names(spectrum) && length(partition_order) > 1) {
    q <- as_numeric_checked(spectrum$q, 'q')
    q0 <- q[seq_len(rows_per_partition)]
    for (i in seq_along(partition_order)[-1L]) {
      inds <- seq.int((i - 1L) * rows_per_partition + 1L, i * rows_per_partition)
      if (!isTRUE(all.equal(q[inds], q0, tolerance = 1e-10)))
        stop('Sequência de q difere na partição ', partition_order[i],
             '; indexação por posição não é segura')
    }
  }
  spectrum$partition <- rep(partition_order, each = rows_per_partition)
  spectrum
}
