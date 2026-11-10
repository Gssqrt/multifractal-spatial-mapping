# Optional downstream analysis; clustering interpolated pixels is not independent replication.
cluster_morphology <- function(matrix_df, descriptors = c('alpha.0','d.alpha','df.alpha'), k = 2L) {
  require_columns(matrix_df, descriptors, 'matriz morfológica')
  ok <- stats::complete.cases(matrix_df[, descriptors, drop = FALSE])
  if (sum(ok) < max(3L, k)) stop('Dados insuficientes')
  X <- scale(matrix_df[ok, descriptors, drop = FALSE])
  if (any(!is.finite(X))) stop('Descritor constante ou não finito')
  hc <- stats::hclust(stats::dist(X), method = 'ward.D2')
  matrix_df$cluster <- NA_integer_
  matrix_df$cluster[ok] <- stats::cutree(hc, k = k)
  list(data = matrix_df, dendrogram = hc)
}
