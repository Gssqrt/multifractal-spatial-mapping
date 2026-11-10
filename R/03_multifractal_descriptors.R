# Extract descriptors from q, a = alpha(q), f.a. = f(alpha), d.q. = D(q).
# Delta f sign follows the original script: f(alpha_max) - f(alpha_min).
extract_descriptors <- function(spectra, q_tolerance = 1e-8) {
  if (is.character(spectra) && length(spectra) == 1L) spectra <- read_semicolon(spectra)
  require_columns(spectra, c('partition', 'q', 'a', 'f.a.', 'd.q.'), 'espectros')
  spectra$partition <- partition_id(spectra$partition)
  for (c in c('q', 'a', 'f.a.', 'd.q.')) spectra[[c]] <- as_numeric_checked(spectra[[c]], c)
  calc <- function(g) {
    if (any(!is.finite(g$a)) || any(!is.finite(g$f.a.))) stop('alpha/f(alpha) não finitos na partição ', g$partition[1])
    ia_min <- which.min(g$a); ia_max <- which.max(g$a); if_max <- which.max(g$f.a.)
    # Refuse ambiguous D(q) lookup; no silent first-match selection.
    d_at <- function(qvalue) {
      where <- which(abs(g$q - qvalue) <= q_tolerance)
      if (length(where) != 1L) return(NA_real_)
      g$d.q.[where]
    }
    a0 <- g$a[if_max]; den <- a0 - min(g$a)
    data.frame(partition = g$partition[1],
      d.alpha = max(g$a) - min(g$a),
      df.alpha = g$f.a.[ia_max] - g$f.a.[ia_min],
      alpha.0 = a0, D0 = d_at(0), D1 = d_at(1), D2 = d_at(2),
      r = if (abs(den) <= .Machine$double.eps^0.5) NA_real_ else (max(g$a) - a0) / den)
  }
  result <- do.call(rbind, lapply(split(spectra, spectra$partition), calc))
  rownames(result) <- NULL
  result
}
