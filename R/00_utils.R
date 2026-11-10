# Shared utilities. All functions are independent of the user's working directory.
read_semicolon <- function(path) {
  if (!file.exists(path)) stop('Arquivo não encontrado: ', path)
  utils::read.csv2(path, check.names = FALSE, stringsAsFactors = FALSE,
                   fileEncoding = 'UTF-8')
}
write_semicolon <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv2(x, path, row.names = FALSE, fileEncoding = 'UTF-8')
  invisible(path)
}
require_columns <- function(df, cols, label = 'tabela') {
  absent <- setdiff(cols, names(df))
  if (length(absent)) stop(label, ': colunas ausentes: ', paste(absent, collapse = ', '))
}
partition_id <- function(x) {
  z <- trimws(as.character(x))
  if (anyNA(z) || any(!nzchar(z))) stop('Identificadores de partição ausentes')
  # For integer identifiers, canonicalize '01', 1, '001' to '1'.
  numeric_only <- grepl('^[0-9]+$', z)
  z[numeric_only] <- as.character(as.integer(z[numeric_only]))
  if (anyNA(z)) stop('Identificador numérico excede intervalo de inteiros')
  z
}
as_numeric_checked <- function(x, name, allow_na = FALSE) {
  if (is.numeric(x)) return(as.numeric(x))
  y <- suppressWarnings(as.numeric(gsub(',', '.', trimws(as.character(x)), fixed = TRUE)))
  bad <- !is.na(x) & nzchar(trimws(as.character(x))) & is.na(y)
  if (any(bad) && !allow_na) stop('Valores não numéricos em ', name, ': ', sum(bad))
  y
}
