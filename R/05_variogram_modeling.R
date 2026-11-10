# Isotropic empirical variograms and model selection.
# Always return cross-validation diagnostics even when selecting by SSErr.
fit_best_variogram <- function(points_sf, parameter, models = c('Sph','Exp','Gau','Mat'),
                               criterion = c('RMSE_CV', 'SSErr')) {
  criterion <- match.arg(criterion)
  if (!requireNamespace('gstat', quietly = TRUE)) stop('Instale gstat')
  if (!requireNamespace('sf', quietly = TRUE)) stop('Instale sf')
  if (!(parameter %in% names(points_sf))) stop('Parâmetro não encontrado: ', parameter)
  x <- points_sf[is.finite(points_sf[[parameter]]), ]
  if (nrow(x) < 10L) stop('Mínimo de 10 observações válidas requerido')
  if (sf::st_is_longlat(x)) stop('Projetar pontos para CRS métrico antes do variograma')
  coords <- sf::st_coordinates(x)
  if (anyDuplicated(data.frame(coords[, 1], coords[, 2]))) stop('Coordenadas duplicadas')
  f <- stats::as.formula(paste0('`', parameter, '` ~ 1'))
  v <- gstat::variogram(f, as(x, 'Spatial'))
  if (nrow(v) < 2L) stop('Pares insuficientes para semivariograma')
  s <- stats::var(x[[parameter]])
  if (!is.finite(s) || s <= 0) stop('Variância nula ou inválida em ', parameter)
  distance <- max(v$dist)
  fits <- list(); summary <- list()
  for (mod in models) {
    init <- gstat::vgm(psill = 0.7 * s, model = mod, range = distance / 3, nugget = 0.3 * s)
    fit <- try(suppressWarnings(gstat::fit.variogram(v, init)), silent = TRUE)
    if (inherits(fit, 'try-error') || any(!is.finite(fit$psill)) ||
        any(fit$psill < 0) || any(!is.finite(fit$range)) ||
        isTRUE(attr(fit, 'singular'))) next
    cv <- try(suppressWarnings(gstat::krige.cv(f, as(x, 'Spatial'), model = fit, nfold = nrow(x))), silent = TRUE)
    rmse <- if (inherits(cv, 'try-error')) NA_real_ else sqrt(mean(cv$residual^2, na.rm = TRUE))
    sse <- as.numeric(attr(fit, 'SSErr'))
    if (!is.finite(sse)) next
    fits[[mod]] <- fit
    summary[[mod]] <- data.frame(model = mod, SSErr = sse, RMSE_CV = rmse,
      nugget = sum(fit$psill[fit$model == 'Nug']),
      partial_sill = sum(fit$psill[fit$model != 'Nug']),
      range_m = fit$range[which(fit$model != 'Nug')[1]])
  }
  if (!length(fits)) stop('Nenhum modelo de semivariograma foi ajustado')
  tab <- do.call(rbind, summary); rownames(tab) <- NULL
  values <- tab[[criterion]]
  if (all(!is.finite(values))) stop('Não há valores válidos para o critério ', criterion)
  chosen <- which.min(ifelse(is.finite(values), values, Inf))
  list(variogram = v, model = fits[[tab$model[chosen]]], selected = tab$model[chosen],
       diagnostics = tab, selection_criterion = criterion)
}
