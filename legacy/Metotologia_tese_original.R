library(dplyr)
library(openxlsx)
library(readxl)
library(sp)
library(terra)
library(gstat)
library(raster)
library(viridis)
library(sf)  

###############################REPROJEÇÃO#######################################
################################################################################
################################################################################

ras_nao_ref <- raster("Band1_class.tif")
ras_ref     <- raster("coord.tif")

# inverter verticalmente (sul↔norte)
#ras_nao_ref <- flip(ras_nao_ref, direction = "z")

# copiar moldura e projeção da referência
extent(ras_nao_ref) <- extent(ras_ref)
crs(ras_nao_ref)    <- crs(ras_ref)

writeRaster(ras_nao_ref, "band1_class_geo.tif", overwrite = TRUE)

########################RECORTE E EXTRAÇÂO DE COORDENADAS#######################
################################################################################
################################################################################

particionar_por_indice_com_centroides <- function(infile, out_dir,
                                                  nx = 17, ny = 17,
                                                  overwrite = FALSE,
                                                  drop_border = TRUE,
                                                  csv_path = NULL,
                                                  imagej_baseline = TRUE,
                                                  binarizar_8bit = TRUE,
                                                  skip_all_zero = TRUE) {
  
  r <- rast(infile)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  base_name <- tools::file_path_sans_ext(basename(infile))
  
  nr <- nrow(r); nc <- ncol(r)
  resx <- res(r)[1]; resy <- res(r)[2]
  
  tile_h <- floor(nr / ny)
  tile_w <- floor(nc / nx)
  
  if (drop_border) {
    use_nr <- tile_h * ny
    use_nc <- tile_w * nx
    x_left  <- xFromCol(r, 1)      - resx/2
    y_top   <- yFromRow(r, 1)      + resy/2
    x_right <- xFromCol(r, use_nc) + resx/2
    y_bot   <- yFromRow(r, use_nr) - resy/2
    r <- crop(r, ext(x_left, x_right, y_bot, y_top))
    nr <- nrow(r); nc <- ncol(r)
  } else {
    tile_h <- ceiling(nr / ny)
    tile_w <- ceiling(nc / nx)
  }
  
  gdal_opts <- if (imagej_baseline) {
    c("PROFILE=BASELINE","COMPRESS=LZW","TILED=NO","TFW=YES",
      "BIGTIFF=IF_NEEDED","INTERLEAVE=BAND")
  } else {
    c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_NEEDED","INTERLEAVE=BAND")
  }
  
  regs <- vector("list", nx * ny)
  idx <- 1L
  
  for (iy in 0:(ny - 1)) {
    r1 <- iy * tile_h + 1L
    r2 <- min(r1 + tile_h - 1L, nr)
    
    y_max <- yFromRow(r, r1) + resy/2
    y_min <- yFromRow(r, r2) - resy/2
    cy    <- (y_max + y_min) / 2
    
    for (ix in 0:(nx - 1)) {
      c1 <- ix * tile_w + 1L
      c2 <- min(c1 + tile_w - 1L, nc)
      
      x_min <- xFromCol(r, c1) - resx/2
      x_max <- xFromCol(r, c2) + resx/2
      cx    <- (x_min + x_max) / 2
      
      etile <- ext(x_min, x_max, y_min, y_max)
      tile  <- crop(r, etile)
      
      suf <- sprintf("%02d", idx)
      outfile <- file.path(out_dir, paste0(base_name, "_", suf, ".tif"))
      
      # ---- BINARIZAÇÃO 8-BIT (0 / 255) ----
      # NA permanece NA; !=0 vira 255; ==0 vira 0
      if (binarizar_8bit) {
        tile <- ifel(is.na(tile), NA, ifel(tile != 0, 255, 0))
      }
      
      # ---- (opcional) NÃO SALVAR tiles all-zero (mas registrar no CSV) ----
      salvar <- TRUE
      motivo <- NA_character_
      
      if (skip_all_zero) {
        v <- values(tile, mat = FALSE)
        v2 <- v[!is.na(v)]
        if (length(v2) == 255) { ## para as de sentinel-3 com sulfixo _geo, o 255 representa os valores ausentes e 0 os valores com informação
          salvar <- FALSE
          motivo <- "all_NA"
        } else if (!any(v2 != 255)) {
          salvar <- FALSE
          motivo <- "all_zero"
        }
      }
      
      if (salvar) {
        writeRaster(tile, outfile,
                    overwrite = overwrite,
                    filetype = "GTiff",
                    datatype = "INT1U",   # 8-bit unsigned
                    gdal = gdal_opts)
      }
      
      regs[[idx]] <- data.frame(
        partition   = suf,
        linha_idx  = iy + 1L,
        coluna_idx = ix + 1L,
        xmin = x_min, xmax = x_max,
        ymin = y_min, ymax = y_max,
        centroid_x = cx, centroid_y = cy,
        arquivo    = if (salvar) outfile else NA_character_,
        salvo      = salvar,
        motivo     = motivo,
        crs        = if (!is.na(crs(r))) as.character(crs(r)) else NA_character_
      )
      
      idx <- idx + 1L
    }
  }
  
  df <- do.call(rbind, regs)
  
  df <- do.call(rbind, regs)
  
  # =========================================================
  # ALTERAÇÃO ADICIONADA:
  # Conta quantas partições foram salvas, ou seja,
  # quantas partições possuem informação de interesse.
  # =========================================================
  n_info <- sum(df$salvo, na.rm = TRUE)
  cat(n_info, "partições têm informação\n")
  # =========================================================
  
  if (is.null(csv_path)) {
    csv_path <- file.path(out_dir, paste0(base_name, "_centroides.csv"))
  }
  
  utils::write.csv2(df, csv_path, row.names = FALSE, fileEncoding = "UTF-8")
  message("OK: centroides salvos em: ", csv_path)
  invisible(df)
}


# Exemplo de uso:
infile  <- "20210530_geo.tif"
out_dir <- "C:/Users/Gssqrt/OneDrive/Área de Trabalho/DIRETOR/Doutorado/Trabalhos/Artigo dissertação/paper/review/new_base/imagens/seccoes/17/partition"

df_centroides <- particionar_por_indice_com_centroides(
  infile, out_dir, nx = 17, ny = 17, overwrite = TRUE,
  imagej_baseline = TRUE
)

head(df_centroides)
View(df_centroides)


#################RECORTE E EXTRAÇÂO DE COORDENADAS CONTINUO#####################
################################################################################
################################################################################

library(terra)

particionar_imagem_original <- function(infile, out_dir,
                                        nx = 17, ny = 17,
                                        overwrite = TRUE,
                                        drop_border = TRUE,
                                        csv_path = NULL,
                                        imagej_baseline = TRUE) {
  
  r <- rast(infile)
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  base_name <- tools::file_path_sans_ext(basename(infile))
  
  nr <- nrow(r)
  nc <- ncol(r)
  resx <- res(r)[1]
  resy <- res(r)[2]
  
  tile_h <- floor(nr / ny)
  tile_w <- floor(nc / nx)
  
  if (drop_border) {
    use_nr <- tile_h * ny
    use_nc <- tile_w * nx
    
    x_left  <- xFromCol(r, 1)      - resx / 2
    y_top   <- yFromRow(r, 1)      + resy / 2
    x_right <- xFromCol(r, use_nc) + resx / 2
    y_bot   <- yFromRow(r, use_nr) - resy / 2
    
    r <- crop(r, ext(x_left, x_right, y_bot, y_top))
    
    nr <- nrow(r)
    nc <- ncol(r)
  }
  
  gdal_opts <- if (imagej_baseline) {
    c(
      "PROFILE=BASELINE",
      "COMPRESS=LZW",
      "TILED=NO",
      "TFW=YES",
      "BIGTIFF=IF_NEEDED",
      "INTERLEAVE=BAND"
    )
  } else {
    c(
      "COMPRESS=LZW",
      "TILED=YES",
      "BIGTIFF=IF_NEEDED",
      "INTERLEAVE=BAND"
    )
  }
  
  regs <- vector("list", nx * ny)
  idx <- 1L
  
  for (iy in 0:(ny - 1)) {
    
    r1 <- iy * tile_h + 1L
    r2 <- min(r1 + tile_h - 1L, nr)
    
    y_max <- yFromRow(r, r1) + resy / 2
    y_min <- yFromRow(r, r2) - resy / 2
    cy <- (y_max + y_min) / 2
    
    for (ix in 0:(nx - 1)) {
      
      c1 <- ix * tile_w + 1L
      c2 <- min(c1 + tile_w - 1L, nc)
      
      x_min <- xFromCol(r, c1) - resx / 2
      x_max <- xFromCol(r, c2) + resx / 2
      cx <- (x_min + x_max) / 2
      
      etile <- ext(x_min, x_max, y_min, y_max)
      tile <- crop(r, etile)
      
      suf <- sprintf("%02d", idx)
      outfile <- file.path(out_dir, paste0(base_name, "_", suf, ".tif"))
      
      writeRaster(
        tile,
        outfile,
        overwrite = overwrite,
        filetype = "GTiff",
        gdal = gdal_opts
      )
      
      regs[[idx]] <- data.frame(
        partition   = suf,
        linha_idx  = iy + 1L,
        coluna_idx = ix + 1L,
        xmin = x_min,
        xmax = x_max,
        ymin = y_min,
        ymax = y_max,
        centroid_x = cx,
        centroid_y = cy,
        arquivo = outfile,
        salvo = TRUE,
        crs = if (!is.na(crs(r))) as.character(crs(r)) else NA_character_
      )
      
      idx <- idx + 1L
    }
  }
  
  df <- do.call(rbind, regs)
  
  if (is.null(csv_path)) {
    csv_path <- file.path(out_dir, paste0(base_name, "_centroides.csv"))
  }
  
  write.csv2(df, csv_path, row.names = FALSE, fileEncoding = "UTF-8")
  
  return(df)
}


infile <- "2023-09-09-00_00_2023-09-09-23_59_Sentinel-2_L2A_B03_(Raw).tiff"

out_dir <- "C:/Users/Gustavo/OneDrive/Área de Trabalho/DIRETOR/Doutorado/Trabalhos/Plumas/Base/São francisco/MF/560"

df_centroides <- particionar_imagem_original(
  infile = infile,
  out_dir = out_dir,
  nx = 17,
  ny = 17,
  overwrite = TRUE,
  drop_border = TRUE,
  imagej_baseline = TRUE
)

head(df_centroides)
View(df_centroides)

######################## REPETIR LINAS PLANILHAS ###############################
################################################################################
################################################################################

# Carregar a planilha
data <- read.csv2("FracLac 2015Sep090313a9330MF DATASET SUMMARY.txt")

# Repetir cada observação 216 vezes consecutivamente
data_repeated <- data[rep(seq_len(nrow(data)), each = 219), ]

# Salvar o resultado em uma nova planilha
write.csv2(data_repeated, "ordem_repetida.csv", row.names = FALSE)

# Verificar as primeiras linhas da nova planilha
head(data_repeated)



###################COMPOSIÇÃO DAS METRICAS MULTIFRACTAIS########################
################################################################################
################################################################################

# Defina o diretório onde estão as planilhas
diretorio <- getwd() 
# data <- read.csv2('mf_30.csv')
# head(data)
# Lista todos os arquivos .xlsx no diretório
arquivos <- list.files(path = diretorio, pattern = "*.csv", full.names = TRUE)

# Função para processar cada arquivo
processa_planilha <- function(arquivo) {
  # Leia a planilha
  sing <- read.csv2(arquivo)
  
  # Converta as colunas relevantes para numérico, se necessário
  sing <- sing |> 
    mutate(
      a = as.numeric(a),
      f.a. = as.numeric(f.a.),
      d.q. = as.numeric(d.q.),
      q = as.numeric(q)
    )
  
  
  metrics_image <- sing |> 
    group_by(partition) |> 
    summarise(
      `d.alpha` = max(a) - min(a),
      `df.alpha` = f.a.[which.max(a)] - f.a.[which.min(a)],
      `alpha.0` = a[which.max(f.a.)],
      `D0` = d.q.[which(q == 0)],
      `D1` = d.q.[which(q == 1)],
      `D2` = d.q.[which(q == 2)],
      r = (max(a) - `alpha.0`) / (`alpha.0` - min(a))
    )
  
  # Retorna o nome da planilha e os resultados
  list(nome = basename(arquivo), dados = metrics_image)
}

# Processar todas as planilhas
resultados <- lapply(arquivos, processa_planilha)

# Salvar cada resultado em um arquivo separado ou combinar
for (resultado in resultados) {
  # Defina o nome do arquivo de saída
  nome_arquivo_saida <- paste0("parameters_", resultado$nome)
  
  # Salvar o arquivo Excel
  write.csv2(resultado$dados, file = file.path(diretorio, nome_arquivo_saida))
}

###########################ADICIONAR COORDENADAS################################
################################################################################
################################################################################



library(dplyr)

diretorio <- getwd()

# Nomes EXATOS das duas planilhas
arquivo1 <- "2023-09-09-00_00_2023-09-09-23_59_Sentinel-2_L2A_B03_(Raw)_centroides.csv"
arquivo2 <- "parameters_ordem_repetida.csv"



path1 <- file.path(diretorio, arquivo1)
path2 <- file.path(diretorio, arquivo2)

df1 <- read.csv2(path1, stringsAsFactors = FALSE)
df2 <- read.csv2(path2, stringsAsFactors = FALSE)

# --- validações para garantir alinhamento ---
if (!("partition" %in% names(df1))) stop("df1 não tem a coluna 'partition'.")
if (!("partition" %in% names(df2))) stop("df2 não tem a coluna 'partition'.")

if (any(duplicated(df1$partition))) stop("df1 tem 'partition' duplicado (isso bagunça o alinhamento).")
if (any(duplicated(df2$partition))) stop("df2 tem 'partition' duplicado (isso bagunça o alinhamento).")

df1 <- df1 %>% arrange(partition)
df2 <- df2 %>% arrange(partition)

# Join por partition + resultado ordenado (alinhado)
resultado <- full_join(df1, df2, by = "partition") %>% arrange(partition)

# (opcional) avisos se estiver faltando partition em algum lado
faltam_em_2 <- setdiff(df1$partition, df2$partition)
faltam_em_1 <- setdiff(df2$partition, df1$partition)

if (length(faltam_em_2) > 0) message("Partitions que estão só no arquivo1: ", paste(faltam_em_2, collapse = ", "))
if (length(faltam_em_1) > 0) message("Partitions que estão só no arquivo2: ", paste(faltam_em_1, collapse = ", "))

# Salvar
saida <- file.path(diretorio, "parametros_coo.csv")
write.csv2(resultado, saida, row.names = FALSE)

###################COMPOSIÇÃO DAS MATRIZAS MORFOLIGICAS#########################
################################################################################
################################################################################

library(sp)
library(gstat)
library(raster)
library(viridis)
library(sf)
library(tools)
library(dplyr)

# ============================================================
# 1) GENERAL SETTINGS
# ============================================================

# Multifractal parameters to interpolate
parametros <- c(
  "d.alpha",
  "df.alpha",
  "alpha.0"
)

# Semivariogram models to test
modelos_testar <- c("Sph", "Exp", "Gau", "Mat")

# Selection criterion:
# "SSErr"   = best fit to the empirical semivariogram
# "RMSE_CV" = best leave-one-out cross-validation RMSE
criterio_escolha <- "SSErr"

# Interpolation grid resolution in UTM
ncol_grid <- 100
nrow_grid <- 100

# Keep zero values because they are real values
remover_zeros <- FALSE

# ============================================================
# 2) SELECT INPUT FILES AND OUTPUT DIRECTORY
# ============================================================

# Select one or more independent CSV files
arquivos_csv <- choose.files(
  caption = "Select the independent CSV files",
  filters = matrix(c("CSV files", "*.csv"), ncol = 2, byrow = TRUE)
)

if (length(arquivos_csv) == 0) {
  stop("No CSV file was selected.")
}

# Select shapefile of the study area
shp_path <- file.choose()

# Select output directory
dir_saida <- choose.dir(caption = "Select the output folder")

if (is.na(dir_saida)) {
  stop("No output folder was selected.")
}

# Create output folders
dir_rasters <- file.path(dir_saida, "ordinary_kriging_rasters_WGS84")
dir_variog <- file.path(dir_saida, "semivariograms")

dir.create(dir_rasters, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_variog, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 3) AUXILIARY FUNCTIONS
# ============================================================

# ------------------------------------------------------------
# Convert columns safely to numeric
# Useful when numbers come with comma decimal separator
# ------------------------------------------------------------

to_num <- function(x) {
  as.numeric(gsub(",", ".", as.character(x)))
}

# ------------------------------------------------------------
# Automatically calculate UTM EPSG from longitude and latitude
# ------------------------------------------------------------

calc_utm_epsg <- function(lon, lat) {
  
  lon <- to_num(lon)
  lat <- to_num(lat)
  
  zona <- floor((mean(lon, na.rm = TRUE) + 180) / 6) + 1
  
  if (mean(lat, na.rm = TRUE) >= 0) {
    epsg <- 32600 + zona
  } else {
    epsg <- 32700 + zona
  }
  
  return(epsg)
}

# ------------------------------------------------------------
# Extract image day from independent CSV file name
# Example:
# parametros_coo_14.csv -> 14
# parameters_mf_coo_28.csv -> 28
# ------------------------------------------------------------

extrair_dia_imagem <- function(nome_base) {
  
  dia <- sub(".*_(\\d+)$", "\\1", nome_base)
  
  if (dia == nome_base) {
    dia <- "unknown"
  }
  
  return(dia)
}

# ------------------------------------------------------------
# Create semivariogram title with mathematical symbols
# ------------------------------------------------------------

titulo_variograma_simbolo <- function(nome_base, param) {
  
  dia <- extrair_dia_imagem(nome_base)
  
  if (param == "alpha.0") {
    return(
      bquote("Image from day " ~ .(dia) ~ "-" ~ alpha[0])
    )
  }
  
  if (param == "d.alpha") {
    return(
      bquote("Image from day " ~ .(dia) ~ "-" ~ Delta * alpha)
    )
  }
  
  if (param == "df.alpha") {
    return(
      bquote("Image from day " ~ .(dia) ~ "-" ~ Delta * f(alpha))
    )
  }
  
  # if (param == "lacunarity") {
  #   return(
  #     bquote("Image from day " ~ .(dia) ~ "-" ~ Lambda)
  #   )
  # }
  
  return(
    paste("Image from day", dia, "-", param)
  )
}

# ------------------------------------------------------------
# Safe names for output files
# ------------------------------------------------------------

nome_parametro_arquivo <- function(param) {
  
  nomes <- c(
    "alpha.0" = "alpha.0",
    "d.alpha" = "d.alpha",
    "df.alpha" = "df.alpha"
  )
  
  if (param %in% names(nomes)) {
    return(unname(nomes[param]))
  } else {
    return(gsub("\\.", "_", param))
  }
}

# ------------------------------------------------------------
# Extract fitted semivariogram parameters
# ------------------------------------------------------------

extrair_parametros_variograma <- function(fit) {
  
  nugget <- fit$psill[fit$model == "Nug"]
  
  if (length(nugget) == 0) {
    nugget <- 0
  }
  
  estrutura <- fit[fit$model != "Nug", ]
  
  partial_sill <- estrutura$psill[1]
  range <- estrutura$range[1]
  modelo <- as.character(estrutura$model[1])
  
  total_sill <- nugget + partial_sill
  
  if ("kappa" %in% names(estrutura)) {
    kappa <- estrutura$kappa[1]
  } else {
    kappa <- NA
  }
  
  return(data.frame(
    model = modelo,
    nugget = nugget,
    partial_sill = partial_sill,
    total_sill = total_sill,
    range = range,
    kappa = kappa
  ))
}

# ------------------------------------------------------------
# Fit several semivariogram models and select the best one
# ------------------------------------------------------------

ajustar_melhor_variograma <- function(dados_sp, parametro, modelos_testar, criterio = "SSErr") {
  
  formula_var <- as.formula(paste(parametro, "~ 1"))
  
  # Empirical semivariogram
  vario_emp <- variogram(formula_var, dados_sp)
  
  resultados <- list()
  
  for (mod in modelos_testar) {
    
    cat("   Testing model:", mod, "\n")
    
    valor_var <- var(dados_sp@data[[parametro]], na.rm = TRUE)
    dist_max <- max(vario_emp$dist, na.rm = TRUE)
    
    modelo_inicial <- vgm(
      psill = valor_var * 0.7,
      model = mod,
      range = dist_max / 3,
      nugget = valor_var * 0.3
    )
    
    fit <- try(
      suppressWarnings(
        fit.variogram(vario_emp, modelo_inicial)
      ),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error")) {
      next
    }
    
    # Avoid problematic models
    if (any(is.na(fit$psill)) || any(fit$psill < 0) || any(is.na(fit$range))) {
      next
    }
    
    sserr <- attr(fit, "SSErr")
    rmse_cv <- NA
    
    # Optional leave-one-out cross-validation
    if (criterio == "RMSE_CV") {
      
      cv <- try(
        suppressWarnings(
          krige.cv(formula_var, dados_sp, model = fit)
        ),
        silent = TRUE
      )
      
      if (!inherits(cv, "try-error")) {
        rmse_cv <- sqrt(mean(cv$residual^2, na.rm = TRUE))
      }
    }
    
    pars <- extrair_parametros_variograma(fit)
    
    resultados[[mod]] <- list(
      model_name = mod,
      empirical_variogram = vario_emp,
      fit = fit,
      SSErr = sserr,
      RMSE_CV = rmse_cv,
      parameters = pars
    )
  }
  
  if (length(resultados) == 0) {
    stop(paste("No semivariogram model could be fitted for", parametro))
  }
  
  # Select best model
  if (criterio == "RMSE_CV") {
    valores <- sapply(resultados, function(x) x$RMSE_CV)
  } else {
    valores <- sapply(resultados, function(x) x$SSErr)
  }
  
  melhor_nome <- names(which.min(valores))
  melhor <- resultados[[melhor_nome]]
  
  return(melhor)
}

# ------------------------------------------------------------
# Save semivariogram plot in English
# ------------------------------------------------------------

salvar_plot_variograma <- function(vario_emp, fit, pars, arquivo_saida, titulo, sserr, rmse_cv) {
  
  png(arquivo_saida, width = 2400, height = 1800, res = 300)
  
  # Theoretical curve of the fitted semivariogram model
  linha_modelo <- variogramLine(
    fit,
    maxdist = max(vario_emp$dist, na.rm = TRUE),
    n = 200
  )
  
  # Y-axis limit
  ylim_max <- max(c(vario_emp$gamma, linha_modelo$gamma), na.rm = TRUE) * 1.15
  
  # Empirical semivariogram
  plot(
    vario_emp$dist,
    vario_emp$gamma,
    pch = 16,
    col = "black",
    xlab = "Distance (m)",
    ylab = "Semivariance",
    main = titulo,
    ylim = c(0, ylim_max),
    cex.main = 1.3,
    cex.lab = 1.2,
    cex.axis = 1.0
  )
  
  # Fitted semivariogram model
  lines(
    linha_modelo$dist,
    linha_modelo$gamma,
    lwd = 2,
    col = "red"
  )
  
  # Model parameters
  legenda <- c(
    paste0("Model: ", pars$model),
    paste0("Nugget: ", round(pars$nugget, 4)),
    paste0("Partial sill: ", round(pars$partial_sill, 4)),
    paste0("Total sill: ", round(pars$total_sill, 4)),
    paste0("Range: ", round(pars$range, 2), " m"),
    paste0("Sum of squared errors: ", signif(sserr, 4))
  )
  
  if (!is.na(rmse_cv)) {
    legenda <- c(
      legenda,
      paste0("Cross-validation RMSE: ", round(rmse_cv, 4))
    )
  }
  
  legend(
    "bottomright",
    legend = legenda,
    bty = "n",
    cex = 0.85
  )
  
  dev.off()
}

# ============================================================
# 4) READ SHAPEFILE AND CREATE INTERPOLATION GRID
# ============================================================

limite_sf <- st_read(shp_path, quiet = TRUE)

# If the shapefile has no CRS, assume WGS84
if (is.na(st_crs(limite_sf))) {
  st_crs(limite_sf) <- 4326
}

# Read first independent CSV to define UTM projection
dados_ref <- read.csv2(
  arquivos_csv[1],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(dados_ref) <- trimws(names(dados_ref))

dados_ref <- dados_ref |>
  mutate(
    longitude = to_num(longitude),
    latitude = to_num(latitude)
  )

epsg_utm <- calc_utm_epsg(
  lon = dados_ref$longitude,
  lat = dados_ref$latitude
)

cat("UTM EPSG used for kriging:", epsg_utm, "\n")

# Transform shapefile to UTM for kriging
limite_utm <- st_transform(limite_sf, epsg_utm)
limite_sp <- as(limite_utm, "Spatial")

# Also create geographic limit in WGS84 for final output masking
limite_geo <- st_transform(limite_sf, 4326)
limite_geo_sp <- as(limite_geo, "Spatial")

# Geographic CRS for final rasters
crs_geo <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Create interpolation grid in UTM
grid_base <- raster(
  extent(limite_sp),
  ncol = ncol_grid,
  nrow = nrow_grid
)

crs(grid_base) <- crs(limite_sp)

# Spatial grid used as newdata for kriging
grid_sp <- as(grid_base, "SpatialPixels")

# ============================================================
# 5) MAIN LOOP: INDEPENDENT FILES x PARAMETERS
# ============================================================

resumo_modelos <- data.frame()

for (arquivo in arquivos_csv) {
  
  nome_base <- file_path_sans_ext(basename(arquivo))
  dia_img <- extrair_dia_imagem(nome_base)
  
  cat("\n====================================================\n")
  cat("Processing independent file:", nome_base, "\n")
  cat("Image day:", dia_img, "\n")
  cat("====================================================\n")
  
  # ----------------------------------------------------------
  # Read independent CSV
  # ----------------------------------------------------------
  
  dados_df <- read.csv2(
    arquivo,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  names(dados_df) <- trimws(names(dados_df))
  
  # ----------------------------------------------------------
  # Convert coordinate columns to numeric
  # ----------------------------------------------------------
  
  dados_df <- dados_df |>
    mutate(
      longitude = to_num(longitude),
      latitude = to_num(latitude)
    )
  
  # ----------------------------------------------------------
  # Check required columns
  # ----------------------------------------------------------
  
  colunas_necessarias <- c("longitude", "latitude", parametros)
  faltantes <- setdiff(colunas_necessarias, names(dados_df))
  
  if (length(faltantes) > 0) {
    warning(
      paste(
        "Skipping file",
        nome_base,
        "because these columns are missing:",
        paste(faltantes, collapse = ", ")
      )
    )
    next
  }
  
  # ----------------------------------------------------------
  # Convert parameter columns to numeric
  # ----------------------------------------------------------
  
  for (p in parametros) {
    dados_df[[p]] <- to_num(dados_df[[p]])
  }
  
  # Remove rows with missing coordinates
  dados_df <- dados_df |>
    filter(
      !is.na(longitude),
      !is.na(latitude)
    )
  
  # ----------------------------------------------------------
  # Loop over parameters
  # ----------------------------------------------------------
  
  for (param in parametros) {
    
    cat("\nParameter:", param, "\n")
    
    dados_param <- dados_df
    
    # Remove NA values only for the selected parameter
    dados_param <- dados_param[!is.na(dados_param[[param]]), ]
    
    # Optional removal of zeros
    # For your case, remover_zeros = FALSE because zero values are real
    if (remover_zeros) {
      dados_param <- dados_param[dados_param[[param]] != 0, ]
    }
    
    if (nrow(dados_param) < 10) {
      warning(
        paste(
          "Too few valid points for",
          nome_base,
          param,
          "- skipped."
        )
      )
      next
    }
    
    # --------------------------------------------------------
    # Convert points to sf and then to UTM
    # --------------------------------------------------------
    
    dados_sf <- st_as_sf(
      dados_param,
      coords = c("longitude", "latitude"),
      crs = 4326,
      remove = FALSE
    )
    
    dados_utm <- st_transform(dados_sf, epsg_utm)
    dados_sp <- as(dados_utm, "Spatial")
    
    # --------------------------------------------------------
    # Fit and select the best semivariogram model
    # --------------------------------------------------------
    
    melhor <- ajustar_melhor_variograma(
      dados_sp = dados_sp,
      parametro = param,
      modelos_testar = modelos_testar,
      criterio = criterio_escolha
    )
    
    pars <- melhor$parameters
    
    # Safe parameter name for output files
    nome_param <- nome_parametro_arquivo(param)
    
    # Semivariogram title with mathematical symbols
    titulo_variograma <- titulo_variograma_simbolo(nome_base, param)
    
    # --------------------------------------------------------
    # Save semivariogram plot
    # --------------------------------------------------------
    
    arq_plot_vario <- file.path(
      dir_variog,
      paste0(nome_base, "_", nome_param, "_semivariogram.png")
    )
    
    salvar_plot_variograma(
      vario_emp = melhor$empirical_variogram,
      fit = melhor$fit,
      pars = pars,
      arquivo_saida = arq_plot_vario,
      titulo = titulo_variograma,
      sserr = melhor$SSErr,
      rmse_cv = melhor$RMSE_CV
    )
    
    # ========================================================
    # Ordinary kriging in UTM
    # ========================================================
    
    formula_krig <- as.formula(paste(param, "~ 1"))
    
    krig_output <- krige(
      formula_krig,
      dados_sp,
      newdata = grid_sp,
      model = melhor$fit
    )
    
    # Prediction raster in UTM
    raster_pred_utm <- raster(krig_output["var1.pred"])
    raster_pred_utm <- mask(raster_pred_utm, limite_sp)
    
    # Kriging variance raster in UTM
    raster_var_utm <- raster(krig_output["var1.var"])
    raster_var_utm <- mask(raster_var_utm, limite_sp)
    
    # ========================================================
    # Reproject final rasters to geographic coordinates - WGS84
    # ========================================================
    
    raster_pred <- projectRaster(
      raster_pred_utm,
      crs = crs_geo,
      method = "bilinear"
    )
    
    raster_pred <- mask(raster_pred, limite_geo_sp)
    
    raster_var <- projectRaster(
      raster_var_utm,
      crs = crs_geo,
      method = "bilinear"
    )
    
    raster_var <- mask(raster_var, limite_geo_sp)
    
    # ========================================================
    # Save final rasters in WGS84
    # ========================================================
    
    arq_raster_pred <- file.path(
      dir_rasters,
      paste0(nome_base, "_", nome_param, "_ordinary_kriging_WGS84.tif")
    )
    
    writeRaster(
      raster_pred,
      filename = arq_raster_pred,
      format = "GTiff",
      overwrite = TRUE
    )
    
    arq_raster_var <- file.path(
      dir_rasters,
      paste0(nome_base, "_", nome_param, "_kriging_variance_WGS84.tif")
    )
    
    writeRaster(
      raster_var,
      filename = arq_raster_var,
      format = "GTiff",
      overwrite = TRUE
    )
    
    # --------------------------------------------------------
    # Save model summary
    # --------------------------------------------------------
    
    linha_resumo <- data.frame(
      file = nome_base,
      image_day = dia_img,
      parameter = param,
      n_points = nrow(dados_param),
      selected_model = pars$model,
      nugget = pars$nugget,
      partial_sill = pars$partial_sill,
      total_sill = pars$total_sill,
      range_m = pars$range,
      kappa = pars$kappa,
      SSErr = melhor$SSErr,
      RMSE_CV = melhor$RMSE_CV,
      selection_criterion = criterio_escolha
    )
    
    resumo_modelos <- rbind(resumo_modelos, linha_resumo)
  }
}

# ============================================================
# 6) SAVE SUMMARY TABLE
# ============================================================

write.csv2(
  resumo_modelos,
  file.path(dir_saida, "semivariogram_model_summary.csv"),
  row.names = FALSE
)


#################################EXTRAÇÃO MATRIZ###################################
###################################################################################
###################################################################################
###################################################################################

library(terra)

rm(list = ls())
graphics.off()
gc()

# Selecionar os 4 rasters
arquivos <- c(
  "parametros_coo_30_alpha_0_ordinary_kriging_WGS84.tif",
  "parametros_coo_30_d_alpha_ordinary_kriging_WGS84.tif",
  "parametros_coo_30_df_alpha_ordinary_kriging_WGS84.tif")

# Ler e empilhar os rasters
r <- rast(arquivos)

# Opcional: renomear as bandas/colunas
names(r) <- c("alpha_0", "d_alpha", "df_alpha")

# Converter para data frame com coordenadas
df <- as.data.frame(
  r,
  xy = TRUE,
  na.rm = FALSE
)

# Renomear coordenadas
names(df)[1:2] <- c("longitude", "latitude")
head(df)


# Salvar planilha
write.csv2(
  df,
  "matriz.csv",
  row.names = FALSE
)

##################################CLUSTERIZING#####################################
###################################################################################
###################################################################################
###################################################################################

library(dplyr)
library(readr)
library(writexl)
library(factoextra)

# ------------------------------------------------------------
# 1. Ler a planilha
# ------------------------------------------------------------

dados <- read.csv2(
  "matriz.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Conferir estrutura
str(dados)

# ------------------------------------------------------------
# 2. Definir as variáveis usadas no agrupamento
# ------------------------------------------------------------

variaveis_cluster <- c("alpha_0", "d_alpha", "df_alpha")

# Criar um ID para manter a posição original das linhas
dados <- dados %>%
  mutate(id_linha = row_number())

# ------------------------------------------------------------
# 3. Opção para remover linhas sem informação
# ------------------------------------------------------------
# No seu arquivo existem muitas linhas com alpha_0 = 0, d_alpha = 0 e df_alpha = 0.
# Se esses zeros forem fundo, ausência de dado ou pixels sem informação, use TRUE.
# Se quiser agrupar tudo, inclusive zeros, use FALSE.

remover_zeros <- FALSE

if (remover_zeros == TRUE) {
  
  dados_cluster <- dados %>%
    filter(!(alpha_0 == 0 & d_alpha == 0 & df_alpha == 0))
  
} else {
  
  dados_cluster <- dados
  
}

# ------------------------------------------------------------
# 4. Remover linhas com NA apenas para a análise
# ------------------------------------------------------------

dados_cluster <- dados_cluster %>%
  filter(if_all(all_of(variaveis_cluster), ~ !is.na(.)))

# ------------------------------------------------------------
# 5. Padronizar as variáveis
# ------------------------------------------------------------
# Isso é importante porque cada variável pode ter escala diferente.

matriz_cluster <- dados_cluster %>%
  dplyr::select(dplyr::all_of(variaveis_cluster)) %>%
  scale()

# ------------------------------------------------------------
# 6. Calcular distância e agrupamento hierárquico
# ------------------------------------------------------------

distancia <- dist(matriz_cluster, method = "euclidean")

agrupamento <- hclust(
  distancia,
  method = "ward.D2"
)

# ------------------------------------------------------------
# 7. Visualizar dendrograma
# ------------------------------------------------------------
# Definir os grupos sem precisar plotar
grupos <- cutree(agrupamento, k = 2)

dados_cluster$grupo_hclust <- paste0("Grupo ", grupos)

write.csv2(
  dados_cluster,
  "matriz_com_grupos_hclust.csv",
  row.names = FALSE
)
# 
# writexl::write_xlsx(
#   dados_cluster,
#   "matriz_com_grupos_hclust.xlsx"
# )

graphics.off()

png(
  filename = "dendrograma_hclust.png",
  width = 2400,
  height = 1600,
  res = 300,
  type = "cairo"
)

par(family = "sans")

plot(
  agrupamento,
  labels = FALSE,
  hang = -1,
  main = "Dendrograma - Agrupamento Hierarquico",
  xlab = "Observacoes",
  ylab = "Distancia"
)

rect.hclust(
  agrupamento,
  k = 4,
  border = 2:5
)

dev.off()




# ------------------------------------------------------------
# 11. Visualizar distribuição dos grupos
# ------------------------------------------------------------

table(dados_cluster$grupo_hclust)

# ------------------------------------------------------------
# 12. Gráfico simples dos grupos no espaço geográfico
# ------------------------------------------------------------

plot(
  dados_cluster$longitude,
  dados_cluster$latitude,
  col = as.factor(dados_cluster$grupo_hclust),
  pch = 16,
  cex = 0.6,
  xlab = "Longitude",
  ylab = "Latitude",
  main = "Distribuição espacial dos grupos hierárquicos"
)

legend(
  "topright",
  legend = unique(dados_cluster$grupo_hclust),
  col = as.factor(unique(dados_cluster$grupo_hclust)),
  pch = 16,
  cex = 0.8
)

 
############################DESCRITIVA COMAPARAçÃO#################################
###################################################################################
###################################################################################
###################################################################################


# ============================================
# Estatísticas descritivas de TODOS os rasters
# Salva CSV2 com o nome do arquivo como ID
# ============================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
})

# --------- CONFIGURE AQUI ---------
dir_rasters <- getwd() 
saida_csv2  <- file.path(dir_rasters, "estatisticas_comparacao.csv")
# ---------------------------------

arquivos <- list.files(
  dir_rasters,
  pattern = "\\.(tif|tiff|grd)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(arquivos) == 0) stop("Nenhum raster .tif/.tiff/.grd encontrado no diretório.")

stats_raster <- function(r) {
  v <- values(r, mat = FALSE)
  v <- v[is.finite(v)]
  
  if (length(v) == 0) {
    return(list(
      n_valid = 0, mean = NA, median = NA, sd = NA, var = NA,
      min = NA, max = NA, q05 = NA, q25 = NA, q50 = NA, q75 = NA, q95 = NA
    ))
  }
  
  qs <- as.numeric(quantile(
    v, probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
    na.rm = TRUE, names = FALSE
  ))
  
  list(
    n_valid = length(v),
    mean   = mean(v),
    median = median(v),
    sd     = sd(v),
    var    = var(v),
    min    = min(v),
    max    = max(v),
    q05    = qs[1],
    q25    = qs[2],
    q50    = qs[3],
    q75    = qs[4],
    q95    = qs[5]
  )
}

df <- lapply(arquivos, function(f) {
  r <- tryCatch(terra::rast(f), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  
  # Se tiver várias bandas, calcula para cada banda separadamente
  out <- lapply(seq_len(nlyr(r)), function(i) {
    ri <- r[[i]]
    st <- stats_raster(ri)
    
    data.frame(
      arquivo   = basename(f),   # nome EXATO do raster
      caminho   = f,
      banda     = i,
      n_valid   = st$n_valid,
      media     = st$mean,
      mediana   = st$median,
      desvio    = st$sd,
      variancia = st$var,
      minimo    = st$min,
      maximo    = st$max,
      q05       = st$q05,
      q25       = st$q25,
      q50       = st$q50,
      q75       = st$q75,
      q95       = st$q95,
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(out)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(arquivo, banda)

readr::write_csv2(df, saida_csv2)

message("OK! CSV2 gerado em: ", saida_csv2)
print(df)


#################################################################################

# ============================================================
# Estabilidade entre escalas (N) para mapas interpolados
# Referência operacional: N = 40 (por parâmetro e data)
# Saídas: CSV2 + gráficos
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

# -------- CONFIGURE AQUI --------
dir_rasters <- getwd()
out_csv2    <- file.path(dir_rasters, "estabilidade_convergencia.csv")
out_dir_fig <- file.path(dir_rasters, "figs_estabilidade")
dir.create(out_dir_fig, showWarnings = FALSE)
# --------------------------------

# 1) listar rasters
files <- list.files(dir_rasters, pattern = "\\.(tif|tiff|grd)$",
                    full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) stop("Nenhum raster (.tif/.tiff/.grd) encontrado.")

# 2) parse do nome: parametro (0/d/df), data (8 dígitos), N (20/25/30/35/40 ou 15 default)
parse_name <- function(path) {
  nm <- tools::file_path_sans_ext(basename(path))
  
  pref <- str_match(nm, "^(df|d|0)_")[,2]
  parametro <- dplyr::case_when(
    pref == "0"  ~ "alpha_0",
    pref == "d"  ~ "d_alpha",
    pref == "df" ~ "df_alpha",
    TRUE ~ NA_character_
  )
  
  data_str <- str_match(nm, "(\\d{8})")[,2]  # primeira ocorrência YYYYMMDD
  
  n_str <- str_match(nm, "_(20|25|30|35|40)$")[,2]
  N <- ifelse(is.na(n_str), 15L, as.integer(n_str))
  
  tibble(
    arquivo = basename(path),
    caminho = path,
    parametro = parametro,
    data = data_str,
    N = N
  )
}

meta <- bind_rows(lapply(files, parse_name)) %>%
  filter(!is.na(parametro)) %>%
  filter(N %in% c(15, 20, 25, 30, 35, 40))

if (nrow(meta) == 0) stop("Nenhum arquivo bateu com o padrão 0_/d_/df_ e N em {15,20,25,30,35,40}.")

# 3) alinhar rasters ao grid da referência
align_to_ref <- function(x, ref) {
  if (!terra::compareGeom(x, ref, stopOnError = FALSE)) {
    x <- terra::resample(x, ref, method = "bilinear")
  }
  x
}

# 4) métricas entre mapas (x vs ref)
map_metrics <- function(x, ref) {
  v_ref <- terra::values(ref, mat = FALSE)
  v_x   <- terra::values(x,   mat = FALSE)
  
  ok <- is.finite(v_ref) & is.finite(v_x)
  if (!any(ok)) {
    return(tibble(
      n_valid = 0,
      cor_pearson = NA_real_,
      cor_spearman = NA_real_,
      mae = NA_real_,
      rmse = NA_real_,
      dq05 = NA_real_,
      dq50 = NA_real_,
      dq95 = NA_real_
    ))
  }
  
  a <- v_x[ok]
  b <- v_ref[ok]
  d <- a - b
  
  qx   <- as.numeric(quantile(a, probs = c(0.05, 0.50, 0.95), names = FALSE, na.rm = TRUE))
  qref <- as.numeric(quantile(b, probs = c(0.05, 0.50, 0.95), names = FALSE, na.rm = TRUE))
  
  tibble(
    n_valid = length(a),
    cor_pearson  = suppressWarnings(cor(a, b, method = "pearson")),
    cor_spearman = suppressWarnings(cor(a, b, method = "spearman")),
    mae  = mean(abs(d)),
    rmse = sqrt(mean(d^2)),
    dq05 = qx[1] - qref[1],
    dq50 = qx[2] - qref[2],
    dq95 = qx[3] - qref[3]
  )
}

# 5) processar por grupos (parametro, data) sem group_modify()
meta$grupo <- paste(meta$parametro, meta$data, sep = "||")
grupos <- split(meta, meta$grupo)

res_list <- lapply(grupos, function(g) {
  g_ref <- g[g$N == 40, , drop = FALSE]
  if (nrow(g_ref) == 0) return(NULL)
  
  ref <- terra::rast(g_ref$caminho[1])
  if (terra::nlyr(ref) > 1) ref <- ref[[1]]
  
  out <- lapply(seq_len(nrow(g)), function(i) {
    row <- g[i, ]
    
    x <- terra::rast(row$caminho)
    if (terra::nlyr(x) > 1) x <- x[[1]]
    x <- align_to_ref(x, ref)
    
    m <- map_metrics(x, ref)
    
    tibble(
      arquivo     = row$arquivo,
      caminho     = row$caminho,
      parametro   = row$parametro,
      data        = row$data,
      N           = row$N,
      ref_arquivo = g_ref$arquivo[1]
    ) %>% bind_cols(m)
  })
  
  bind_rows(out)
})

res <- bind_rows(res_list) %>%
  arrange(parametro, data, N)

if (nrow(res) == 0) stop("Nenhum grupo tinha N=40 como referência. Verifique seus arquivos _40.")

# 6) salvar CSV2
readr::write_csv2(res, out_csv2)
message("OK! CSV2 salvo em: ", out_csv2)

# 7) gráficos (métrica vs N), excluindo N=40 (é a referência)
plot_metric <- function(df, y, ylab, fname) {
  p <- ggplot(df %>% filter(N != 40), aes(x = N, y = .data[[y]])) +
    geom_line() +
    geom_point() +
    facet_grid(data ~ parametro, scales = "free_y") +
    labs(x = "N (partições por eixo)", y = ylab) +
    theme_bw()
  ggsave(file.path(out_dir_fig, fname), p, width = 10, height = 6, dpi = 200)
}

plot_metric(res, "cor_pearson",  "Correlação (Pearson) vs ref N=40",  "cor_pearson_vsN.png")
plot_metric(res, "cor_spearman", "Correlação (Spearman) vs ref N=40", "cor_spearman_vsN.png")
plot_metric(res, "rmse",         "RMSE (distância entre mapas) vs ref N=40", "rmse_vsN.png")
plot_metric(res, "mae",          "MAE (distância entre mapas) vs ref N=40",  "mae_vsN.png")
plot_metric(res, "dq50",         "Δmediana (q50) vs ref N=40", "dq50_vsN.png")
plot_metric(res, "dq05",         "Δq05 vs ref N=40",           "dq05_vsN.png")
plot_metric(res, "dq95",         "Δq95 vs ref N=40",           "dq95_vsN.png")

message("OK! Figuras salvas em: ", out_dir_fig)
