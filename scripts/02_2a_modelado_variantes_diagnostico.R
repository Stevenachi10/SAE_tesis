# ==============================================================================
# 02_2a_modelado_variantes_diagnostico.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.2a — Variantes del modelo Fay-Herriot (DIAGNOSTICO)
#
# Requiere output/modelado_variantes.rds, producido por
# 02_2a_modelado_variantes.R.
#
# ------------------------------------------------------------------------------
# CRITERIOS
# ------------------------------------------------------------------------------
# COEFICIENTE DE VARIACION. Se calcula únicamente sobre dominios con estimación
# admisible: predicción finita en el intervalo (0, 1] y error cuadrático medio
# finito y no negativo. Los dominios excluidos se cuentan por causa y se
# reportan.
#
# El conteo de predicciones fuera del intervalo unitario es un diagnóstico de
# validez de la predicción, no una medida de precisión. Ambos se reportan por
# separado y no deben leerse como si midieran lo mismo.
#
# COMPARACION ENTRE ESCALAS. Solo el predictor y el error cuadrático medio son
# comparables entre escalas, porque fh() los devuelve retrotransformados. La
# componente de varianza y el factor de contracción no lo son. Además, la
# escala logit emplea un estimador del error cuadrático medio distinto del de
# las otras dos, de modo que esa comparación enfrenta estimadores distintos.
#
# NORMALIDAD. Los valores p de Shapiro-Wilk se reportan como diagnóstico
# complementario y no constituyen criterio de descarte. El efecto aleatorio de
# área es una cantidad latente: su predicción está contraída y su varianza
# depende del dominio, de modo que no constituye una muestra directa de los
# efectos aleatorios del modelo.
#
# FACTOR DE CONTRACCION. Se toma del objeto ajustado cuando está disponible. La
# reconstrucción a partir de sigma2_u y de la varianza de muestreo solo es
# válida en la escala original: bajo una transformación ambas cantidades viven
# en escalas distintas y el cociente carece de sentido. Se registra la fuente
# empleada en cada caso.
#
# ------------------------------------------------------------------------------
# SALIDAS (output/)
# ------------------------------------------------------------------------------
#   mod_cv.csv               Precisión y validez de la predicción.
#   mod_normalidad.csv       Pruebas de normalidad.
#   mod_gamma.csv            Factor de contracción, escala original.
#   modelado_diagnostico.rds
# ==============================================================================

library(emdi)
library(dplyr)
library(here)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")

mv     <- readRDS(file.path(ruta_out, "modelado_variantes.rds"))
matriz <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

datos            <- mv$datos
ESPECIFICACIONES <- mv$especificaciones
registro         <- mv$registro
n_dom            <- nrow(datos)

stopifnot(
  "El orden de dominios guardado no coincide con la matriz de entrada" =
    identical(as.character(datos$cod_mun), as.character(matriz$cod_mun)))

# Umbral descriptivo para el conteo de dominios. No constituye criterio formal
# de publicacion mientras no se adopte y se cite una referencia.
CV_UMBRAL <- 20

escalar <- function(x) {
  v <- tryCatch(as.numeric(x), error = function(e) NA_real_)
  if (length(v) != 1 || !is.finite(v)) return(NA_real_)
  v
}


# ==============================================================================
# FUNCIONES
# ==============================================================================

# ------------------------------------------------------------------------------
# Extraccion de predicciones y errores cuadraticos medios
# ------------------------------------------------------------------------------
extraer <- function(aj, etiqueta = "") {
  
  col_de <- function(df, patron) {
    if (is.null(df)) return(NULL)
    i <- grep(patron, names(df))
    if (length(i) == 0) return(NULL)
    as.numeric(df[[i[1]]])
  }
  
  dom <- tryCatch(as.character(aj$ind$Domain), error = function(e) NULL)
  
  if (is.null(dom) || !identical(dom, as.character(datos$cod_mun))) {
    stop("El orden de dominios del ajuste ", etiqueta,
         " no coincide con 'datos'.")
  }
  
  pred <- col_de(aj$ind, "^FH$")
  
  if (is.null(pred) || length(pred) != n_dom) {
    stop("No se pudo extraer la prediccion del ajuste ", etiqueta, ".")
  }
  
  mse <- col_de(aj$MSE, "^FH$")
  if (is.null(mse) || length(mse) != n_dom) mse <- rep(NA_real_, n_dom)
  
  list(pred = pred, mse = mse)
}

# ------------------------------------------------------------------------------
# Coeficiente de variacion con exclusion explicita
# ------------------------------------------------------------------------------
cv_de <- function(pred, mse) {
  
  admisible <- is.finite(pred) & pred > 0 & pred <= 1 &
    is.finite(mse) & mse >= 0
  
  cv <- rep(NA_real_, length(pred))
  cv[admisible] <- sqrt(mse[admisible]) / pred[admisible] * 100
  
  list(cv = cv,
       n_admisibles   = sum(admisible),
       n_pred_na      = sum(!is.finite(pred)),
       n_no_positiva  = sum(is.finite(pred) & pred <= 0),
       n_fuera_uno    = sum(is.finite(pred) & pred > 1),
       n_mse_invalido = sum(!is.finite(mse) | (is.finite(mse) & mse < 0)))
}

# ------------------------------------------------------------------------------
# Pruebas de normalidad
# ------------------------------------------------------------------------------
# Se exige la columna del valor p de manera explicita: un patron que solo
# busque "shapiro" devuelve primero el estadistico, cuyo rango habitual entre
# 0,9 y 1 puede confundirse con un valor p no significativo.
extraer_normalidad <- function(aj) {
  
  m <- tryCatch(summary(aj)$normality, error = function(e) NULL)
  
  if (is.null(m) || is.null(rownames(m)) || is.null(colnames(m))) {
    return(c(efectos = NA_real_, residuos = NA_real_))
  }
  
  col <- if ("Shapiro_p" %in% colnames(m)) {
    "Shapiro_p"
  } else {
    cand <- grep("_p$|p_value|pvalue|p\\.value", colnames(m),
                 ignore.case = TRUE, value = TRUE)
    if (length(cand) > 0) cand[1] else NA_character_
  }
  
  if (is.na(col)) {
    warning("No se identifico la columna del valor p en la matriz de ",
            "normalidad. Columnas: ", paste(colnames(m), collapse = ", "))
    return(c(efectos = NA_real_, residuos = NA_real_))
  }
  
  f_ef  <- grep("random", rownames(m), ignore.case = TRUE)
  f_res <- grep("resid",  rownames(m), ignore.case = TRUE)
  
  c(efectos  = if (length(f_ef)  > 0) escalar(m[f_ef[1],  col]) else NA_real_,
    residuos = if (length(f_res) > 0) escalar(m[f_res[1], col]) else NA_real_)
}

# ------------------------------------------------------------------------------
# Factor de contraccion
# ------------------------------------------------------------------------------
extraer_gamma <- function(aj, escala) {
  
  g <- tryCatch(aj$model$gamma, error = function(e) NULL)
  
  if (!is.null(g)) {
    if (is.data.frame(g)) {
      col <- grep("gamma", names(g), ignore.case = TRUE, value = TRUE)
      if (length(col) > 0 && nrow(g) == n_dom) {
        return(list(valores = as.numeric(g[[col[1]]]), fuente = "objeto"))
      }
    } else if (is.numeric(g) && length(g) == n_dom) {
      return(list(valores = as.numeric(g), fuente = "objeto"))
    }
  }
  
  if (escala == "original") {
    s2u <- escalar(tryCatch(aj$model$variance, error = function(e) NA))
    if (!is.na(s2u)) {
      return(list(valores = s2u / (s2u + datos$varianza_pobreza),
                  fuente  = "reconstruido"))
    }
  }
  
  list(valores = rep(NA_real_, n_dom), fuente = "no disponible")
}


# ==============================================================================
# PARTE 1 — PRECISION Y VALIDEZ DE LA PREDICCION
# ==============================================================================

cat("\n============================================================\n")
cat("PRECISION Y VALIDEZ DE LA PREDICCION\n")
cat("============================================================\n")

filas <- list(); filas_n <- list(); filas_g <- list(); cv_lista <- list()

estimados <- registro %>% filter(estimado)

for (i in seq_len(nrow(estimados))) {
  
  r   <- estimados[i, ]
  aj  <- mv$ajustes[[r$modelo]]
  if (is.null(aj)) next
  
  ex <- extraer(aj, r$modelo)
  cc <- cv_de(ex$pred, ex$mse)
  v  <- cc$cv[is.finite(cc$cv)]
  
  cv_lista[[r$modelo]] <- cc$cv
  
  filas[[length(filas) + 1]] <- data.frame(
    modelo = r$modelo, especificacion = r$especificacion,
    metodo = r$metodo, escala = r$escala, papel = r$papel,
    mse_type = r$mse_type,
    n = length(v),
    cv_medio   = if (length(v)) round(mean(v), 3)   else NA_real_,
    cv_mediana = if (length(v)) round(median(v), 3) else NA_real_,
    cv_max     = if (length(v)) round(max(v), 3)    else NA_real_,
    n_cv_mayor = if (length(v)) sum(v > CV_UMBRAL)  else NA_integer_,
    n_excluidos    = n_dom - cc$n_admisibles,
    n_pred_na      = cc$n_pred_na,
    n_no_positiva  = cc$n_no_positiva,
    n_fuera_uno    = cc$n_fuera_uno,
    n_mse_invalido = cc$n_mse_invalido,
    stringsAsFactors = FALSE)
  
  sh <- extraer_normalidad(aj)
  
  filas_n[[length(filas_n) + 1]] <- data.frame(
    modelo = r$modelo, especificacion = r$especificacion,
    metodo = r$metodo, escala = r$escala,
    shapiro_ef_p  = signif(unname(sh["efectos"]),  4),
    shapiro_res_p = signif(unname(sh["residuos"]), 4),
    stringsAsFactors = FALSE)
  
  if (r$escala == "original") {
    gm <- extraer_gamma(aj, r$escala)
    filas_g[[length(filas_g) + 1]] <- data.frame(
      modelo = r$modelo, especificacion = r$especificacion,
      metodo = r$metodo, fuente = gm$fuente,
      gamma_medio = round(mean(gm$valores, na.rm = TRUE), 4),
      gamma_min   = round(min(gm$valores,  na.rm = TRUE), 4),
      gamma_max   = round(max(gm$valores,  na.rm = TRUE), 4),
      stringsAsFactors = FALSE)
  }
}

tabla_cv   <- do.call(rbind, filas)
normalidad <- do.call(rbind, filas_n)
gammas     <- do.call(rbind, filas_g)

cat("\n-- Escala original, por metodo --\n\n")
print(tabla_cv %>%
        filter(escala == "original") %>%
        select(especificacion, metodo, papel, cv_medio, cv_mediana, cv_max,
               n_cv_mayor, n_fuera_uno, n_excluidos) %>%
        arrange(especificacion, metodo),
      row.names = FALSE)

cat("\n-- Comparacion entre escalas (metodo principal) --\n")
cat("Advertencia: la escala logit emplea un estimador del MSE distinto.\n\n")
print(tabla_cv %>%
        filter(metodo == mv$metodo_principal) %>%
        select(especificacion, escala, mse_type, cv_medio, cv_mediana,
               cv_max, n_cv_mayor, n_fuera_uno, n_excluidos) %>%
        arrange(especificacion, escala),
      row.names = FALSE)

if (any(tabla_cv$n_excluidos > 0)) {
  cat("\n-- Dominios excluidos del calculo, por causa --\n")
  print(tabla_cv %>%
          filter(n_excluidos > 0) %>%
          select(modelo, n_excluidos, n_pred_na, n_no_positiva,
                 n_fuera_uno, n_mse_invalido),
        row.names = FALSE)
}

write.csv(tabla_cv, file.path(ruta_out, "mod_cv.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# PARTE 2 — FACTOR DE CONTRACCION
# ==============================================================================

if (!is.null(gammas) && nrow(gammas) > 0) {
  
  cat("\n\n============================================================\n")
  cat("FACTOR DE CONTRACCION (ESCALA ORIGINAL)\n")
  cat("============================================================\n\n")
  
  print(gammas %>% arrange(especificacion, metodo), row.names = FALSE)
  
  if (any(gammas$fuente != "objeto")) {
    cat("\nAviso: hay factores reconstruidos a partir de sigma2_u.\n")
  }
  
  write.csv(gammas, file.path(ruta_out, "mod_gamma.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}


# ==============================================================================
# PARTE 3 — NORMALIDAD
# ==============================================================================

cat("\n\n============================================================\n")
cat("PRUEBAS DE NORMALIDAD (DIAGNOSTICO COMPLEMENTARIO)\n")
cat("============================================================\n")
cat("Valores p de Shapiro-Wilk. No constituyen criterio de descarte: el\n")
cat("efecto aleatorio de area es una cantidad latente y su prediccion no\n")
cat("es una muestra directa de los efectos aleatorios del modelo.\n\n")

print(normalidad %>% arrange(especificacion, escala, metodo),
      row.names = FALSE)

write.csv(normalidad, file.path(ruta_out, "mod_normalidad.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# GUARDADO
# ==============================================================================

saveRDS(
  list(
    tabla_cv       = tabla_cv,
    normalidad     = normalidad,
    gammas         = gammas,
    cv_por_dominio = cv_lista,
    cv_umbral      = CV_UMBRAL,
    session        = sessionInfo()
  ),
  file.path(ruta_out, "modelado_diagnostico.rds"))

cat("\n============================================================\n")
cat("DIAGNOSTICO COMPLETADO\n")
cat("============================================================\n")