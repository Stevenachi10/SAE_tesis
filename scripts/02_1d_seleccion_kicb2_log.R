# ==============================================================================
# 02_1d_seleccion_kicb2_log.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.1d: Sensibilidad de la selección a la escala de la respuesta
#
# ------------------------------------------------------------------------------
# OBJETIVO
# ------------------------------------------------------------------------------
# M1 se seleccionó con KICb2 en la escala original y después se ajustó en
# escala logarítmica. Este script repite la selección con KICb2 directamente
# sobre la respuesta en escala logarítmica, con el mismo universo de
# candidatas, la misma búsqueda en ambas direcciones, el mismo filtro de VIF,
# B = 500 y semilla 2906. Responde a la pregunta: ¿se mantiene M1 si la
# selección se hace en la escala del modelo final?
#
# ------------------------------------------------------------------------------
# ESCALA LOGARITMICA
# ------------------------------------------------------------------------------
# Un modelo Fay-Herriot en escala logarítmica es un Fay-Herriot estándar sobre
# la respuesta transformada, con varianzas aproximadas por el método delta:
#
#     y_log   = log(pobreza_monetaria)
#     var_log = varianza_pobreza / pobreza_monetaria^2
#
# La retrotransformación no interviene en la selección, porque el criterio se
# calcula con la verosimilitud en la escala transformada. Antes de la búsqueda
# se verifica que este ajuste coincide con el de fh(transformation = "log").
#
# Los valores de KICb2 en escala logarítmica no son comparables con los de la
# escala original, porque la respuesta es distinta. Solo se comparan entre
# especificaciones ajustadas en la misma escala.
#
# ------------------------------------------------------------------------------
# ENTRADAS
# ------------------------------------------------------------------------------
#   output/matriz_sae_transformada_v2.rds
#   output/seleccion_stepwise_completo_both.rds   (para M1 y el universo)
#
# ------------------------------------------------------------------------------
# SALIDAS (output/)
# ------------------------------------------------------------------------------
#   seleccion_kicb2_log.rds
#   seleccion_kicb2_log_historial.csv
#   seleccion_kicb2_log_comparacion.csv
# ==============================================================================


# ==============================================================================
# PAQUETES
# ==============================================================================

library(emdi)
library(car)
library(dplyr)
library(here)

select <- dplyr::select
filter <- dplyr::filter


# ==============================================================================
# RUTAS
# ==============================================================================

ruta_out <- here("output")
matriz   <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))
corrida  <- readRDS(file.path(ruta_out, "seleccion_stepwise_completo_both.rds"))


# ==============================================================================
# CONFIGURACION (identica a 02_1c, salvo el criterio)
# ==============================================================================

CRITERIO         <- "KICb2"
DIRECCION        <- "both"
METODO_SELECCION <- "ml"
B_BOOT           <- 500
SEMILLA          <- 2906
TOL_MEJORA       <- 1e-6
VIF_MAX          <- 10
COR_MAX          <- 0.99

NO_COVARIABLES <- c(
  "cod_mun", "Municipio", "cod_dpto", "Departamento",
  "pobreza_monetaria", "varianza_pobreza",
  "ee_pobreza", "cv_directo", "n_eff", "n_muestra",
  "error_estandar", "cvlog"
)


# ==============================================================================
# UNIVERSO DE CANDIDATAS
# ------------------------------------------------------------------------------
# Se reconstruye con los mismos filtros de 02_1c y se exige que coincida con el
# universo guardado en la corrida original. Si no coincide, la comparacion no
# seria valida y el script se detiene.
# ==============================================================================

numericas  <- names(matriz)[sapply(matriz, is.numeric)]
CANDIDATAS <- setdiff(numericas, NO_COVARIABLES)

constantes <- CANDIDATAS[sapply(CANDIDATAS, function(v) {
  x <- matriz[[v]]
  length(unique(x[!is.na(x)])) < 2
})]
CANDIDATAS <- setdiff(CANDIDATAS, constantes)

con_na <- CANDIDATAS[sapply(CANDIDATAS, function(v) any(is.na(matriz[[v]])))]
CANDIDATAS <- setdiff(CANDIDATAS, con_na)

if (length(CANDIDATAS) > 1) {
  M  <- as.matrix(matriz[, CANDIDATAS, drop = FALSE])
  cm <- abs(cor(M, use = "pairwise.complete.obs"))
  cm[!upper.tri(cm)] <- 0
  redundantes <- character(0)
  for (j in seq_along(CANDIDATAS)) {
    if (CANDIDATAS[j] %in% redundantes) next
    pares <- setdiff(CANDIDATAS[which(cm[j, ] > COR_MAX)], redundantes)
    if (length(pares) > 0) redundantes <- c(redundantes, pares)
  }
  CANDIDATAS <- setdiff(CANDIDATAS, redundantes)
}

X   <- cbind(1, as.matrix(matriz[, CANDIDATAS, drop = FALSE]))
colnames(X)[1] <- "(Intercept)"
qrX <- qr(X)
if (qrX$rank < ncol(X)) {
  CANDIDATAS <- intersect(CANDIDATAS, colnames(X)[qrX$pivot[seq_len(qrX$rank)]])
}

if (!setequal(CANDIDATAS, corrida$candidatas)) {
  stop("El universo reconstruido no coincide con el de la corrida original.\n",
       "Solo aqui: ", paste(setdiff(CANDIDATAS, corrida$candidatas),
                            collapse = ", "), "\n",
       "Solo en la original: ", paste(setdiff(corrida$candidatas, CANDIDATAS),
                                      collapse = ", "))
}

# Mismo orden que en la corrida original, para que los empates se resuelvan
# igual.
CANDIDATAS <- corrida$candidatas

VARS_M1 <- strsplit(
  corrida$tabla$covariables[corrida$tabla$criterio == "KICb2"],
  " + ", fixed = TRUE)[[1]]

cat("\n============================================================\n")
cat("SELECCION KICb2 EN ESCALA LOGARITMICA\n")
cat("============================================================\n")
cat("Candidatas:  ", length(CANDIDATAS), " (igual a la corrida original)\n",
    sep = "")
cat("M1 original: ", paste(VARS_M1, collapse = ", "), "\n", sep = "")
cat("B = ", B_BOOT, " | semilla = ", SEMILLA, " | VIF < ", VIF_MAX, "\n",
    sep = "")
cat("============================================================\n")


# ==============================================================================
# DATOS EN ESCALA LOGARITMICA
# ==============================================================================

stopifnot(all(matriz$pobreza_monetaria > 0),
          all(matriz$varianza_pobreza > 0))

datos <- matriz %>%
  select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza,
         all_of(CANDIDATAS)) %>%
  mutate(y_log   = log(pobreza_monetaria),
         var_log = varianza_pobreza / pobreza_monetaria^2) %>%
  as.data.frame()


# ==============================================================================
# FUNCIONES
# ==============================================================================

formula_log <- function(vars) {
  if (length(vars) == 0) return(as.formula("y_log ~ 1"))
  as.formula(paste("y_log ~", paste(vars, collapse = " + ")))
}

formula_orig <- function(vars) {
  as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
}

escalar <- function(x) {
  v <- tryCatch(as.numeric(x), error = function(e) NA_real_)
  if (length(v) != 1 || !is.finite(v)) return(NA_real_)
  v
}

clave_de <- function(vars) {
  if (length(vars) == 0) return("<intercepto>")
  paste(sort(vars), collapse = "|")
}

# El VIF depende solo de las covariables, no de la respuesta.
vif_de <- function(vars) {
  if (length(vars) < 2) return(1)
  tryCatch({
    v <- car::vif(lm(formula_log(vars), data = datos))
    if (is.matrix(v)) v <- v[, ncol(v)]
    escalar(max(v, na.rm = TRUE))
  }, error = function(e) NA_real_)
}

coefs_de <- function(aj) {
  cf <- aj$model$coefficients
  if (is.data.frame(cf) || is.matrix(cf)) cf <- cf[, 1]
  as.numeric(cf)
}

ajustar_log <- function(vars, B = c(0, B_BOOT)) {
  tryCatch(
    withCallingHandlers(
      fh(formula_log(vars),
         vardir        = "var_log",
         combined_data = datos,
         domains       = "cod_mun",
         method        = METODO_SELECCION,
         MSE           = FALSE,
         B             = B,
         seed          = SEMILLA),
      warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) NULL)
}

kicb2_de <- function(aj) {
  if (is.null(aj)) return(NA_real_)
  escalar(aj$model$model_select[[CRITERIO]])
}


# ==============================================================================
# VERIFICACION DE EQUIVALENCIA CON fh(transformation = "log")
# ------------------------------------------------------------------------------
# Se ajusta M1 de las dos formas bajo ML. Si la varianza del efecto de area o
# los coeficientes difieren, la transformacion manual no reproduce la de emdi
# y el script se detiene.
# ==============================================================================

aj_emdi <- withCallingHandlers(
  fh(formula_orig(VARS_M1),
     vardir             = "varianza_pobreza",
     combined_data      = datos,
     domains            = "cod_mun",
     method             = METODO_SELECCION,
     transformation     = "log",
     backtransformation = "bc_crude",
     MSE                = FALSE,
     B                  = c(0, 0),
     seed               = SEMILLA),
  warning = function(w) invokeRestart("muffleWarning"))

aj_manual <- ajustar_log(VARS_M1, B = c(0, 0))

if (is.null(aj_manual)) stop("No se pudo ajustar M1 con la transformacion manual.")

dif_s2u <- abs(escalar(aj_emdi$model$variance) -
                 escalar(aj_manual$model$variance))
dif_bet <- max(abs(coefs_de(aj_emdi) - coefs_de(aj_manual)))

cat("\nVerificacion de equivalencia (M1, ML):\n")
cat("  Diferencia en sigma2_u:      ", signif(dif_s2u, 3), "\n", sep = "")
cat("  Diferencia maxima en betas:  ", signif(dif_bet, 3), "\n", sep = "")

if (!is.finite(dif_s2u) || !is.finite(dif_bet) ||
    dif_s2u > 1e-6 || dif_bet > 1e-6) {
  stop("La transformacion manual no coincide con fh(transformation = 'log'). ",
       "Revisar como transforma emdi las varianzas antes de continuar.")
}
cat("  Coinciden: la seleccion en escala logaritmica es valida.\n")


# ==============================================================================
# BUSQUEDA SECUENCIAL (misma logica que 02_1c, un solo criterio)
# ==============================================================================

memoria   <- new.env(parent = emptyenv())
n_ajustes <- 0L

evaluar <- function(vars) {
  k <- clave_de(vars)
  if (exists(k, envir = memoria, inherits = FALSE)) {
    return(get(k, envir = memoria, inherits = FALSE))
  }
  v <- kicb2_de(ajustar_log(vars))
  n_ajustes <<- n_ajustes + 1L
  assign(k, v, envir = memoria)
  v
}

seleccionadas   <- character(0)
criterio_actual <- evaluar(character(0))
historial       <- list()
movs_max        <- 3 * length(CANDIDATAS)

if (is.na(criterio_actual)) stop("No se pudo calcular KICb2 del modelo nulo.")

cat("\nPaso 0: KICb2 = ", round(criterio_actual, 4), "\n", sep = "")

t_inicio <- Sys.time()

for (paso in seq_len(movs_max)) {
  
  t_paso     <- Sys.time()
  candidatos <- list()
  
  # ---- Incorporaciones -------------------------------------------------------
  for (variable in setdiff(CANDIDATAS, seleccionadas)) {
    vars_c <- c(seleccionadas, variable)
    vm <- vif_de(vars_c)
    if (is.na(vm) || vm >= VIF_MAX) next
    v <- evaluar(vars_c)
    if (is.na(v)) next
    candidatos[[paste0("+", variable)]] <- list(vars = vars_c, valor = v,
                                                mov = "+", variable = variable)
  }
  
  # ---- Eliminaciones ---------------------------------------------------------
  if (DIRECCION == "both" && length(seleccionadas) >= 2) {
    for (variable in seleccionadas) {
      vars_c <- setdiff(seleccionadas, variable)
      v <- evaluar(vars_c)
      if (is.na(v)) next
      candidatos[[paste0("-", variable)]] <- list(vars = vars_c, valor = v,
                                                  mov = "-", variable = variable)
    }
  }
  
  if (length(candidatos) == 0) break
  
  valores <- sapply(candidatos, function(x) x$valor)
  mejor   <- candidatos[[which.min(valores)]]
  
  if ((criterio_actual - mejor$valor) <= TOL_MEJORA) {
    cat("\nNingun movimiento reduce el criterio. Fin de la busqueda.\n")
    break
  }
  
  mejora          <- criterio_actual - mejor$valor
  seleccionadas   <- mejor$vars
  criterio_actual <- mejor$valor
  
  historial[[length(historial) + 1]] <- data.frame(
    paso        = paso,
    movimiento  = mejor$mov,
    variable    = mejor$variable,
    n_vars      = length(seleccionadas),
    kicb2_log   = round(criterio_actual, 4),
    mejora      = round(mejora, 4),
    vif_max     = round(vif_de(seleccionadas), 2),
    tiempo_min  = round(as.numeric(difftime(Sys.time(), t_paso,
                                            units = "mins")), 2),
    covariables = paste(seleccionadas, collapse = " + "),
    stringsAsFactors = FALSE)
  
  cat(sprintf("Paso %2d: %s %-30s KICb2 = %.4f  (%.1f min)\n",
              paso, mejor$mov, mejor$variable, criterio_actual,
              as.numeric(difftime(Sys.time(), t_paso, units = "mins"))))
}

historial_df <- if (length(historial) > 0) do.call(rbind, historial) else
  data.frame()

VARS_LOG <- seleccionadas


# ==============================================================================
# COMPARACION CON M1
# ------------------------------------------------------------------------------
# 1. Composicion: variables comunes y propias de cada especificacion.
# 2. KICb2 en escala logaritmica de ambas (mismo B y semilla, comparables).
# 3. Ajuste final como M1 log de 02_2a: REML, bc_crude, MSE analitico.
# ==============================================================================

kicb2_m1_log  <- evaluar(VARS_M1)
kicb2_sel_log <- criterio_actual

ajustar_final <- function(vars) {
  withCallingHandlers(
    fh(formula_orig(vars),
       vardir             = "varianza_pobreza",
       combined_data      = datos,
       domains            = "cod_mun",
       method             = "reml",
       transformation     = "log",
       backtransformation = "bc_crude",
       MSE                = TRUE,
       mse_type           = "analytical",
       B                  = c(0, 0),
       seed               = SEMILLA),
    warning = function(w) invokeRestart("muffleWarning"))
}

resumen_final <- function(vars, etiqueta, kicb2) {
  aj  <- ajustar_final(vars)
  cvs <- sqrt(aj$MSE$FH) / aj$ind$FH * 100
  data.frame(
    especificacion = etiqueta,
    n_vars         = length(vars),
    kicb2_log      = round(kicb2, 4),
    vif_max        = round(vif_de(vars), 2),
    sigma2_u_log   = signif(escalar(aj$model$variance), 6),
    cv_medio       = round(mean(cvs), 3),
    cv_mediana     = round(median(cvs), 3),
    cv_max         = round(max(cvs), 3),
    covariables    = paste(vars, collapse = " + "),
    stringsAsFactors = FALSE)
}

comparacion <- rbind(
  resumen_final(VARS_M1,  "M1 (seleccion en escala original)", kicb2_m1_log),
  resumen_final(VARS_LOG, "Seleccion KICb2 en escala log",      kicb2_sel_log))

comunes   <- intersect(VARS_M1, VARS_LOG)
solo_m1   <- setdiff(VARS_M1, VARS_LOG)
solo_log  <- setdiff(VARS_LOG, VARS_M1)
identicas <- setequal(VARS_M1, VARS_LOG)

cat("\n\n============================================================\n")
cat("RESULTADO\n")
cat("============================================================\n")
cat("Seleccion en log (", length(VARS_LOG), "): ",
    paste(VARS_LOG, collapse = ", "), "\n", sep = "")
cat("M1 original     (", length(VARS_M1), "): ",
    paste(VARS_M1, collapse = ", "), "\n\n", sep = "")

if (identicas) {
  cat(">>> La seleccion en escala logaritmica reproduce M1.\n")
} else {
  cat("Comunes:    ", if (length(comunes)  == 0) "ninguna" else
    paste(comunes,  collapse = ", "), "\n", sep = "")
  cat("Solo en M1: ", if (length(solo_m1)  == 0) "ninguna" else
    paste(solo_m1,  collapse = ", "), "\n", sep = "")
  cat("Solo en log:", if (length(solo_log) == 0) " ninguna" else
    paste0(" ", paste(solo_log, collapse = ", ")), "\n", sep = "")
}

cat("\nComparacion en escala logaritmica (KICb2 comparable entre filas):\n\n")
print(comparacion %>% select(-covariables), row.names = FALSE)

cat("\nAjustes efectivos: ", n_ajustes, " | tiempo total: ",
    round(as.numeric(difftime(Sys.time(), t_inicio, units = "mins")), 1),
    " min\n", sep = "")


# ==============================================================================
# GUARDADO
# ==============================================================================

if (nrow(historial_df) > 0) {
  write.csv(historial_df,
            file.path(ruta_out, "seleccion_kicb2_log_historial.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}

write.csv(comparacion,
          file.path(ruta_out, "seleccion_kicb2_log_comparacion.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(
  list(vars_log       = VARS_LOG,
       vars_m1        = VARS_M1,
       identicas      = identicas,
       comunes        = comunes,
       solo_m1        = solo_m1,
       solo_log       = solo_log,
       kicb2_m1_log   = kicb2_m1_log,
       kicb2_sel_log  = kicb2_sel_log,
       comparacion    = comparacion,
       historial      = historial_df,
       verificacion   = c(dif_sigma2_u = dif_s2u, dif_betas = dif_bet),
       candidatas     = CANDIDATAS,
       B_boot         = B_BOOT,
       semilla        = SEMILLA,
       n_ajustes      = n_ajustes,
       session        = sessionInfo()),
  file.path(ruta_out, "seleccion_kicb2_log.rds"))

cat("\nGuardado en output/seleccion_kicb2_log*.\n")