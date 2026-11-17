# ==============================================================================
# 03_9_evaluacion_seleccion_log.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA III (anexo): Evaluación de la especificación seleccionada en escala
#                    logarítmica frente a M1 log
#
# ------------------------------------------------------------------------------
# OBJETIVO
# ------------------------------------------------------------------------------
# Aplicar a la especificación obtenida en 02_1d (selección KICb2 en escala
# logarítmica) la misma batería de la Etapa III que se aplicó a M1 log, y
# presentar ambas lado a lado:
#
#   1. Diagnóstico distribucional (summary()$normality de emdi), como 03_0.
#   2. Brown: regresión de la referencia sobre la estimación modelada,
#      H0 conjunta b0 = 0, b1 = 1, con covarianzas convencionales y HC3, y
#      Goldfeld-Quandt ordenado por psi, como 03_1 bloque 5b.
#   3. Formulación inversa (modelada sobre referencia), como 03_1 bloque 11b.
#   4. Estadístico W con m y m - 1 grados de libertad, como 03_1 bloque 6.
#   5. Concordancia de intervalos de Brown (recálculo manual de 03_1 bloque 8b)
#      y número de municipios cuya referencia cae en el intervalo al 95 %.
#   6. Precisión: CV medio, mediano y máximo, y municipios con CV menor que el
#      de la referencia.
#   7. Agregados departamentales frente a la cifra oficial, como 03_2.
#
# M1 log se toma del objeto ya guardado en modelado_variantes.rds, de modo que
# sus valores deben coincidir con los reportados en la tesis. El bloque final
# compara algunos de ellos con las cifras del documento como verificación.
#
# ------------------------------------------------------------------------------
# ENTRADAS
# ------------------------------------------------------------------------------
#   output/matriz_sae_transformada_v2.rds
#   output/modelado_variantes.rds
#   output/seleccion_kicb2_log.rds
#   data/raw/Agregados.xlsx
#
# ------------------------------------------------------------------------------
# SALIDAS (output/)
# ------------------------------------------------------------------------------
#   sens_escala_evaluacion.csv
#   sens_escala_agregados.csv
#   sens_escala_coeficientes.csv
#   sens_escala_evaluacion.rds
# ==============================================================================


# ==============================================================================
# PAQUETES
# ==============================================================================

library(emdi)
library(dplyr)
library(here)
library(car)
library(lmtest)
library(sandwich)
library(readxl)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")
SEMILLA  <- 2906
B_MSE    <- 500


# ==============================================================================
# INSUMOS
# ==============================================================================

matriz    <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))
variantes <- readRDS(file.path(ruta_out, "modelado_variantes.rds"))
sel_log   <- readRDS(file.path(ruta_out, "seleccion_kicb2_log.rds"))

datos <- variantes$datos
datos$cod_mun <- as.character(datos$cod_mun)

DIRECTO <- datos$pobreza_monetaria
PSI     <- datos$varianza_pobreza
n_dom   <- nrow(datos)

VARS_LOG <- sel_log$vars_log

cat("\n============================================================\n")
cat("EVALUACION DE LA SELECCION EN ESCALA LOGARITMICA\n")
cat("============================================================\n")
cat("Seleccion en log: ", paste(VARS_LOG, collapse = ", "), "\n", sep = "")
cat("M1:               ", paste(sel_log$vars_m1, collapse = ", "), "\n",
    sep = "")


# ==============================================================================
# AJUSTES
# ------------------------------------------------------------------------------
# M1 log: el objeto de 02_2a.
# Alternativa: mismos argumentos que M1 log en 02_2a (REML, log, bc_crude,
# MSE analitico, B = c(500, 0), semilla 2906).
# ==============================================================================

aj_m1 <- variantes$ajustes[["M1_reml_log"]]
if (is.null(aj_m1)) stop("No se encontro M1_reml_log en modelado_variantes.rds.")

datos_alt <- matriz %>%
  as.data.frame() %>%
  mutate(cod_mun = as.character(cod_mun)) %>%
  select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza,
         all_of(VARS_LOG))

datos_alt <- datos_alt[match(datos$cod_mun, datos_alt$cod_mun), ]

stopifnot(
  "Dominios distintos entre la matriz y los datos de 02_2a" =
    identical(datos_alt$cod_mun, datos$cod_mun),
  "La referencia difiere entre fuentes" =
    isTRUE(all.equal(datos_alt$pobreza_monetaria, DIRECTO)),
  "La varianza difiere entre fuentes" =
    isTRUE(all.equal(datos_alt$varianza_pobreza, PSI)))

aj_alt <- withCallingHandlers(
  fh(as.formula(paste("pobreza_monetaria ~", paste(VARS_LOG, collapse = " + "))),
     vardir             = "varianza_pobreza",
     combined_data      = datos_alt,
     domains            = "cod_mun",
     method             = "reml",
     transformation     = "log",
     backtransformation = "bc_crude",
     MSE                = TRUE,
     mse_type           = "analytical",
     B                  = c(B_MSE, 0),
     seed               = SEMILLA),
  warning = function(w) invokeRestart("muffleWarning"))

AJUSTES <- list("M1 log" = aj_m1, "Seleccion log" = aj_alt)


# ==============================================================================
# FUNCIONES
# ==============================================================================

escalar <- function(x) {
  v <- suppressWarnings(tryCatch(as.numeric(x), error = function(e) NA_real_))
  if (length(v) < 1 || !is.finite(v[1])) NA_real_ else v[1]
}

# Estimacion y MSE alineados con datos$cod_mun
extraer <- function(aj) {
  dom <- as.character(aj$ind$Domain)
  idx <- match(datos$cod_mun, dom)
  if (anyNA(idx)) stop("Dominios del ajuste no alineables con los datos.")
  list(est = as.numeric(aj$ind$FH)[idx],
       mse = as.numeric(aj$MSE$FH)[idx])
}

# summary()$normality, como en 03_0
via_summary <- function(aj) {
  m <- tryCatch(summary(aj)$normality, error = function(e) NULL)
  if (is.null(m)) return(NULL)
  col_de <- function(patron) {
    i <- grep(patron, colnames(m), ignore.case = TRUE)
    if (length(i) == 0) NA_character_ else colnames(m)[i[1]]
  }
  c_skew <- col_de("skew"); c_kurt <- col_de("kurt")
  c_p    <- if ("Shapiro_p" %in% colnames(m)) "Shapiro_p" else
    col_de("_p$|p.?value")
  f_ef  <- grep("random", rownames(m), ignore.case = TRUE)
  f_res <- grep("resid",  rownames(m), ignore.case = TRUE)
  saca  <- function(f, c) if (length(f) == 0 || is.na(c)) NA_real_ else
    escalar(m[f[1], c])
  list(res_asim = saca(f_res, c_skew), res_curt = saca(f_res, c_kurt),
       res_p    = saca(f_res, c_p),
       ef_asim  = saca(f_ef,  c_skew), ef_curt  = saca(f_ef,  c_kurt),
       ef_p     = saca(f_ef,  c_p))
}

# Contraste conjunto b0 = 0, b1 = 1 (convencional y HC3) y Goldfeld-Quandt
contraste <- function(y, x) {
  fit <- lm(y ~ x)
  p   <- car::linearHypothesis(fit, diag(2), c(0, 1))[["Pr(>F)"]][2]
  p3  <- car::linearHypothesis(fit, diag(2), c(0, 1),
                               vcov. = sandwich::vcovHC(fit, type = "HC3")
  )[["Pr(>F)"]][2]
  gq  <- tryCatch(lmtest::gqtest(fit, order.by = PSI,
                                 alternative = "two.sided")$p.value,
                  error = function(e) NA_real_)
  c(b0 = unname(coef(fit)[1]), b1 = unname(coef(fit)[2]),
    p = p, p_hc3 = p3, gq = gq)
}


# ==============================================================================
# CIFRAS DEPARTAMENTALES OFICIALES (como 03_2)
# ==============================================================================

normalizar <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- chartr("áéíóúüñ", "aeiouun", x)
  gsub("[[:space:]]+", " ", x)
}

crudo <- read_excel(here("data", "raw", "Agregados.xlsx"), sheet = 1,
                    col_names = TRUE, .name_repair = "minimal")
stopifnot(ncol(crudo) == 6)
names(crudo) <- c("departamento_fuente", "estimacion", "error_estandar",
                  "limite_inferior", "limite_superior", "cv_pct")

oficial <- crudo %>%
  mutate(clave = normalizar(departamento_fuente),
         across(c(estimacion, error_estandar, limite_inferior,
                  limite_superior), ~ suppressWarnings(as.numeric(.x)) / 100)) %>%
  filter(clave %in% c("cauca", "valle del cauca")) %>%
  mutate(cod_dpto = ifelse(clave == "cauca", "19", "76")) %>%
  select(cod_dpto, estimacion, error_estandar, limite_inferior,
         limite_superior)

stopifnot("Faltan departamentos en Agregados.xlsx" = nrow(oficial) == 2)

PESO <- exp(matriz$log_pob_total[match(datos$cod_mun,
                                       as.character(matriz$cod_mun))])
DPTO <- substr(datos$cod_mun, 1, 2)

stopifnot(all(is.finite(PESO)), all(PESO > 0),
          all(DPTO %in% c("19", "76")))


# ==============================================================================
# EVALUACION
# ==============================================================================

z95     <- qnorm(0.975)
cv_ref  <- sqrt(PSI) / DIRECTO * 100
filas   <- list()
agrs    <- list()

for (nm in names(AJUSTES)) {
  
  aj <- AJUSTES[[nm]]
  ex <- extraer(aj)
  est <- ex$est; mse <- ex$mse
  
  # 1. Normalidad
  nor <- via_summary(aj)
  
  # 2-3. Brown y formulacion inversa
  br  <- contraste(DIRECTO, est)
  inv <- contraste(est, DIRECTO)
  
  # 4. W
  W <- sum((DIRECTO - est)^2 / (PSI + mse))
  
  # 5. Intervalos (recalculo de 03_1, bloque 8b)
  sy <- sqrt(PSI); sx <- sqrt(mse)
  incluido <- (DIRECTO > est - z95 * sx) & (DIRECTO < est + z95 * sx)
  zb <- z95 * sqrt(sx^2 + sy^2) / (sx + sy)
  solapa <- (est - zb * sx < DIRECTO + zb * sy) &
    (DIRECTO - zb * sy < est + zb * sx)
  bt <- binom.test(sum(solapa), length(solapa), p = 0.95)
  
  # 6. Precision
  cv <- sqrt(mse) / est * 100
  
  # 7. Agregados
  ag <- data.frame(cod_dpto = DPTO, est = est, peso = PESO) %>%
    group_by(cod_dpto) %>%
    summarise(agregado = weighted.mean(est, peso), .groups = "drop") %>%
    left_join(oficial, by = "cod_dpto") %>%
    mutate(especificacion = nm,
           brecha_pp = (agregado - estimacion) * 100,
           brecha_ee = (agregado - estimacion) / error_estandar,
           dentro_ic = agregado >= limite_inferior &
             agregado <= limite_superior) %>%
    as.data.frame()
  agrs[[nm]] <- ag
  
  filas[[nm]] <- data.frame(
    especificacion      = nm,
    res_asimetria       = round(nor$res_asim, 3),
    res_curtosis        = round(nor$res_curt, 3),
    res_shapiro_p       = round(nor$res_p, 4),
    ef_asimetria        = round(nor$ef_asim, 3),
    ef_curtosis         = round(nor$ef_curt, 3),
    ef_shapiro_p        = round(nor$ef_p, 4),
    brown_b0            = round(br[["b0"]], 4),
    brown_b1            = round(br[["b1"]], 4),
    brown_p             = round(br[["p"]], 4),
    brown_p_hc3         = round(br[["p_hc3"]], 4),
    brown_gq_p          = round(br[["gq"]], 4),
    inversa_b1          = round(inv[["b1"]], 4),
    inversa_p           = round(inv[["p"]], 4),
    inversa_p_hc3       = round(inv[["p_hc3"]], 4),
    W                   = round(W, 2),
    W_crit_m            = round(qchisq(0.95, n_dom), 2),
    W_crit_m1           = round(qchisq(0.95, n_dom - 1), 2),
    solapados_brown     = sum(solapa),
    intervalos_p        = round(bt$p.value, 4),
    referencia_en_ic95  = sum(incluido),
    cv_medio            = round(mean(cv), 2),
    cv_mediana          = round(median(cv), 2),
    cv_max              = round(max(cv), 2),
    cv_menor_referencia = sum(cv < cv_ref),
    agregado_cauca      = round(ag$agregado[ag$cod_dpto == "19"] * 100, 2),
    agregado_valle      = round(ag$agregado[ag$cod_dpto == "76"] * 100, 2),
    brecha_cauca_pp     = round(ag$brecha_pp[ag$cod_dpto == "19"], 2),
    brecha_valle_pp     = round(ag$brecha_pp[ag$cod_dpto == "76"], 2),
    brecha_media_pp     = round(mean(abs(ag$brecha_pp)), 2),
    agregados_en_ic     = sum(ag$dentro_ic),
    stringsAsFactors = FALSE)
}

evaluacion <- do.call(rbind, filas)
agregados  <- do.call(rbind, agrs)


# ==============================================================================
# COEFICIENTES DE LA ALTERNATIVA (signos, caracter descriptivo)
# ==============================================================================

coefs <- function(aj, nm) {
  cf <- as.data.frame(aj$model$coefficients)
  cf$termino <- rownames(cf)
  cf$especificacion <- nm
  cf
}

coeficientes <- rbind(coefs(aj_m1, "M1 log"), coefs(aj_alt, "Seleccion log"))


# ==============================================================================
# RESULTADOS
# ==============================================================================

cat("\n\n============================================================\n")
cat("COMPARACION (filas = indicadores)\n")
cat("============================================================\n\n")

tabla <- as.data.frame(t(evaluacion[, -1]))
names(tabla) <- evaluacion$especificacion
print(tabla)

cat("\n\n-- Agregados departamentales --\n\n")
print(agregados %>%
        mutate(across(c(agregado, estimacion, limite_inferior,
                        limite_superior), ~ round(.x * 100, 2)),
               across(c(brecha_pp, brecha_ee), ~ round(.x, 2))) %>%
        select(especificacion, cod_dpto, agregado, estimacion,
               limite_inferior, limite_superior, brecha_pp, brecha_ee,
               dentro_ic),
      row.names = FALSE)

cat("\n\n-- Coeficientes (escala log, caracter descriptivo) --\n\n")
print(coeficientes, row.names = FALSE)


# ==============================================================================
# VERIFICACION CONTRA LAS CIFRAS DE LA TESIS (M1 log)
# ------------------------------------------------------------------------------
# Si alguna fila sale con diferencia, revisar antes de usar la tabla.
# ==============================================================================

m1 <- evaluacion[evaluacion$especificacion == "M1 log", ]

verif <- data.frame(
  indicador = c("CV medio", "CV maximo", "W", "Referencia en IC95",
                "Agregado Cauca", "Agregado Valle", "Brecha media",
                "Inversa p", "Inversa p HC3"),
  tesis     = c(11.57, 14.58, 22.83, 83, 42.59, 25.65, 0.28, 0.275, 0.221),
  script    = c(m1$cv_medio, m1$cv_max, m1$W, m1$referencia_en_ic95,
                m1$agregado_cauca, m1$agregado_valle, m1$brecha_media_pp,
                round(m1$inversa_p, 3), round(m1$inversa_p_hc3, 3)))
verif$coincide <- abs(verif$tesis - verif$script) <= 0.011

cat("\n\n-- Verificacion de M1 log contra la tesis --\n\n")
print(verif, row.names = FALSE)
cat(if (all(verif$coincide)) "\nTodo coincide.\n" else
  "\nAVISO: hay diferencias con la tesis. Revisar antes de usar la tabla.\n")


# ==============================================================================
# GUARDADO
# ==============================================================================

write.csv(evaluacion, file.path(ruta_out, "sens_escala_evaluacion.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(agregados, file.path(ruta_out, "sens_escala_agregados.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(coeficientes, file.path(ruta_out, "sens_escala_coeficientes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(evaluacion = evaluacion, agregados = agregados,
             coeficientes = coeficientes, verificacion = verif,
             vars_log = VARS_LOG, ajuste_alternativa = aj_alt,
             session = sessionInfo()),
        file.path(ruta_out, "sens_escala_evaluacion.rds"))

cat("\nGuardado en output/sens_escala_*.\n")