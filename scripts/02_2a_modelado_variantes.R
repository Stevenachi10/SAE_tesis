# ==============================================================================
# 02_2a_modelado_variantes.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.2a — Estimación de variantes del modelo Fay-Herriot (ESTIMACION)
#
# Este script solo estima y guarda. El cálculo de coeficientes de variación, la
# comparación entre métodos y escalas y las pruebas de normalidad corresponden
# a 02_2a_modelado_variantes_diagnostico.R.
#
# ------------------------------------------------------------------------------
# OBJETIVO
# ------------------------------------------------------------------------------
# Ajustar las tres especificaciones seleccionadas en la subetapa de selección de
# variables bajo distintas variantes de estimación, sin decidir entre ellas. La
# comparación de supuestos, la evaluación del error cuadrático medio y la
# selección del modelo final corresponden a la Etapa III.
#
#   Bloque A — escala original, variando el método de estimación de la
#              componente de varianza.
#   Bloque B — transformaciones logarítmica y logit bajo verosimilitud
#              restringida.
#
# ------------------------------------------------------------------------------
# METODO PRINCIPAL Y ANALISIS DE SENSIBILIDAD
# ------------------------------------------------------------------------------
# La componente de varianza se estima por máxima verosimilitud restringida
# (REML), que constituye el método principal. Los procedimientos de
# verosimilitud ajustada se estiman como análisis de sensibilidad al método, no
# como especificaciones candidatas.
#
# Esos procedimientos incorporan un factor de ajuste a la función de
# verosimilitud que garantiza una solución estrictamente positiva para la
# componente de varianza, situación relevante cuando la estimación se aproxima
# a cero y el predictor colapsa al estimador sintético. Li y Lahiri (2010)
# introducen ese factor; Yoshimori y Lahiri (2014) proponen variantes que
# modifican sus propiedades. Verificar la formulación exacta en las fuentes
# antes de citarlas en el documento.
#
# ------------------------------------------------------------------------------
# ERROR CUADRATICO MEDIO
# ------------------------------------------------------------------------------
# El tipo de estimador se fija por escala y no se sustituye de manera
# automática: si el ajuste falla, se registra el error y la variante queda como
# no estimada. Un cambio silencioso de procedimiento produciría coeficientes de
# variación no comparables entre variantes.
#
# La escala logit no admite el estimador analítico, de modo que la comparación
# entre las dos transformaciones enfrenta estimadores distintos y debe
# declararse como tal.
#
# ------------------------------------------------------------------------------
# SOBRE LA TRANSFORMACION ARCOSENO
# ------------------------------------------------------------------------------
# No se evalúa. Requiere el tamaño de muestra efectivo por dominio, que no es
# recuperable a partir de la información disponible: se trabaja con las
# estimaciones y varianzas publicadas por el DANE y no con los microdatos, cuyo
# identificador municipal no se difunde por reserva estadística. La
# transformación logit cumple el mismo propósito de acotar las predicciones al
# intervalo unitario sin requerir esa información.
#
# ------------------------------------------------------------------------------
# SOBRE LA ESCALA
# ------------------------------------------------------------------------------
# La selección de covariables se realizó sobre la escala original de la
# proporción. Las transformaciones se evalúan aquí como variantes de estimación
# sobre las especificaciones ya seleccionadas, no como alternativas de
# selección: los criterios de información calculados en escalas distintas no
# son comparables y por ese motivo no se reportan.
#
# La componente de varianza tampoco es comparable entre escalas. El predictor y
# el error cuadrático medio sí lo son, porque fh() los devuelve
# retrotransformados a la escala original.
#
# ------------------------------------------------------------------------------
# ENTRADAS
# ------------------------------------------------------------------------------
#   output/matriz_sae_transformada_v2.rds
#   output/seleccion_stepwise_completo_both.rds
#
# ------------------------------------------------------------------------------
# SALIDAS (output/)
# ------------------------------------------------------------------------------
#   modelado_variantes.rds       Objetos de ajuste y registro de la estimación.
#   mod_estimacion.csv           Registro de qué se estimó y qué no.
# ==============================================================================

library(emdi)
library(dplyr)
library(here)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")

matriz  <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))
corrida <- readRDS(file.path(ruta_out, "seleccion_stepwise_completo_both.rds"))


# ==============================================================================
# CONFIGURACION
# ==============================================================================

SEMILLA <- 2906

# ---- Metodo principal -------------------------------------------------------
METODO_PRINCIPAL <- "reml"

# ---- Metodos de sensibilidad al procedimiento de estimacion -----------------
METODOS_SENSIBILIDAD <- c("ml", "amrl", "ampl", "amrl_yl", "ampl_yl")

# ---- Estimador del error cuadratico medio, fijo por escala ------------------
MSE_ORIGINAL <- "analytical"
MSE_LOG      <- "analytical"
MSE_LOGIT    <- "boot"

# ---- Replicas para el MSE por bootstrap -------------------------------------
B_MSE <- 500

# ---- Transformaciones -------------------------------------------------------
# vardir se pasa siempre en la escala original: fh() transforma la variable
# dependiente y las varianzas de muestreo de manera interna.
TRANSFORMACIONES <- list(
  list(nombre = "log",
       transformation     = "log",
       backtransformation = "bc_crude",
       mse = MSE_LOG),
  list(nombre = "logit",
       transformation     = "logit",
       backtransformation = "bc",
       mse = MSE_LOGIT)
)

formula_de <- function(vars) {
  as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
}

escalar <- function(x) {
  v <- tryCatch(as.numeric(x), error = function(e) NA_real_)
  if (length(v) != 1 || !is.finite(v)) return(NA_real_)
  v
}


# ==============================================================================
# ESPECIFICACIONES Y DATOS
# ------------------------------------------------------------------------------
# Se leen del objeto de la etapa anterior en lugar de escribirse a mano, para
# que cualquier cambio en la seleccion se propague sin intervencion.
# ==============================================================================

tabla_sel <- corrida$tabla
etiquetas <- c(KICb2 = "M1", KICc = "M2", KIC = "M3")

ESPECIFICACIONES <- lapply(seq_len(nrow(tabla_sel)), function(i) {
  list(etiqueta = unname(etiquetas[tabla_sel$criterio[i]]),
       criterio = tabla_sel$criterio[i],
       vars     = strsplit(tabla_sel$covariables[i], " + ", fixed = TRUE)[[1]])
})

names(ESPECIFICACIONES) <- vapply(ESPECIFICACIONES, function(x) x$etiqueta, "")
ESPECIFICACIONES <- ESPECIFICACIONES[order(names(ESPECIFICACIONES))]

cat("\n== ESPECIFICACIONES ==\n")
for (esp in ESPECIFICACIONES) {
  cat(esp$etiqueta, " (", esp$criterio, ", ", length(esp$vars), "): ",
      paste(esp$vars, collapse = ", "), "\n", sep = "")
}

vars_todas <- unique(unlist(lapply(ESPECIFICACIONES, function(x) x$vars)))

faltan <- setdiff(vars_todas, names(matriz))

if (length(faltan) > 0) {
  stop("Covariables ausentes en la matriz: ", paste(faltan, collapse = ", "))
}

datos <- matriz %>%
  select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza,
         all_of(vars_todas)) %>%
  as.data.frame()

n_dom <- nrow(datos)

stopifnot(
  "Hay varianzas de muestreo no positivas" = all(datos$varianza_pobreza > 0),
  "Hay valores faltantes en las covariables" =
    !any(is.na(datos[, vars_todas, drop = FALSE])),
  "Hay codigos municipales duplicados" = !any(duplicated(datos$cod_mun)))

# Las transformaciones requieren una proporcion en el intervalo abierto (0, 1).
rango <- range(datos$pobreza_monetaria, na.rm = TRUE)
cat("\nRecorrido de la variable dependiente: [", signif(rango[1], 4), ", ",
    signif(rango[2], 4), "]\n", sep = "")

if (rango[1] <= 0 || rango[2] >= 1) {
  warning("La variable dependiente alcanza los extremos del intervalo ",
          "unitario: las transformaciones pueden no ser estimables.")
}


# ==============================================================================
# FUNCIONES
# ==============================================================================

# ------------------------------------------------------------------------------
# Alineacion de dominios
# ------------------------------------------------------------------------------
# Predicciones, errores cuadraticos medios y residuos se combinan por posicion
# en la etapa de diagnostico. Una reordenacion de la salida de fh() invalidaria
# toda comparacion sin producir error.
alineado <- function(aj) {
  if (is.null(aj)) return(NA)
  dom <- tryCatch(as.character(aj$ind$Domain), error = function(e) NULL)
  if (is.null(dom)) return(NA)
  identical(dom, as.character(datos$cod_mun))
}

# ------------------------------------------------------------------------------
# Ajuste con registro de errores y avisos
# ------------------------------------------------------------------------------
# El tipo de MSE viene fijado por la escala. No hay sustitucion automatica: si
# el ajuste falla, se conserva el mensaje del error.
#
# La bandera var_limite se obtiene del texto del aviso, mecanismo que depende
# de la version del paquete. Se conserva el texto completo de los avisos para
# que la deteccion pueda revisarse.
ajustar <- function(argumentos, mse_type) {
  
  avisos <- character(0)
  error  <- NA_character_
  limite <- FALSE
  
  args <- c(argumentos, list(MSE = TRUE, mse_type = mse_type))
  
  aj <- withCallingHandlers(
    tryCatch(do.call(fh, args),
             error = function(e) { error <<- conditionMessage(e); NULL }),
    warning = function(w) {
      msg <- conditionMessage(w)
      avisos <<- c(avisos, msg)
      if (grepl("interval limit|singular", msg, ignore.case = TRUE)) {
        limite <<- TRUE
      }
      invokeRestart("muffleWarning")
    })
  
  list(ajuste = aj, mse_type = mse_type, error = error,
       avisos = unique(avisos), limite = limite)
}

# ------------------------------------------------------------------------------
# Registro de un ajuste
# ------------------------------------------------------------------------------
# Solo se registra lo ocurrido durante la estimacion. Las cantidades derivadas
# (coeficiente de variacion, normalidad, factor de contraccion) se calculan en
# el script de diagnostico.
registrar <- function(res, esp, metodo, escala, papel, segundos) {
  
  aj <- res$ajuste
  
  data.frame(
    modelo         = paste(esp$etiqueta, metodo, escala, sep = "_"),
    especificacion = esp$etiqueta,
    criterio       = esp$criterio,
    n_vars         = length(esp$vars),
    metodo         = metodo,
    escala         = escala,
    papel          = papel,
    mse_type       = res$mse_type,
    estimado       = !is.null(aj),
    alineado       = alineado(aj),
    sigma2_u       = signif(escalar(tryCatch(aj$model$variance,
                                             error = function(e) NA)), 6),
    var_limite     = res$limite,
    n_avisos       = length(res$avisos),
    avisos         = if (length(res$avisos) == 0) NA_character_
    else paste(res$avisos, collapse = " | "),
    error          = res$error,
    segundos       = round(segundos, 1),
    stringsAsFactors = FALSE)
}


# ==============================================================================
# BLOQUE A — ESCALA ORIGINAL
# ==============================================================================

cat("\n\n============================================================\n")
cat("BLOQUE A — ESCALA ORIGINAL\n")
cat("============================================================\n")
cat("Metodo principal:  ", METODO_PRINCIPAL, "\n", sep = "")
cat("Sensibilidad:      ", paste(METODOS_SENSIBILIDAD, collapse = ", "),
    "\n", sep = "")
cat("Estimador del MSE: ", MSE_ORIGINAL, "\n\n", sep = "")

filas <- list(); objetos <- list()

metodos_A <- c(METODO_PRINCIPAL, METODOS_SENSIBILIDAD)

for (esp in ESPECIFICACIONES) {
  for (metodo in metodos_A) {
    
    papel <- if (metodo == METODO_PRINCIPAL) "principal" else "sensibilidad"
    
    t0  <- Sys.time()
    res <- ajustar(list(formula_de(esp$vars),
                        vardir        = "varianza_pobreza",
                        combined_data = datos,
                        domains       = "cod_mun",
                        method        = metodo,
                        B             = c(0, 0),
                        seed          = SEMILLA),
                   MSE_ORIGINAL)
    tt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    fila <- registrar(res, esp, metodo, "original", papel, tt)
    filas[[length(filas) + 1]] <- fila
    
    if (!is.null(res$ajuste)) objetos[[fila$modelo]] <- res$ajuste
    
    cat(sprintf("  %-4s %-9s %-13s %s  alineado: %s%s\n",
                esp$etiqueta, metodo, papel,
                if (is.null(res$ajuste)) "NO ESTIMADO" else "estimado",
                as.character(fila$alineado),
                if (is.na(res$error)) "" else paste0("  ERROR: ", res$error)))
  }
}


# ==============================================================================
# BLOQUE B — TRANSFORMACIONES
# ------------------------------------------------------------------------------
# fh() transforma internamente tanto la variable dependiente como las varianzas
# de muestreo, de modo que vardir se pasa en la escala original.
# ==============================================================================

cat("\n\n============================================================\n")
cat("BLOQUE B — TRANSFORMACIONES\n")
cat("============================================================\n")
cat("Metodo:            ", METODO_PRINCIPAL, "\n", sep = "")
cat("Estimador del MSE: log = ", MSE_LOG, " | logit = ", MSE_LOGIT,
    "\n", sep = "")
cat("Replicas bootstrap:", B_MSE, "\n\n")

for (esp in ESPECIFICACIONES) {
  for (tr in TRANSFORMACIONES) {
    
    t0  <- Sys.time()
    res <- ajustar(list(formula_de(esp$vars),
                        vardir             = "varianza_pobreza",
                        combined_data      = datos,
                        domains            = "cod_mun",
                        method             = METODO_PRINCIPAL,
                        transformation     = tr$transformation,
                        backtransformation = tr$backtransformation,
                        B                  = c(B_MSE, 0),
                        seed               = SEMILLA),
                   tr$mse)
    tt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    fila <- registrar(res, esp, METODO_PRINCIPAL, tr$nombre, "principal", tt)
    filas[[length(filas) + 1]] <- fila
    
    if (!is.null(res$ajuste)) objetos[[fila$modelo]] <- res$ajuste
    
    cat(sprintf("  %-4s %-9s %-13s %s  alineado: %s  (%.1f s)%s\n",
                esp$etiqueta, tr$nombre, tr$mse,
                if (is.null(res$ajuste)) "NO ESTIMADO" else "estimado",
                as.character(fila$alineado), tt,
                if (is.na(res$error)) "" else paste0("  ERROR: ", res$error)))
  }
}


# ==============================================================================
# REGISTRO
# ==============================================================================

registro <- do.call(rbind, filas)

cat("\n\n============================================================\n")
cat("REGISTRO DE LA ESTIMACION\n")
cat("============================================================\n\n")

print(registro %>%
        select(especificacion, metodo, escala, papel, mse_type, estimado,
               alineado, sigma2_u, var_limite, n_avisos, segundos),
      row.names = FALSE)

fallidos <- registro %>% filter(!estimado)

if (nrow(fallidos) > 0) {
  cat("\n-- Ajustes no estimados --\n")
  print(fallidos %>% select(modelo, mse_type, error), row.names = FALSE)
}

con_avisos <- registro %>% filter(n_avisos > 0)

if (nrow(con_avisos) > 0) {
  cat("\n-- Ajustes con avisos --\n")
  print(con_avisos %>% select(modelo, avisos), row.names = FALSE)
}

desalineados <- registro %>% filter(estimado & (is.na(alineado) | !alineado))

if (nrow(desalineados) > 0) {
  stop("Hay ajustes cuya salida no coincide con el orden de 'datos': ",
       paste(desalineados$modelo, collapse = ", "))
}

# ---- Sensibilidad de la componente de varianza al metodo --------------------
# Solo en la escala original: sigma2_u no es comparable entre escalas.

comp <- registro %>%
  filter(estimado, escala == "original") %>%
  select(especificacion, metodo, papel, sigma2_u)

ref <- comp %>% filter(metodo == METODO_PRINCIPAL) %>%
  select(especificacion, sigma2_ref = sigma2_u)

cat("\n-- Desviacion porcentual de sigma2_u respecto de ",
    METODO_PRINCIPAL, " --\n", sep = "")

print(comp %>%
        left_join(ref, by = "especificacion") %>%
        filter(metodo != METODO_PRINCIPAL) %>%
        mutate(desvio_rel = round(100 * (sigma2_u / sigma2_ref - 1), 2)) %>%
        select(especificacion, metodo, sigma2_u, sigma2_ref, desvio_rel) %>%
        as.data.frame(),
      row.names = FALSE)

write.csv(registro, file.path(ruta_out, "mod_estimacion.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# GUARDADO
# ------------------------------------------------------------------------------
# Se conservan los objetos fh completos: el diagnostico y la Etapa III
# requieren residuos, efectos aleatorios y errores cuadraticos medios, y
# reajustar con bootstrap resultaria costoso.
# ==============================================================================

saveRDS(
  list(
    registro             = registro,
    ajustes              = objetos,
    especificaciones     = ESPECIFICACIONES,
    datos                = datos,
    metodo_principal     = METODO_PRINCIPAL,
    metodos_sensibilidad = METODOS_SENSIBILIDAD,
    transformaciones     = TRANSFORMACIONES,
    mse_original         = MSE_ORIGINAL,
    mse_log              = MSE_LOG,
    mse_logit            = MSE_LOGIT,
    B_mse                = B_MSE,
    semilla              = SEMILLA,
    session              = sessionInfo()
  ),
  file.path(ruta_out, "modelado_variantes.rds"))

cat("\n============================================================\n")
cat("ESTIMACION COMPLETADA\n")
cat("Variantes estimadas: ", sum(registro$estimado), " de ", nrow(registro),
    "\n", sep = "")
cat("============================================================\n")