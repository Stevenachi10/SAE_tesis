# ==============================================================================
# 02_4_robusto.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.4 — Extensión robusta del modelo Fay-Herriot (ESTIMACION)
#
# Este script solo estima y guarda. El cálculo de coeficientes de variación, la
# comparación de coeficientes, el desplazamiento por dominio y las figuras
# corresponden a 02_4_robusto_diagnostico.R.
#
# ------------------------------------------------------------------------------
# OBJETIVO
# ------------------------------------------------------------------------------
# Ajustar las variantes robustas sobre las tres especificaciones seleccionadas
# en la subetapa de selección de variables, y establecer en qué medida las
# estimaciones municipales dependen de dominios influyentes.
#
# No se realiza un diagnóstico previo de dominios atípicos. La variante robusta
# no requiere que se incumpla un supuesto para ser aplicable: se comporta de
# manera equivalente al estimador estándar en ausencia de observaciones
# influyentes, de modo que la comparación entre ambos ajustes constituye por sí
# misma la evidencia disponible sobre la presencia de tales dominios. Las
# herramientas de influencia específicas para estimación en áreas pequeñas
# (Marcis et al., 2023) no cuentan con implementación en los paquetes de uso
# corriente.
#
# ------------------------------------------------------------------------------
# DECISIONES DE ESTIMACION
# ------------------------------------------------------------------------------
# CONSTANTE DE SINTONIZACION. Se fija k = 1,345, valor convencional de la
# función de Huber por corresponder a una eficiencia asintótica del 95 % bajo
# normalidad, y valor por defecto en el paquete. Los valores 1,0 y 2,0 se
# ajustan por separado como análisis de sensibilidad, no como especificaciones
# candidatas.
#
# METODOS. reblup y reblupbc son dos procedimientos distintos: el segundo
# incorpora una corrección de sesgo con constante multiplicativa. Se ajustan
# ambos y se comparan por separado con el modelo convencional, sin agregarlos.
#
# CONSTANTE MULTIPLICATIVA. mult_constant = 1 es una elección, no un valor
# neutro. Se fija en 1 y se declara como tal.
#
# ERROR CUADRATICO MEDIO. Se emplea la aproximación por linealización
# ("pseudo") como estimador principal. No se aplica ningún mecanismo de
# sustitución automática: si el ajuste falla, el error se registra y se reporta.
# La comparación de coeficientes de variación frente al modelo convencional
# enfrenta estimadores distintos del error cuadrático medio (linealización
# frente a aproximación analítica de Prasad-Rao), circunstancia que debe
# declararse y que impide que esa comparación sustente por sí sola una
# conclusión.
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
#   robusto.rds        Objetos de ajuste y registro de la estimación.
#   rob_estimacion.csv Registro de qué se estimó y qué no.
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

# Metodos robustos. reblup: verosimilitud robustificada con prediccion robusta.
# reblupbc: la anterior con correccion de sesgo.
METODOS <- c("reblup", "reblupbc")

# Constante de sintonizacion principal y valores de sensibilidad.
K_PRINCIPAL   <- 1.345
K_SENSIBILIDAD <- c(1.0, 2.0)

# Constante multiplicativa de la correccion de sesgo (valor por defecto).
MULT_CONSTANT <- 1

# Estimador del error cuadratico medio. Sin sustitucion automatica.
MSE_ROBUSTO <- "pseudo"

# Metodo y estimador de referencia.
METODO_BASE <- "reml"
MSE_BASE    <- "analytical"

formula_de <- function(vars) {
  as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
}


# ==============================================================================
# ESPECIFICACIONES Y DATOS
# ==============================================================================

tabla_sel <- corrida$tabla
etiquetas <- c(KICb2 = "M1", KICc = "M2", KIC = "M3")

ESPECIFICACIONES <- list()

for (i in seq_len(nrow(tabla_sel))) {
  et <- unname(etiquetas[tabla_sel$criterio[i]])
  ESPECIFICACIONES[[et]] <- strsplit(tabla_sel$covariables[i],
                                     " + ", fixed = TRUE)[[1]]
}

ESPECIFICACIONES <- ESPECIFICACIONES[order(names(ESPECIFICACIONES))]

cat("\n== ESPECIFICACIONES ==\n")
for (nm in names(ESPECIFICACIONES)) {
  cat(nm, " (", length(ESPECIFICACIONES[[nm]]), "): ",
      paste(ESPECIFICACIONES[[nm]], collapse = ", "), "\n", sep = "")
}

vars_todas <- unique(unlist(ESPECIFICACIONES))

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


# ==============================================================================
# FUNCIONES
# ==============================================================================

# ------------------------------------------------------------------------------
# Alineacion de dominios
# ------------------------------------------------------------------------------
# Las estimaciones, los errores cuadraticos medios y los residuos se combinan
# por posicion en la etapa de diagnostico. Una reordenacion de la salida de
# fh() invalidaria toda comparacion sin producir error.
alineado <- function(aj) {
  if (is.null(aj)) return(NA)
  dom <- tryCatch(as.character(aj$ind$Domain), error = function(e) NULL)
  if (is.null(dom)) return(NA)
  identical(dom, as.character(datos$cod_mun))
}

# ------------------------------------------------------------------------------
# Ajuste con registro de errores y avisos
# ------------------------------------------------------------------------------
# Los errores no se silencian: se devuelve el mensaje para que quede en la
# tabla de registro. Un ajuste que no se estima debe ser visible.
ajustar <- function(argumentos) {
  
  avisos <- character(0)
  error  <- NA_character_
  
  aj <- withCallingHandlers(
    tryCatch(do.call(fh, argumentos),
             error = function(e) { error <<- conditionMessage(e); NULL }),
    warning = function(w) {
      avisos <<- c(avisos, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  
  list(ajuste = aj, error = error, avisos = unique(avisos))
}

ajustar_base <- function(vars) {
  ajustar(list(formula_de(vars),
               vardir        = "varianza_pobreza",
               combined_data = datos,
               domains       = "cod_mun",
               method        = METODO_BASE,
               MSE           = TRUE,
               mse_type      = MSE_BASE,
               B             = c(0, 0),
               seed          = SEMILLA))
}

ajustar_robusto <- function(vars, metodo, k) {
  ajustar(list(formula_de(vars),
               vardir        = "varianza_pobreza",
               combined_data = datos,
               domains       = "cod_mun",
               method        = metodo,
               k             = k,
               mult_constant = MULT_CONSTANT,
               MSE           = TRUE,
               mse_type      = MSE_ROBUSTO,
               seed          = SEMILLA))
}

# ------------------------------------------------------------------------------
# Registro de un ajuste
# ------------------------------------------------------------------------------
# Solo se registra lo que ocurrio durante la estimacion. Las cantidades
# derivadas se calculan en el script de diagnostico.
registrar <- function(res, etiqueta, metodo, k, papel, segundos) {
  
  aj <- res$ajuste
  
  s2u <- tryCatch({
    v <- as.numeric(aj$model$variance)
    if (length(v) == 1 && is.finite(v)) v else NA_real_
  }, error = function(e) NA_real_)
  
  data.frame(
    modelo         = paste(etiqueta, metodo, format(k), sep = "_"),
    especificacion = etiqueta,
    n_vars         = length(ESPECIFICACIONES[[etiqueta]]),
    metodo         = metodo,
    k              = k,
    papel          = papel,
    mse_type       = if (metodo == METODO_BASE) MSE_BASE else MSE_ROBUSTO,
    estimado       = !is.null(aj),
    alineado       = alineado(aj),
    sigma2_u       = signif(s2u, 6),
    n_avisos       = length(res$avisos),
    avisos         = if (length(res$avisos) == 0) NA_character_
    else paste(res$avisos, collapse = " | "),
    error          = res$error,
    segundos       = round(segundos, 1),
    stringsAsFactors = FALSE)
}


# ==============================================================================
# AJUSTE DE REFERENCIA
# ==============================================================================

cat("\n== AJUSTE DE REFERENCIA (", METODO_BASE, ", ", MSE_BASE, ") ==\n",
    sep = "")

filas <- list(); base_ajustes <- list()

for (nm in names(ESPECIFICACIONES)) {
  
  t0  <- Sys.time()
  res <- ajustar_base(ESPECIFICACIONES[[nm]])
  tt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  base_ajustes[[nm]] <- res$ajuste
  
  filas[[length(filas) + 1]] <-
    registrar(res, nm, METODO_BASE, NA_real_, "referencia", tt)
  
  cat(sprintf("  %-4s  %s  alineado: %s%s\n", nm,
              if (is.null(res$ajuste)) "NO ESTIMADO" else "estimado",
              as.character(alineado(res$ajuste)),
              if (is.na(res$error)) "" else paste0("  ERROR: ", res$error)))
}


# ==============================================================================
# AJUSTE ROBUSTO — CONSTANTE PRINCIPAL
# ==============================================================================

cat("\n== AJUSTE ROBUSTO, k = ", K_PRINCIPAL, " ==\n", sep = "")

rob_ajustes <- list()

for (nm in names(ESPECIFICACIONES)) {
  for (metodo in METODOS) {
    
    t0  <- Sys.time()
    res <- ajustar_robusto(ESPECIFICACIONES[[nm]], metodo, K_PRINCIPAL)
    tt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    fila <- registrar(res, nm, metodo, K_PRINCIPAL, "principal", tt)
    filas[[length(filas) + 1]] <- fila
    
    if (!is.null(res$ajuste)) rob_ajustes[[fila$modelo]] <- res$ajuste
    
    cat(sprintf("  %-4s %-9s  %s  alineado: %s%s\n", nm, metodo,
                if (is.null(res$ajuste)) "NO ESTIMADO" else "estimado",
                as.character(fila$alineado),
                if (is.na(res$error)) "" else paste0("  ERROR: ", res$error)))
  }
}


# ==============================================================================
# AJUSTE ROBUSTO — SENSIBILIDAD A LA CONSTANTE
# ==============================================================================

cat("\n== SENSIBILIDAD A LA CONSTANTE DE SINTONIZACION ==\n")
cat("Valores adicionales:", paste(K_SENSIBILIDAD, collapse = ", "), "\n")

for (nm in names(ESPECIFICACIONES)) {
  for (metodo in METODOS) {
    for (k in K_SENSIBILIDAD) {
      
      t0  <- Sys.time()
      res <- ajustar_robusto(ESPECIFICACIONES[[nm]], metodo, k)
      tt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      
      fila <- registrar(res, nm, metodo, k, "sensibilidad", tt)
      filas[[length(filas) + 1]] <- fila
      
      if (!is.null(res$ajuste)) rob_ajustes[[fila$modelo]] <- res$ajuste
      
      cat(sprintf("  %-4s %-9s k = %-5s  %s\n", nm, metodo, format(k),
                  if (is.null(res$ajuste)) "NO ESTIMADO" else "estimado"))
    }
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
        select(especificacion, metodo, k, papel, estimado, alineado,
               sigma2_u, n_avisos, segundos),
      row.names = FALSE)

fallidos <- registro %>% filter(!estimado)

if (nrow(fallidos) > 0) {
  cat("\n-- Ajustes no estimados --\n")
  print(fallidos %>% select(modelo, error), row.names = FALSE)
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

write.csv(registro, file.path(ruta_out, "rob_estimacion.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(
  list(
    registro         = registro,
    ajustes_base     = base_ajustes,
    ajustes_robustos = rob_ajustes,
    especificaciones = ESPECIFICACIONES,
    datos            = datos,
    metodos          = METODOS,
    k_principal      = K_PRINCIPAL,
    k_sensibilidad   = K_SENSIBILIDAD,
    mult_constant    = MULT_CONSTANT,
    mse_robusto      = MSE_ROBUSTO,
    metodo_base      = METODO_BASE,
    mse_base         = MSE_BASE,
    semilla          = SEMILLA,
    session          = sessionInfo()
  ),
  file.path(ruta_out, "robusto.rds"))

cat("\n============================================================\n")
cat("ESTIMACION COMPLETADA\n")
cat("Ajustes estimados: ", sum(registro$estimado), " de ", nrow(registro),
    "\n", sep = "")
cat("============================================================\n")
