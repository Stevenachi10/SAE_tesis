# ==============================================================================
# 02_4_robusto.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.4. Extensión robusta del modelo Fay-Herriot (ESTIMACION)
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
# neutro. Se fija en 1 y se declara como tal. En saeRobust, del que emdi toma
# la estimación robusta, reblupbc es un estimador de traslación limitada: si el
# REBLUP se aleja del directo más de mult_constant errores estándar, se desplaza
# hasta esa distancia.
#
# ERROR CUADRATICO MEDIO. El estimador principal es el bootstrap paramétrico
# (Warnholz, 2016), con B = 500 réplicas. La aproximación por
# pseudolinealización (Warnholz, 2016, basada en la representación pseudolineal
# de Chambers, Chandra y Tzavidis, 2011) se estima como análisis de
# sensibilidad con la constante principal. La elección responde a dos
# propiedades de la pseudolinealización verificadas en 02_4d:
#
#   1. Trata como fijos los pesos del predictor y omite la variabilidad debida
#      a la estimación de beta y sigma2_u. Es el único de los estimadores del
#      error cuadrático medio empleados en la tesis que no la incorpora.
#
#   2. En reblupbc, la representación pseudolineal de los dominios recortados
#      asigna un único peso al estimador directo y atribuye un término de sesgo
#      del orden de un error estándar, que en estos datos representa cerca de
#      la mitad del error cuadrático medio de esos dominios.
#
# El bootstrap paramétrico simula bajo el modelo normal de trabajo, sin efectos
# de área atípicos, de modo que puede subestimar el error cuadrático medio en
# dominios con efectos realmente atípicos. La limitación se declara.
#
# No se aplica ningún mecanismo de sustitución automática: si un ajuste falla,
# el error se registra y se reporta.
#
# CONVERGENCIA. Se emplean tol = 1e-8 y maxit = 1000, los mismos valores de
# 02_3_espacial.R. Con los valores por defecto de fh() (1e-4 y 100), la
# representación pseudolineal W y difería del REBLUP hasta en 0,0035 (02_4d),
# señal de que el ajuste no alcanzaba el punto fijo con la precisión necesaria.
# El bloque final de este script verifica la convergencia.
#
# ------------------------------------------------------------------------------
# GUARDADO PARCIAL
# ------------------------------------------------------------------------------
# Los ajustes con bootstrap son costosos. Cada uno se guarda al terminar en
# output/robusto_parcial.rds, y una ejecución interrumpida retoma solo los
# pendientes. El archivo parcial se descarta si la configuración o los datos
# cambiaron, y se elimina al completar la estimación.
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
#   robusto.rds                  Objetos de ajuste y registro de la estimación.
#   rob_estimacion.csv           Registro de qué se estimó y qué no.
#   rob_verificacion_convergencia.csv
# ==============================================================================

library(emdi)
library(saeRobust)
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
K_PRINCIPAL    <- 1.345
K_SENSIBILIDAD <- c(1.0, 2.0)

# Constante multiplicativa de la correccion de sesgo.
MULT_CONSTANT <- 1

# Estimador principal del error cuadratico medio y replicas.
MSE_ROBUSTO <- "boot"
B_BOOT      <- 500

# Estimador de sensibilidad, solo con la constante principal.
MSE_SENSIBILIDAD <- "pseudo"

# Convergencia, como en 02_3_espacial.R.
TOL   <- 1e-8
MAXIT <- 1000

# Metodo y estimador de referencia.
METODO_BASE <- "reml"
MSE_BASE    <- "analytical"

ARCHIVO_PARCIAL <- file.path(ruta_out, "robusto_parcial.rds")

formula_de <- function(vars) {
  as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
}

clave_de <- function(nm, metodo, k) paste(nm, metodo, format(k), sep = "_")


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
               tol           = TOL,
               maxit         = MAXIT,
               seed          = SEMILLA))
}

ajustar_robusto <- function(vars, metodo, k, mse_type) {
  ajustar(list(formula_de(vars),
               vardir        = "varianza_pobreza",
               combined_data = datos,
               domains       = "cod_mun",
               method        = metodo,
               k             = k,
               mult_constant = MULT_CONSTANT,
               MSE           = TRUE,
               mse_type      = mse_type,
               B             = if (mse_type == "boot") c(B_BOOT, 0) else c(0, 0),
               tol           = TOL,
               maxit         = MAXIT,
               seed          = SEMILLA))
}

# ------------------------------------------------------------------------------
# Registro de un ajuste
# ------------------------------------------------------------------------------
# Solo se registra lo que ocurrio durante la estimacion. Las cantidades
# derivadas se calculan en el script de diagnostico. El metodo del ECM se toma
# del propio objeto, no de la configuracion.
registrar <- function(res, etiqueta, metodo, k, papel, mse_type, segundos) {
  
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
    mse_type       = mse_type,
    mse_objeto     = if (is.null(aj)) NA_character_ else aj$method$MSE_method,
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

filas_base <- list(); base_ajustes <- list()

for (nm in names(ESPECIFICACIONES)) {
  
  t0  <- Sys.time()
  res <- ajustar_base(ESPECIFICACIONES[[nm]])
  tt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  base_ajustes[[nm]] <- res$ajuste
  
  filas_base[[nm]] <-
    registrar(res, nm, METODO_BASE, NA_real_, "referencia", MSE_BASE, tt)
  
  cat(sprintf("  %-4s  %s  alineado: %s%s\n", nm,
              if (is.null(res$ajuste)) "NO ESTIMADO" else "estimado",
              as.character(alineado(res$ajuste)),
              if (is.na(res$error)) "" else paste0("  ERROR: ", res$error)))
}


# ==============================================================================
# AJUSTE ROBUSTO CON BOOTSTRAP: CONSTANTE PRINCIPAL Y SENSIBILIDAD A k
# ==============================================================================

CONFIG <- list(tol = TOL, maxit = MAXIT, B = B_BOOT, semilla = SEMILLA,
               mult = MULT_CONSTANT, mse = MSE_ROBUSTO,
               especificaciones = ESPECIFICACIONES, datos = datos)

parcial <- if (file.exists(ARCHIVO_PARCIAL)) readRDS(ARCHIVO_PARCIAL) else NULL

if (!is.null(parcial) && !identical(parcial$config, CONFIG)) {
  cat("\nEl archivo parcial corresponde a otra configuración o a otros datos;",
      "se descarta.\n")
  parcial <- NULL
}

if (is.null(parcial)) {
  parcial <- list(config = CONFIG, ajustes = list(), filas = list())
} else {
  cat("\nSe retoman", length(parcial$ajustes), "ajustes ya guardados.\n")
}

plan <- rbind(
  expand.grid(nm = names(ESPECIFICACIONES), metodo = METODOS,
              k = K_PRINCIPAL, papel = "principal",
              stringsAsFactors = FALSE),
  expand.grid(nm = names(ESPECIFICACIONES), metodo = METODOS,
              k = K_SENSIBILIDAD, papel = "sensibilidad_k",
              stringsAsFactors = FALSE))

plan <- plan[order(plan$papel != "principal", plan$nm, plan$metodo, plan$k), ]

cat("\n== AJUSTE ROBUSTO, ECM POR BOOTSTRAP (B = ", B_BOOT, ") ==\n", sep = "")
cat("Ajustes previstos:", nrow(plan), "\n")

for (i in seq_len(nrow(plan))) {
  
  nm     <- plan$nm[i]
  metodo <- plan$metodo[i]
  k      <- plan$k[i]
  papel  <- plan$papel[i]
  clave  <- clave_de(nm, metodo, k)
  
  if (!is.null(parcial$ajustes[[clave]])) {
    cat(sprintf("  %-22s ya estimado, se omite\n", clave))
    next
  }
  
  t0  <- Sys.time()
  res <- ajustar_robusto(ESPECIFICACIONES[[nm]], metodo, k, MSE_ROBUSTO)
  tt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  parcial$filas[[clave]] <-
    registrar(res, nm, metodo, k, papel, MSE_ROBUSTO, tt)
  
  if (!is.null(res$ajuste)) parcial$ajustes[[clave]] <- res$ajuste
  
  saveRDS(parcial, ARCHIVO_PARCIAL)
  
  cat(sprintf("  %-22s %-14s %s  (%.1f min)%s\n", clave, papel,
              if (is.null(res$ajuste)) "NO ESTIMADO" else "estimado", tt / 60,
              if (is.na(res$error)) "" else paste0("  ERROR: ", res$error)))
}

rob_ajustes <- parcial$ajustes


# ==============================================================================
# SENSIBILIDAD AL ESTIMADOR DEL ECM: PSEUDOLINEALIZACION, CONSTANTE PRINCIPAL
# ==============================================================================

cat("\n== SENSIBILIDAD: ECM POR PSEUDOLINEALIZACION, k = ", K_PRINCIPAL,
    " ==\n", sep = "")

rob_pseudo <- list(); filas_pseudo <- list()

for (nm in names(ESPECIFICACIONES)) {
  for (metodo in METODOS) {
    
    clave <- clave_de(nm, metodo, K_PRINCIPAL)
    
    t0  <- Sys.time()
    res <- ajustar_robusto(ESPECIFICACIONES[[nm]], metodo, K_PRINCIPAL,
                           MSE_SENSIBILIDAD)
    tt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    fila <- registrar(res, nm, metodo, K_PRINCIPAL, "sensibilidad_mse",
                      MSE_SENSIBILIDAD, tt)
    filas_pseudo[[clave]] <- fila
    
    if (!is.null(res$ajuste)) rob_pseudo[[clave]] <- res$ajuste
    
    cat(sprintf("  %-22s %s  alineado: %s%s\n", clave,
                if (is.null(res$ajuste)) "NO ESTIMADO" else "estimado",
                as.character(fila$alineado),
                if (is.na(res$error)) "" else paste0("  ERROR: ", res$error)))
  }
}


# ==============================================================================
# REGISTRO
# ==============================================================================

registro <- do.call(rbind, c(filas_base, parcial$filas, filas_pseudo))
rownames(registro) <- NULL

cat("\n\n============================================================\n")
cat("REGISTRO DE LA ESTIMACION\n")
cat("============================================================\n\n")

print(registro %>%
        select(especificacion, metodo, k, papel, mse_objeto, estimado,
               alineado, sigma2_u, n_avisos, segundos),
      row.names = FALSE)

fallidos <- registro %>% filter(!estimado)

if (nrow(fallidos) > 0) {
  cat("\n-- Ajustes no estimados --\n")
  print(fallidos %>% select(modelo, papel, error), row.names = FALSE)
}

con_avisos <- registro %>% filter(n_avisos > 0)

if (nrow(con_avisos) > 0) {
  cat("\n-- Ajustes con avisos --\n")
  print(con_avisos %>% select(modelo, papel, avisos), row.names = FALSE)
}

desalineados <- registro %>% filter(estimado & (is.na(alineado) | !alineado))

if (nrow(desalineados) > 0) {
  stop("Hay ajustes cuya salida no coincide con el orden de 'datos': ",
       paste(desalineados$modelo, collapse = ", "))
}

# Las versiones bootstrap y pseudolinealizacion deben compartir las
# estimaciones puntuales: solo cambia el estimador del ECM.
for (clave in names(rob_pseudo)) {
  if (!is.null(rob_ajustes[[clave]])) {
    d <- max(abs(rob_ajustes[[clave]]$ind$FH - rob_pseudo[[clave]]$ind$FH))
    if (d > 1e-10) {
      stop("Las estimaciones puntuales de ", clave,
           " difieren entre bootstrap y pseudolinealizacion (", signif(d, 3), ").")
    }
  }
}

cat("\nOK  bootstrap y pseudolinealización comparten las estimaciones puntuales.\n")


# ==============================================================================
# VERIFICACION DE CONVERGENCIA
# ------------------------------------------------------------------------------
# fh() no conserva el objeto interno de saeRobust. Se reajusta con la misma
# llamada que emdi, con la misma tolerancia, y se comprueba:
#   - que el reajuste reproduce el REBLUP de emdi;
#   - que la representación pseudolineal W y reproduce el REBLUP. Una
#     diferencia apreciable indica que el ajuste no alcanzó el punto fijo.
# ==============================================================================

cat("\n== VERIFICACION DE CONVERGENCIA ==\n\n")

verif_conv <- do.call(rbind, lapply(names(ESPECIFICACIONES), function(nm) {
  
  clave <- clave_de(nm, "reblup", K_PRINCIPAL)
  aj    <- rob_ajustes[[clave]]
  
  if (is.null(aj)) {
    return(data.frame(especificacion = nm, dif_reajuste = NA_real_,
                      dif_Wy_reblup = NA_real_, stringsAsFactors = FALSE))
  }
  
  fit <- saeRobust::rfh(formula_de(ESPECIFICACIONES[[nm]]),
                        data        = datos,
                        samplingVar = "varianza_pobreza",
                        k           = K_PRINCIPAL,
                        tol         = TOL,
                        maxIter     = MAXIT)
  
  W <- as.matrix(weights(fit, c = MULT_CONSTANT)$W)
  
  data.frame(
    especificacion = nm,
    dif_reajuste   = signif(max(abs(as.numeric(fit$reblup) - aj$ind$FH)), 3),
    dif_Wy_reblup  = signif(max(abs(as.numeric(W %*% datos$pobreza_monetaria) -
                                      as.numeric(fit$reblup))), 3),
    stringsAsFactors = FALSE)
}))

print(verif_conv, row.names = FALSE)

if (any(verif_conv$dif_Wy_reblup > 1e-6, na.rm = TRUE)) {
  cat("\nAVISO  W y difiere del REBLUP en más de 1e-6: revisar tol y maxit.\n")
} else {
  cat("\nOK  la representación pseudolineal reproduce el REBLUP.\n")
}

write.csv(verif_conv, file.path(ruta_out, "rob_verificacion_convergencia.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# GUARDADO
# ==============================================================================

write.csv(registro, file.path(ruta_out, "rob_estimacion.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(
  list(
    registro                 = registro,
    ajustes_base             = base_ajustes,
    ajustes_robustos         = rob_ajustes,
    ajustes_robustos_pseudo  = rob_pseudo,
    especificaciones         = ESPECIFICACIONES,
    datos                    = datos,
    metodos                  = METODOS,
    k_principal              = K_PRINCIPAL,
    k_sensibilidad           = K_SENSIBILIDAD,
    mult_constant            = MULT_CONSTANT,
    mse_robusto              = MSE_ROBUSTO,
    mse_sensibilidad         = MSE_SENSIBILIDAD,
    b_boot                   = B_BOOT,
    tol                      = TOL,
    maxit                    = MAXIT,
    metodo_base              = METODO_BASE,
    mse_base                 = MSE_BASE,
    semilla                  = SEMILLA,
    verificacion_convergencia = verif_conv,
    session                  = sessionInfo()
  ),
  file.path(ruta_out, "robusto.rds"))

# La estimación está completa y guardada: el archivo parcial ya no hace falta.
if (file.exists(ARCHIVO_PARCIAL)) file.remove(ARCHIVO_PARCIAL)

cat("\n============================================================\n")
cat("ESTIMACION COMPLETADA\n")
cat("Ajustes estimados: ", sum(registro$estimado), " de ", nrow(registro),
    "\n", sep = "")
cat("============================================================\n")