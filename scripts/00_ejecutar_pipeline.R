# ==============================================================================
# 00_ejecutar_pipeline.R
# ------------------------------------------------------------------------------
# Ejecuta los scripts del análisis en orden, cada uno en un proceso de R
# separado y no interactivo (callr::rscript). Así las pausas que emdi introduce
# entre figuras no detienen la ejecución, y cada script parte de un entorno
# limpio, como ocurriría al reproducir el trabajo desde cero.
#
# Si un script falla, la ejecución se detiene y se muestra el error.
#
# Uso: abrir el proyecto SAEtesis.Rproj y ejecutar este archivo completo
# (source o Ctrl+Shift+Enter).
# ==============================================================================

library(here)
library(callr)

CARPETA <- here("scripts")

# Orden de ejecución. Poner FALSE para omitir un paso cuyas salidas ya existen.
#
# Pasos lentos, desactivados por defecto porque sus resultados ya están en
# output/:
#   02_1c  selección de covariables (cerca de tres horas);
#   02_1d  selección en escala logarítmica, sensibilidad del anexo;
#   02_1e  estabilidad de la selección (57 corridas, ejecutadas en una máquina
#          virtual).
#
# 02_4 estima los ajustes robustos con bootstrap (B = 500) y también tarda;
# guarda un archivo parcial y retoma lo pendiente si se interrumpe.
PASOS <- c(
  # Etapa I: base integrada
  "01_etapa1_construccion_base.R"            = TRUE,
  "01_1_verificacion_imputacion.R"           = TRUE,

  # Etapa II: selección y ajuste de los modelos
  "02_1c_seleccion_stepwise.R"               = FALSE,
  "02_1d_seleccion_kicb2_log.R"              = FALSE,
  "02_1e_estabilidad_seleccion_kicb2.R"      = FALSE,
  "02_2a_modelado_variantes.R"               = TRUE,
  "02_2a_modelado_variantes_diagnostico.R"   = TRUE,
  "02_3_espacial.R"                          = TRUE,
  "02_3b_ganancia_espacial.R"                = TRUE,
  "02_4_robusto.R"                           = TRUE,
  "02_4_robusto_diagnostico.R"               = TRUE,
  "02_4d_verificacion_warnholz.R"            = TRUE,
  "02_4e_comparacion_mse_pseudo_boot.R"      = TRUE,
  "02_4f_cifras_robusto_mse.R"               = TRUE,
  "02_4g_figura_coeficientes.R"              = TRUE,

  # Etapa III: supuestos, precisión, consistencia y benchmarking
  "03_0_normalidad.R"                        = TRUE,
  "03_1_precision_diagnosticos.R"            = TRUE,
  "03_1d_incertidumbre_referencia.R"         = TRUE,
  "03_2_agregado_departamental.R"            = TRUE,
  "03_9_evaluacion_seleccion_log.R"          = TRUE,
  "04_1_benchmarking_m1log.R"                = TRUE,

  # Estimaciones finales de M1 log. Corre antes del Shapley, que verifica que
  # la coalición completa coincide con la estimación final.
  "06_0_estimaciones_coef_cv_m1_log.R"       = TRUE,

  # Etapa IV: contribución de las variables auxiliares (Shapley)
  "05_1_shapley_etapa4_completo.R"             = TRUE,

  # Mapas y análisis exploratorio de las especificaciones
  "06_1_mapas_m1_log.R"                      = TRUE,
  "07_2_exploratorio_especificaciones.R"                                   = TRUE,

  # Tablero interactivo
  "08_exportar_datos_dashboard.R"            = TRUE,
  "09_dashboard_construir.R"                 = TRUE
)

# Solo se exige que existan los scripts que se van a ejecutar.
activos <- names(PASOS)[PASOS]
faltan  <- activos[!file.exists(file.path(CARPETA, activos))]

if (length(faltan) > 0) {
  stop("No se encontraron en ", CARPETA, ": ", paste(faltan, collapse = ", "))
}

if (!PASOS[["02_1c_seleccion_stepwise.R"]] &&
    !file.exists(here("output", "seleccion_stepwise_completo_both.rds"))) {
  stop("La selección está desactivada pero no existe ",
       "output/seleccion_stepwise_completo_both.rds.")
}

if (!PASOS[["02_1d_seleccion_kicb2_log.R"]] &&
    !file.exists(here("output", "seleccion_kicb2_log.rds"))) {
  stop("La selección en escala logarítmica está desactivada pero no existe ",
       "output/seleccion_kicb2_log.rds.")
}

t_total <- Sys.time()

for (s in activos) {

  cat("\n############################################################\n")
  cat("# ", s, "\n", sep = "")
  cat("############################################################\n")

  t0 <- Sys.time()

  callr::rscript(file.path(CARPETA, s), wd = here(), show = TRUE,
                 fail_on_status = TRUE)

  cat("\n>>> ", s, " completado en ",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
      " minutos\n", sep = "")
}

cat("\nCadena completada en ",
    round(as.numeric(difftime(Sys.time(), t_total, units = "mins")), 2),
    " minutos\n", sep = "")
