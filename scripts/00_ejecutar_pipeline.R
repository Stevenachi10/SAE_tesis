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
# 02_1c tarda cerca de tres horas: si output/seleccion_stepwise_completo_both.rds
# ya existe, dejarlo en FALSE.
#
# 02_4 estima los ajustes robustos con bootstrap (B = 500) y también tarda;
# guarda un archivo parcial y retoma lo pendiente si se interrumpe.
PASOS <- c(
  # Etapa I: base integrada
  "01_etapa1_construccion_base.R"            = TRUE,

  # Etapa II: selección y ajuste de los modelos
  "02_1c_seleccion_stepwise.R"               = FALSE,
  "02_2a_modelado_variantes.R"               = TRUE,
  "02_2a_modelado_variantes_diagnostico.R"   = TRUE,
  "02_3_espacial.R"                          = TRUE,
  "02_4_robusto.R"                           = TRUE,
  "02_4_robusto_diagnostico.R"               = TRUE,
  "02_4d_verificacion_warnholz.R"            = TRUE,
  "02_4e_comparacion_mse_pseudo_boot.R"      = TRUE,
  "02_4f_cifras_robusto_mse.R"               = TRUE,

  # Etapa III: supuestos, precisión, consistencia y benchmarking
  "03_0_normalidad.R"                        = TRUE,
  "03_1_precision_diagnosticos.R"            = TRUE,
  "03_2_agregado_departamental.R"            = TRUE,
  "03_3_distribucion_gamma.R"                = TRUE,
  "04_1_benchmarking_m1log.R"                = TRUE,

  # Etapa IV: contribución de las variables auxiliares (Shapley)
  "05_1_shapley_kicb2_m1_original.R"         = TRUE,
  "05_2_shapley_mse_m1_log.R"                = TRUE,
  "05_3_shapley_cv_m1_log.R"                 = TRUE,

  # Estimaciones finales y mapas
  "06_0_estimaciones_coef_cv_m1_log.R"       = TRUE,
  "06_1_mapas_m1_log.R"                      = TRUE
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
