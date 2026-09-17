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
# Uso: abrir el proyecto de RStudio y ejecutar este archivo completo
# (source o Ctrl+Shift+Enter).
# ==============================================================================

library(here)
library(callr)

# Carpeta donde estan los scripts, relativa a la raiz del proyecto.
CARPETA <- here("scripts")

# Orden de ejecucion. Poner FALSE para omitir un paso cuyas salidas ya existen.
# 02_1c tarda cerca de tres horas: si output/seleccion_stepwise_completo_both.rds
# ya existe, dejarlo en FALSE.
PASOS <- c(
  "01_etapa1_construccion_base.R"            = TRUE,
  "02_1c_seleccion_stepwise.R"               = FALSE,
  "02_2a_modelado_variantes.R"               = TRUE,
  "02_2a_modelado_variantes_diagnostico.R"   = TRUE,
  "02_3_espacial.R"                          = TRUE,
  "02_4_robusto.R"                           = TRUE,
  "02_4_robusto_diagnostico.R"               = TRUE,
  "03_1_normalidad.R"                        = TRUE
)

faltan <- names(PASOS)[!file.exists(file.path(CARPETA, names(PASOS)))]
if (length(faltan) > 0) {
  stop("No se encontraron en ", CARPETA, ": ", paste(faltan, collapse = ", "))
}

if (!PASOS[["02_1c_seleccion_stepwise.R"]] &&
    !file.exists(here("output", "seleccion_stepwise_completo_both.rds"))) {
  stop("La seleccion esta desactivada pero no existe ",
       "output/seleccion_stepwise_completo_both.rds.")
}

t_total <- Sys.time()

for (s in names(PASOS)[PASOS]) {

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
