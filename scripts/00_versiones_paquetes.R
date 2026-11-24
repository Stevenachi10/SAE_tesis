# ==============================================================================
# 00_versiones_paquetes.R
# ------------------------------------------------------------------------------
# Registra la versión de R y las versiones de los paquetes empleados por los
# scripts del análisis. Escribe versiones_paquetes.txt en la raíz del proyecto.
#
# Uso: abrir SAEtesis.Rproj y ejecutar este archivo completo. Debe correrse en
# el mismo equipo y con la misma instalación de R con que se obtuvieron los
# resultados de la tesis.
# ==============================================================================

paquetes <- c(
  "callr", "car", "dplyr", "emdi", "ggplot2", "ggspatial", "here", "jsonlite",
  "knitr", "lmtest", "Matrix", "patchwork", "readxl", "robustbase", "sae",
  "saeRobust", "SAEval", "sandwich", "scales", "sf", "shadowtext", "spdep",
  "stringr", "tidyr", "TUvalues", "VIM"
)

version_de <- function(p) {
  v <- tryCatch(as.character(utils::packageVersion(p)),
                error = function(e) "NO INSTALADO")
  v
}

tabla <- data.frame(
  paquete = paquetes,
  version = vapply(paquetes, version_de, character(1)),
  row.names = NULL
)

salida <- here::here("versiones_paquetes.txt")

lineas <- c(
  "Versiones con que se obtuvieron los resultados de la tesis",
  paste("Fecha:", format(Sys.Date())),
  paste("R:", R.version.string),
  paste("Sistema:", utils::sessionInfo()$running),
  "",
  sprintf("%-12s %s", "Paquete", "Versión"),
  sprintf("%-12s %s", tabla$paquete, tabla$version)
)

writeLines(lineas, salida, useBytes = TRUE)
cat(lineas, sep = "\n")
cat("\nGuardado en:", salida, "\n")

if (any(tabla$version == "NO INSTALADO")) {
  warning("Hay paquetes sin instalar: ",
          paste(tabla$paquete[tabla$version == "NO INSTALADO"], collapse = ", "))
}
