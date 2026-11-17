# ==============================================================================
# 02_3b_ganancia_espacial.R
# ------------------------------------------------------------------------------
# Ganancia de precision del modelo Fay-Herriot espacial frente al convencional.
# Lee output/espacial.rds (generado por 02_3_espacial.R) y no reestima nada.
#
# Para cada especificacion compara, sobre los dominios con estimacion positiva en
# ambos modelos, el MSE y el CV del modelo espacial con los del convencional.
# La reduccion relativa se define como (base - espacial) / base * 100: un valor
# positivo indica que el componente espacial mejora la precision.
#
# SALIDAS (output/)
#   esp_ganancia_mse.csv   Una fila por especificacion.
#   esp_ganancia_dominio.csv   Detalle por municipio y especificacion.
# ==============================================================================

library(dplyr)
library(here)

ruta_out <- here("output")
cat("02_3b_ganancia_espacial.R: ganancia de MSE y CV del modelo espacial\n")

esp <- readRDS(file.path(ruta_out, "espacial.rds"))

extraer <- function(r) {
  est <- as.numeric(r$ind$FH)
  mse <- as.numeric(r$MSE$FH)
  list(est = est, mse = mse,
       cv = ifelse(is.finite(est) & est > 0 & is.finite(mse),
                   sqrt(pmax(mse, 0)) / est * 100, NA_real_))
}

filas <- list(); detalle <- list()

for (nm in names(esp$ajustes_esp)) {
  
  ab <- esp$ajustes_base[[nm]]
  ae <- esp$ajustes_esp[[nm]]
  if (is.null(ab) || is.null(ae)) next
  
  b <- extraer(ab); s <- extraer(ae)
  comunes <- is.finite(b$cv) & is.finite(s$cv) &
    is.finite(b$mse) & is.finite(s$mse)
  
  red_mse_dom <- (b$mse - s$mse) / b$mse * 100
  
  rho <- esp$ajuste_queen$rho[esp$ajuste_queen$especificacion == nm][1]
  
  filas[[nm]] <- data.frame(
    especificacion = nm,
    rho = rho,
    n_comunes = sum(comunes),
    mse_medio_fh  = mean(b$mse[comunes]),
    mse_medio_sfh = mean(s$mse[comunes]),
    reduccion_mse_medio_pct = round((1 - mean(s$mse[comunes]) /
                                       mean(b$mse[comunes])) * 100, 2),
    mediana_reduccion_mse_dom_pct = round(median(red_mse_dom[comunes]), 2),
    n_dom_mse_menor_sfh = sum(s$mse[comunes] < b$mse[comunes]),
    cv_medio_fh  = round(mean(b$cv[comunes]), 3),
    cv_medio_sfh = round(mean(s$cv[comunes]), 3),
    reduccion_cv_medio_pct = round((1 - mean(s$cv[comunes]) /
                                      mean(b$cv[comunes])) * 100, 2),
    stringsAsFactors = FALSE)
  
  detalle[[nm]] <- data.frame(
    especificacion = nm, cod_mun = esp$cod_mun,
    mse_fh = b$mse, mse_sfh = s$mse, reduccion_mse_pct = round(red_mse_dom, 2),
    cv_fh = b$cv, cv_sfh = s$cv, stringsAsFactors = FALSE)
}

resumen <- do.call(rbind, filas)
print(resumen %>%
        mutate(across(c(mse_medio_fh, mse_medio_sfh), ~ signif(.x, 4))),
      row.names = FALSE)

write.csv(resumen, file.path(ruta_out, "esp_ganancia_mse.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(do.call(rbind, detalle), file.path(ruta_out, "esp_ganancia_dominio.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")