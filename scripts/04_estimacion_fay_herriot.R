# ============================================================
# 04 - ESTIMACIÓN FAY-HERRIOT Y SELECCIÓN FINAL
# Tesis SAE - Pobreza Monetaria Municipal (Cauca y Valle)
# ============================================================
# Entrada: matriz_sae_transformada_v2.rds
# Salida: comparacion_modelos_fh.csv, predicciones_eblup_M09.rds
# ============================================================

library(sae)
library(dplyr)
library(car)
library(here) # Para manejo dinámico de rutas

# ============================================================
# CONFIGURACIÓN DE RUTAS 
# ============================================================
ruta_out <- here("output")
matriz <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

# Núcleo de 8 (en los 3 métodos, datos corregidos)
nucleo8 <- c("extension","saber11_lectura_2022","partos_calificados_2020",
             "dist_popayan_km","dist_cali_km","pct_contributivo",
             "rural_disperso","tiene_coca")

# ============================================================
# 1. DEFINIR MODELOS A COMPARAR
# ============================================================
modelos <- list(
  # --- Base ---
  "M01_nucleo8"               = nucleo8,
  
  # --- Núcleo + una candidata de 2 métodos ---
  "M02_+desplaz"              = c(nucleo8, "desplaz_tasa_1000"),
  "M03_+categoria617"         = c(nucleo8, "categoria_617"),
  "M04_+terciario"            = c(nucleo8, "pct_terciario_2023"),
  "M05_+iica"                 = c(nucleo8, "iica_2023"),
  "M06_+mort_infantil"        = c(nucleo8, "mort_infantil_1a_2020"),
  "M07_+victimas"             = c(nucleo8, "victimas_tasa_1000"),
  "M08_+viol_intra"           = c(nucleo8, "tasa_violencia_intra_2019"),
  
  # --- Núcleo + dos candidatas ---
  "M09_+desplaz+cat617"       = c(nucleo8, "desplaz_tasa_1000", "categoria_617"),
  "M10_+desplaz+terciario"    = c(nucleo8, "desplaz_tasa_1000", "pct_terciario_2023"),
  "M11_+cat617+terciario"     = c(nucleo8, "categoria_617", "pct_terciario_2023"),
  "M12_+desplaz+iica"         = c(nucleo8, "desplaz_tasa_1000", "iica_2023"),
  "M13_+cat617+iica"          = c(nucleo8, "categoria_617", "iica_2023"),
  
  # --- Núcleo + tres candidatas ---
  "M14_+desplaz+terc+cat617"  = c(nucleo8, "desplaz_tasa_1000", "pct_terciario_2023", "categoria_617"),
  "M15_+desplaz+cat617+iica"  = c(nucleo8, "desplaz_tasa_1000", "categoria_617", "iica_2023"),
  "M16_nucleo9+cat617"        = c(nucleo8, "desplaz_tasa_1000", "categoria_617")
)

# ============================================================
# 2. FUNCIÓN DE EVALUACIÓN
# ============================================================
evaluar <- function(vars, nombre) {
  datos <- matriz %>% dplyr::select(pobreza_monetaria, varianza_pobreza, all_of(vars))
  formula_fh <- as.formula(paste("pobreza_monetaria ~", paste(vars, collapse=" + ")))
  
  modelo <- tryCatch(mseFH(formula_fh, vardir = varianza_pobreza, data = datos),
                     error = function(e) NULL)
  
  if (is.null(modelo)) return(data.frame(modelo=nombre, n_vars=length(vars),
                                         conv="ERROR", cv_fh=NA, mse=NA, n_signif=NA, vif_max=NA))
  
  cv_fh <- mean(sqrt(modelo$mse)/modelo$est$eblup*100)
  pvals <- modelo$est$fit$estcoef$pvalue[-1]
  
  vif_max <- tryCatch(max(vif(lm(as.formula(paste("pobreza_monetaria ~", paste(vars, collapse="+"))),
                                 data = matriz %>% dplyr::select(pobreza_monetaria, all_of(vars))))),
                      error=function(e) NA)
  
  data.frame(modelo=nombre, n_vars=length(vars), conv=modelo$est$fit$convergence,
             cv_fh=round(cv_fh,3), mse=round(mean(modelo$mse),6),
             n_signif=paste0(sum(pvals<0.05),"/",length(vars)),
             vif_max=round(vif_max,2))
}

# ============================================================
# 3. TABLA COMPARATIVA
# ============================================================
comparacion <- do.call(rbind, lapply(names(modelos), function(n) evaluar(modelos[[n]], n)))
comparacion <- comparacion[order(comparacion$cv_fh), ]

cat("\n========== COMPARACIÓN DE MODELOS (ordenado por CV) ==========\n")
print(comparacion, row.names = FALSE)

# Guardar la tabla comparativa para el documento de tesis
write.csv(comparacion, file.path(ruta_out, "comparacion_modelos_fh.csv"), row.names = FALSE)

# ============================================================
# 4. SUMMARY DE MODELOS PROMETEDORES
# ============================================================
prometedores <- c("M01_nucleo8","M02_+desplaz","M03_+categoria617",
                  "M09_+desplaz+cat617","M10_+desplaz+terciario","M14_+desplaz+terc+cat617")

cat("\n\n========== COEFICIENTES DE MODELOS CLAVE ==========\n")
for (nombre in prometedores) {
  vars <- modelos[[nombre]]
  datos <- matriz %>% dplyr::select(pobreza_monetaria, varianza_pobreza, all_of(vars))
  formula_fh <- as.formula(paste("pobreza_monetaria ~", paste(vars, collapse=" + ")))
  modelo <- tryCatch(mseFH(formula_fh, vardir = varianza_pobreza, data = datos), error=function(e) NULL)
  
  cat("\n=================================================\n")
  cat("MODELO:", nombre, "| CV promedio:", round(mean(sqrt(modelo$mse)/modelo$est$eblup*100),3), "%\n")
  cat("=================================================\n")
  if (!is.null(modelo)) print(round(modelo$est$fit$estcoef, 5))
}

# ============================================================
# 5. ESTIMACIÓN FINAL Y EXTRACCIÓN DEL EBLUP (MODELO M09)
# ============================================================
cat("\n\n========== EXTRACCIÓN DEL MODELO FINAL (M09) ==========\n")

vars_finales <- modelos[["M09_+desplaz+cat617"]]
formula_final <- as.formula(paste("pobreza_monetaria ~", paste(vars_finales, collapse=" + ")))
datos_finales <- matriz %>% dplyr::select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza, cvlog, all_of(vars_finales))

modelo_final <- mseFH(formula_final, vardir = varianza_pobreza, data = datos_finales)


# Ensamblar el dataset final con las estimaciones
resultados_eblup <- datos_finales %>%
  mutate(
    eblup_pobreza = modelo_final$est$eblup,
    mse_eblup     = modelo_final$mse,
    cv_eblup      = (sqrt(mse_eblup) / eblup_pobreza) * 100,
    gamma_fh      = modelo_final$est$fit$refvar / (modelo_final$est$fit$refvar + varianza_pobreza)
  ) %>%
  dplyr::select(cod_mun, Municipio, 
                directo_pobreza = pobreza_monetaria, cv_directo = cvlog,
                eblup_pobreza, cv_eblup, gamma_fh)

# Mostrar resumen de ganancia en precisión
cat("CV Promedio Directo:", round(mean(resultados_eblup$cv_directo * 100, na.rm=TRUE), 2), "%\n")
cat("CV Promedio EBLUP:  ", round(mean(resultados_eblup$cv_eblup, na.rm=TRUE), 2), "%\n")

# Guardar resultados finales
saveRDS(resultados_eblup, file.path(ruta_out, "predicciones_eblup_M09.rds"))
write.csv(resultados_eblup, file.path(ruta_out, "predicciones_eblup_M09.csv"), row.names = FALSE)

cat("\nEstimaciones finales extraídas y guardadas exitosamente en 'output/'.\n")