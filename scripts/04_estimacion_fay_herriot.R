# ============================================================
# 04 - ESTIMACIÓN FAY-HERRIOT, SELECCIÓN FINAL Y ROBUSTEZ
# Tesis SAE - Pobreza Monetaria Municipal (Cauca y Valle del Cauca)
# ============================================================
# Entrada : output/matriz_sae_transformada_v2.rds
# Salidas : output/comparacion_modelos_fh.csv
#           output/predicciones_eblup_M09.rds (.csv)
#           output/robustez_codificacion_ley617.csv
#           output/robustez_transformaciones.csv
#           output/robustez_coeficientes_transf.csv
# ------------------------------------------------------------
# Estructura:
#   A. Selección de especificaciones        (sae::mseFH)
#   B. Modelo final M09 y extracción EBLUP   (sae::mseFH)
#   C. Robustez 1: codificación de Cat.617   (sae::mseFH)
#   D. Robustez 2: transformaciones respuesta (emdi::fh)
# ============================================================

library(sae)
library(dplyr)
library(car)
library(here)
library(emdi)   # install.packages("emdi") una sola vez si no está instalado

# ------------------------------------------------------------
# Configuración de rutas y carga de datos
# ------------------------------------------------------------
ruta_out <- here("output")
matriz   <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

# Núcleo de 8 covariables (estable, datos corregidos)
nucleo8 <- c("extension", "saber11_lectura_2022", "partos_calificados_2020",
             "dist_popayan_km", "dist_cali_km", "pct_contributivo",
             "rural_disperso", "tiene_coca")

# ============================================================
# A. SELECCIÓN DE ESPECIFICACIONES (Fay-Herriot, paquete sae)
# ============================================================

# A.1 Modelos candidatos: núcleo + 1, 2 y 3 covariables adicionales
#     (se eliminó el antiguo M16, idéntico a M09)
modelos <- list(
  # Base
  "M01_nucleo8"              = nucleo8,
  # Núcleo + 1 candidata
  "M02_+desplaz"             = c(nucleo8, "desplaz_tasa_1000"),
  "M03_+cat617"              = c(nucleo8, "categoria_617"),
  "M04_+terciario"           = c(nucleo8, "pct_terciario_2023"),
  "M05_+iica"                = c(nucleo8, "iica_2023"),
  "M06_+mort_infantil"       = c(nucleo8, "mort_infantil_1a_2020"),
  "M07_+victimas"            = c(nucleo8, "victimas_tasa_1000"),
  "M08_+viol_intra"          = c(nucleo8, "tasa_violencia_intra_2019"),
  # Núcleo + 2 candidatas
  "M09_+desplaz+cat617"      = c(nucleo8, "desplaz_tasa_1000", "categoria_617"),
  "M10_+desplaz+terciario"   = c(nucleo8, "desplaz_tasa_1000", "pct_terciario_2023"),
  "M11_+cat617+terciario"    = c(nucleo8, "categoria_617", "pct_terciario_2023"),
  "M12_+desplaz+iica"        = c(nucleo8, "desplaz_tasa_1000", "iica_2023"),
  "M13_+cat617+iica"         = c(nucleo8, "categoria_617", "iica_2023"),
  # Núcleo + 3 candidatas
  "M14_+desplaz+terc+cat617" = c(nucleo8, "desplaz_tasa_1000", "pct_terciario_2023", "categoria_617"),
  "M15_+desplaz+cat617+iica" = c(nucleo8, "desplaz_tasa_1000", "categoria_617", "iica_2023")
)

# A.2 Función de evaluación: CV del EBLUP, significancia y VIF máximo
evaluar <- function(vars, nombre) {
  datos      <- matriz %>% dplyr::select(pobreza_monetaria, varianza_pobreza, all_of(vars))
  formula_fh <- as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
  
  modelo <- tryCatch(mseFH(formula_fh, vardir = varianza_pobreza, data = datos),
                     error = function(e) NULL)
  if (is.null(modelo)) {
    return(data.frame(modelo = nombre, n_vars = length(vars), conv = "ERROR",
                      cv_fh = NA, mse = NA, n_signif = NA, vif_max = NA))
  }
  
  cv_fh   <- mean(sqrt(modelo$mse) / modelo$est$eblup * 100)
  pvals   <- modelo$est$fit$estcoef$pvalue[-1]   # se excluye el intercepto
  vif_max <- tryCatch(
    max(vif(lm(formula_fh, data = matriz %>% dplyr::select(pobreza_monetaria, all_of(vars))))),
    error = function(e) NA)
  
  data.frame(modelo = nombre, n_vars = length(vars),
             conv = modelo$est$fit$convergence,
             cv_fh = round(cv_fh, 3), mse = round(mean(modelo$mse), 6),
             n_signif = paste0(sum(pvals < 0.05), "/", length(vars)),
             vif_max = round(vif_max, 2))
}

# A.3 Tabla comparativa ordenada por CV
comparacion <- do.call(rbind, lapply(names(modelos), \(n) evaluar(modelos[[n]], n)))
comparacion <- comparacion[order(comparacion$cv_fh), ]
cat("\n===== COMPARACIÓN DE MODELOS (ordenado por CV) =====\n")
print(comparacion, row.names = FALSE)
write.csv(comparacion, file.path(ruta_out, "comparacion_modelos_fh.csv"), row.names = FALSE)

# A.4 Coeficientes de modelos clave (para el cuerpo del documento)
clave <- c("M01_nucleo8", "M02_+desplaz", "M03_+cat617",
           "M09_+desplaz+cat617", "M10_+desplaz+terciario", "M14_+desplaz+terc+cat617")
cat("\n===== COEFICIENTES DE MODELOS CLAVE =====\n")
for (nombre in clave) {
  vars       <- modelos[[nombre]]
  datos      <- matriz %>% dplyr::select(pobreza_monetaria, varianza_pobreza, all_of(vars))
  formula_fh <- as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
  modelo     <- mseFH(formula_fh, vardir = varianza_pobreza, data = datos)
  cat("\n----- ", nombre, " | CV:",
      round(mean(sqrt(modelo$mse) / modelo$est$eblup * 100), 3), "% -----\n")
  print(round(modelo$est$fit$estcoef, 5))
}

# ============================================================
# B. MODELO FINAL M09 Y EXTRACCIÓN DEL EBLUP (sae)
# ============================================================
vars_finales  <- modelos[["M09_+desplaz+cat617"]]
formula_final <- as.formula(paste("pobreza_monetaria ~", paste(vars_finales, collapse = " + ")))
datos_finales <- matriz %>%
  dplyr::select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza, cvlog,
                all_of(vars_finales))

modelo_final <- mseFH(formula_final, vardir = varianza_pobreza, data = datos_finales)

resultados_eblup <- datos_finales %>%
  mutate(
    eblup_pobreza = modelo_final$est$eblup,
    mse_eblup     = modelo_final$mse,
    cv_eblup      = (sqrt(mse_eblup) / eblup_pobreza) * 100,
    gamma_fh      = modelo_final$est$fit$refvar /
      (modelo_final$est$fit$refvar + varianza_pobreza)
  ) %>%
  dplyr::select(cod_mun, Municipio,
                directo_pobreza = pobreza_monetaria, cv_directo = cvlog,
                eblup_pobreza, cv_eblup, gamma_fh)

cat("\nCV promedio directo:", round(mean(resultados_eblup$cv_directo * 100, na.rm = TRUE), 2), "%\n")
cat("CV promedio EBLUP  :", round(mean(resultados_eblup$cv_eblup, na.rm = TRUE), 2), "%\n")

saveRDS(resultados_eblup, file.path(ruta_out, "predicciones_eblup_M09.rds"))
write.csv(resultados_eblup, file.path(ruta_out, "predicciones_eblup_M09.csv"), row.names = FALSE)

# ============================================================
# C. ROBUSTEZ 1 - CODIFICACIÓN DE LA CATEGORÍA LEY 617 (sae)
#    Se mantiene fijo el resto de M09 y se cambia solo Cat.617
# ============================================================
matriz <- matriz %>%
  mutate(cat617_d6 = ifelse(categoria_617 == 6, 1L, 0L))   # 1 = sexta, 0 = resto

base_m09 <- c(nucleo8, "desplaz_tasa_1000")
modelos_cod <- list(
  "cat617_numerica"    = c(base_m09, "categoria_617"),   # codificación usada en M09
  "cat617_dummy_sexta" = c(base_m09, "cat617_d6")
)
tabla_cod <- do.call(rbind, lapply(names(modelos_cod), \(n) evaluar(modelos_cod[[n]], n)))
tabla_cod <- tabla_cod[order(tabla_cod$cv_fh), ]
cat("\n===== ROBUSTEZ: codificación de Cat.617 =====\n")
print(tabla_cod, row.names = FALSE)
write.csv(tabla_cod, file.path(ruta_out, "robustez_codificacion_ley617.csv"), row.names = FALSE)

# ============================================================
# D. ROBUSTEZ 2 - TRANSFORMACIONES DE LA RESPUESTA (emdi)
#    Se reajusta TODO en emdi para que los CV sean comparables:
#    cada modelo retro-transforma con corrección de sesgo y el CV
#    queda expresado en la escala original de la pobreza.
# ============================================================

# D.1 emdi exige el tamaño de muestra efectivo para la transf. arcoseno.
#     Se deriva de la varianza de diseño: Var(p) ~= p(1-p)/n_eff
datos_emdi <- datos_finales %>%
  mutate(n_eff = pobreza_monetaria * (1 - pobreza_monetaria) / varianza_pobreza)

f_m09 <- pobreza_monetaria ~ extension + saber11_lectura_2022 +
  partos_calificados_2020 + dist_popayan_km + dist_cali_km +
  pct_contributivo + rural_disperso + tiene_coca +
  desplaz_tasa_1000 + categoria_617

# D.2 Sin transformación (línea base, mismo marco que las otras dos; REML)
fh_no <- fh(f_m09, vardir = "varianza_pobreza", combined_data = datos_emdi,
            domains = "cod_mun", transformation = "no",
            MSE = TRUE, mse_type = "analytical")

# D.3 Logarítmica. Combinación válida en emdi:
#     backtransformation = "bc_sm" (Slud-Maiti) -> exige method = "ml" y mse_type = "analytical".
#     La diferencia ML/REML es de segundo orden y no altera el ordenamiento de los CV.
fh_log <- fh(f_m09, vardir = "varianza_pobreza", combined_data = datos_emdi,
             domains = "cod_mun", transformation = "log",
             backtransformation = "bc_sm", method = "ml",
             MSE = TRUE, mse_type = "analytical")

# D.4 Arcoseno (estabilizadora de varianza para proporciones). Usa n_eff; MSE por bootstrap.
fh_arc <- fh(f_m09, vardir = "varianza_pobreza", combined_data = datos_emdi,
             domains = "cod_mun", transformation = "arcsin",
             backtransformation = "bc", eff_smpsize = "n_eff",
             MSE = TRUE, mse_type = "boot", B = c(200, 0))

# D.5 CV promedio en ESCALA NATURAL (emdi ya retro-transformó)
cv_prom <- function(m) mean(sqrt(m$MSE$FH) / m$ind$FH, na.rm = TRUE) * 100
cv_transf <- round(c(sin_transf = cv_prom(fh_no),
                     log        = cv_prom(fh_log),
                     arcoseno   = cv_prom(fh_arc)), 3)
cat("\n===== ROBUSTEZ: CV por transformación (escala de pobreza) =====\n")
print(cv_transf)
write.csv(data.frame(transformacion = names(cv_transf), cv = as.numeric(cv_transf)),
          file.path(ruta_out, "robustez_transformaciones.csv"), row.names = FALSE)

# D.6 Tabla comparativa de coeficientes (beta + estrellas), en R base.
#     OJO: los coef. de log/arcoseno están en su escala transformada;
#     se compara SIGNO y SIGNIFICANCIA entre columnas, NO la magnitud.
#     Convención de estrellas: *** p<0,01   ** p<0,05   * p<0,10
coef_celda <- function(m) {
  cf  <- as.data.frame(m$model$coefficients)
  sig <- as.character(cut(cf$p.value, c(-Inf, .01, .05, .1, Inf),
                          labels = c("***", "**", "*", "")))
  data.frame(variable = rownames(cf),
             celda = paste0(formatC(cf$coefficients, format = "f", digits = 5), sig),
             stringsAsFactors = FALSE)
}
t_no  <- coef_celda(fh_no);  names(t_no)[2]  <- "Sin_transf"
t_log <- coef_celda(fh_log); names(t_log)[2] <- "Log"
t_arc <- coef_celda(fh_arc); names(t_arc)[2] <- "Arcoseno"

tabla_coefs <- Reduce(function(a, b) merge(a, b, by = "variable", all = TRUE),
                      list(t_no, t_log, t_arc))
orden       <- coef_celda(fh_no)$variable        # orden natural (Intercept primero)
tabla_coefs <- tabla_coefs[match(orden, tabla_coefs$variable), ]
cat("\n===== ROBUSTEZ: coeficientes por transformación =====\n")
print(tabla_coefs, row.names = FALSE)
write.csv(tabla_coefs, file.path(ruta_out, "robustez_coeficientes_transf.csv"), row.names = FALSE)

cat("\nListo. Resultados y robustez guardados en 'output/'.\n")


# ============================================================
# ROBUSTEZ 3 - Sensibilidad a observaciones influyentes (M09)
# ============================================================

# 1. Conjuntos a excluir
cod_capitales <- c(76001, 19001)                 # Cali, Popayán

# Atípicos por código DANE (recomputa la Mahalanobis -> cod_mun, que es único)
vars_continuas <- c("extension","saber11_lectura_2022","partos_calificados_2020",
                    "dist_popayan_km","dist_cali_km","desplaz_tasa_1000","pct_contributivo")
Xc  <- matriz %>% select(all_of(vars_continuas)) %>% as.matrix()
dM  <- mahalanobis(Xc, colMeans(Xc), cov(Xc))
p   <- length(vars_continuas)
atip95 <- matriz$cod_mun[dM > qchisq(0.95, df = p)]   # 7 municipios
atip99 <- matriz$cod_mun[dM > qchisq(0.99, df = p)]   # 3 municipios

datos_conf <- list(
  "Completo"        = datos_finales,
  "Sin capitales"   = datos_finales %>% filter(!cod_mun %in% cod_capitales),
  "Sin atipicos 95" = datos_finales %>% filter(!cod_mun %in% atip95),
  "Sin atipicos 99" = datos_finales %>% filter(!cod_mun %in% atip99)
)
sapply(datos_conf, nrow)            # ahora sí: 84 / 82 / 77 / 81

# 2. Ajuste de M09 en cada configuración
fit_conf <- function(datos) {
  m <- mseFH(formula_final, vardir = varianza_pobreza, data = datos)
  list(coef = m$est$fit$estcoef,
       resumen = data.frame(n = nrow(datos),
                            sigma2_u = round(m$est$fit$refvar, 6),
                            cv_fh    = round(mean(sqrt(m$mse) / m$est$eblup * 100), 3)))
}
ajustes <- lapply(datos_conf, fit_conf)

# 3. Resumen por configuración
resumen <- do.call(rbind, lapply(names(ajustes),
                                 \(k) cbind(config = k, ajustes[[k]]$resumen)))
cat("\n===== Sensibilidad: resumen =====\n"); print(resumen, row.names = FALSE)

# 4. Coeficientes comparados (*** p<0,01  ** p<0,05  * p<0,10)
coef_celda2 <- function(cf) {
  cf  <- as.data.frame(cf)
  sig <- as.character(cut(cf$pvalue, c(-Inf, .01, .05, .1, Inf),
                          labels = c("***", "**", "*", "")))
  out <- paste0(formatC(cf$beta, format = "f", digits = 5), sig)
  names(out) <- sub("^X", "", rownames(cf))   # mseFH antepone "X" al nombre
  out
}
tab_coef <- sapply(ajustes, \(a) coef_celda2(a$coef))
cat("\n===== Sensibilidad: coeficientes =====\n"); print(tab_coef, quote = FALSE)
