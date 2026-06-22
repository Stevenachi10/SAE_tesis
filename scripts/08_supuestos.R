# ============================================================
# 08 DIAGNÓSTICOS DEL MODELO FAY-HERRIOT (M09):
#       SUPUESTOS, AUTOCORRELACIÓN ESPACIAL Y MODELO SFH
# Tesis SAE - Pobreza Monetaria Municipal (Cauca y Valle del Cauca)
# ============================================================
# Entradas: output/matriz_sae_transformada_v2.rds
#           data/pivoteadas/municipios_cauca_valle_84_detallado.geojson
# Salidas : output/fig_supuestos.png
#           (consola) test de Moran, ajuste SFH y comparación de modelos
# ------------------------------------------------------------
# Flujo:
#   1. Reajuste de M09 en emdi -> residuos y efectos aleatorios
#   2. Panel de supuestos (4 diagnósticos)
#   3. Matriz de vecindad W ALINEADA al orden de los datos
#   4. I de Moran sobre los residuos (autocorrelación espacial)
#   5. Fay-Herriot espacial (SFH) y comparación directo / FH / SFH
# ============================================================

library(emdi)
library(sae)
library(dplyr)
library(ggplot2)
library(patchwork)
library(sf)
library(spdep)
library(here)

ruta_out <- here("output")
matriz   <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

COL_PUNTO <- "#b5482e"; COL_BANDA  <- "#c9c2b3"
COL_MIN   <- "#2c6e6b"; COL_OSCURO <- "#1a1814"

# ------------------------------------------------------------
# 1. Reajuste de M09 en emdi (sin transformación = mismo modelo)
# ------------------------------------------------------------
datos_emdi <- matriz %>% mutate(cod_mun = as.character(cod_mun))

f_m09 <- pobreza_monetaria ~ extension + saber11_lectura_2022 +
  partos_calificados_2020 + dist_popayan_km + dist_cali_km +
  pct_contributivo + rural_disperso + tiene_coca +
  desplaz_tasa_1000 + categoria_617

fh_m09 <- fh(f_m09, vardir = "varianza_pobreza", combined_data = datos_emdi,
             domains = "cod_mun", transformation = "no",
             MSE = TRUE, mse_type = "analytical")

res_std <- as.numeric(fh_m09$model$std_real_residuals)   
re      <- as.numeric(fh_m09$model$random_effects)        
ajust   <- as.numeric(fh_m09$ind$FH)                     

p_res <- shapiro.test(res_std)$p.value
p_re  <- shapiro.test(re)$p.value
cat(sprintf("Shapiro residuos: p = %.4f | Shapiro efectos aleatorios: p = %.2e\n",
            p_res, p_re))

# ------------------------------------------------------------
# 2. Panel de supuestos
# ------------------------------------------------------------
tema <- theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 11, color = COL_OSCURO),
        plot.subtitle = element_text(size = 8.5, color = "#6b6b6b"),
        panel.grid.minor = element_blank())

# (a) QQ de los residuos estandarizados
p1 <- ggplot(data.frame(r = res_std), aes(sample = r)) +
  stat_qq(color = COL_MIN, size = 1.4, alpha = 0.85) +
  stat_qq_line(color = COL_OSCURO, linewidth = 0.5) +
  labs(title = "Normalidad de los residuos",
       subtitle = sprintf("Shapiro-Wilk: p = %.3f", p_res),
       x = "Cuantiles teóricos", y = "Residuos estandarizados") + tema

# (b) QQ de los efectos aleatorios
p2 <- ggplot(data.frame(u = re), aes(sample = u)) +
  stat_qq(color = COL_PUNTO, size = 1.4, alpha = 0.85) +
  stat_qq_line(color = COL_OSCURO, linewidth = 0.5) +
  labs(title = "Normalidad de los efectos aleatorios",
       subtitle = sprintf("Shapiro-Wilk: p = %.2e", p_re),
       x = "Cuantiles teóricos", y = "Efectos aleatorios") + tema

# (c) Residuos vs. valores ajustados (homocedasticidad / patrón)
p3 <- ggplot(data.frame(aj = ajust, r = res_std), aes(aj, r)) +
  geom_hline(yintercept = 0, color = COL_OSCURO, linewidth = 0.4) +
  geom_point(color = COL_MIN, size = 1.4, alpha = 0.85) +
  geom_smooth(method = "loess", se = FALSE, color = COL_PUNTO, linewidth = 0.6) +
  labs(title = "Residuos vs. valores ajustados",
       subtitle = "Sin patrón ni embudo = homocedasticidad",
       x = "EBLUP (ajustado)", y = "Residuos estandarizados") + tema

# (d) Densidad de los efectos aleatorios vs. normal teórica
p4 <- ggplot(data.frame(u = re), aes(u)) +
  geom_density(fill = COL_PUNTO, alpha = 0.25, color = COL_PUNTO, linewidth = 0.6) +
  stat_function(fun = dnorm, args = list(mean = mean(re), sd = sd(re)),
                color = COL_OSCURO, linewidth = 0.5, linetype = "dashed") +
  labs(title = "Distribución de los efectos aleatorios",
       subtitle = "Línea discontinua: normal teórica",
       x = "Efectos aleatorios", y = "Densidad") + tema

g_supuestos <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title   = "Diagnóstico de supuestos del modelo Fay-Herriot (M09)",
    caption = "Fuente: elaboración propia.",
    theme = theme(plot.title   = element_text(face = "bold", size = 13, color = COL_OSCURO),
                  plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0)))

ggsave(file.path(ruta_out, "fig_supuestos.png"), g_supuestos,
       width = 9, height = 8, dpi = 300, bg = "white")
cat("Figura guardada: output/fig_supuestos.png\n")

# ------------------------------------------------------------
# 3. Matriz de vecindad W ALINEADA al orden de datos_emdi
#    OJO: eblupSFH empareja 'data' y 'proxmat' por POSICIÓN de fila,
#    así que el shapefile debe ir en el mismo orden que los datos.
# ------------------------------------------------------------
shp <- st_read(here("data", "pivoteadas",
                    "municipios_cauca_valle_84_detallado.geojson"), quiet = TRUE) %>%
  st_make_valid() %>%
  mutate(cod_mun = as.character(cod_mun))


shp <- shp[match(datos_emdi$cod_mun, shp$cod_mun), ]
stopifnot(all(shp$cod_mun == datos_emdi$cod_mun))  
shp$res_std <- res_std                              
# Vecindad por contigüidad reina.
# 
nb <- poly2nb(shp, queen = TRUE, snap = 0.001)      
print(nb)  
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
W  <- nb2mat(nb,  style = "W", zero.policy = TRUE) 

# ------------------------------------------------------------
# 4. I de Moran sobre los residuos estandarizados
# ------------------------------------------------------------
moran <- moran.test(shp$res_std, lw, zero.policy = TRUE)
print(moran)

# ------------------------------------------------------------
# 5. Fay-Herriot espacial (SFH) y comparación de modelos
#
# ------------------------------------------------------------
f_sfh <- pobreza_monetaria ~ extension + saber11_lectura_2022 +
  partos_calificados_2020 + dist_popayan_km + dist_cali_km +
  pct_contributivo + rural_disperso + tiene_coca +
  desplaz_tasa_1000 + categoria_617

sfh <- eblupSFH(f_sfh, vardir = varianza_pobreza, proxmat = W,
                method = "REML", data = datos_emdi)

cat("\n--- Fay-Herriot espacial (SFH) ---\n")
cat(sprintf("rho (correlacion espacial): %.4f\n", sfh$fit$spatialcorr))
cat(sprintf("Convergencia: %s | Iteraciones: %s\n",
            sfh$fit$convergence, sfh$fit$iterations))
cat("Coeficientes:\n");  print(round(sfh$fit$estcoef, 5))
cat(sprintf("refvar (sigma2_u SFH): %.6f\n", sfh$fit$refvar))
cat("Bondad de ajuste:\n"); print(sfh$fit$goodness)

# MSE / CV analítico del SFH (paralelo al mseFH analítico)
mse_sfh <- mseSFH(f_sfh, vardir = varianza_pobreza, proxmat = W,
                  method = "REML", data = datos_emdi)
cv_sfh  <- mean(sqrt(mse_sfh$mse) / mse_sfh$est$eblup) * 100

# CV del directo y del FH no espacial (para la tabla comparativa)
cv_dir <- mean(sqrt(matriz$varianza_pobreza) / matriz$pobreza_monetaria) * 100
cv_fh  <- mean(sqrt(fh_m09$MSE$FH) / fh_m09$ind$FH) * 100
sigma2u_fh <- tryCatch(fh_m09$model$variance, error = function(e) NA_real_)

tabla_modelos <- data.frame(
  Modelo  = c("Directo", "FH no espacial (M09)", "FH espacial (SFH)"),
  CV_pct  = round(c(cv_dir, cv_fh, cv_sfh), 2),
  sigma2u = round(c(NA, sigma2u_fh, sfh$fit$refvar), 6),
  rho     = round(c(NA, NA, sfh$fit$spatialcorr), 4)
)
cat("\n--- Comparación de modelos ---\n")
print(tabla_modelos, row.names = FALSE)