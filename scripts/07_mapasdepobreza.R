# ============================================================
# 07 - ESTIMACIÓN M09, EBLUP Y FIGURAS DE RESULTADOS
# Tesis SAE - Pobreza Monetaria Municipal (Cauca y Valle del Cauca)
# ============================================================
# Entrada : output/matriz_sae_transformada_v2.rds
#           data/pivoteadas/municipios_cauca_valle_84_detallado.geojson
# Salidas : output/mapa_eblup_pobreza.png
#           output/mapa_cv_eblup.png
#           output/mapa_gamma.png
#           output/fig_cv_directo_vs_eblup.png
# ------------------------------------------------------------
# Reserva estadística: el estimador directo NO se cartografía;
# solo se usa de forma agregada (boxplot de CV).
# ============================================================

library(sae)
library(dplyr)
library(ggplot2)
library(sf)
library(ggspatial)   # install.packages("ggspatial") una vez si no está
library(here)

ruta_out <- here("output")
matriz   <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

# Paleta (consistente con la figura de Mahalanobis)
COL_PUNTO <- "#b5482e"; COL_BANDA  <- "#c9c2b3"
COL_MIN   <- "#2c6e6b"; COL_OSCURO <- "#1a1814"
COL_CLARO <- "#ece6da"

# ============================================================
# 1. AJUSTE DE M09 Y EXTRACCIÓN DEL EBLUP
# ============================================================
vars_m09 <- c("extension", "saber11_lectura_2022", "partos_calificados_2020",
              "dist_popayan_km", "dist_cali_km", "pct_contributivo",
              "rural_disperso", "tiene_coca", "desplaz_tasa_1000", "categoria_617")

formula_m09 <- as.formula(paste("pobreza_monetaria ~", paste(vars_m09, collapse = " + ")))

datos <- matriz %>%
  select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza, cvlog, all_of(vars_m09))

modelo <- mseFH(formula_m09, vardir = varianza_pobreza, data = datos)

resultados <- datos %>%
  mutate(
    eblup    = modelo$est$eblup,
    mse      = modelo$mse,
    cv_eblup = sqrt(mse) / eblup * 100,                          # %
    cv_dir   = cvlog * 100,                                      # cvlog fracción -> %
    gamma    = modelo$est$fit$refvar / (modelo$est$fit$refvar + varianza_pobreza)
  ) %>%
  select(cod_mun, Municipio, directo = pobreza_monetaria, eblup,
         cv_dir, cv_eblup, gamma) %>%
  mutate(cod_mun = as.character(cod_mun))

cat("CV directo:", round(mean(resultados$cv_dir, na.rm = TRUE), 2),
    "% | CV EBLUP:", round(mean(resultados$cv_eblup, na.rm = TRUE), 2), "%\n")

# ============================================================
# 2. UNIR AL SHAPEFILE DETALLADO (84 municipios, con cod_mun)
# ============================================================
ruta_shp <- here("data", "pivoteadas", "municipios_cauca_valle_84_detallado.geojson")
shp <- st_read(ruta_shp, quiet = TRUE) %>%
  st_make_valid() %>%
  mutate(cod_mun = as.character(cod_mun))

mapa <- shp %>% left_join(resultados, by = "cod_mun")
cat("Municipios sin estimación (debe ser 0):", sum(is.na(mapa$eblup)), "\n")

# Bordes departamentales (disolver municipios por departamento)
deptos <- mapa %>% group_by(departamento) %>% summarise(.groups = "drop")

# Tema base y capas comunes a los tres mapas
tema_mapa <- theme_void(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13, color = COL_OSCURO),
        plot.caption  = element_text(size = 8, color = "#6b6b6b", hjust = 0,
                                     margin = margin(t = 8)),
        legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"),
        plot.title.position = "plot")

capas_base <- list(
  geom_sf(data = deptos, fill = NA, color = COL_OSCURO, linewidth = 0.5),
  annotation_scale(location = "bl", height = unit(0.15, "cm"), text_cex = 0.6),
  annotation_north_arrow(location = "tr",
                         height = unit(0.9, "cm"), width = unit(0.9, "cm"),
                         style = north_arrow_minimal()),
  tema_mapa
)

# ============================================================
# 3. MAPA A - Pobreza monetaria estimada (EBLUP)   [CUERPO]
# ============================================================
g_eblup <- ggplot(mapa) +
  geom_sf(aes(fill = eblup), color = "white", linewidth = 0.08) +
  scale_fill_stepsn(colours = c(COL_CLARO, "#e3b9a3", "#d8956f", COL_PUNTO, "#7d2f1c"),
                    n.breaks = 6, labels = scales::percent_format(accuracy = 1),
                    name = "Pobreza (EBLUP)") +
  labs(title = "Pobreza monetaria municipal estimada (EBLUP Fay-Herriot)",
       caption = "Fuente: elaboración propia. Base cartográfica: DANE, MGN 2018.") +
  capas_base
ggsave(file.path(ruta_out, "mapa_eblup_pobreza.png"), g_eblup,
       width = 8, height = 8.5, dpi = 300, bg = "white")

# ============================================================
# 4. MAPA B - Precisión: CV del EBLUP   [CUERPO]
# ============================================================
g_cv <- ggplot(mapa) +
  geom_sf(aes(fill = cv_eblup), color = "white", linewidth = 0.08) +
  scale_fill_stepsn(colours = c(COL_MIN, "#7fa6a0", COL_BANDA, "#d8956f", COL_PUNTO),
                    n.breaks = 6, name = "CV (%)") +
  labs(title = "Precisión de las estimaciones: CV del EBLUP",
       caption = "Fuente: elaboración propia. Valores menores indican mayor precisión.") +
  capas_base
ggsave(file.path(ruta_out, "mapa_cv_eblup.png"), g_cv,
       width = 8, height = 8.5, dpi = 300, bg = "white")

# ============================================================
# 5. MAPA C - Factor de encogimiento gamma   [ANEXO / alternativa segura]
#    gamma bajo -> el área se apoya en el modelo (área pequeña)
#    gamma alto -> confía en su estimación directa (muestra grande)
# ============================================================
g_gamma <- ggplot(mapa) +
  geom_sf(aes(fill = gamma), color = "white", linewidth = 0.08) +
  scale_fill_stepsn(colours = c(COL_PUNTO, "#d8956f", COL_BANDA, "#7fa6a0", COL_MIN),
                    n.breaks = 6, name = expression(gamma)) +
  labs(title = expression(paste("Factor de encogimiento ", gamma, " por municipio")),
       caption = "Fuente: elaboración propia. Valores bajos: mayor apoyo en el modelo.") +
  capas_base
ggsave(file.path(ruta_out, "mapa_gamma.png"), g_gamma,
       width = 8, height = 8.5, dpi = 300, bg = "white")

# ============================================================
# 6. FIGURA - Reducción del CV (agregado, sin exponer el directo)
# ============================================================
df_cv <- data.frame(
  tipo = rep(c("Directo", "EBLUP"), each = nrow(resultados)),
  cv   = c(resultados$cv_dir, resultados$cv_eblup)
)

g_dist <- ggplot(df_cv, aes(x = tipo, y = cv, fill = tipo)) +
  geom_boxplot(width = 0.5, alpha = 0.85, outlier.color = COL_OSCURO) +
  scale_fill_manual(values = c("Directo" = COL_BANDA, "EBLUP" = COL_MIN), guide = "none") +
  labs(x = NULL, y = "Coeficiente de variación (%)",
       title = "Reducción del CV: estimador directo vs. EBLUP",
       caption = "Fuente: elaboración propia.") +
  theme_minimal(base_size = 11) +
  theme(plot.title   = element_text(face = "bold", size = 13, color = COL_OSCURO),
        plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0),
        panel.grid.minor = element_blank())
ggsave(file.path(ruta_out, "fig_cv_directo_vs_eblup.png"), g_dist,
       width = 6, height = 5, dpi = 300, bg = "white")

cat("\nFiguras guardadas en 'output/'.\n")


# ============================================================
# TABLA POR MUNICIPIO PARA EL ANEXO (EBLUP + CV)
# ============================================================
tabla_anexo <- resultados %>%
  mutate(
    departamento  = ifelse(substr(cod_mun, 1, 2) == "19", "Cauca", "Valle del Cauca"),
    pobreza_eblup = round(eblup * 100, 1),   # %
    cv            = round(cv_eblup, 1),       # %
    gamma         = round(gamma, 3)
  ) %>%
  arrange(departamento, Municipio) %>%
  select(departamento, Municipio, pobreza_eblup, cv, gamma)

print(tabla_anexo, row.names = FALSE)
write.csv(tabla_anexo, file.path(ruta_out, "tabla_eblup_anexo.csv"), row.names = FALSE)
