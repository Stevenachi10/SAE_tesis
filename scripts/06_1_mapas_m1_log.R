# ==============================================================================
# 06_1_mapas_m1_log.R
# ------------------------------------------------------------------------------
# Etapa IV: estimaciones municipales finales de M1 log.
#
#   1. Coeficientes estimados de M1 (sin valores p).
#   2. Estimación puntual y coeficiente de variación por municipio.
#   3. Mapas lado a lado: pobreza estimada (a) y CV (b).
#
# ENTRADAS
#   output/modelado_variantes.rds          (ajuste M1_reml_log)
#   output/matriz_sae_transformada_v2.rds
#   data/pivoteadas/MGN2025_MPIO_GRAFICO/MGN_ADM_MPIO_GRAFICO.shp
#
# SALIDAS (output/)
#   final_coeficientes_m1_log.csv
#   final_estimaciones_m1_log.csv
#   fig_mapas_m1_log.png
# ==============================================================================


# ==============================================================================
# 0. PAQUETES Y PARÁMETROS
# ==============================================================================

library(emdi)
library(sf)
library(dplyr)
library(ggplot2)
library(here)

stopifnot(
  "Falta el paquete patchwork: install.packages(\"patchwork\")" =
    requireNamespace("patchwork", quietly = TRUE)
)

library(patchwork)

# Escala gráfica y flecha de norte (opcional pero recomendado)
tiene_ggspatial <- requireNamespace("ggspatial", quietly = TRUE)

if (!tiene_ggspatial) {
  message("ggspatial no está instalado: el mapa sale sin escala ni norte. ",
          "Instale con install.packages(\"ggspatial\").")
}

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")
ruta_car <- here("data", "pivoteadas", "MGN2025_MPIO_GRAFICO")

AJUSTE_FINAL   <- "M1_reml_log"
CRS_PROYECTADO <- 9377

# Etiquetas con coma decimal
coma <- function(x) gsub(".", ",", format(x, trim = TRUE), fixed = TRUE)

tema_mapa <- theme_void(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 9),
    plot.tag        = element_text(size = 12),
    plot.margin     = margin(4, 8, 4, 8)
  )

# Capitales departamentales. Coordenadas aproximadas de la cabecera urbana
# (WGS84); solo sirven para ubicar el punto en el mapa.
capitales <- data.frame(
  nombre = c("Cali", "Popayán"),
  lon    = c(-76.532, -76.606),
  lat    = c(3.452, 2.444)
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(CRS_PROYECTADO)


# ==============================================================================
# 1. CARGA DEL AJUSTE FINAL
# ==============================================================================

variantes <- readRDS(file.path(ruta_out, "modelado_variantes.rds"))
matriz    <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

aj <- variantes$ajustes[[AJUSTE_FINAL]]

stopifnot("No se encuentra el ajuste final." = !is.null(aj))

cat("\n============================================================\n")
cat("ESTIMACIONES FINALES:", AJUSTE_FINAL, "\n")
cat("============================================================\n")
cat("Transformación:",
    paste(unlist(aj$transformation), collapse = " / "), "\n")
cat("Método de estimación:",
    paste(unlist(aj$method$method), collapse = " / "), "\n")


# ==============================================================================
# 2. COEFICIENTES (SIN VALORES P)
# ==============================================================================

coefs <- as.data.frame(aj$model$coefficients)

stopifnot(
  "La tabla de coeficientes no tiene la columna esperada." =
    "coefficients" %in% names(coefs)
)

tabla_coef <- data.frame(
  variable    = rownames(coefs),
  coeficiente = coefs$coefficients,
  row.names   = NULL,
  stringsAsFactors = FALSE
)

cat("\nCoeficientes (escala logarítmica):\n")
print(tabla_coef, row.names = FALSE, digits = 6)

cat("\nsigma_u^2 (escala logarítmica):",
    format(aj$model$variance, digits = 6), "\n")

write.csv(
  tabla_coef,
  file.path(ruta_out, "final_coeficientes_m1_log.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ==============================================================================
# 3. ESTIMACIÓN PUNTUAL Y CV POR MUNICIPIO
# ==============================================================================

est <- data.frame(
  cod_mun = as.character(aj$ind$Domain),
  directo = as.numeric(aj$ind$Direct),
  fh      = as.numeric(aj$ind$FH),
  mse_fh  = as.numeric(aj$MSE$FH),
  stringsAsFactors = FALSE
)

stopifnot(
  "El orden de dominios de ind y MSE no coincide." =
    identical(as.character(aj$MSE$Domain), est$cod_mun)
)

nombres <- matriz %>%
  as.data.frame() %>%
  select(cod_mun, Municipio) %>%
  mutate(cod_mun = as.character(cod_mun))

est <- est %>%
  left_join(nombres, by = "cod_mun") %>%
  mutate(
    departamento = ifelse(substr(cod_mun, 1, 2) == "19",
                          "Cauca", "Valle del Cauca"),
    directo_pct  = 100 * directo,
    pobreza_pct  = 100 * fh,
    cv_pct       = 100 * sqrt(mse_fh) / fh
  )

stopifnot(
  "No hay 84 municipios." = nrow(est) == 84,
  "Hay municipios sin nombre." = !anyNA(est$Municipio),
  "Hay CV no finitos." = all(is.finite(est$cv_pct))
)


# ==============================================================================
# 4. DESCRIPCIÓN
# ==============================================================================

cat("\n============================================================\n")
cat("POBREZA MONETARIA ESTIMADA (%)\n")
cat("============================================================\n")
print(summary(est$pobreza_pct))

cat("\nPor departamento (promedio simple de municipios):\n")
print(
  est %>%
    group_by(departamento) %>%
    summarise(
      n       = n(),
      minimo  = min(pobreza_pct),
      mediana = median(pobreza_pct),
      media   = mean(pobreza_pct),
      maximo  = max(pobreza_pct),
      .groups = "drop"
    ),
  n = Inf
)

cat("\nCinco municipios con mayor pobreza estimada:\n")
print(
  est %>%
    arrange(desc(pobreza_pct)) %>%
    select(Municipio, departamento, pobreza_pct, cv_pct) %>%
    head(5),
  row.names = FALSE
)

cat("\nCinco municipios con menor pobreza estimada:\n")
print(
  est %>%
    arrange(pobreza_pct) %>%
    select(Municipio, departamento, pobreza_pct, cv_pct) %>%
    head(5),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("COEFICIENTE DE VARIACIÓN (%)\n")
cat("============================================================\n")
print(summary(est$cv_pct))
cat("Percentil 90:", round(quantile(est$cv_pct, 0.90), 2), "\n")

cat("\nCinco municipios con mayor CV:\n")
print(
  est %>%
    arrange(desc(cv_pct)) %>%
    select(Municipio, departamento, pobreza_pct, cv_pct) %>%
    head(5),
  row.names = FALSE
)

cat("\nCinco municipios con menor CV:\n")
print(
  est %>%
    arrange(cv_pct) %>%
    select(Municipio, departamento, pobreza_pct, cv_pct) %>%
    head(5),
  row.names = FALSE
)

write.csv(
  est %>%
    select(cod_mun, Municipio, departamento,
           directo_pct, pobreza_pct, cv_pct),
  file.path(ruta_out, "final_estimaciones_m1_log.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ==============================================================================
# 5. CARTOGRAFÍA
# ==============================================================================

poligonos <- st_read(
  file.path(ruta_car, "MGN_ADM_MPIO_GRAFICO.shp"),
  quiet = TRUE
) %>%
  filter(dpto_ccdgo %in% c("19", "76")) %>%
  mutate(cod_mun = as.character(mpio_cdpmp)) %>%
  select(cod_mun, geometry) %>%
  st_make_valid() %>%
  st_transform(CRS_PROYECTADO)

faltan <- setdiff(est$cod_mun, poligonos$cod_mun)

stopifnot(
  "Faltan municipios del estudio en la cartografía." = length(faltan) == 0
)

mapa <- poligonos %>%
  inner_join(est, by = "cod_mun")

stopifnot("El mapa no tiene 84 municipios." = nrow(mapa) == 84)

deptos <- mapa %>%
  group_by(departamento) %>%
  summarise(.groups = "drop")


# ==============================================================================
# 6. MAPAS
# ==============================================================================

# Clases con cortes redondos; una escala secuencial de un solo tono por mapa.
br_pob <- pretty(mapa$pobreza_pct, n = 6)
br_cv  <- pretty(mapa$cv_pct, n = 6)

guia <- function() {
  guide_colorsteps(
    title.position = "top",
    show.limits    = TRUE,
    barwidth       = unit(5.5, "cm"),
    barheight      = unit(0.35, "cm")
  )
}

# Capas comunes: límite departamental y capitales con su nombre.
# El nombre lleva un halo blanco para leerse sobre cualquier color; si no está
# instalado shadowtext, sale sin halo. Cali se rotula hacia el oeste para no
# tapar los municipios pequeños al oriente de la ciudad.
cap_xy <- data.frame(
  nombre = capitales$nombre,
  st_coordinates(capitales),
  hjust  = c(1, 0),
  desp   = c(-15000, 15000)
) %>%
  mutate(x_lab = X + desp)

capa_nombres <- if (requireNamespace("shadowtext", quietly = TRUE)) {
  shadowtext::geom_shadowtext(
    data = cap_xy,
    aes(x = x_lab, y = Y, label = nombre, hjust = hjust),
    colour = "black", bg.colour = "white", bg.r = 0.15,
    size = 3.4, fontface = "bold", inherit.aes = FALSE
  )
} else {
  message("shadowtext no está instalado: los nombres salen sin halo.")
  geom_text(
    data = cap_xy,
    aes(x = x_lab, y = Y, label = nombre, hjust = hjust),
    colour = "black", size = 3.4, fontface = "bold", inherit.aes = FALSE
  )
}

capas_base <- list(
  geom_sf(data = deptos, fill = NA, color = "grey25", linewidth = 0.5),
  geom_sf(data = capitales, shape = 21, fill = "white", color = "black",
          size = 2.2, stroke = 0.6),
  capa_nombres,
  coord_sf(crs = CRS_PROYECTADO, datum = NA, expand = FALSE)
)

# Rampas secuenciales de un solo tono. Se omiten los dos pasos más claros de
# Brewer para que la primera clase no se confunda con el fondo blanco.
rampa <- function(paleta) scales::brewer_pal(palette = paleta)(9)[3:9]

# Sin título dentro de la imagen: el título va en el caption de LaTeX.
p_pob <- ggplot(mapa) +
  geom_sf(aes(fill = pobreza_pct), color = "white", linewidth = 0.15) +
  capas_base +
  scale_fill_stepsn(
    colours = rampa("Blues"),
    breaks  = br_pob,
    limits  = range(br_pob),
    labels  = coma,
    name    = "Pobreza monetaria (%)",
    guide   = guia()
  ) +
  tema_mapa

p_cv <- ggplot(mapa) +
  geom_sf(aes(fill = cv_pct), color = "white", linewidth = 0.15) +
  capas_base +
  scale_fill_stepsn(
    colours = rampa("Oranges"),
    breaks  = br_cv,
    limits  = range(br_cv),
    labels  = coma,
    name    = "Coeficiente de variación (%)",
    guide   = guia()
  ) +
  tema_mapa

# Escala y norte solo en el panel (a): ambos mapas tienen la misma extensión.
if (tiene_ggspatial) {
  p_pob <- p_pob +
    ggspatial::annotation_scale(
      location   = "bl",
      width_hint = 0.3,
      text_cex   = 0.7,
      pad_x      = unit(0.1, "cm"),
      pad_y      = unit(0.1, "cm")
    ) +
    # Esquina superior izquierda: queda libre porque el norte del Valle es
    # angosto y está hacia el oriente.
    ggspatial::annotation_north_arrow(
      location    = "tl",
      which_north = "true",
      height      = unit(0.9, "cm"),
      width       = unit(0.7, "cm"),
      pad_x       = unit(0.1, "cm"),
      pad_y       = unit(0.1, "cm"),
      style       = ggspatial::north_arrow_orienteering(text_size = 7)
    )
}

figura <- p_pob + p_cv +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

# Tamaño cercano al ancho de texto de la tesis: al insertarla con
# width=\textwidth las letras casi no se reducen.
ggsave(
  file.path(ruta_out, "fig_mapas_m1_log.png"),
  figura,
  width  = 7.5,
  height = 7.8,
  dpi    = 400,
  bg     = "white"
)

cat("\nArchivos generados en", ruta_out, ":\n")
cat("  final_coeficientes_m1_log.csv\n")
cat("  final_estimaciones_m1_log.csv\n")
cat("  fig_mapas_m1_log.png\n")
cat("\nProceso terminado.\n")