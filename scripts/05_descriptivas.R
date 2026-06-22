# ============================================================
# 05 - ANÁLISIS DESCRIPTIVO
# Tesis SAE - Pobreza Monetaria Municipal (Cauca y Valle)
# ============================================================
# Entrada: matriz_sae_transformada_v2.rds, cov_descripciongeneral.rds,
#          municipios_geo.geojson
# Salida:  figuras .png en output/
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(sf)
library(patchwork)
library(here)

# Paleta del proyecto
COL_PUNTO  <- "#b5482e"   # teja
COL_BANDA  <- "#c9c2b3"   # crema
COL_MIN    <- "#2c6e6b"   # verde azulado
COL_OSCURO <- "#1a1814"   # casi negro

ruta_out <- here("output")
fix_cod  <- function(x) str_pad(as.character(as.integer(x)), 5, pad = "0")

matriz <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

# Recuperar cat_ruralidad desde el cov original
descripcion <- readRDS(here("data", "pivoteadas", "cov_descripciongeneral.rds"))
rural_cat <- descripcion %>%
  mutate(cod_mun = fix_cod(`Código Entidad`)) %>%
  select(cod_mun, cat_ruralidad = `Categoría de ruralidad_2000`)
matriz <- matriz %>% left_join(rural_cat, by = "cod_mun")

# ============================================================
# 1. HISTOGRAMAS DE VARIABLES CONTINUAS  [ANEXO]
# ============================================================
vars_cont <- c("extension", "saber11_lectura_2022", "partos_calificados_2020",
               "dist_popayan_km", "dist_cali_km", "desplaz_tasa_1000",
               "pct_contributivo")

etiquetas <- c(
  extension = "Extensión territorial (km²)",
  saber11_lectura_2022 = "Puntaje Saber 11 lectura crítica",
  partos_calificados_2020 = "Partos con atención calificada (%)",
  dist_popayan_km = "Distancia a Popayán (km)",
  dist_cali_km = "Distancia a Cali (km)",
  desplaz_tasa_1000 = "Tasa de desplazamiento (por 1.000 hab.)",
  pct_contributivo = "Afiliación al régimen contributivo (%)"
)

df_hist <- matriz %>%
  select(all_of(vars_cont)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "valor") %>%
  mutate(variable = factor(etiquetas[variable], levels = etiquetas))

g_hist <- ggplot(df_hist, aes(x = valor)) +
  geom_histogram(bins = 14, fill = COL_PUNTO, color = "white", linewidth = 0.3) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(x = NULL, y = "Número de municipios",
       title = "Distribución de las variables auxiliares continuas",
       caption = "Fuente: elaboración propia con datos del DANE, ICFES,\nMinSalud, UARIV y ADRES.") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, color = COL_OSCURO),
        plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0, margin = margin(t = 10)),
        strip.text = element_text(face = "bold", size = 8.5, color = COL_OSCURO),
        strip.background = element_rect(fill = COL_BANDA, color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 8, color = "#5b5b5b"),
        plot.title.position = "plot", plot.caption.position = "plot")

ggsave(file.path(ruta_out, "fig_histogramas.png"), g_hist,
       width = 10, height = 6.5, dpi = 300, bg = "white")

# ============================================================
# 2. BARRAS DE CATEGÓRICAS  [ANEXO]
# ============================================================

# --- 2a. Categoría de ruralidad ---
df_rural <- matriz %>%
  filter(!is.na(cat_ruralidad)) %>%
  count(cat_ruralidad) %>%
  mutate(cat_ruralidad = factor(cat_ruralidad,
                                levels = c("Ciudades y aglomeraciones", "Intermedio", "Rural", "Rural disperso")))

g_rural <- ggplot(df_rural, aes(x = cat_ruralidad, y = n)) +
  geom_col(fill = COL_PUNTO, color = "white", linewidth = 0.3, width = 0.68) +
  geom_text(aes(label = n), vjust = -0.4, size = 4, color = COL_OSCURO, fontface = "bold") +
  labs(x = NULL, y = "Número de municipios",
       title = "Distribución de municipios por categoría de ruralidad",
       caption = "Fuente: elaboración propia con datos del DNP (TerriData).") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, color = COL_OSCURO),
        plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0, margin = margin(t = 10)),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 12, hjust = 1, size = 10),
        plot.title.position = "plot", plot.caption.position = "plot")

ggsave(file.path(ruta_out, "fig_ruralidad.png"), g_rural, width = 8, height = 5, dpi = 300, bg = "white")

# --- 2b. Categoría Ley 617 ---
df_cat617 <- matriz %>% count(categoria_617) %>% mutate(categoria_617 = factor(categoria_617))

g_cat617 <- ggplot(df_cat617, aes(x = categoria_617, y = n)) +
  geom_col(fill = COL_MIN, color = "white", linewidth = 0.3, width = 0.68) +
  geom_text(aes(label = n), vjust = -0.4, size = 4, color = COL_OSCURO, fontface = "bold") +
  labs(x = "Categoría Ley 617", y = "Número de municipios",
       title = "Distribución de municipios por categoría de la Ley 617",
       caption = "Fuente: elaboración propia con datos del DNP (TerriData).") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, color = COL_OSCURO),
        plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0, margin = margin(t = 10)),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title.position = "plot", plot.caption.position = "plot")

ggsave(file.path(ruta_out, "fig_cat617.png"), g_cat617, width = 8, height = 5, dpi = 300, bg = "white")

# ============================================================
# 3. MAPAS COROPLÉTICOS  [CUERPO]
# ============================================================

# --- Cargar geometría y unir por nombre + departamento ---
geo <- st_read(here("data", "pivoteadas", "municipios_geo.geojson"))

norm_nombre <- function(x) trimws(tolower(iconv(x, "UTF-8", "ASCII//TRANSLIT")))
corregir_nombre <- function(x) {
  case_when(
    x == "piendamo"         ~ "piendamo - tunia",
    x == "sotara"           ~ "sotara paispamba",
    x == "santiago de cali" ~ "cali",
    TRUE ~ x
  )
}

geo <- geo %>%
  mutate(municipio_norm = corregir_nombre(norm_nombre(municipio)),
         depto = departamento)

matriz_geo <- matriz %>%
  mutate(municipio_norm = corregir_nombre(norm_nombre(Municipio)),
         depto = ifelse(substr(cod_mun, 1, 2) == "19", "Cauca", "Valle del Cauca"))

mapa <- geo %>% left_join(matriz_geo, by = c("municipio_norm", "depto"))

# Verificar el join
cat("Municipios en geojson:", nrow(geo),
    "| unidos:", sum(!is.na(mapa$pobreza_monetaria)),
    "| sin unir:", sum(is.na(mapa$pobreza_monetaria)), "\n")

# --- Mapa 1: Pobreza monetaria ---
g_mapa_pobreza <- ggplot(mapa) +
  geom_sf(aes(fill = pobreza_monetaria), color = "white", linewidth = 0.12) +
  scale_fill_gradient(low = COL_BANDA, high = COL_PUNTO, name = "Incidencia",
                      na.value = "grey85", labels = scales::percent) +
  labs(title = "Incidencia de pobreza monetaria",
       caption = "Fuente: elaboración propia con datos del DANE.") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = COL_OSCURO),
        plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0.5),
        legend.position = "bottom")
ggsave(file.path(ruta_out, "fig_mapa_pobreza.png"), g_mapa_pobreza, width = 7, height = 8, dpi = 300, bg = "white")

# --- Mapa 2: Presencia de coca (dummy) ---
mapa <- mapa %>% mutate(coca_lab = case_when(
  tiene_coca == 1 ~ "Con coca", tiene_coca == 0 ~ "Sin coca", TRUE ~ NA_character_))

g_mapa_coca <- ggplot(mapa) +
  geom_sf(aes(fill = coca_lab), color = "white", linewidth = 0.12) +
  scale_fill_manual(values = c("Sin coca" = COL_BANDA, "Con coca" = COL_PUNTO),
                    name = NULL, na.value = "grey85") +
  labs(title = "Presencia de cultivos de coca",
       caption = "Fuente: elaboración propia con datos de UNODC-SIMCI.") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = COL_OSCURO),
        plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0.5),
        legend.position = "bottom")
ggsave(file.path(ruta_out, "fig_mapa_coca.png"), g_mapa_coca, width = 7, height = 8, dpi = 300, bg = "white")

# --- Mapa 3: Tasa de desplazamiento (escala raíz para el sesgo) ---
g_mapa_desplaz <- ggplot(mapa) +
  geom_sf(aes(fill = desplaz_tasa_1000), color = "white", linewidth = 0.12) +
  scale_fill_gradient(low = COL_BANDA, high = COL_PUNTO, name = "Por 1.000 hab.",
                      na.value = "grey85", trans = "sqrt") +
  labs(title = "Tasa de desplazamiento forzado",
       caption = "Fuente: elaboración propia con datos de la UARIV.") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = COL_OSCURO),
        plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0.5),
        legend.position = "bottom")
ggsave(file.path(ruta_out, "fig_mapa_desplaz.png"), g_mapa_desplaz, width = 7, height = 8, dpi = 300, bg = "white")

cat("\nGráficos guardados en output/\n")
###############################


library(patchwork)

# (Reutiliza el objeto 'mapa' del join que ya hiciste arriba)

# Etiqueta de coca
mapa <- mapa %>% mutate(coca_lab = case_when(
  tiene_coca == 1 ~ "Con coca", tiene_coca == 0 ~ "Sin coca", TRUE ~ NA_character_))

# Tema común para los 3 (compacto, para que quepan)
tema_mapa <- theme_void(base_size = 10) +
  theme(plot.subtitle = element_text(face = "bold", size = 10, hjust = 0.5, color = COL_OSCURO),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))

# --- Mapa 1: Pobreza ---
m_pobreza <- ggplot(mapa) +
  geom_sf(aes(fill = pobreza_monetaria), color = "white", linewidth = 0.1) +
  scale_fill_gradient(low = COL_BANDA, high = COL_PUNTO, name = "Incidencia",
                      na.value = "grey85", labels = scales::percent) +
  labs(subtitle = "Pobreza monetaria") + tema_mapa

# --- Mapa 2: Coca ---
m_coca <- ggplot(mapa) +
  geom_sf(aes(fill = coca_lab), color = "white", linewidth = 0.1) +
  scale_fill_manual(values = c("Sin coca" = COL_BANDA, "Con coca" = COL_PUNTO),
                    name = NULL, na.value = "grey85") +
  labs(subtitle = "Presencia de coca") + tema_mapa

# --- Mapa 3: Desplazamiento ---
m_desplaz <- ggplot(mapa) +
  geom_sf(aes(fill = desplaz_tasa_1000), color = "white", linewidth = 0.1) +
  scale_fill_gradient(low = COL_BANDA, high = COL_PUNTO, name = "Por 1.000 hab.",
                      na.value = "grey85", trans = "sqrt") +
  labs(subtitle = "Tasa de desplazamiento") + tema_mapa

# --- Combinar los 3 lado a lado ---
mapa_panel <- (m_pobreza | m_coca | m_desplaz) +
  plot_annotation(
    title = "Distribución territorial de la pobreza, la coca y el desplazamiento forzado",
    caption = "Fuente: elaboración propia con datos del DANE, UNODC-SIMCI y la UARIV.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = COL_OSCURO),
      plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0.5)
    )
  )

print(mapa_panel)
ggsave(file.path(ruta_out, "fig_mapa_panel.png"), mapa_panel,
       width = 13, height = 6, dpi = 300, bg = "white")
##############################3

# --- Combinar: 2 arriba (coca, desplazamiento) + 1 abajo (pobreza) ---
mapa_panel <- (m_coca | m_desplaz) / m_pobreza +
  plot_annotation(
    title = "Distribución territorial de la pobreza, la coca y el desplazamiento forzado",
    caption = "Fuente: elaboración propia con datos del DANE, UNODC-SIMCI y la UARIV.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = COL_OSCURO),
      plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0.5)
    )
  )

print(mapa_panel)
ggsave(file.path(ruta_out, "fig_mapa_panel.png"), mapa_panel,
       width = 10, height = 11, dpi = 300, bg = "white")
#############################################

# ============================================================
# 1. Municipios CON coca
# ============================================================
cat("=== MUNICIPIOS CON PRESENCIA DE COCA ===\n")
con_coca <- matriz %>%
  filter(tiene_coca == 1) %>%
  mutate(depto = ifelse(substr(cod_mun,1,2)=="19","Cauca","Valle")) %>%
  select(Municipio, depto, pobreza_monetaria) %>%
  arrange(depto, desc(pobreza_monetaria))

cat("Total con coca:", nrow(con_coca), "de", nrow(matriz), "\n\n")
print(as.data.frame(con_coca), row.names = FALSE)

# Cuántos por departamento
cat("\nPor departamento:\n")
print(table(con_coca$depto))

# ============================================================
# 2. Municipios con MAYOR tasa de desplazamiento
# ============================================================
cat("\n\n=== TOP 10 MAYOR TASA DE DESPLAZAMIENTO ===\n")
top_desplaz <- matriz %>%
  mutate(depto = ifelse(substr(cod_mun,1,2)=="19","Cauca","Valle")) %>%
  select(Municipio, depto, desplaz_tasa_1000, pobreza_monetaria, tiene_coca) %>%
  arrange(desc(desplaz_tasa_1000)) %>%
  head(10)

print(as.data.frame(top_desplaz), row.names = FALSE)

cat("\n\n=== TOP 10 MAYOR POBREZA ===\n")
top_pobreza <- matriz %>%
  mutate(depto = ifelse(substr(cod_mun,1,2)=="19","Cauca","Valle")) %>%
  select(Municipio, depto, pobreza_monetaria, tiene_coca, desplaz_tasa_1000) %>%
  arrange(desc(pobreza_monetaria)) %>%
  head(10)

print(as.data.frame(top_pobreza), row.names = FALSE)

cat("\n\n=== MUNICIPIOS CRÍTICOS (coca + desplaz alto + pobreza alta) ===\n")
criticos <- matriz %>%
  mutate(depto = ifelse(substr(cod_mun,1,2)=="19","Cauca","Valle")) %>%
  filter(tiene_coca == 1 &
           desplaz_tasa_1000 > median(matriz$desplaz_tasa_1000) &
           pobreza_monetaria > median(matriz$pobreza_monetaria)) %>%
  select(Municipio, depto, pobreza_monetaria, desplaz_tasa_1000) %>%
  arrange(desc(pobreza_monetaria))

print(as.data.frame(criticos), row.names = FALSE)


###########################

library(reshape2)

ruta_out <- here("output")
matriz <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

COL_PUNTO <- "#b5482e"; COL_BANDA <- "#c9c2b3"
COL_MIN   <- "#2c6e6b"; COL_OSCURO <- "#1a1814"

# Variables del modelo (continuas + categoria_617 ordinal)
# Se excluyen las dummies puras (tiene_coca, rural_disperso) porque la
# correlación de Spearman con binarias es poco informativa.
vars_corr <- c("pobreza_monetaria",
               "extension", "saber11_lectura_2022", "partos_calificados_2020",
               "dist_popayan_km", "dist_cali_km", "desplaz_tasa_1000",
               "pct_contributivo", "categoria_617")

# Etiquetas cortas para que quepan
etiq <- c(
  pobreza_monetaria = "Pobreza",
  extension = "Extensión",
  saber11_lectura_2022 = "Saber 11",
  partos_calificados_2020 = "Partos calif.",
  dist_popayan_km = "Dist. Popayán",
  dist_cali_km = "Dist. Cali",
  desplaz_tasa_1000 = "Desplazamiento",
  pct_contributivo = "% Contributivo",
  categoria_617 = "Categoría 617"
)

# Matriz de correlación de Spearman
M <- matriz %>% select(all_of(vars_corr))
cor_sp <- cor(M, method = "spearman", use = "complete.obs")

# Renombrar con etiquetas
rownames(cor_sp) <- etiq[rownames(cor_sp)]
colnames(cor_sp) <- etiq[colnames(cor_sp)]

# Quedarnos con el triángulo inferior (evita redundancia)
cor_sp[upper.tri(cor_sp)] <- NA

# Derretir a formato largo
df_cor <- melt(cor_sp, na.rm = TRUE)

# Heatmap
g_cor <- ggplot(df_cor, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3, color = COL_OSCURO) +
  scale_fill_gradient2(low = COL_MIN, mid = "white", high = COL_PUNTO,
                       midpoint = 0, limit = c(-1, 1), name = "Spearman") +
  labs(title = "Matriz de correlación de Spearman",
       caption = "Fuente: elaboración propia.") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, color = COL_OSCURO),
        plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0, margin = margin(t = 8)),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.title.position = "plot", plot.caption.position = "plot")

print(g_cor)
ggsave(file.path(ruta_out, "fig_correlacion_spearman.png"), g_cor,
       width = 8, height = 7, dpi = 300, bg = "white")
