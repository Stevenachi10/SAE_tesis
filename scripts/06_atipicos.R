library(dplyr)
library(ggplot2)
library(here)

ruta_out <- here("output")
matriz <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

COL_PUNTO <- "#b5482e"; COL_BANDA <- "#c9c2b3"
COL_MIN   <- "#2c6e6b"; COL_OSCURO <- "#1a1814"

# ============================================================
# ATÍPICOS MULTIDIMENSIONALES - Distancia de Mahalanobis
# Solo variables CONTINUAS (NO dummies, NO ordinales, NO la pobreza)
# ============================================================
vars_continuas <- c("extension", "saber11_lectura_2022", "partos_calificados_2020",
                    "dist_popayan_km", "dist_cali_km", "desplaz_tasa_1000",
                    "pct_contributivo")

X <- matriz %>% select(all_of(vars_continuas)) %>% as.matrix()

# Distancia de Mahalanobis (invariante a escala)
centro   <- colMeans(X)
cov_mat  <- cov(X)
d_mahal  <- mahalanobis(X, center = centro, cov = cov_mat)

# Umbral chi-cuadrado (p = nº de variables continuas)
p <- length(vars_continuas)
umbral_95 <- qchisq(0.95, df = p)
umbral_99 <- qchisq(0.99, df = p)

# Tabla de resultados (con nombres únicos para duplicados)
depto <- ifelse(substr(matriz$cod_mun, 1, 2) == "19", "Cauca", "Valle")
nombre_mun <- matriz$Municipio
# Diferenciar solo los duplicados (Argelia, Bolívar)
dups <- nombre_mun[duplicated(nombre_mun) | duplicated(nombre_mun, fromLast = TRUE)]
nombre_final <- ifelse(nombre_mun %in% dups,
                       paste0(nombre_mun, " (", depto, ")"),
                       nombre_mun)

resultado <- data.frame(
  municipio = nombre_final,
  d_mahalanobis = round(d_mahal, 2),
  atipico_95 = d_mahal > umbral_95,
  atipico_99 = d_mahal > umbral_99
) %>% arrange(desc(d_mahalanobis))

cat("Umbral chi2 95%:", round(umbral_95, 2), "| 99%:", round(umbral_99, 2), "\n")
cat("Atípicos al 95%:", sum(resultado$atipico_95),
    "| al 99%:", sum(resultado$atipico_99), "\n\n")
print(head(resultado, 12), row.names = FALSE)

# ============================================================
# GRÁFICO: distancia de Mahalanobis por municipio
# ============================================================
df_plot <- resultado %>%
  mutate(municipio = factor(municipio, levels = rev(municipio)),
         categoria = case_when(
           atipico_99 ~ "Atípico (99%)",
           atipico_95 ~ "Atípico (95%)",
           TRUE ~ "No atípico"))

# Solo etiquetar/mostrar los más altos para que sea legible
df_top <- df_plot %>% arrange(desc(d_mahalanobis)) %>% head(20)

g_mahal <- ggplot(df_top, aes(x = d_mahalanobis, y = municipio, fill = categoria)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = umbral_95, linetype = "dashed", color = COL_MIN, linewidth = 0.6) +
  geom_vline(xintercept = umbral_99, linetype = "dashed", color = COL_OSCURO, linewidth = 0.6) +
  annotate("text", x = umbral_95, y = 1, label = "95%", hjust = -0.2, size = 3, color = COL_MIN) +
  annotate("text", x = umbral_99, y = 1, label = "99%", hjust = -0.2, size = 3, color = COL_OSCURO) +
  scale_fill_manual(values = c("Atípico (99%)" = COL_PUNTO,
                               "Atípico (95%)" = "#d8956f",
                               "No atípico" = COL_BANDA), name = NULL) +
  labs(x = "Distancia de Mahalanobis", y = NULL,
       title = "Atípicos multidimensionales (distancia de Mahalanobis)",
       caption = "Fuente: elaboración propia. Umbrales según distribución chi-cuadrado.") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, color = COL_OSCURO),
        plot.caption = element_text(size = 8, color = "#6b6b6b", hjust = 0, margin = margin(t = 8)),
        axis.text.y = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        plot.title.position = "plot", plot.caption.position = "plot")

print(g_mahal)
ggsave(file.path(ruta_out, "fig_mahalanobis.png"), g_mahal,
       width = 8, height = 7, dpi = 300, bg = "white")
############################

# Ver los valores de los atípicos en cada variable continua
atipicos_nombres <- c("Argelia", "Buenaventura", "Jambaló", "López de Micay",
                      "Caldono", "Páez", "Suárez")

matriz %>%
  filter(Municipio %in% atipicos_nombres & substr(cod_mun,1,2) %in% c("19","76")) %>%
  mutate(depto = ifelse(substr(cod_mun,1,2)=="19","Cauca","Valle")) %>%
  select(Municipio, depto, extension, desplaz_tasa_1000, dist_popayan_km,
         dist_cali_km, pct_contributivo, partos_calificados_2020, saber11_lectura_2022) %>%
  arrange(desc(extension)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Comparar con las medias de la muestra
cat("\n=== MEDIAS DE LA MUESTRA (referencia) ===\n")
matriz %>%
  summarise(extension = mean(extension),
            desplaz = mean(desplaz_tasa_1000),
            dist_pop = mean(dist_popayan_km),
            dist_cali = mean(dist_cali_km),
            contrib = mean(pct_contributivo),
            partos = mean(partos_calificados_2020),
            saber11 = mean(saber11_lectura_2022)) %>%
  round(2) %>% as.data.frame() %>% print(row.names = FALSE)
