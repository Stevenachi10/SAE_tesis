# ==============================================================================
# 02_4g_figura_coeficientes.R
# ------------------------------------------------------------------------------
# Figura de la estabilidad de los coeficientes del modelo robusto respecto del
# modelo convencional (EBLUP). Lee output/rob_coeficientes.csv, generado por
# 02_4_robusto_diagnostico.R, y no reestima nada.
#
# Se grafica, para cada covariable, la diferencia (robusto - convencional) en
# unidades del error estándar del ajuste convencional. Es una medida descriptiva
# de magnitud relativa: no es el error estándar de la diferencia ni sustenta un
# contraste de significancia (ver 02_4_robusto_diagnostico.R).
#
# SALIDAS (output/)
#   fig_rob_coeficientes.png
#   rob_coeficientes_resumen.csv   Resumen por especificación y método.
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

ruta_out <- here("output")
cat("02_4g_figura_coeficientes.R: version v4, nombres académicos de las covariables\n")


# Nombres académicos de las covariables, como en la tabla de coeficientes.
ETIQUETAS <- c(
  alcantarillado_2018       = "Cobertura de alcantarillado (2018)",
  pct_rural                 = "Población rural (%, 2024)",
  partos_calificados_2020   = "Partos con atención calificada (2020)",
  transferencias_pc         = "Transferencias corrientes per cápita (2024)",
  saber11_lectura_2022      = "Saber 11, lectura crítica (2022)",
  tasa_violencia_intra_2019 = "Tasa de violencia intrafamiliar (2019)",
  pct_secundario_2023       = "Participación del sector secundario (2023)",
  tiene_empresas_formales   = "Presencia de empresas formales",
  dist_cali_km              = "Distancia a Cali (km)"
)
etiqueta <- function(x) ifelse(x %in% names(ETIQUETAS), ETIQUETAS[x], x)

COL_PUNTO <- "#b5482e"; COL_LINEA <- "#2c6e6b"; COL_BANDA <- "#c9c2b3"

tema_tesis <- theme_minimal(base_size = 13) +
  theme(plot.title = element_blank(), panel.grid.minor = element_blank())

cf <- read.csv(file.path(ruta_out, "rob_coeficientes.csv"),
               stringsAsFactors = FALSE, fileEncoding = "UTF-8")

stopifnot(all(c("especificacion", "metodo", "covariable",
                "desplazamiento_ee_base", "cambia_signo") %in% names(cf)))

# El intercepto se omite de la figura (se conserva en el resumen y en la tabla).
cf_fig <- cf %>% filter(covariable != "(Intercept)")
sin_nombre <- setdiff(unique(cf_fig$covariable), names(ETIQUETAS))
if (length(sin_nombre)) warning("Covariables sin nombre académico: ", paste(sin_nombre, collapse = ", "))

# reblup y reblupbc comparten el ajuste de los coeficientes. Si los coeficientes
# son idénticos, se grafica una sola serie para no duplicar puntos.
ancho <- cf_fig %>%
  select(especificacion, covariable, metodo, beta_robusto) %>%
  tidyr::pivot_wider(names_from = metodo, values_from = beta_robusto)
metodos_cf <- setdiff(names(ancho), c("especificacion", "covariable"))
idénticos <- length(metodos_cf) < 2 ||
  all(abs(ancho[[metodos_cf[1]]] - ancho[[metodos_cf[2]]]) < 1e-8)

if (idénticos) {
  cf_fig <- cf_fig %>% filter(metodo == metodos_cf[1]) %>%
    mutate(metodo = "Robusto (REBLUP y REBLUPbc)")
  colores <- c("Robusto (REBLUP y REBLUPbc)" = COL_LINEA)
} else {
  colores <- c(reblup = COL_LINEA, reblupbc = COL_PUNTO)
  cf_fig$metodo <- factor(cf_fig$metodo, levels = c("reblup", "reblupbc"))
}

# Orden de las covariables dentro de cada especificación: mayor desplazamiento
# absoluto arriba.
cf_fig <- cf_fig %>%
  group_by(especificacion, covariable) %>%
  mutate(orden = max(abs(desplazamiento_ee_base), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(covariable_f = paste(especificacion, covariable, sep = "___"),
         covariable_f = reorder(covariable_f, orden))

# Límite del eje: deja ver la banda de ±1 error estándar.
lim <- max(1.5, ceiling(max(abs(cf_fig$desplazamiento_ee_base), na.rm = TRUE) * 2) / 2)

# La leyenda de la forma se muestra siempre, aunque ningún coeficiente cambie de
# signo, para que el lector sepa interpretar la marca.
cf_fig <- cf_fig %>%
  mutate(cambia_signo = factor(as.logical(cambia_signo), levels = c(FALSE, TRUE)))

g <- ggplot(cf_fig, aes(x = desplazamiento_ee_base, y = covariable_f)) +
  annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf,
           fill = COL_BANDA, alpha = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", color = COL_BANDA) +
  geom_point(aes(color = metodo, shape = cambia_signo),
             size = 2.4, alpha = 0.9,
             position = position_dodge(width = 0.5),
             show.legend = TRUE) +
  scale_color_manual(values = colores, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                     labels = c(`FALSE` = "Mismo signo",
                                `TRUE`  = "Signo diferente"),
                     limits = c("FALSE", "TRUE"),
                     name = NULL, drop = FALSE) +
  guides(shape = guide_legend(order = 2,
                              override.aes = list(colour = "black", alpha = 1,
                                                  size = 2.8, stroke = 0.9)),
         colour = guide_legend(order = 1)) +
  scale_y_discrete(labels = function(x) etiqueta(sub("^.*___", "", x))) +
  coord_cartesian(xlim = c(-lim, lim)) +
  facet_grid(especificacion ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Diferencia (robusto menos convencional)\nen errores estándar del ajuste convencional",
       y = NULL) +
  tema_tesis +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 11))

ggsave(file.path(ruta_out, "fig_rob_coeficientes.png"), g,
       width = 9, height = 7, dpi = 400)

# ---- Resumen numérico --------------------------------------------------------
resumen <- cf %>%
  group_by(especificacion, metodo) %>%
  summarise(n_coef = n(),
            max_abs_ee = round(max(abs(desplazamiento_ee_base), na.rm = TRUE), 3),
            media_abs_ee = round(mean(abs(desplazamiento_ee_base), na.rm = TRUE), 3),
            cambios_signo = sum(cambia_signo),
            .groups = "drop")

print(as.data.frame(resumen), row.names = FALSE)
write.csv(resumen, file.path(ruta_out, "rob_coeficientes_resumen.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")