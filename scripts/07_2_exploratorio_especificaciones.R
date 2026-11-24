# ==============================================================================
# 07_2  EDA FOCALIZADO EN M1 Y SIGNO DE LA RURALIDAD
# ==============================================================================
# Responde a los puntos 17 (analisis descriptivo focalizado en las variables
# finales) y 18 (signo negativo de ruralidad) de la revision.
#
# ------------------------------------------------------------------------------
# PARTE A. EDA DE M1 (punto 17)
# ------------------------------------------------------------------------------
#   - distribucion de la respuesta operacional y de las seis covariables;
#   - relaciones bivariadas con la respuesta, con ajuste lineal y loess para
#     revisar la forma funcional;
#   - asociaciones entre predictores (Pearson y Spearman).
#
# ------------------------------------------------------------------------------
# PARTE B. SIGNO DE pct_rural (punto 18)
# ------------------------------------------------------------------------------
# La revision pide explicar por que el coeficiente de pct_rural es negativo en
# M1-log cuando, descriptivamente, los municipios mas rurales son mas pobres.
# La asociacion condicional puede diferir de la marginal; el script documenta
# esa diferencia con cuatro piezas:
#
#   1. Asociacion marginal: correlacion y pendiente simple con la respuesta.
#   2. Coeficiente ajustado: pendiente de pct_rural dentro de M1.
#   3. Correlaciones de pct_rural con las otras cinco covariables.
#   4. Trayectoria del coeficiente: (a) al ir agregando las covariables de
#      M1 una a una y (b) al retirar cada una de las otras cinco. Identifica
#      que covariable hace cambiar el signo (supresion o confusion).
#
# Cada coeficiente se reporta de dos formas: MCO sobre log(pobreza) y modelo
# Fay-Herriot con transformacion logaritmica (emdi::fh, REML, bc_crude, igual
# que M1-log en 02_2a). MCO no usa las varianzas de muestreo y sirve para
# entender la estructura de asociacion; el coeficiente de FH es el del modelo.
# Si emdi no esta instalado, solo se calcula la parte MCO.
#
# Todo es descriptivo: no se interpreta el modelo como prueba causal de que la
# ruralidad aumente o reduzca la pobreza.
#
# ------------------------------------------------------------------------------
# ENTRADAS
# ------------------------------------------------------------------------------
#   output/matriz_sae_transformada_v2.rds
#   output/seleccion_stepwise_completo_both.rds
#
# ------------------------------------------------------------------------------
# SALIDAS (output/)
# ------------------------------------------------------------------------------
#   eda_m1_descriptivo.csv           Resumen de la respuesta y las seis covariables.
#   eda_m1_correlaciones_pearson.csv
#   eda_m1_correlaciones_spearman.csv
#   ruralidad_marginal_vs_ajustado.csv
#   ruralidad_trayectoria_agregando.csv
#   ruralidad_trayectoria_retirando.csv
#   figuras/fig_eda_m1_distribuciones.png
#   figuras/fig_eda_m1_bivariadas.png
#   figuras/fig_eda_m1_bivariadas_log.png
#   figuras/fig_eda_m1_correlaciones.png
#   figuras/fig_ruralidad_variable_agregada.png
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")
ruta_fig <- file.path(ruta_out, "figuras")
dir.create(ruta_fig, showWarnings = FALSE, recursive = TRUE)

SEMILLA <- 2906
HAY_EMDI <- requireNamespace("emdi", quietly = TRUE)

matriz  <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))
corrida <- readRDS(file.path(ruta_out, "seleccion_stepwise_completo_both.rds"))

# M1 es la especificacion elegida por KICb2.
fila_m1 <- corrida$tabla[corrida$tabla$criterio == "KICb2", ]
VARS_M1 <- strsplit(fila_m1$covariables, " + ", fixed = TRUE)[[1]]

if (length(VARS_M1) != 6 || !"pct_rural" %in% VARS_M1) {
  stop("M1 no tiene la composicion esperada (seis variables con pct_rural): ",
       paste(VARS_M1, collapse = ", "))
}

faltan <- setdiff(VARS_M1, names(matriz))
if (length(faltan) > 0) {
  stop("Covariables ausentes en la matriz: ", paste(faltan, collapse = ", "))
}

datos <- matriz %>%
  select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza,
         all_of(VARS_M1)) %>%
  mutate(
    depto = ifelse(substr(cod_mun, 1, 2) == "19", "Cauca", "Valle del Cauca"),
    log_pobreza = log(pobreza_monetaria)
  ) %>%
  as.data.frame()

cat("Dominios:", nrow(datos), "| Cauca:", sum(datos$depto == "Cauca"),
    "| Valle:", sum(datos$depto != "Cauca"), "\n")
cat("M1:", paste(VARS_M1, collapse = ", "), "\n")

COL_DEPTO <- c("Cauca" = "#B03A2E", "Valle del Cauca" = "#1F618D")

asimetria <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x)
  mean(((x - m) / s)^3)
}


# ==============================================================================
# PARTE A. EDA DE M1
# ==============================================================================

variables_eda <- c("pobreza_monetaria", VARS_M1)

# ---- A1. Resumen ------------------------------------------------------------
descriptivo <- do.call(rbind, lapply(variables_eda, function(v) {
  x <- datos[[v]]
  data.frame(
    variable = v,
    n        = sum(is.finite(x)),
    media    = mean(x), de = sd(x),
    minimo   = min(x), p25 = unname(quantile(x, 0.25)),
    mediana  = median(x), p75 = unname(quantile(x, 0.75)),
    maximo   = max(x), asimetria = asimetria(x),
    stringsAsFactors = FALSE)
}))
descriptivo[, -1] <- round(descriptivo[, -1], 3)

cat("\n== RESUMEN DE LA RESPUESTA Y LAS SEIS COVARIABLES ==\n")
print(descriptivo, row.names = FALSE)
write.csv(descriptivo, file.path(ruta_out, "eda_m1_descriptivo.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- A2. Distribuciones -----------------------------------------------------
largo <- datos %>%
  select(cod_mun, depto, all_of(variables_eda)) %>%
  pivot_longer(all_of(variables_eda), names_to = "variable",
               values_to = "valor") %>%
  mutate(variable = factor(variable, levels = variables_eda))

g_dist <- ggplot(largo, aes(valor, fill = depto)) +
  geom_histogram(bins = 20, colour = "white", linewidth = 0.2) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  scale_fill_manual(values = COL_DEPTO, name = NULL) +
  labs(x = NULL, y = "Municipios",
       title = "Distribucion de la respuesta operacional y de las covariables de M1") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(ruta_fig, "fig_eda_m1_distribuciones.png"), g_dist,
       width = 10, height = 8, dpi = 200)

# ---- A3. Relaciones bivariadas con la respuesta -----------------------------
# Nombres cortos de los paneles. El título y la explicación de las líneas van
# en el pie de la figura en la tesis, no dentro de la imagen.
PANELES <- c(
  alcantarillado_2018       = "Alcantarillado (%)",
  pct_rural                 = "Población rural (%)",
  partos_calificados_2020   = "Partos calificados (%)",
  transferencias_pc         = "Transferencias per cápita",
  saber11_lectura_2022      = "Saber 11, lectura crítica",
  tasa_violencia_intra_2019 = "Violencia intrafamiliar"
)
etiqueta_panel <- function(x) ifelse(x %in% names(PANELES), PANELES[x], x)
bivariada <- function(respuesta, etiqueta_y, archivo) {
  b <- datos %>%
    select(cod_mun, depto, y = all_of(respuesta), all_of(VARS_M1)) %>%
    pivot_longer(all_of(VARS_M1), names_to = "variable", values_to = "x") %>%
    mutate(variable = factor(variable, levels = VARS_M1))
  
  g <- ggplot(b, aes(x, y)) +
    geom_point(aes(colour = depto), alpha = 0.75, size = 1.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey25", linewidth = 0.6) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                colour = "darkorange", linewidth = 0.6, linetype = "dashed") +
    facet_wrap(~ variable, scales = "free_x", ncol = 3,
               labeller = as_labeller(etiqueta_panel)) +
    scale_colour_manual(values = COL_DEPTO, name = NULL) +
    labs(x = NULL, y = etiqueta_y) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(ruta_fig, archivo), g, width = 10, height = 7, dpi = 200)
  invisible(g)
}

bivariada("pobreza_monetaria", "Pobreza monetaria (referencia DANE)",
          "fig_eda_m1_bivariadas.png")
bivariada("log_pobreza", "log(pobreza monetaria)",
          "fig_eda_m1_bivariadas_log.png")

# ---- A4. Asociaciones entre predictores -------------------------------------
cor_p <- round(cor(datos[, VARS_M1], method = "pearson"), 3)
cor_s <- round(cor(datos[, VARS_M1], method = "spearman"), 3)

cat("\n== CORRELACION DE PEARSON ENTRE PREDICTORES ==\n")
print(cor_p)

write.csv(cor_p, file.path(ruta_out, "eda_m1_correlaciones_pearson.csv"),
          fileEncoding = "UTF-8")
write.csv(cor_s, file.path(ruta_out, "eda_m1_correlaciones_spearman.csv"),
          fileEncoding = "UTF-8")

cor_larga <- as.data.frame(as.table(cor_p)) %>%
  rename(x = Var1, y = Var2, r = Freq)

g_cor <- ggplot(cor_larga, aes(x, y, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "#1F618D", mid = "white", high = "#B03A2E",
                       limits = c(-1, 1), name = "Pearson") +
  labs(x = NULL, y = NULL, title = "Correlacion entre las covariables de M1") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(ruta_fig, "fig_eda_m1_correlaciones.png"), g_cor,
       width = 7.5, height = 6, dpi = 200)


# ==============================================================================
# PARTE B. SIGNO DE pct_rural
# ==============================================================================

# ---- Ajustes ----------------------------------------------------------------
ols_coef <- function(vars) {
  f <- as.formula(paste("log_pobreza ~", paste(vars, collapse = " + ")))
  s <- summary(lm(f, data = datos))$coefficients
  data.frame(b_mco = s["pct_rural", 1], p_mco = s["pct_rural", 4])
}

# Misma llamada que M1-log en 02_2a (REML, log, bc_crude, MSE analitico).
# Los coeficientes se leen de fit$model$coefficients, con la misma estructura
# que usa 06_0; si la version de emdi la entregara distinta, el script se
# detiene con un mensaje y se ajusta aqui.
fh_coef <- function(vars) {
  if (!HAY_EMDI) return(data.frame(b_fh = NA_real_, p_fh = NA_real_))
  f <- as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
  fit <- suppressWarnings(emdi::fh(
    f, vardir = "varianza_pobreza", combined_data = datos,
    domains = "cod_mun", method = "reml",
    transformation = "log", backtransformation = "bc_crude",
    MSE = TRUE, mse_type = "analytical", B = c(0, 0), seed = SEMILLA))
  cf <- as.data.frame(fit$model$coefficients)
  # Misma estructura que usa 06_0: coefficients, std.error, t.value, p.value.
  if (!all(c("coefficients", "p.value") %in% names(cf)) ||
      !"pct_rural" %in% rownames(cf)) {
    stop("Estructura inesperada de fit$model$coefficients. Revisar ",
         "str(fit$model$coefficients) y ajustar fh_coef().")
  }
  data.frame(b_fh = as.numeric(cf["pct_rural", "coefficients"]),
             p_fh = as.numeric(cf["pct_rural", "p.value"]))
}

coef_rural <- function(vars) cbind(ols_coef(vars), fh_coef(vars))

signo <- function(b) ifelse(is.na(b), NA_character_, ifelse(b < 0, "-", "+"))

# ---- B1 y B2. Marginal frente a ajustado ------------------------------------
marg <- coef_rural("pct_rural")
ajus <- coef_rural(VARS_M1)

r_marg  <- cor(datos$pct_rural, datos$pobreza_monetaria)
rs_marg <- cor(datos$pct_rural, datos$pobreza_monetaria, method = "spearman")
r_marg_log <- cor(datos$pct_rural, datos$log_pobreza)

tabla_marg <- rbind(
  cbind(especificacion = "Marginal: solo pct_rural", marg),
  cbind(especificacion = "Ajustada: M1 (seis covariables)", ajus))
tabla_marg$signo_mco <- signo(tabla_marg$b_mco)
tabla_marg$signo_fh  <- signo(tabla_marg$b_fh)

cat("\n== pct_rural: ASOCIACION MARGINAL FRENTE A AJUSTADA (respuesta en log) ==\n")
print(tabla_marg, row.names = FALSE, digits = 4)
cat(sprintf(paste0("Correlacion marginal pct_rural con pobreza: Pearson %.3f ",
                   "| Spearman %.3f | con log(pobreza): %.3f\n"),
            r_marg, rs_marg, r_marg_log))

write.csv(tabla_marg, file.path(ruta_out, "ruralidad_marginal_vs_ajustado.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- B3. Correlacion de pct_rural con las otras covariables -----------------
otras <- setdiff(VARS_M1, "pct_rural")
cat("\n== CORRELACION DE pct_rural CON LAS OTRAS CINCO COVARIABLES ==\n")
print(round(cor_p["pct_rural", otras], 3))

# ---- B4a. Trayectoria agregando covariables ---------------------------------
# Se agregan en el orden en que aparecen en M1, empezando por pct_rural para
# que el coeficiente exista desde el primer paso.
orden_agregar <- c("pct_rural", setdiff(VARS_M1, "pct_rural"))

trayectoria_ag <- do.call(rbind, lapply(seq_along(orden_agregar), function(k) {
  vars <- orden_agregar[seq_len(k)]
  cbind(paso = k,
        agregada = orden_agregar[k],
        n_covariables = k,
        coef_rural(vars))
}))
trayectoria_ag$signo_mco <- signo(trayectoria_ag$b_mco)
trayectoria_ag$signo_fh  <- signo(trayectoria_ag$b_fh)

cat("\n== COEFICIENTE DE pct_rural AL AGREGAR COVARIABLES (log) ==\n")
print(trayectoria_ag, row.names = FALSE, digits = 4)
write.csv(trayectoria_ag,
          file.path(ruta_out, "ruralidad_trayectoria_agregando.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- B4b. Trayectoria retirando una covariable ------------------------------
trayectoria_ret <- do.call(rbind, lapply(otras, function(v) {
  cbind(retirada = v, coef_rural(setdiff(VARS_M1, v)))
}))
trayectoria_ret$signo_mco <- signo(trayectoria_ret$b_mco)
trayectoria_ret$signo_fh  <- signo(trayectoria_ret$b_fh)

cat("\n== COEFICIENTE DE pct_rural AL RETIRAR UNA COVARIABLE DE M1 (log) ==\n")
print(trayectoria_ret, row.names = FALSE, digits = 4)
write.csv(trayectoria_ret,
          file.path(ruta_out, "ruralidad_trayectoria_retirando.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- B5. Grafico de variable agregada (MCO sobre log) -----------------------
# Residuos de log(pobreza) y de pct_rural tras quitar el efecto de las otras
# cinco. La pendiente de esta nube es el coeficiente ajustado de MCO.
res_y <- resid(lm(as.formula(paste("log_pobreza ~",
                                   paste(otras, collapse = " + "))),
                  data = datos))
res_x <- resid(lm(as.formula(paste("pct_rural ~",
                                   paste(otras, collapse = " + "))),
                  data = datos))

va <- data.frame(res_x = res_x, res_y = res_y, depto = datos$depto)

g_va <- ggplot(va, aes(res_x, res_y)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(aes(colour = depto), alpha = 0.8, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "grey25", linewidth = 0.7) +
  scale_colour_manual(values = COL_DEPTO, name = NULL) +
  labs(x = "pct_rural, descontado el efecto de las otras cinco covariables",
       y = "log(pobreza), descontado el efecto de las otras cinco covariables",
       title = "Asociacion condicional de pct_rural (variable agregada)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(ruta_fig, "fig_ruralidad_variable_agregada.png"), g_va,
       width = 7, height = 5.5, dpi = 200)

cat("\nCorrelacion parcial de pct_rural con log(pobreza), dadas las otras ",
    "cinco: ", round(cor(res_x, res_y), 3), "\n", sep = "")

cat("\nListo. Salidas en output/ y output/figuras/\n")