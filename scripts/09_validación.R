# ============================================================
# 07 - VALIDACIÓN DE LA PRECISIÓN (MSE) DEL MODELO M09
#      Analítico (Prasad-Rao)  vs  Bootstrap paramétrico
# Tesis SAE - Pobreza Monetaria Municipal (Cauca y Valle del Cauca)
# ============================================================
# Entrada: output/matriz_sae_transformada_v2.rds
# Salidas: (consola) tabla comparativa + correlación de los MSE
#          output/fig_mse_validacion.png   (opcional, para anexo)


library(sae)
library(dplyr)
library(ggplot2)
library(here)

matriz <- readRDS(here("output", "matriz_sae_transformada_v2.rds")) %>%
  mutate(cod_mun = as.character(cod_mun))

f_m09 <- pobreza_monetaria ~ extension + saber11_lectura_2022 +
  partos_calificados_2020 + dist_popayan_km + dist_cali_km +
  pct_contributivo + rural_disperso + tiene_coca +
  desplaz_tasa_1000 + categoria_617

# ------------------------------------------------------------
# 1. MSE analítico (Prasad-Rao) -> medida OFICIAL
# ------------------------------------------------------------
fh_an  <- mseFH(f_m09, vardir = varianza_pobreza, method = "REML", data = matriz)
eblup  <- fh_an$est$eblup
mse_an <- fh_an$mse
cv_an  <- sqrt(mse_an) / eblup * 100

# ------------------------------------------------------------
# 2. MSE por bootstrap paramétrico (loop manual, optimizado)
# ------------------------------------------------------------
set.seed(2906)   # reproducibilidad
B <- 1000         # sube a B = 1000 para la versión final
D <- nrow(matriz)

# Parámetros del modelo ajustado en la Sección 1
sigma2_u_hat <- fh_an$est$fit$refvar
beta_hat     <- fh_an$est$fit$estcoef$beta
Xbeta        <- as.vector(model.matrix(f_m09, data = matriz) %*% beta_hat)  # parte fija

# Objetos invariantes, fuera del loop
matriz_boot  <- matriz
f_boot       <- update(f_m09, pob_boot ~ .)
mse_boot_sim <- matrix(NA_real_, nrow = D, ncol = B)

cat("Iniciando bootstrap paramétrico (B =", B, ")...\n")
for (b in 1:B) {
 
  mu_boot <- Xbeta + rnorm(D, 0, sqrt(sigma2_u_hat))             # verdadero
  y_boot  <- mu_boot + rnorm(D, 0, sqrt(matriz$varianza_pobreza)) # observado
  
  matriz_boot$pob_boot <- y_boot
  fh_sim <- try(eblupFH(f_boot, vardir = varianza_pobreza, method = "REML",
                        data = matriz_boot), silent = TRUE)
  
  if (!inherits(fh_sim, "try-error"))
    mse_boot_sim[, b] <- (fh_sim$eblup - mu_boot)^2   
  
  if (b %% 100 == 0) cat("  Iteración", b, "completada...\n")
}

mse_pb <- rowMeans(mse_boot_sim, na.rm = TRUE)
cv_pb  <- sqrt(mse_pb) / eblup * 100   # CV con el EBLUP original
cat(sprintf("Réplicas válidas: %d de %d\n",
            sum(colSums(!is.na(mse_boot_sim)) > 0), B))

# ------------------------------------------------------------
# 3. Comparación
# ------------------------------------------------------------
comp <- data.frame(
  Metodo   = c("Analitico (Prasad-Rao)", "Bootstrap parametrico (B=500)"),
  CV_medio = round(c(mean(cv_an), mean(cv_pb)), 2),
  CV_min   = round(c(min(cv_an),  min(cv_pb)),  2),
  CV_max   = round(c(max(cv_an),  max(cv_pb)),  2)
)
cat("\n--- Validación del MSE: analítico vs bootstrap ---\n")
print(comp, row.names = FALSE)

cat(sprintf("\nCorrelación MSE analítico vs bootstrap: %.4f\n", cor(mse_an, mse_pb)))
cat(sprintf("Diferencia relativa media en el MSE: %.2f %%\n",
            mean(abs(mse_pb - mse_an) / mse_an) * 100))

# ------------------------------------------------------------
# 4. Figura de validación (opcional, para el anexo)
#    Cada punto = un municipio; la diagonal es la identidad.
# ------------------------------------------------------------
g_mse <- ggplot(data.frame(cv_an = cv_an, cv_pb = cv_pb), aes(cv_an, cv_pb)) +
  geom_abline(slope = 1, intercept = 0, color = "#1a1814", linetype = "dashed") +
  geom_point(color = "#2c6e6b", size = 1.6, alpha = 0.85) +
  labs(title = "Validación del MSE: analítico vs. bootstrap paramétrico",
       subtitle = "Cada punto es un municipio; la línea discontinua es la identidad",
       x = "CV analítico (%)", y = "CV bootstrap (%)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", color = "#1a1814"))

ggsave(here("output", "fig_mse_validacion.png"), g_mse,
       width = 6, height = 5, dpi = 300, bg = "white")
cat("\nFigura guardada: output/fig_mse_validacion.png\n")
mean(mse_an)                 # MSE medio (diminuto, NO presentar crudo)
mean(sqrt(mse_an)) * 100     # RMSE medio en puntos porcentuales (= error estándar) 

# Intervalos de confianza al 95% (si pobreza está en proporción 0-1):
li <- (eblup - 1.96 * sqrt(mse_an)) * 100   # límite inferior en %
ls <- (eblup + 1.96 * sqrt(mse_an)) * 100   # límite superior en %
li
ls
mean(sqrt(matriz$varianza_pobreza)) * 100
