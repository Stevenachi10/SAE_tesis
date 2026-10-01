# ==============================================================================
# 06_0_estimaciones_coef_cv_m1_log.R
# ------------------------------------------------------------------------------
# Modelo final M1 log: coeficientes, estimaciones puntuales por municipio,
# coeficiente de variación e intervalos.
#
# 1. Coeficientes con error estándar, z y valor p de Wald (sin marcas de
#    significancia). No incorporan la incertidumbre de la selección de
#    variables. Coeficientes en escala logarítmica.
# 2. Estimación puntual por municipio (directo y M1), MSE, CV e intervalo
#    aproximado al 95 % (estimación ± 1,96 raíz del MSE, escala original).
# 3. Cifras para el texto de la sección de estimaciones finales.
#
# SALIDAS (output/)
#   final_coeficientes_m1_log.csv
#   final_estimaciones_municipales_m1_log.csv
# ==============================================================================


# ==============================================================================
# 0. PAQUETES Y PARÁMETROS
# ==============================================================================

library(emdi)
library(dplyr)
library(here)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")

AJUSTE_FINAL <- "M1_reml_log"

# Códigos DIVIPOLA de los municipios que se citan en el texto
COD_CALI         <- "76001"
COD_POPAYAN      <- "19001"
COD_BUENAVENTURA <- "76109"

z_975 <- qnorm(0.975)


# ==============================================================================
# 1. CARGA
# ==============================================================================

variantes <- readRDS(file.path(ruta_out, "modelado_variantes.rds"))
matriz    <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

aj <- variantes$ajustes[[AJUSTE_FINAL]]

stopifnot("No se encuentra el ajuste final." = !is.null(aj))

cat("\n============================================================\n")
cat("MODELO FINAL:", AJUSTE_FINAL, "\n")
cat("============================================================\n")
cat("Transformación:",
    paste(unlist(aj$transformation), collapse = " / "), "\n")
cat("Método de estimación:",
    paste(unlist(aj$method$method), collapse = " / "), "\n")
cat("Método del MSE:",
    paste(unlist(aj$method$MSE_method), collapse = " / "), "\n")
cat("sigma_u^2 (escala logarítmica):",
    format(aj$model$variance, digits = 6), "\n")


# ==============================================================================
# 2. COEFICIENTES
# ==============================================================================

coefs <- as.data.frame(aj$model$coefficients)

stopifnot(
  "La tabla de coeficientes no tiene las columnas esperadas." =
    all(c("coefficients", "std.error", "t.value", "p.value") %in% names(coefs))
)

tabla_coef <- data.frame(
  variable = rownames(coefs),
  coef     = coefs$coefficients,
  se       = coefs$std.error,
  z        = coefs$t.value,
  p        = coefs$p.value,
  row.names = NULL,
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("COEFICIENTES (escala logarítmica)\n")
cat("============================================================\n")
print(tabla_coef, row.names = FALSE, digits = 6)


# ------------------------------------------------------------------------------
# Filas listas para LaTeX (coma decimal; notación científica si |coef| < 0,001)
# ------------------------------------------------------------------------------

coma <- function(x) gsub(".", "{,}", x, fixed = TRUE)

coef_tex <- function(x) {
  if (abs(x) < 0.001) {
    e <- floor(log10(abs(x)))
    m <- x / 10^e
    paste0("$", coma(sprintf("%.2f", m)), " \\times 10^{", e, "}$")
  } else {
    paste0("$", coma(sprintf("%.5f", x)), "$")
  }
}

p_tex <- function(p) {
  if (p < 0.001) "$<0{,}001$" else paste0("$", coma(sprintf("%.3f", p)), "$")
}

nombre_tex <- function(v) {
  if (v == "(Intercept)") "Intercepto"
  else paste0("\\texttt{", gsub("_", "\\\\_", v), "}")
}

cat("\nFilas para la tabla de LaTeX:\n\n")
for (i in seq_len(nrow(tabla_coef))) {
  cat(sprintf("%-40s & %-26s & %s \\\\\n",
              nombre_tex(tabla_coef$variable[i]),
              coef_tex(tabla_coef$coef[i]),
              p_tex(tabla_coef$p[i])))
}
cat(sprintf("\n$\\hat{\\sigma}^2_u$ & %s & \\\\\n",
            coef_tex(as.numeric(aj$model$variance))))

write.csv(
  tabla_coef,
  file.path(ruta_out, "final_coeficientes_m1_log.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ==============================================================================
# 3. ESTIMACIONES MUNICIPALES
# ==============================================================================

est <- data.frame(
  cod_mun = as.character(aj$ind$Domain),
  fh      = as.numeric(aj$ind$FH),
  mse_fh  = as.numeric(aj$MSE$FH),
  stringsAsFactors = FALSE
)

stopifnot(
  "El orden de dominios de ind y MSE no coincide." =
    identical(as.character(aj$MSE$Domain), est$cod_mun)
)

base <- matriz %>%
  as.data.frame() %>%
  select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza) %>%
  mutate(cod_mun = as.character(cod_mun))

est <- est %>%
  left_join(base, by = "cod_mun") %>%
  mutate(
    departamento   = ifelse(substr(cod_mun, 1, 2) == "19",
                            "Cauca", "Valle del Cauca"),
    directo_pct    = 100 * pobreza_monetaria,
    cv_directo_pct = 100 * sqrt(varianza_pobreza) / pobreza_monetaria,
    fh_pct         = 100 * fh,
    cv_fh_pct      = 100 * sqrt(mse_fh) / fh,
    ic95_inf_pct   = 100 * (fh - z_975 * sqrt(mse_fh)),
    ic95_sup_pct   = 100 * (fh + z_975 * sqrt(mse_fh)),
    amplitud_pct   = ic95_sup_pct - ic95_inf_pct
  )

stopifnot(
  "No hay 84 municipios." = nrow(est) == 84,
  "Hay municipios sin nombre." = !anyNA(est$Municipio),
  "Hay CV no finitos." = all(is.finite(est$cv_fh_pct))
)

# Control: el directo que guarda emdi debe ser el mismo de la matriz.
dif_directo <- max(abs(as.numeric(aj$ind$Direct) - est$pobreza_monetaria))
cat("\nMáxima diferencia directo emdi vs matriz:",
    format(dif_directo, scientific = TRUE), "\n")
if (dif_directo > 1e-8) {
  warning("El directo de emdi no coincide con la matriz. Revise el orden.")
}

tabla_mun <- est %>%
  select(cod_mun, Municipio, departamento,
         directo_pct, cv_directo_pct,
         fh_pct, mse_fh, cv_fh_pct,
         ic95_inf_pct, ic95_sup_pct, amplitud_pct) %>%
  arrange(departamento, desc(fh_pct))

write.csv(
  tabla_mun,
  file.path(ruta_out, "final_estimaciones_municipales_m1_log.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ==============================================================================
# 4. CIFRAS PARA EL TEXTO
# ==============================================================================

r2 <- function(x) format(round(x, 2), nsmall = 2, decimal.mark = ",")

fila <- function(cod) {
  f <- est[est$cod_mun == cod, ]
  if (nrow(f) != 1) stop("No se encuentra el municipio ", cod)
  f
}

mas_pobre  <- est[which.max(est$fh_pct), ]
menos_pobre <- est[which.min(est$fh_pct), ]
cv_max     <- est[which.max(est$cv_fh_pct), ]
cv_min     <- est[which.min(est$cv_fh_pct), ]

cali  <- fila(COD_CALI)
popa  <- fila(COD_POPAYAN)
buena <- fila(COD_BUENAVENTURA)

med_dpto <- est %>%
  group_by(departamento) %>%
  summarise(mediana = median(fh_pct), .groups = "drop")

cat("\n============================================================\n")
cat("CIFRAS PARA EL TEXTO\n")
cat("============================================================\n")

cat("\nPobreza estimada (%):\n")
cat("  Mínimo: ", r2(menos_pobre$fh_pct), " en ", menos_pobre$Municipio,
    " (", menos_pobre$departamento, ")\n", sep = "")
cat("  Máximo: ", r2(mas_pobre$fh_pct), " en ", mas_pobre$Municipio,
    " (", mas_pobre$departamento, ")\n", sep = "")
for (i in seq_len(nrow(med_dpto))) {
  cat("  Mediana ", med_dpto$departamento[i], ": ",
      r2(med_dpto$mediana[i]), "\n", sep = "")
}
cat("  Buenaventura: ", r2(buena$fh_pct), " (CV ", r2(buena$cv_fh_pct),
    ")\n", sep = "")

cat("\nCoeficiente de variación (%):\n")
cat("  Media:       ", r2(mean(est$cv_fh_pct)), "\n", sep = "")
cat("  Mediana:     ", r2(median(est$cv_fh_pct)), "\n", sep = "")
cat("  Percentil 90:", r2(quantile(est$cv_fh_pct, 0.90)), "\n")
cat("  Máximo:      ", r2(cv_max$cv_fh_pct), " en ", cv_max$Municipio,
    " (", cv_max$departamento, ")\n", sep = "")
cat("  Mínimo:      ", r2(cv_min$cv_fh_pct), " en ", cv_min$Municipio,
    " (", cv_min$departamento, ")\n", sep = "")
cat("  Cali:        ", r2(cali$cv_fh_pct), " (pobreza ",
    r2(cali$fh_pct), ")\n", sep = "")
cat("  Popayán:     ", r2(popa$cv_fh_pct), " (pobreza ",
    r2(popa$fh_pct), ")\n", sep = "")

cat("\nCinco municipios con menor CV:\n")
print(
  est %>%
    arrange(cv_fh_pct) %>%
    select(Municipio, departamento, fh_pct, cv_fh_pct) %>%
    head(5),
  row.names = FALSE, digits = 4
)

cat("\nCinco municipios con mayor CV:\n")
print(
  est %>%
    arrange(desc(cv_fh_pct)) %>%
    select(Municipio, departamento, fh_pct, cv_fh_pct) %>%
    head(5),
  row.names = FALSE, digits = 4
)

cat("\nCinco municipios con mayor pobreza estimada:\n")
print(
  est %>%
    arrange(desc(fh_pct)) %>%
    select(Municipio, departamento, fh_pct, cv_fh_pct) %>%
    head(5),
  row.names = FALSE, digits = 4
)

cat("\nCinco municipios con menor pobreza estimada:\n")
print(
  est %>%
    arrange(fh_pct) %>%
    select(Municipio, departamento, fh_pct, cv_fh_pct) %>%
    head(5),
  row.names = FALSE, digits = 4
)

cat("\nArchivos generados en", ruta_out, ":\n")
cat("  final_coeficientes_m1_log.csv\n")
cat("  final_estimaciones_municipales_m1_log.csv\n")
cat("\nProceso terminado.\n")

library(here)
variantes <- readRDS(here("output", "modelado_variantes.rds"))
format(variantes$ajustes[["M1_reml_log"]]$model$variance, digits = 6)
