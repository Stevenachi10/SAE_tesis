# ==============================================================================
# 02_4f_cifras_robusto_mse.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# Cifras para el anexo de sensibilidad del MSE de los ajustes robustos
# (pseudolinealización frente a bootstrap).
#
# No estima nada nuevo. Lee las salidas de 02_4d y 02_4e e imprime en consola
# las cifras ya redondeadas.
#
# Requiere en output/:
#   mse_pseudo_boot_detalle.csv              (02_4e)
#   mse_pseudo_boot_recorte_resumen.csv      (02_4e)
#   mse_pseudo_boot_recorte_municipios.csv   (02_4e)
#   warnholz_verificacion_resumen.csv        (02_4d)
#
# Si no existen, correr antes 02_4d y 02_4e (ambos leen output/robusto.rds).
# ==============================================================================

library(dplyr)
library(here)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")

leer <- function(archivo) {
  ruta <- file.path(ruta_out, archivo)
  if (!file.exists(ruta)) stop("No existe ", ruta, ". Correr antes 02_4d y 02_4e.")
  read.csv(ruta, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
}

detalle <- leer("mse_pseudo_boot_detalle.csv")
rec     <- leer("mse_pseudo_boot_recorte_resumen.csv")
munis   <- leer("mse_pseudo_boot_recorte_municipios.csv")
warn    <- leer("warnholz_verificacion_resumen.csv")

options(width = 200)


# ==============================================================================
# 1. MSE Y CV POR ESPECIFICACIÓN Y MÉTODO (84 municipios)
# ==============================================================================

tab1 <- detalle %>%
  group_by(especificacion, metodo) %>%
  summarise(
    razon_mse_mediana   = round(median(razon_mse), 2),
    n_boot_mayor_pseudo = sum(razon_mse > 1),
    cv_medio_pseudo     = round(100 * mean(cv_pseudo), 2),
    cv_medio_boot       = round(100 * mean(cv_boot), 2),
    cv_p90_pseudo       = round(100 * quantile(cv_pseudo, 0.90, names = FALSE), 2),
    cv_p90_boot         = round(100 * quantile(cv_boot, 0.90, names = FALSE), 2),
    cv_max_pseudo       = round(100 * max(cv_pseudo), 2),
    cv_max_boot         = round(100 * max(cv_boot), 2),
    n_pseudo_menor_dir  = sum(as.logical(pseudo_menor_directo)),
    n_boot_menor_dir    = sum(as.logical(boot_menor_directo)),
    .groups = "drop"
  )


# ==============================================================================
# 2. REBLUP-BC EN LOS MUNICIPIOS CORREGIDOS
# ==============================================================================

tab2 <- warn %>%
  select(especificacion, n_corregidos = n_recortados,
         peso_sesgo_pct = cuota_sesgo_bc_recortados) %>%
  mutate(peso_sesgo_pct = round(100 * peso_sesgo_pct, 1)) %>%
  left_join(
    rec %>%
      filter(as.logical(recortado)) %>%
      select(especificacion, n,
             razon_cv_bc_rb_pseudo = razon_pseudo_med,
             razon_cv_bc_rb_boot   = razon_boot_med),
    by = "especificacion"
  ) %>%
  mutate(across(starts_with("razon"), ~ round(.x, 2)))

stopifnot(
  "El número de municipios corregidos no coincide entre 02_4d y 02_4e" =
    all(tab2$n_corregidos == tab2$n)
)
tab2 <- tab2 %>% select(-n)


# ==============================================================================
# 3. DETALLE POR MUNICIPIO CORREGIDO
# ==============================================================================

tab3 <- munis %>%
  arrange(especificacion, desc(razon_bc_rb_pseudo)) %>%
  transmute(
    especificacion,
    Municipio,
    distancia_ee       = round(distancia_ee, 2),
    cv_reblup_pseudo,
    cv_reblupbc_pseudo,
    cv_reblup_boot,
    cv_reblupbc_boot
  )

# Municipios corregidos en las tres especificaciones
comunes <- munis %>%
  distinct(especificacion, Municipio) %>%
  count(Municipio, name = "n_especificaciones") %>%
  arrange(desc(n_especificaciones), Municipio)


# ==============================================================================
# 4. IMPRESIÓN
# ==============================================================================

cat("\n============================================================\n")
cat("1. MSE Y CV POR ESPECIFICACIÓN Y MÉTODO (84 municipios)\n")
cat("============================================================\n")
cat("razon_mse_mediana: MSE boot / MSE pseudo, mediana\n")
cat("n_boot_mayor_pseudo: municipios con MSE boot > MSE pseudo\n")
cat("n_*_menor_dir: municipios con CV menor al del directo\n\n")
print(as.data.frame(tab1), row.names = FALSE)

cat("\n============================================================\n")
cat("2. REBLUP-BC EN LOS MUNICIPIOS CORREGIDOS\n")
cat("============================================================\n")
cat("peso_sesgo_pct: mediana del sesgo^2 / MSE pseudo de REBLUP-BC (%)\n")
cat("razon_cv_bc_rb: CV REBLUP-BC / CV REBLUP, mediana\n\n")
print(as.data.frame(tab2), row.names = FALSE)

cat("\n============================================================\n")
cat("3. DETALLE POR MUNICIPIO CORREGIDO (CV en %)\n")
cat("============================================================\n\n")
print(as.data.frame(tab3), row.names = FALSE)

cat("\nMunicipios corregidos y en cuántas especificaciones:\n\n")
print(as.data.frame(comunes), row.names = FALSE)

munis %>% filter(Municipio == "Argelia") %>% select(especificacion, cod_mun, distancia_ee)
