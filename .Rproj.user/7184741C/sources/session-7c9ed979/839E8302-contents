# ============================================================
# 02 - TRANSFORMACIONES Y CODIFICACION DE VARIABLES
# Tesis SAE - Pobreza Monetaria Municipal (Cauca y Valle)
# ============================================================
# Entrada:  matriz_sae_final_v2.rds   (84 mun., normalizada, varianza cvlog)
# Salida:   matriz_sae_transformada_v2.rds
# ============================================================

library(dplyr)
library(robustbase)   # funcion mc() para el medcouple
library(ggplot2)
library(tidyr)
library(here)         # Para manejo dinámico de rutas en repositorios

# ============================================================
# CONFIGURACIÓN DE RUTAS (Adaptado para GitHub)
# ============================================================
# Leemos desde output (resultado del script 01) y guardamos en output
ruta_out <- here("output")

# ---- CARGAR MATRIZ NORMALIZADA (v2, con varianza cvlog) ----
matriz_sae <- readRDS(file.path(ruta_out, "matriz_sae_final_v2.rds"))

# ============================================================
# 1. DIAGNOSTICO DE ASIMETRIA CON MEDCOUPLE (ANTES de transformar)
# ============================================================
vars_excluir <- c("cod_mun", "Municipio", "pobreza_monetaria",
                  "varianza_pobreza", "error_estandar", "cvlog",
                  "pdet", "aislado_vial", "snap_km",
                  "categoria_617", "cat_ruralidad")

candidatas <- setdiff(names(matriz_sae), vars_excluir)
candidatas <- candidatas[sapply(matriz_sae[candidatas], is.numeric)]

mc_antes <- sapply(matriz_sae[candidatas], function(x) robustbase::mc(x, na.rm = TRUE))
mc_df_antes <- data.frame(variable = names(mc_antes),
                          mc_antes = round(mc_antes, 3)) %>%
  arrange(desc(abs(mc_antes)))

cat("=== MEDCOUPLE ANTES DE TRANSFORMAR ===\n")
print(mc_df_antes, row.names = FALSE)

# ============================================================
# 2. MANEJO DE VARIABLES CON EXCESO DE CEROS (MC = 1.000)
# ============================================================
vars_sospechosas <- c("coca_pct_territorio", "mort_materna_2020", "mort_eda_2020",
                      "mort_ira_2020", "mort_desnutricion_2020",
                      "pagos_educacion_pc", "pagos_salud_pc",
                      "pagos_serv_publicos_pc", "empresas_formales_2016")

cat("\n=== NUMERO DE CEROS (de 84 municipios) ===\n")
ceros <- sapply(matriz_sae[vars_sospechosas], function(x) sum(x == 0, na.rm = TRUE))
print(ceros)

# ---- 2a. GRAFICO: distribucion de las 9 variables con exceso de ceros ----
df_ceros <- matriz_sae %>%
  select(all_of(vars_sospechosas)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "valor")

g_ceros <- ggplot(df_ceros, aes(x = valor)) +
  geom_histogram(bins = 30, fill = "#2c7fb8", color = "white") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(title = "Distribución de las variables con exceso de ceros (MC = 1.000)",
       x = "Valor", y = "Frecuencia (municipios)") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 8))

print(g_ceros)

ggsave(file.path(ruta_out, "fig_exceso_ceros.png"), g_ceros,
       width = 10, height = 7, dpi = 300)

# Guardar el MC de las variables que se transforman a log (ANTES), para comparar
vars_log <- c("valor_agregado_pc", "avaluo_promedio", "ingresos_tributarios_pc",
              "densidad_2024", "pob_total", "ingresos_pc")
mc_log_antes <- round(sapply(matriz_sae[vars_log], function(x) mc(x, na.rm = TRUE)), 3)

# ---- 2b. Crear dummies de presencia / eliminar las de varianza casi nula ----
matriz_sae <- matriz_sae %>%
  mutate(
    tiene_coca               = ifelse(coca_pct_territorio > 0, 1, 0),
    tiene_mort_materna       = ifelse(mort_materna_2020 > 0, 1, 0),
    tiene_mort_ira           = ifelse(mort_ira_2020 > 0, 1, 0),
    tiene_pagos_salud        = ifelse(pagos_salud_pc > 0, 1, 0),
    tiene_empresas_formales  = ifelse(empresas_formales_2016 > 0, 1, 0)
  ) %>%
  select(
    -pagos_serv_publicos_pc,   # 92% ceros
    -pagos_educacion_pc,       # 88% ceros
    -mort_eda_2020,            # 89% ceros
    -mort_desnutricion_2020,   # 89% ceros
    -mort_materna_2020,        # -> tiene_mort_materna
    -mort_ira_2020,            # -> tiene_mort_ira
    -pagos_salud_pc,           # -> tiene_pagos_salud
    -empresas_formales_2016    # -> tiene_empresas_formales
  )

# ---- 2c. GRAFICO: balance de las dummies creadas ----
vars_dummy <- c("tiene_coca", "tiene_mort_materna", "tiene_mort_ira",
                "tiene_pagos_salud", "tiene_empresas_formales")

df_dummy <- matriz_sae %>%
  select(all_of(vars_dummy)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "valor") %>%
  mutate(valor = factor(valor, levels = c(0, 1),
                        labels = c("Ausencia (0)", "Presencia (1)"))) %>%
  count(variable, valor)

g_dummy <- ggplot(df_dummy, aes(x = variable, y = n, fill = valor)) +
  geom_col(position = "stack") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            color = "white", size = 3.5) +
  scale_fill_manual(values = c("Ausencia (0)" = "#d95f0e",
                               "Presencia (1)" = "#2c7fb8")) +
  labs(title = "Balance de las variables indicadoras (dummy) creadas",
       x = NULL, y = "Número de municipios", fill = NULL) +
  coord_flip() +
  theme_minimal(base_size = 11)

print(g_dummy)

ggsave(file.path(ruta_out, "fig_dummies_balance.png"), g_dummy,
       width = 9, height = 5, dpi = 300)

# ============================================================
# 3. TRANSFORMACION LOGARITMICA (sesgo continuo real, sin ceros)
# ============================================================
cat("\n=== MINIMOS ANTES DEL LOG (deben ser > 0) ===\n")
print(sapply(matriz_sae[vars_log], min, na.rm = TRUE))

matriz_sae <- matriz_sae %>%
  mutate(across(all_of(vars_log), ~ log(.), .names = "log_{.col}")) %>%
  select(-all_of(vars_log))

# ============================================================
# 4. COMPARACION MEDCOUPLE: ANTES vs DESPUES DEL LOG
# ============================================================
mc_log_despues <- round(sapply(matriz_sae[paste0("log_", vars_log)],
                               function(x) mc(x, na.rm = TRUE)), 3)

comparacion_log <- data.frame(
  variable   = vars_log,
  mc_antes   = mc_log_antes,
  mc_despues = as.numeric(mc_log_despues),
  row.names  = NULL
)

cat("\n=== MEDCOUPLE ANTES vs DESPUES DEL LOG ===\n")
print(comparacion_log, row.names = FALSE)

write.csv(comparacion_log, file.path(ruta_out, "comparacion_medcouple_log.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ============================================================
# 5. CODIFICACION DE VARIABLES CATEGORICAS
# ============================================================

# ---- 5a. categoria_617: ordinal (Especial -> 0, luego 1 a 6) ----
matriz_sae <- matriz_sae %>%
  mutate(
    categoria_617 = case_when(
      categoria_617 == "Especial" ~ "0",
      TRUE ~ as.character(categoria_617)
    ),
    categoria_617 = as.numeric(categoria_617)
  )

cat("\n=== categoria_617 (NA esperado = 0) ===\n")
cat("NA:", sum(is.na(matriz_sae$categoria_617)),
    "| Valores:", paste(sort(unique(matriz_sae$categoria_617)), collapse = ", "), "\n")

# ---- 5b. cat_ruralidad: dummies, referencia = "Ciudades y aglomeraciones" ----
matriz_sae <- matriz_sae %>%
  mutate(
    rural_intermedio = ifelse(cat_ruralidad == "Intermedio", 1, 0),
    rural            = ifelse(cat_ruralidad == "Rural", 1, 0),
    rural_disperso   = ifelse(cat_ruralidad == "Rural disperso", 1, 0)
  ) %>%
  select(-cat_ruralidad)

# ============================================================
# 6. VERIFICACION Y GUARDADO
# ============================================================
cat("\n=== VERIFICACION FINAL ===\n")
cat("Dimensiones:", dim(matriz_sae), "\n")
cat("NA totales:", sum(is.na(matriz_sae)), "\n")

saveRDS(matriz_sae, file.path(ruta_out, "matriz_sae_transformada_v2.rds"))
write.csv(matriz_sae, file.path(ruta_out, "matriz_sae_transformada_v2.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\nMatriz transformada (v2) guardada correctamente.\n")