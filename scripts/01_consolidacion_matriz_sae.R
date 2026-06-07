# ============================================================
# PIPELINE COMPLETO — SAE POBREZA CAUCA Y VALLE
# Versión con varianza real del estimador directo (cvlog)
# ============================================================

library(dplyr)
library(stringr)
library(sf)
library(geodata)
library(readxl)
library(VIM)
library(here) # Para manejo dinámico de rutas en repositorios

# ============================================================
# CONFIGURACIÓN DE RUTAS (Adaptado para GitHub)
# ============================================================
# Se utilizan rutas relativas basadas en la raíz del proyecto.

ruta_piv <- here("data", "pivoteadas")
ruta_raw <- here("data", "raw")
ruta_out <- here("output")

fix_cod <- function(x) str_pad(as.character(as.integer(x)), width = 5, pad = "0")

# ============================================================
# PARTE 0 — VARIABLE OBJETIVO CON VARIANZA REAL (cvlog, método delta)
# ============================================================

varianza_estimador <- read_excel(file.path(ruta_raw, "varianza estimador.xlsx"))

pobreza <- varianza_estimador %>%
  mutate(cod_mun = fix_cod(`Código Municipio`)) %>%
  filter(`Código Departamento` %in% c("19", "76")) %>%
  filter(!grepl("000$", cod_mun)) %>%
  mutate(
    varianza_pobreza = (`Estimación Pobreza Monetaria` * cvlog)^2   # método delta
  ) %>%
  select(cod_mun, Municipio,
         pobreza_monetaria = `Estimación Pobreza Monetaria`,
         error_estandar    = `Error estandar`,
         varianza_pobreza,
         cvlog)

# ============================================================
# PARTE 1 — BASES CON CÓDIGO DANE
# ============================================================

afiliados <- readRDS(file.path(ruta_piv, "afiliados_regimen_2022.rds")) %>%
  mutate(cod_mun = fix_cod(CodMunicipio)) %>%
  select(cod_mun, C_2022, E_2022, S_2022, I_2022)

demo_etnia <- readRDS(file.path(ruta_piv, "cov_demografia_etnia.rds")) %>%
  mutate(cod_mun = fix_cod(`Código Entidad`)) %>%
  select(cod_mun, pct_indigena_2024, pct_afro_2024,
         pct_menores15_2024, indice_envejecimiento_2024, sisben_grupoA_2024)

demografica <- readRDS(file.path(ruta_piv, "cov_demografica_cauca_valle.rds")) %>%
  mutate(cod_mun = fix_cod(MPIO)) %>%          # MPIO (no DPMP)
  select(cod_mun, pob_total, pct_rural, pct_urbano, tasa_dependencia)

descripcion <- readRDS(file.path(ruta_piv, "cov_descripciongeneral.rds")) %>%
  mutate(cod_mun = fix_cod(`Código Entidad`)) %>%
  select(cod_mun,
         categoria_617 = `Categoría ley 617 de 2000_2024`,
         cat_ruralidad = `Categoría de ruralidad_2000`,
         extension = `Extensión_2017`,
         densidad_2024 = `Densidad poblacional_2024`,
         var_intercensal = `Variación porcentual intercensal 2005 - 2018_2018`)

desempeno <- readRDS(file.path(ruta_piv, "cov_desempeñofiscal.rds")) %>%
  mutate(cod_mun = fix_cod(`Código.Entidad`)) %>%
  select(cod_mun, idf_2024 = `Indice de Desempeño Fiscal 2024`)

educacion <- readRDS(file.path(ruta_piv, "cov_educacion_cauca_valle.rds")) %>%
  mutate(cod_mun = fix_cod(`Código.Entidad`)) %>%
  select(cod_mun,
         cobertura_educ_2024 = `Cobertura neta en educación - Total_2024`,
         desercion_2024 = `Tasa de deserción intra-anual del sector oficial en educación básica y media (Desde transición hasta once)_2024`,
         analfabetismo_2018 = `Tasa de Analfabetismo (Censo)_2018`,
         analfabetismo_rural_2018 = `Tasa de Analfabetismo Rural (Censo)_2018`,
         saber11_mate_2022 = `Puntaje promedio Pruebas Saber 11 - Matemáticas_2022`,
         saber11_lectura_2022 = `Puntaje promedio Pruebas Saber 11 - Lectura crítica_2022`)

ingresos <- readRDS(file.path(ruta_piv, "cov_ejecucion_ingresos.rds")) %>%
  mutate(cod_mun = fix_cod(`Código DANE`)) %>%
  select(cod_mun, INGRESOS, `INGRESOS CORRIENTES`, `INGRESOS TRIBUTARIOS`,
         `INGRESOS NO TRIBUTARIOS`, `TRANSFERENCIAS CORRIENTES`,
         `SISTEMA GENERAL DE PARTICIPACIONES`, `RECURSOS DE CAPITAL`)

fiscal <- readRDS(file.path(ruta_piv, "cov_fiscal.rds")) %>%
  mutate(cod_mun = fix_cod(`Código.DANE`)) %>%
  select(cod_mun, pagos_admin_2024, pagos_educacion_2024,
         pagos_salud_2024, pagos_serv_publicos_2024, pagos_total_2024)

laboral <- readRDS(file.path(ruta_piv, "cov_laboral.rds")) %>%
  mutate(cod_mun = fix_cod(`Código Entidad`)) %>%
  select(cod_mun, cotizantes_2016, empresas_formales_2016, pct_ocupados_form_2016)

pib <- readRDS(file.path(ruta_piv, "cov_pib.rds")) %>%
  mutate(cod_mun = fix_cod(cod_mun)) %>%
  select(cod_mun, pct_primario_2023, pct_secundario_2023,
         pct_terciario_2023, valor_agregado_2023, peso_relativo_2023)

salud <- readRDS(file.path(ruta_piv, "cov_salud.rds")) %>%
  mutate(cod_mun = fix_cod(`Código.Entidad`)) %>%
  select(cod_mun, tasa_mortalidad_2020, mort_materna_2020, mort_infantil_1a_2020,
         cobertura_pentavalente_2020, bajo_peso_2020, control_prenatal_4_2020,
         promedio_prenatales_2020, mort_infantil_5a_2020, mort_eda_2020,
         mort_ira_2020, mort_desnutricion_2020, fecundidad_15_19_2020,
         mort_neonatal_2020, partos_calificados_2020)

seguridad <- readRDS(file.path(ruta_piv, "cov_seguridad.rds")) %>%
  mutate(cod_mun = fix_cod(`Código Entidad`)) %>%
  select(cod_mun, tasa_hurto_2019, tasa_homicidios_2019,
         tasa_extorsion_2019, tasa_violencia_intra_2019)

territorial <- readRDS(file.path(ruta_piv, "cov_territorial_cauca_valle.rds")) %>%
  mutate(cod_mun = fix_cod(`Código Entidad`)) %>%
  select(cod_mun,
         predial_2020 = `Predial_2020`,
         recaudo_predial_2020 = `Recaudo efectivo por impuesto predial_2020`,
         total_predios_2022 = `Total de predios_2022`,
         avaluo_2020 = `Avalúo catastral total_2020`,
         pct_predios_urb_2022 = `Porcentaje de predios urbanos_2022`,
         pct_predios_rur_2022 = `Porcentaje de predios rurales_2022`)

vivienda <- readRDS(file.path(ruta_piv, "cov_vivienda.rds")) %>%
  mutate(cod_mun = fix_cod(`Código Entidad`)) %>%
  select(cod_mun, acueducto_2018, alcantarillado_2018, energia_2018, internet_2018)

conflicto <- readRDS(file.path(ruta_piv, "terridata_conflicto_2023.rds")) %>%
  mutate(cod_mun = fix_cod(`Código Entidad`)) %>%
  select(cod_mun, coca_2023, iica_2023, desplaz_2023, victimas_2023)

# ============================================================
# PARTE 2 — GEOGRÁFICAS (clave doble nombre+depto)
# ============================================================

distancias <- readRDS(file.path(ruta_piv, "cov_distancias.rds"))
carretera  <- readRDS(file.path(ruta_piv, "cov_carretera.rds"))

norm_nombre <- function(x) trimws(tolower(iconv(x, "UTF-8", "ASCII//TRANSLIT")))

corregir_nombre <- function(x) {
  case_when(
    x == "piendamo"         ~ "piendamo - tunia",
    x == "sotara"           ~ "sotara paispamba",
    x == "santiago de cali" ~ "cali",
    TRUE ~ x
  )
}

puente <- pobreza %>%
  mutate(
    municipio_norm = corregir_nombre(norm_nombre(Municipio)),
    depto = ifelse(substr(cod_mun, 1, 2) == "19", "Cauca", "Valle del Cauca")
  ) %>%
  select(cod_mun, municipio_norm, depto)

distancias <- distancias %>%
  mutate(municipio_norm = corregir_nombre(norm_nombre(municipio)),
         depto = departamento)

carretera <- carretera %>%
  mutate(municipio_norm = corregir_nombre(norm_nombre(municipio)))

distancias_cod <- distancias %>%
  left_join(puente, by = c("municipio_norm", "depto")) %>%
  select(cod_mun, dist_popayan_km, dist_cali_km, dist_capital_km)

carretera_cod <- carretera %>%
  left_join(distancias %>% select(municipio_norm, depto), by = "municipio_norm") %>%
  left_join(puente, by = c("municipio_norm", "depto")) %>%
  distinct(cod_mun, .keep_all = TRUE) %>%
  select(cod_mun, dist_carretera_km, snap_km, aislado_vial)

# ============================================================
# PARTE 3 — MERGE SECUENCIAL
# ============================================================

matriz_sae <- pobreza %>%
  left_join(afiliados,      by = "cod_mun") %>%
  left_join(demo_etnia,     by = "cod_mun") %>%
  left_join(demografica,    by = "cod_mun") %>%
  left_join(descripcion,    by = "cod_mun") %>%
  left_join(desempeno,      by = "cod_mun") %>%
  left_join(educacion,      by = "cod_mun") %>%
  left_join(ingresos,       by = "cod_mun") %>%
  left_join(fiscal,         by = "cod_mun") %>%
  left_join(laboral,        by = "cod_mun") %>%
  left_join(pib,            by = "cod_mun") %>%
  left_join(salud,          by = "cod_mun") %>%
  left_join(seguridad,      by = "cod_mun") %>%
  left_join(territorial,    by = "cod_mun") %>%
  left_join(vivienda,       by = "cod_mun") %>%
  left_join(conflicto,      by = "cod_mun") %>%
  left_join(distancias_cod, by = "cod_mun") %>%
  left_join(carretera_cod,  by = "cod_mun")

# ============================================================
# PARTE 4 — CONVERTIR TEXTO A NUMÉRICO
# ============================================================

cols_continuas <- c("extension", "densidad_2024", "var_intercensal",
                    "recaudo_predial_2020", "total_predios_2022", "avaluo_2020",
                    "pct_predios_urb_2022", "pct_predios_rur_2022")

matriz_sae <- matriz_sae %>%
  mutate(across(all_of(cols_continuas),
                ~ as.numeric(gsub(",", ".", gsub("\\.", "", as.character(.)))))) %>%
  select(-predial_2020)   # categórica problemática, se elimina

# ============================================================
# PARTE 5 — DUMMY PDET
# ============================================================

codigos_pdet <- c("19050","19075","19110","19137","19142","19212","19256",
                  "19364","19450","19455","19473","19532","19548","19698",
                  "19780","19821","76109","76275","76563","76233")

matriz_sae <- matriz_sae %>%
  mutate(pdet = ifelse(cod_mun %in% codigos_pdet, 1, 0))

# ============================================================
# PARTE 6 — IMPUTACIONES
# ============================================================

# 6a. Analfabetismo capitales (DANE Censo 2018)
matriz_sae <- matriz_sae %>%
  mutate(
    analfabetismo_2018 = case_when(
      cod_mun == "19001" ~ 3.9, cod_mun == "76001" ~ 3.1,
      TRUE ~ analfabetismo_2018),
    analfabetismo_rural_2018 = case_when(
      cod_mun == "19001" ~ 8.0, cod_mun == "76001" ~ 6.0,
      TRUE ~ analfabetismo_rural_2018)
  )

# 6b. SISBÉN Jambaló con KNN
matriz_sae <- kNN(
  matriz_sae,
  variable = "sisben_grupoA_2024",
  dist_var = c("pct_indigena_2024", "pct_rural", "pob_total",
               "tasa_dependencia", "acueducto_2018", "internet_2018"),
  k = 5
) %>% select(-ends_with("_imp"))

# 6c. dist_carretera litoral (valor alto)
max_dist <- max(matriz_sae$dist_carretera_km, na.rm = TRUE)
matriz_sae <- matriz_sae %>%
  mutate(dist_carretera_km = ifelse(is.na(dist_carretera_km),
                                    max_dist * 1.5, dist_carretera_km))

# ============================================================
# PARTE 7 — NORMALIZACIÓN (per cápita / tasas / proporciones)
# ============================================================

matriz_sae <- matriz_sae %>%
  mutate(
    coca_pct_territorio  = (coca_2023 / (extension * 100)) * 100,
    desplaz_tasa_1000    = (desplaz_2023  / pob_total) * 1000,
    victimas_tasa_1000   = (victimas_2023 / pob_total) * 1000,
    cotizantes_tasa_1000 = (cotizantes_2016 / pob_total) * 1000,
    pagos_admin_pc         = pagos_admin_2024         / pob_total,
    pagos_educacion_pc     = pagos_educacion_2024     / pob_total,
    pagos_salud_pc         = pagos_salud_2024         / pob_total,
    pagos_serv_publicos_pc = pagos_serv_publicos_2024 / pob_total,
    pagos_total_pc         = pagos_total_2024         / pob_total,
    ingresos_pc            = INGRESOS                 / pob_total,
    ingresos_corrientes_pc = `INGRESOS CORRIENTES`    / pob_total,
    ingresos_tributarios_pc= `INGRESOS TRIBUTARIOS`   / pob_total,
    ingresos_no_trib_pc    = `INGRESOS NO TRIBUTARIOS`/ pob_total,
    transferencias_pc      = `TRANSFERENCIAS CORRIENTES`/ pob_total,
    sgp_pc                 = `SISTEMA GENERAL DE PARTICIPACIONES`/ pob_total,
    recursos_capital_pc    = `RECURSOS DE CAPITAL`    / pob_total,
    valor_agregado_pc      = valor_agregado_2023 / pob_total,
    predios_pc             = total_predios_2022 / pob_total,
    avaluo_promedio        = avaluo_2020 / total_predios_2022,
    total_afiliados   = C_2022 + E_2022 + S_2022 + I_2022,
    pct_contributivo  = (C_2022 / total_afiliados) * 100,
    pct_subsidiado    = (S_2022 / total_afiliados) * 100,
    pct_excepcion     = (E_2022 / total_afiliados) * 100,
    pct_sisben_A      = (sisben_grupoA_2024 / pob_total) * 100
  ) %>%
  select(
    -coca_2023, -desplaz_2023, -victimas_2023, -cotizantes_2016,
    -pagos_admin_2024, -pagos_educacion_2024, -pagos_salud_2024,
    -pagos_serv_publicos_2024, -pagos_total_2024,
    -INGRESOS, -`INGRESOS CORRIENTES`, -`INGRESOS TRIBUTARIOS`,
    -`INGRESOS NO TRIBUTARIOS`, -`TRANSFERENCIAS CORRIENTES`,
    -`SISTEMA GENERAL DE PARTICIPACIONES`, -`RECURSOS DE CAPITAL`,
    -valor_agregado_2023, -total_predios_2022, -avaluo_2020,
    -C_2022, -E_2022, -S_2022, -I_2022, -total_afiliados,
    -sisben_grupoA_2024
  )

# ============================================================
# PARTE 8 — VERIFICACIÓN Y GUARDADO
# ============================================================

cat("Dimensiones:", dim(matriz_sae), "\n")
cat("NA totales:", sum(is.na(matriz_sae)), "\n")
na_final <- colSums(is.na(matriz_sae))
print(na_final[na_final > 0])

cat("\nResumen varianza (cvlog, método delta):\n")
print(summary(matriz_sae$varianza_pobreza))

# Guardado usando rutas relativas
saveRDS(matriz_sae, file.path(ruta_out, "matriz_sae_final_v2.rds"))
write.csv(matriz_sae, file.path(ruta_out, "matriz_sae_final_v2.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")