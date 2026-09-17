# ==============================================================================
# 01_etapa1_construccion_base.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas (Fay-Herriot). Departamentos de Cauca y Valle del Cauca, 2024.
#
# ETAPA I — Construcción de la base de datos integrada.
#
# Consolida la matriz de 84 dominios municipales que alimenta la selección de
# variables (Etapa II), aplica las transformaciones derivadas del diagnóstico
# de asimetría y produce los insumos de trazabilidad que se reportan en el
# capítulo de resultados.
#
# ------------------------------------------------------------------------------
# CÓMO EJECUTARLO
# ------------------------------------------------------------------------------
# 1. Abrir el proyecto de RStudio en la raíz del repositorio; las rutas se
#    resuelven con here() y no requieren ajuste manual.
# 2. Verificar que existan las carpetas data/raw, data/pivoteadas y output.
# 3. Ejecutar el script completo. Tiempo aproximado: dos minutos.
# 4. Cerrar cualquier archivo de output/ abierto en Excel antes de ejecutar:
#    Windows bloquea la escritura sobre archivos abiertos.
#
# Los archivos .rds de data/pivoteadas corresponden a las descargas de las
# fuentes oficiales ya depuradas y pivoteadas a formato ancho (una fila por
# municipio). Los scripts de ese procesamiento previo no forman parte del
# repositorio; el inventario de cada archivo (fuente, dimensiones y columnas)
# se documenta en data/pivoteadas/README.md. Este script no descarga datos:
# asume que los insumos están disponibles.
#
# ------------------------------------------------------------------------------
# ENTRADAS
# ------------------------------------------------------------------------------
#   data/raw/varianza estimador.xlsx
#       Estimaciones de pobreza monetaria municipal y errores estándar
#       (DANE, metodología SAE, 2024). Variable dependiente del estudio.
#
#   data/pivoteadas/*.rds
#       Covariables auxiliares por dimensión temática. Fuentes: TerriData
#       (DNP), CUIPO, Censo 2018 (DANE), Sisbén IV, UNODC, ADRES, IGAC, GADM
#       y OSRM.
#
# ------------------------------------------------------------------------------
# SALIDAS (output/)
# ------------------------------------------------------------------------------
#   matriz_sae_final_v2.rds / .csv         Matriz en escala original, base del
#                                          análisis descriptivo.
#   matriz_sae_transformada_v2.rds / .csv  Matriz con las transformaciones
#                                          aplicadas, insumo de la Etapa II.
#   listado_covariables.txt                Nombres de las covariables candidatas.
#   diagnostico_faltantes_pre.csv          Faltantes antes de imputar.
#   descriptivo_cuantitativas.csv/.tex     Estadísticas de variables continuas.
#   descriptivo_cualitativas.csv/.tex      Frecuencias de categóricas.
#   descriptivo_departamentos.csv          Contraste Cauca y Valle del Cauca.
#   diagnostico_asimetria.csv              Medcouple y proporción de ceros.
#   comparacion_medcouple_log.csv          Asimetría antes y después del log.
#   fig_exceso_ceros.png                   Distribución de variables con ceros.
#   fig_dummies_balance.png                Balance de las indicadoras creadas.
#
# ------------------------------------------------------------------------------
# DEPENDENCIAS
# ------------------------------------------------------------------------------
#   dplyr, stringr, readxl, VIM, here, robustbase, tidyr, ggplot2, knitr
#
# ------------------------------------------------------------------------------
# DECISIONES METODOLÓGICAS IMPLEMENTADAS
# ------------------------------------------------------------------------------
# Las decisiones que siguen están documentadas en el capítulo de metodología.
# Se enuncian aquí para que el código sea legible sin consultar el documento.
#
#   Varianza de muestreo. La fuente publica la estimación y su coeficiente de
#   variación en escala logarítmica; la varianza se reconstruye por el método
#   delta. El script imprime su diferencia con el error estándar publicado
#   como verificación de consistencia.
#
#   Clave de unión. Código Divipola de cinco dígitos. Las fuentes
#   cartográficas no incluyen ese código, por lo que se unen mediante el
#   nombre normalizado del municipio combinado con el departamento; la clave
#   doble es necesaria porque Argelia y Bolívar existen en ambos
#   departamentos.
#
#   Municipios sin conexión vial. Guapi, López de Micay y Timbiquí carecen de
#   ruta terrestre hacia la capital departamental. No se imputa una distancia
#   inexistente: la variable se excluye del conjunto de candidatas y el
#   aislamiento se captura mediante la indicadora `aislado_vial`, definida por
#   una distancia de enganche (snap) de OSRM superior a 20 km.
#
#   Analfabetismo de las capitales. TerriData no reporta el dato para Popayán
#   ni Cali. Se recupera de los planes de desarrollo municipales: 2,97 % para
#   Popayán (PDM 2024-2027) y 3,2 % para Cali (PDM 2016-2019). Ninguna de las
#   dos fuentes desagrega el dato rural, por lo que la tasa total se emplea
#   como aproximación.
#
#   Colinealidad exacta. `pct_urbano` y `pct_predios_rur_2022` son
#   complementos aritméticos de `pct_rural` y `pct_predios_urb_2022`. Se
#   conserva una variable de cada par.
#
#   Valores atípicos. No se aplica una detección previa ni se elimina ningún
#   municipio: los valores extremos reflejan condiciones reales y su
#   influencia se trata en la Etapa II con el modelo Fay-Herriot robusto. En
#   esta etapa solo se verifica que los valores sean coherentes con la
#   naturaleza de cada variable (control de rangos, Parte 7).
#
#   Escala del análisis descriptivo. Las estadísticas descriptivas se calculan
#   sobre las variables en escala original, antes de aplicar las
#   transformaciones logarítmicas. La asimetría se refleja mediante la mediana
#   y el medcouple.
#
#   Criterios de transformación. Un medcouple igual o muy cercano a la unidad
#   no refleja asimetría continua sino concentración de ceros. Las variables
#   con más del 88 % de ceros se eliminan por su escasa variabilidad; las que
#   presentan entre 75 % y 83 % se recodifican como indicadoras de presencia;
#   cuando la relevancia territorial del fenómeno lo justifica se conservan
#   tanto la indicadora como la medida de intensidad. Las variables con
#   asimetría positiva moderada y sin ceros reciben transformación
#   logarítmica.
#
# ------------------------------------------------------------------------------
# LIMITACIÓN CONOCIDA
# ------------------------------------------------------------------------------
#   `cov_carretera.rds` no incluye la columna de departamento, por lo que la
#   unión de esa fuente recae en el nombre del municipio y los homónimos
#   requieren verificación manual, que el script imprime en consola. Si se
#   agrega la columna `departamento` a ese archivo, el script la detecta y
#   aplica la clave doble de forma automática.
# ==============================================================================

library(dplyr)
library(stringr)
library(readxl)
library(VIM)
library(here)
library(robustbase)
library(tidyr)
library(ggplot2)

ruta_piv <- here("data", "pivoteadas")
ruta_raw <- here("data", "raw")
ruta_out <- here("output")

# Normaliza el código Divipola a cinco dígitos con ceros a la izquierda.
fix_cod <- function(x) str_pad(as.character(as.integer(x)), width = 5, pad = "0")

# Tratamiento de los municipios sin conexión vial.
#   FALSE (configuración del estudio): se excluye la distancia por carretera.
#   TRUE: imputa un valor centinela de 1,5 veces el máximo observado.
#         Alternativa evaluada y descartada; se conserva para replicación.
usar_centinela <- FALSE

# Tema común de las figuras. Sin título: la descripción de cada figura se
# incorpora en el caption del documento. Tipografía amplia para que las
# etiquetas resulten legibles con el documento al 100 %.
tema_tesis <- theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_blank(),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12),
    legend.text      = element_text(size = 12),
    strip.text       = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  )

COL_PRINCIPAL <- "#2c7fb8"
COL_CONTRASTE <- "#d95f0e"


# ==============================================================================
# PARTE 0 — VARIABLE DEPENDIENTE
# ------------------------------------------------------------------------------
# La fuente entrega la estimación de pobreza monetaria y su coeficiente de
# variación en escala logarítmica. La varianza se reconstruye por el método
# delta: Var(theta) = (theta * cvlog)^2.
# ==============================================================================

varianza_estimador <- read_excel(file.path(ruta_raw, "varianza estimador.xlsx"))

n_total_fuente <- nrow(varianza_estimador)

pobreza_pre <- varianza_estimador %>%
  mutate(cod_mun = fix_cod(`Código Municipio`)) %>%
  filter(`Código Departamento` %in% c("19", "76"))

n_cauca_valle   <- nrow(pobreza_pre)
n_agregados_dep <- sum(grepl("000$", pobreza_pre$cod_mun))

pobreza <- pobreza_pre %>%
  filter(!grepl("000$", cod_mun)) %>%
  mutate(varianza_pobreza = (`Estimación Pobreza Monetaria` * cvlog)^2) %>%
  select(cod_mun, Municipio,
         pobreza_monetaria = `Estimación Pobreza Monetaria`,
         error_estandar    = `Error estandar`,
         varianza_pobreza,
         cvlog)

cat("== DEPURACIÓN DE LA VARIABLE DEPENDIENTE ==\n")
cat("Registros en la fuente:", n_total_fuente, "\n")
cat("Registros de Cauca y Valle:", n_cauca_valle, "\n")
cat("Agregados departamentales eliminados:", n_agregados_dep, "\n")
cat("Dominios municipales:", nrow(pobreza), "\n\n")

# Verificación de consistencia. El error estándar implícito en la varianza
# reconstruida debería coincidir con el publicado por la fuente; diferencias
# apreciables indicarían que cvlog no corresponde al error estándar del
# logaritmo. Solo se imprime: no modifica ningún objeto.
dif_ee <- sqrt(pobreza$varianza_pobreza) -
  suppressWarnings(as.numeric(pobreza$error_estandar))

cat("== CONSISTENCIA DE LA VARIANZA DE MUESTREO ==\n")
cat("Diferencia entre el error estándar reconstruido y el publicado:\n")
print(summary(dif_ee))
cat("\n")


# ==============================================================================
# PARTE 1 — COVARIABLES CON CÓDIGO DIVIPOLA
# ------------------------------------------------------------------------------
# Cada fuente se carga desde su archivo pivoteado, se estandariza la clave y
# se seleccionan las columnas de interés. La correspondencia entre estas
# variables y las dimensiones temáticas se detalla en la sección de datos del
# documento.
# ==============================================================================

afiliados <- readRDS(file.path(ruta_piv, "afiliados_regimen_2022.rds")) %>%
  mutate(cod_mun = fix_cod(CodMunicipio)) %>%
  select(cod_mun, C_2022, E_2022, S_2022, I_2022)

demo_etnia <- readRDS(file.path(ruta_piv, "cov_demografia_etnia.rds")) %>%
  mutate(cod_mun = fix_cod(`Código Entidad`)) %>%
  select(cod_mun, pct_indigena_2024, pct_afro_2024,
         pct_menores15_2024, indice_envejecimiento_2024, sisben_grupoA_2024)

demografica <- readRDS(file.path(ruta_piv, "cov_demografica_cauca_valle.rds")) %>%
  mutate(cod_mun = fix_cod(MPIO)) %>%
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


# ==============================================================================
# PARTE 2 — VARIABLES GEOGRÁFICAS
# ------------------------------------------------------------------------------
# Las distancias se construyeron a partir de la cartografía municipal de GADM
# (centroides) y del motor de enrutamiento OSRM. Estas fuentes identifican los
# municipios por nombre y no por código, de modo que la unión requiere
# normalizar el nombre y combinarlo con el departamento.
#
# La normalización corrige además las diferencias de denominación detectadas
# entre fuentes: Santiago de Cali frente a Cali, Piendamó frente a Piendamó -
# Tunía, y Sotará frente a Sotará Paispamba.
# ==============================================================================

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

# Tabla puente entre el código Divipola y el par (nombre normalizado, depto).
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

# Si el archivo de carretera incluye el departamento se aplica la clave doble.
# En caso contrario se une por nombre y se imprimen los municipios homónimos
# para su revisión: la distancia por carretera nunca puede ser inferior a la
# distancia geodésica al mismo destino.
if ("departamento" %in% names(carretera)) {
  
  carretera_cod <- carretera %>%
    mutate(depto = departamento) %>%
    left_join(puente, by = c("municipio_norm", "depto")) %>%
    distinct(cod_mun, .keep_all = TRUE) %>%
    select(cod_mun, dist_carretera_km, snap_km, aislado_vial)
  
  cat("Variables de carretera unidas con clave doble (municipio y depto).\n\n")
  
} else {
  
  carretera_cod <- carretera %>%
    left_join(distancias %>% select(municipio_norm, depto), by = "municipio_norm") %>%
    left_join(puente, by = c("municipio_norm", "depto")) %>%
    distinct(cod_mun, .keep_all = TRUE) %>%
    select(cod_mun, dist_carretera_km, snap_km, aislado_vial)
  
  cat("AVISO: el archivo de carretera no incluye departamento.\n")
  cat("Revisar la asignación de los municipios homónimos:\n")
  print(
    carretera_cod %>%
      filter(cod_mun %in% c("19050", "76054", "19100", "76100")) %>%
      left_join(distancias_cod, by = "cod_mun")
  )
  cat("\n")
}


# ==============================================================================
# PARTE 3 — INTEGRACIÓN DE LAS FUENTES
# ------------------------------------------------------------------------------
# La tabla de pobreza fija el universo de estudio. Cada fuente se incorpora
# mediante left_join sobre el código Divipola, de modo que la matriz conserva
# exactamente los 84 dominios.
# ==============================================================================

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


# ==============================================================================
# PARTE 4 — ARMONIZACIÓN DE FORMATOS NUMÉRICOS
# ------------------------------------------------------------------------------
# Varias variables de las fuentes catastral y territorial llegan con
# convención hispanohablante (punto de miles, coma decimal) y se convierten al
# estándar del lenguaje R.
# ==============================================================================

cols_texto <- c("extension", "densidad_2024", "var_intercensal",
                "recaudo_predial_2020", "total_predios_2022", "avaluo_2020",
                "pct_predios_urb_2022", "pct_predios_rur_2022")

a_numerico <- function(x) {
  if (is.numeric(x)) return(x)
  as.numeric(gsub(",", ".", gsub("\\.", "", as.character(x))))
}

matriz_sae <- matriz_sae %>%
  mutate(across(all_of(cols_texto), a_numerico))


# ==============================================================================
# PARTE 5 — INDICADORA DE MUNICIPIO PDET
# ------------------------------------------------------------------------------
# Programas de Desarrollo con Enfoque Territorial, según el listado oficial de
# la Agencia de Renovación del Territorio.
# ==============================================================================

codigos_pdet <- c("19050","19075","19110","19137","19142","19212","19256",
                  "19364","19450","19455","19473","19532","19548","19698",
                  "19780","19821","76109","76275","76563","76233")

matriz_sae <- matriz_sae %>%
  mutate(pdet = ifelse(cod_mun %in% codigos_pdet, 1, 0))


# ==============================================================================
# PARTE 6 — DATOS FALTANTES
# ------------------------------------------------------------------------------
# El diagnóstico se ejecuta antes de imputar, de modo que quede registro de
# qué variables y qué municipios presentaban ausencias.
#
# Se aplican dos estrategias según la naturaleza del dato ausente:
#   (a) recuperación desde una fuente oficial alternativa cuando el dato
#       existe y es verificable;
#   (b) imputación por k vecinos más cercanos cuando el dato no está
#       disponible en ninguna fuente.
# Los municipios sin conexión vial no se imputan.
# ==============================================================================

diag_na <- matriz_sae %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_na") %>%
  filter(n_na > 0) %>%
  mutate(pct_na = round(100 * n_na / nrow(matriz_sae), 1)) %>%
  arrange(desc(n_na))

cat("== DATOS FALTANTES ANTES DE IMPUTAR ==\n")
print(as.data.frame(diag_na))

cat("\nMunicipios afectados:\n")
print(as.data.frame(
  matriz_sae %>% filter(if_any(everything(), is.na)) %>% select(cod_mun, Municipio)
))
cat("\n")

write.csv(diag_na, file.path(ruta_out, "diagnostico_faltantes_pre.csv"),
          row.names = FALSE)

# (a) Analfabetismo de las capitales, recuperado de los planes de desarrollo.
matriz_sae <- matriz_sae %>%
  mutate(
    analfabetismo_2018 = case_when(
      cod_mun == "19001" ~ 2.97,   # Popayán, PDM 2024-2027
      cod_mun == "76001" ~ 3.20,   # Cali, PDM 2016-2019
      TRUE ~ analfabetismo_2018),
    analfabetismo_rural_2018 = case_when(
      cod_mun == "19001" ~ 2.97,   # la tasa total aproxima la rural
      cod_mun == "76001" ~ 3.20,
      TRUE ~ analfabetismo_rural_2018)
  )

# (b) Clasificación Sisbén de Jambaló, imputada por KNN con predictoras
#     asociadas a la focalización del instrumento.
matriz_sae <- kNN(
  matriz_sae,
  variable = "sisben_grupoA_2024",
  dist_var = c("pct_indigena_2024", "pct_rural", "pob_total",
               "tasa_dependencia", "acueducto_2018", "internet_2018"),
  k = 5
) %>% select(-ends_with("_imp"))

# (c) Municipios sin conexión vial.
if (usar_centinela) {
  
  max_dist <- max(matriz_sae$dist_carretera_km, na.rm = TRUE)
  matriz_sae <- matriz_sae %>%
    mutate(dist_carretera_km = ifelse(is.na(dist_carretera_km),
                                      max_dist * 1.5, dist_carretera_km))
  cat("Distancia por carretera: centinela aplicado.\n\n")
  
} else {
  
  matriz_sae <- matriz_sae %>% select(-dist_carretera_km)
  cat("Distancia por carretera excluida del conjunto de candidatas.\n\n")
  
}


# ==============================================================================
# PARTE 7 — CONSTRUCCIÓN DE MEDIDAS RELATIVAS
# ------------------------------------------------------------------------------
# Las variables en valores absolutos se transforman para hacerlas comparables
# entre municipios de distinto tamaño:
#   tasas por mil habitantes para los conteos de personas;
#   valores per cápita para las magnitudes monetarias y los predios;
#   proporciones respecto a un total para las variables de composición.
#
# El área de coca se expresa como porcentaje del territorio municipal; la
# extensión viene en km2 y el cultivo en hectáreas, de ahí el factor 100.
#
# Las columnas en niveles se descartan una vez construida su versión relativa.
# Al final se verifica que los valores estén dentro de su rango admisible.
# ==============================================================================

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
    
    ingresos_pc             = INGRESOS                          / pob_total,
    ingresos_corrientes_pc  = `INGRESOS CORRIENTES`             / pob_total,
    ingresos_tributarios_pc = `INGRESOS TRIBUTARIOS`            / pob_total,
    ingresos_no_trib_pc     = `INGRESOS NO TRIBUTARIOS`         / pob_total,
    transferencias_pc       = `TRANSFERENCIAS CORRIENTES`       / pob_total,
    sgp_pc                  = `SISTEMA GENERAL DE PARTICIPACIONES` / pob_total,
    recursos_capital_pc     = `RECURSOS DE CAPITAL`             / pob_total,
    
    valor_agregado_pc = valor_agregado_2023 / pob_total,
    predios_pc        = total_predios_2022  / pob_total,
    avaluo_promedio   = avaluo_2020 / total_predios_2022,
    
    total_afiliados  = C_2022 + E_2022 + S_2022 + I_2022,
    pct_contributivo = (C_2022 / total_afiliados) * 100,
    pct_subsidiado   = (S_2022 / total_afiliados) * 100,
    pct_excepcion    = (E_2022 / total_afiliados) * 100,
    pct_sisben_A     = (sisben_grupoA_2024 / pob_total) * 100
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

# Complementos aritméticos: se conserva una variable de cada par.
matriz_sae <- matriz_sae %>%
  select(-pct_urbano, -pct_predios_rur_2022)

# Control de rangos admisibles. Un valor fuera de rango indica un problema de
# escala en la fuente y debe resolverse antes de continuar.
cat("== CONTROL DE RANGOS ==\n")
print(list(
  "coca_pct_territorio (0-100)" = range(matriz_sae$coca_pct_territorio),
  "pct_contributivo (0-100)"    = range(matriz_sae$pct_contributivo),
  "pct_subsidiado (0-100)"      = range(matriz_sae$pct_subsidiado),
  "pct_sisben_A (0-100)"        = range(matriz_sae$pct_sisben_A),
  "razón pagos/ingresos"        = range(matriz_sae$pagos_total_pc /
                                          matriz_sae$ingresos_pc)
))
cat("\n")

stopifnot(
  "El área de coca supera la extensión municipal" =
    max(matriz_sae$coca_pct_territorio) <= 100,
  "El gasto supera ampliamente los ingresos: revisar la fuente fiscal" =
    max(matriz_sae$pagos_total_pc / matriz_sae$ingresos_pc) < 2
)


# ==============================================================================
# PARTE 8 — MATRIZ EN ESCALA ORIGINAL
# ------------------------------------------------------------------------------
# Esta versión conserva las unidades interpretables y es la base del análisis
# descriptivo. Las transformaciones se aplican más adelante, sobre una copia.
# ==============================================================================

cols_id_dep <- c("cod_mun", "Municipio", "pobreza_monetaria",
                 "error_estandar", "varianza_pobreza", "cvlog")

cat("== MATRIZ EN ESCALA ORIGINAL ==\n")
cat("Dimensiones:", dim(matriz_sae)[1], "x", dim(matriz_sae)[2], "\n")
cat("Covariables candidatas:", ncol(matriz_sae) - length(cols_id_dep), "\n")
cat("Datos ausentes:", sum(is.na(matriz_sae)), "\n\n")

stopifnot(
  "El universo no contiene 84 dominios"  = nrow(matriz_sae) == 84,
  "La matriz contiene datos ausentes"    = sum(is.na(matriz_sae)) == 0
)

writeLines(setdiff(names(matriz_sae), cols_id_dep),
           file.path(ruta_out, "listado_covariables.txt"))

saveRDS(matriz_sae, file.path(ruta_out, "matriz_sae_final_v2.rds"))
write.csv(matriz_sae, file.path(ruta_out, "matriz_sae_final_v2.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# PARTE 9 — ANÁLISIS DESCRIPTIVO
# ------------------------------------------------------------------------------
# Se calcula sobre las variables en escala original, con unidades
# interpretables. Produce los tres insumos del análisis descriptivo del
# capítulo de resultados, en CSV y en LaTeX.
#
# Las covariables continuas se identifican aquí y se reutilizan en el
# diagnóstico de asimetría de la Parte 10. Las indicadoras quedan fuera
# porque las medidas de forma no aplican sobre ellas.
# ==============================================================================

es_continua <- function(v) is.numeric(v) && length(unique(na.omit(v))) >= 3

covs_cont <- matriz_sae %>%
  select(where(is.numeric), -pobreza_monetaria, -error_estandar,
         -varianza_pobreza, -cvlog) %>%
  select(where(es_continua))

vars_cuanti <- bind_cols(matriz_sae %>% select(pobreza_monetaria), covs_cont)

resumen_cuanti <- lapply(names(vars_cuanti), function(nm) {
  v <- vars_cuanti[[nm]]
  data.frame(variable = nm,
             media   = mean(v),   de      = sd(v),
             minimo  = min(v),    mediana = median(v),
             maximo  = max(v),    mc      = mc(v))
}) %>%
  bind_rows() %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

write.csv(resumen_cuanti, file.path(ruta_out, "descriptivo_cuantitativas.csv"),
          row.names = FALSE)

if (requireNamespace("knitr", quietly = TRUE)) {
  writeLines(
    knitr::kable(resumen_cuanti, format = "latex", booktabs = TRUE,
                 longtable = TRUE,
                 col.names = c("Variable", "Media", "DE", "Mín.",
                               "Mediana", "Máx.", "MC"),
                 caption = "Estadísticas descriptivas de las covariables candidatas",
                 label = "descriptivo-cuanti"),
    file.path(ruta_out, "descriptivo_cuantitativas.tex")
  )
}

frec_tabla <- function(v, nombre) {
  tb <- table(v, useNA = "ifany")
  data.frame(variable = nombre, categoria = names(tb),
             n = as.integer(tb),
             porcentaje = round(100 * as.integer(tb) / sum(tb), 1))
}

resumen_cuali <- bind_rows(
  frec_tabla(matriz_sae$cat_ruralidad, "Categoría de ruralidad"),
  frec_tabla(matriz_sae$categoria_617, "Categoría Ley 617"),
  frec_tabla(matriz_sae$pdet,          "Municipio PDET"),
  frec_tabla(matriz_sae$aislado_vial,  "Aislamiento vial")
)

cat("== COVARIABLES CATEGÓRICAS ==\n")
print(resumen_cuali)
cat("\n")

write.csv(resumen_cuali, file.path(ruta_out, "descriptivo_cualitativas.csv"),
          row.names = FALSE)

if (requireNamespace("knitr", quietly = TRUE)) {
  writeLines(
    knitr::kable(resumen_cuali, format = "latex", booktabs = TRUE,
                 col.names = c("Variable", "Categoría", "Municipios", "\\%"),
                 caption = "Distribución de las covariables categóricas e indicadoras",
                 label = "descriptivo-cuali", escape = FALSE),
    file.path(ruta_out, "descriptivo_cualitativas.tex")
  )
}

contraste_depto <- matriz_sae %>%
  mutate(depto = ifelse(substr(cod_mun, 1, 2) == "19",
                        "Cauca", "Valle del Cauca")) %>%
  group_by(depto) %>%
  summarise(
    municipios      = n(),
    pobreza_media   = round(mean(pobreza_monetaria), 3),
    pobreza_mediana = round(median(pobreza_monetaria), 3),
    pobreza_min     = round(min(pobreza_monetaria), 3),
    pobreza_max     = round(max(pobreza_monetaria), 3),
    pct_rural_medio = round(mean(pct_rural), 1),
    pdet_n          = sum(pdet),
    .groups = "drop"
  )

cat("== CONTRASTE DEPARTAMENTAL ==\n")
print(as.data.frame(contraste_depto))

cat("\nDominios con mayor y menor tasa de pobreza:\n")
print(matriz_sae %>%
        arrange(desc(pobreza_monetaria)) %>%
        select(Municipio, pobreza_monetaria) %>%
        slice(c(1:5, (n() - 4):n())) %>%
        as.data.frame())

write.csv(contraste_depto, file.path(ruta_out, "descriptivo_departamentos.csv"),
          row.names = FALSE)
cat("\n")


# ==============================================================================
# PARTE 10 — DIAGNÓSTICO DE ASIMETRÍA
# ------------------------------------------------------------------------------
# El medcouple mide la asimetría de forma robusta y toma valores en [-1, 1].
# Un valor igual o muy cercano a la unidad no refleja asimetría continua sino
# una distribución degenerada por concentración de ceros.
#
# El diagnóstico se reporta completo y se verifica que las variables tratadas
# en las secciones siguientes cumplan efectivamente los criterios sobre los
# datos vigentes, de modo que un cambio en las fuentes no pase inadvertido.
# ==============================================================================

diagnostico <- data.frame(
  variable  = names(covs_cont),
  mc        = sapply(covs_cont, function(x) mc(x, na.rm = TRUE)),
  pct_ceros = sapply(covs_cont,
                     function(x) 100 * sum(x == 0, na.rm = TRUE) / length(x)),
  minimo    = sapply(covs_cont, min, na.rm = TRUE),
  row.names = NULL
) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  arrange(desc(mc))

cat("== DIAGNÓSTICO DE ASIMETRÍA ==\n")
cat("Variables con medcouple igual o superior a 0,95:\n")
print(diagnostico %>% filter(mc >= 0.95), row.names = FALSE)

cat("\nVariables con asimetría moderada (entre 0,20 y 0,95):\n")
print(diagnostico %>% filter(mc > 0.20, mc < 0.95), row.names = FALSE)
cat("\n")

write.csv(diagnostico, file.path(ruta_out, "diagnostico_asimetria.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# PARTE 11 — TRATAMIENTO DE LAS VARIABLES CON EXCESO DE CEROS
# ==============================================================================

vars_exceso_ceros <- c("coca_pct_territorio", "mort_materna_2020", "mort_eda_2020",
                       "mort_ira_2020", "mort_desnutricion_2020",
                       "pagos_educacion_pc", "pagos_salud_pc",
                       "pagos_serv_publicos_pc", "empresas_formales_2016")

cat("== VERIFICACIÓN DEL GRUPO CON EXCESO DE CEROS ==\n")
print(diagnostico %>%
        filter(variable %in% vars_exceso_ceros) %>%
        select(variable, mc, pct_ceros) %>%
        arrange(desc(pct_ceros)),
      row.names = FALSE)

no_listadas <- setdiff(diagnostico$variable[diagnostico$mc >= 0.95],
                       vars_exceso_ceros)
if (length(no_listadas) > 0) {
  warning("Variables con medcouple mayor o igual a 0,95 fuera de la lista: ",
          paste(no_listadas, collapse = ", "))
}

g_ceros <- matriz_sae %>%
  select(all_of(vars_exceso_ceros)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "valor") %>%
  ggplot(aes(x = valor)) +
  geom_histogram(bins = 30, fill = COL_PRINCIPAL, color = "white") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(x = "Valor", y = "Número de municipios") +
  tema_tesis

ggsave(file.path(ruta_out, "fig_exceso_ceros.png"), g_ceros,
       width = 11, height = 8, dpi = 400)

# Medcouple previo de las variables que reciben transformación logarítmica,
# necesario para documentar el efecto de la transformación.
vars_log <- c("valor_agregado_pc", "avaluo_promedio", "ingresos_tributarios_pc",
              "densidad_2024", "pob_total", "ingresos_pc")

mc_log_antes <- sapply(matriz_sae[vars_log], function(x) mc(x, na.rm = TRUE))

# Recodificación: indicadoras de presencia y eliminación de las variables cuya
# proporción de ceros supera el 88 %.
matriz_trans <- matriz_sae %>%
  mutate(
    tiene_coca              = ifelse(coca_pct_territorio > 0, 1, 0),
    tiene_mort_materna      = ifelse(mort_materna_2020 > 0, 1, 0),
    tiene_mort_ira          = ifelse(mort_ira_2020 > 0, 1, 0),
    tiene_pagos_salud       = ifelse(pagos_salud_pc > 0, 1, 0),
    tiene_empresas_formales = ifelse(empresas_formales_2016 > 0, 1, 0)
  ) %>%
  select(
    -pagos_serv_publicos_pc,
    -pagos_educacion_pc,
    -mort_eda_2020,
    -mort_desnutricion_2020,
    -mort_materna_2020,
    -mort_ira_2020,
    -pagos_salud_pc,
    -empresas_formales_2016
  )

vars_dummy <- c("tiene_coca", "tiene_mort_materna", "tiene_mort_ira",
                "tiene_pagos_salud", "tiene_empresas_formales")

g_dummy <- matriz_trans %>%
  select(all_of(vars_dummy)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "valor") %>%
  mutate(valor = factor(valor, levels = c(0, 1),
                        labels = c("Ausencia", "Presencia"))) %>%
  count(variable, valor) %>%
  ggplot(aes(x = variable, y = n, fill = valor)) +
  geom_col(position = "stack") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            color = "white", size = 4.5) +
  scale_fill_manual(values = c("Ausencia" = COL_CONTRASTE,
                               "Presencia" = COL_PRINCIPAL)) +
  labs(x = NULL, y = "Número de municipios", fill = NULL) +
  coord_flip() +
  tema_tesis

ggsave(file.path(ruta_out, "fig_dummies_balance.png"), g_dummy,
       width = 10, height = 6, dpi = 400)


# ==============================================================================
# PARTE 12 — TRANSFORMACIÓN LOGARÍTMICA
# ------------------------------------------------------------------------------
# Se aplica a las variables con asimetría positiva moderada, de naturaleza
# monetaria o poblacional, en las que la concentración de valores altos en unos
# pocos municipios responde al fenómeno medido y no a un exceso de ceros.
# ==============================================================================

minimos <- sapply(matriz_trans[vars_log], min, na.rm = TRUE)

cat("\n== MÍNIMOS ANTES DE LA TRANSFORMACIÓN LOGARÍTMICA ==\n")
print(round(minimos, 4))

stopifnot("Alguna variable a transformar contiene valores no positivos" =
            all(minimos > 0))

matriz_trans <- matriz_trans %>%
  mutate(across(all_of(vars_log), ~ log(.), .names = "log_{.col}")) %>%
  select(-all_of(vars_log))

comparacion_log <- data.frame(
  variable   = vars_log,
  mc_antes   = round(as.numeric(mc_log_antes), 3),
  mc_despues = round(sapply(matriz_trans[paste0("log_", vars_log)],
                            function(x) mc(x, na.rm = TRUE)), 3),
  row.names  = NULL
)

cat("\n== MEDCOUPLE ANTES Y DESPUÉS DE LA TRANSFORMACIÓN ==\n")
print(comparacion_log, row.names = FALSE)

write.csv(comparacion_log, file.path(ruta_out, "comparacion_medcouple_log.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# PARTE 13 — CODIFICACIÓN DE VARIABLES CATEGÓRICAS
# ------------------------------------------------------------------------------
# La categoría de la Ley 617 es ordinal: establece una jerarquía legal de mayor
# a menor capacidad fiscal, por lo que se traduce a una escala numérica que
# preserva ese orden, asignando el valor cero a la categoría especial.
#
# La categoría de ruralidad es nominal y se recodifica mediante indicadoras,
# tomando como referencia la categoría más urbana, ciudades y aglomeraciones.
# ==============================================================================

matriz_trans <- matriz_trans %>%
  mutate(
    categoria_617 = as.numeric(ifelse(categoria_617 == "Especial", "0",
                                      as.character(categoria_617)))
  )

cat("\n== CODIFICACIÓN DE LA CATEGORÍA LEY 617 ==\n")
cat("Valores:", paste(sort(unique(matriz_trans$categoria_617)), collapse = ", "), "\n")

stopifnot("La codificación de categoria_617 generó valores ausentes" =
            sum(is.na(matriz_trans$categoria_617)) == 0)

matriz_trans <- matriz_trans %>%
  mutate(
    rural_intermedio = ifelse(cat_ruralidad == "Intermedio", 1, 0),
    rural            = ifelse(cat_ruralidad == "Rural", 1, 0),
    rural_disperso   = ifelse(cat_ruralidad == "Rural disperso", 1, 0)
  ) %>%
  select(-cat_ruralidad)


# ==============================================================================
# PARTE 14 — MATRIZ TRANSFORMADA
# ------------------------------------------------------------------------------
# Insumo de la Etapa II. Conserva la variable dependiente y su varianza junto
# con las covariables ya transformadas y codificadas.
# ==============================================================================

cat("\n== MATRIZ TRANSFORMADA ==\n")
cat("Dimensiones:", dim(matriz_trans)[1], "x", dim(matriz_trans)[2], "\n")
cat("Covariables:", ncol(matriz_trans) - length(cols_id_dep), "\n")
cat("Datos ausentes:", sum(is.na(matriz_trans)), "\n")

stopifnot("La matriz transformada contiene datos ausentes" =
            sum(is.na(matriz_trans)) == 0)

saveRDS(matriz_trans, file.path(ruta_out, "matriz_sae_transformada_v2.rds"))
write.csv(matriz_trans, file.path(ruta_out, "matriz_sae_transformada_v2.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\nEtapa I completada. Resultados en output/.\n")