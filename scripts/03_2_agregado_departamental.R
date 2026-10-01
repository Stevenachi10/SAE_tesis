# ==============================================================================
# 03_2_agregado_departamental.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA III.3 — Contraste externo frente al agregado departamental
#
# Objetivo
# --------
# Agregar las estimaciones municipales al nivel departamental, utilizando
# población municipal como ponderador, y contrastarlas con la estimación
# departamental oficial publicada para Cauca y Valle del Cauca.
#
# El contraste se realiza para:
#   1. estimador directo municipal;
#   2. modelos FH base;
#   3. modelos FH espaciales;
#   4. modelos robustos disponibles de la etapa anterior.
#
# La comparación con la cifra oficial es un contraste externo de consistencia.
# No se utiliza para establecer por sí sola la validez del modelo ni para
# seleccionar automáticamente una especificación.
#
# ------------------------------------------------------------------------------
# AGREGACIÓN
# ------------------------------------------------------------------------------
#
# Para cada departamento d se calcula:
#
#                 sum_i N_i * P_i
#   P_d = ------------------------------
#                sum_i N_i
#
# donde:
#   N_i = población del municipio i
#   P_i = incidencia de pobreza estimada para el municipio i.
#
# Los pesos normalizados son:
#
#                  N_i
#   w_i = -------------------------
#             sum_i N_i
#
# ------------------------------------------------------------------------------
# INCERTIDUMBRE DEL AGREGADO
# ------------------------------------------------------------------------------
#
# Para un estimador agregado:
#
#   P_d = sum_i w_i * P_i
#
# la incertidumbre completa del agregado depende de la covarianza conjunta
# entre las estimaciones municipales.
#
# Por ello, este script NO construye una varianza ni un intervalo del agregado
# del modelo a partir de la suma de los MSE municipales.
#
# El contraste externo utiliza el intervalo publicado por la fuente oficial
# como referencia descriptiva de compatibilidad.
#
# ------------------------------------------------------------------------------
# INTERPRETACIÓN
# ------------------------------------------------------------------------------
#
# dentro_ic_oficial:
#   TRUE  -> el agregado municipal de la especificación se encuentra dentro del
#            intervalo reportado para la estimación oficial.
#   FALSE -> el agregado queda fuera de dicho intervalo.
#
# Esto se interpreta como evidencia descriptiva de compatibilidad con la
# referencia oficial, no como una prueba formal de igualdad entre estimadores.
#
# brecha_en_ee:
#   distancia entre el agregado y la estimación oficial medida en unidades del
#   error estándar reportado por la fuente oficial. Es una medida descriptiva;
#   no constituye un estadístico de prueba porque no incorpora la incertidumbre
#   completa del agregado del modelo ni su covarianza con la referencia oficial.
#
# factor_benchmark:
#   razón oficial/agregado. Se reporta únicamente con carácter informativo y no
#   se utiliza como criterio de selección.
#
# ------------------------------------------------------------------------------
# ALCANCE DEL RESUMEN DE COMPATIBILIDAD
# ------------------------------------------------------------------------------
#
# El resumen del bloque 12 se calcula únicamente sobre las especificaciones
# candidatas. Quedan fuera:
#
#   - el estimador directo, que es el insumo del modelo y no una especificación
#     que compita con las demás;
#
#   - M1 en escala original, que se conserva como referencia de contraste y no
#     como candidata.
#
# Ambos siguen apareciendo en la tabla general por departamento.
#
# ------------------------------------------------------------------------------
# ENTRADAS
# ------------------------------------------------------------------------------
#   data/raw/Agregados.xlsx
#       Cifras departamentales oficiales.
#
#   output/matriz_sae_transformada_v2.rds
#       Matriz de insumos municipales.
#
#   output/modelado_variantes.rds
#       Estimador directo y modelos FH base.
#
#   output/espacial.rds
#       Modelos FH espaciales.
#
#   output/robusto.rds
#       Modelos robustos.
#
# ------------------------------------------------------------------------------
# SALIDAS
# ------------------------------------------------------------------------------
#   output/agr_referencia_departamental.csv
#       Referencia oficial de los departamentos estudiados.
#
#   output/agr_comparacion.csv
#       Agregados municipales, brechas y contraste descriptivo frente a la
#       referencia oficial.
#
#   output/agr_resumen_compatibilidad.csv
#       Recuento de especificaciones candidatas compatibles por departamento.
# ==============================================================================


# ==============================================================================
# 0. PAQUETES Y PARÁMETROS
# ==============================================================================

library(readxl)
library(dplyr)
library(here)

# Evita conflictos con funciones de otros paquetes
select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")
ruta_agr <- here("data", "raw", "Agregados.xlsx")


# ------------------------------------------------------------------------------
# Departamentos del estudio
# ------------------------------------------------------------------------------

DEPARTAMENTOS <- data.frame(
  cod_dpto = c("19", "76"),
  nombre_fuente = c("cauca", "valle del cauca"),
  departamento = c("Cauca", "Valle del Cauca"),
  municipios_esperados = c(42L, 42L),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# Ponderador de población
# ------------------------------------------------------------------------------

COL_PESO <- "log_pob_total"

PESO_EN_LOGARITMO <- TRUE


# ------------------------------------------------------------------------------
# Especificaciones espaciales que avanzan de la etapa anterior
# ------------------------------------------------------------------------------

ESPACIALES <- c("M1", "M2")


# ------------------------------------------------------------------------------
# Constante de sintonización del estimador robusto
# ------------------------------------------------------------------------------
# Se declara como cadena para que la clave de búsqueda coincida exactamente
# con la utilizada al guardar los ajustes robustos.
# ------------------------------------------------------------------------------

K_ROBUSTO <- "1.345"


# ==============================================================================
# 1. FUNCIONES AUXILIARES
# ==============================================================================


# ------------------------------------------------------------------------------
# Normalización de nombres
# ------------------------------------------------------------------------------

normalizar <- function(x) {
  
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- chartr("áéíóúüñ", "aeiouun", x)
  x <- gsub("[[:space:]]+", " ", x)
  
  x
}


# ------------------------------------------------------------------------------
# Agregación departamental
# ------------------------------------------------------------------------------
#
# El agregado departamental es:
#
#       sum_i N_i * P_i
#   -----------------------
#          sum_i N_i
#
# Se utiliza weighted.mean(), que implementa exactamente esta expresión.
#
# No se construye aquí una varianza del agregado porque la incertidumbre
# completa requiere considerar la covarianza conjunta entre estimaciones
# municipales.
# ------------------------------------------------------------------------------

agregar_por_dpto <- function(datos, est) {
  
  # .env$est fuerza a tomar el vector recibido como argumento y no una
  # columna de 'datos' que pudiera llamarse igual.
  d <- datos %>%
    mutate(
      est = as.numeric(.env$est)
    )
  
  d %>%
    group_by(cod_dpto) %>%
    summarise(
      
      agregado = weighted.mean(
        x = est,
        w = peso,
        na.rm = FALSE
      ),
      
      .groups = "drop"
      
    ) %>%
    as.data.frame()
}


# ==============================================================================
# 2. CIFRAS DEPARTAMENTALES OFICIALES
# ==============================================================================


# ------------------------------------------------------------------------------
# Comprobar existencia del archivo oficial
# ------------------------------------------------------------------------------

stopifnot(
  "No se encuentra Agregados.xlsx" =
    file.exists(ruta_agr)
)


# ------------------------------------------------------------------------------
# Lectura de la fuente oficial
# ------------------------------------------------------------------------------

crudo <- read_excel(
  ruta_agr,
  sheet = 1,
  col_names = TRUE,
  .name_repair = "minimal"
)


cat(
  "\n============================================================\n"
)

cat(
  "FUENTE DE LAS CIFRAS DEPARTAMENTALES\n"
)

cat(
  "============================================================\n"
)


cat(
  "Filas:",
  nrow(crudo),
  "| Columnas:",
  ncol(crudo),
  "\n"
)


cat(
  "Encabezado de origen:\n"
)

print(
  names(crudo)
)


# ------------------------------------------------------------------------------
# Validación estructural
# ------------------------------------------------------------------------------

stopifnot(
  "La hoja no tiene las seis columnas esperadas" =
    ncol(crudo) == 6
)


# ------------------------------------------------------------------------------
# Estandarización de nombres
# ------------------------------------------------------------------------------

names(crudo) <- c(
  "departamento_fuente",
  "estimacion",
  "error_estandar",
  "limite_inferior",
  "limite_superior",
  "cv_pct"
)


# ------------------------------------------------------------------------------
# Conversión numérica
# ------------------------------------------------------------------------------

numericas <- c(
  "estimacion",
  "error_estandar",
  "limite_inferior",
  "limite_superior",
  "cv_pct"
)


for (v in numericas) {
  
  crudo[[v]] <- suppressWarnings(
    as.numeric(crudo[[v]])
  )
  
}


# ------------------------------------------------------------------------------
# Normalización de nombres departamentales
# ------------------------------------------------------------------------------

crudo$clave <- normalizar(
  crudo$departamento_fuente
)


# ------------------------------------------------------------------------------
# Comprobar que cada departamento del estudio aparece una sola vez
# ------------------------------------------------------------------------------

coincidencias <- crudo %>%
  filter(
    clave %in% DEPARTAMENTOS$nombre_fuente
  )


stopifnot(
  "Hay departamentos del estudio repetidos en la fuente oficial" =
    !anyDuplicated(coincidencias$clave)
)


# ------------------------------------------------------------------------------
# Extracción de los dos departamentos
# ------------------------------------------------------------------------------

referencia <- DEPARTAMENTOS %>%
  left_join(
    crudo %>%
      select(
        clave,
        all_of(numericas)
      ),
    by = c(
      "nombre_fuente" = "clave"
    )
  )


# ------------------------------------------------------------------------------
# Comprobar faltantes
# ------------------------------------------------------------------------------

faltantes <- referencia$departamento[
  is.na(referencia$estimacion)
]


if (length(faltantes) > 0) {
  
  cat(
    "\nNo se encontró la cifra de:",
    paste(faltantes, collapse = ", "),
    "\n"
  )
  
  cat(
    "\nDepartamentos disponibles en la fuente:\n"
  )
  
  print(
    sort(crudo$departamento_fuente)
  )
  
  stop(
    "Falta al menos un departamento del estudio en Agregados.xlsx."
  )
  
}


# ------------------------------------------------------------------------------
# Validación de valores oficiales
# ------------------------------------------------------------------------------

stopifnot(
  
  "Estimación oficial no válida" =
    all(
      is.finite(
        referencia$estimacion
      )
    ),
  
  "Error estándar oficial no válido" =
    all(
      is.finite(
        referencia$error_estandar
      )
    ),
  
  "Límite inferior no válido" =
    all(
      is.finite(
        referencia$limite_inferior
      )
    ),
  
  "Límite superior no válido" =
    all(
      is.finite(
        referencia$limite_superior
      )
    )
  
)


# ------------------------------------------------------------------------------
# Conversión de puntos porcentuales a proporciones
# ------------------------------------------------------------------------------

referencia <- referencia %>%
  mutate(
    
    estimacion =
      estimacion / 100,
    
    error_estandar =
      error_estandar / 100,
    
    limite_inferior =
      limite_inferior / 100,
    
    limite_superior =
      limite_superior / 100
    
  )


cat(
  "\n-- Referencia oficial, en proporción --\n\n"
)


print(
  
  referencia %>%
    
    select(
      cod_dpto,
      departamento,
      estimacion,
      error_estandar,
      limite_inferior,
      limite_superior,
      cv_pct
    ),
  
  row.names = FALSE
  
)


# ==============================================================================
# 3. INSUMOS MUNICIPALES
# ==============================================================================


matriz <- readRDS(
  file.path(
    ruta_out,
    "matriz_sae_transformada_v2.rds"
  )
)


variantes <- readRDS(
  file.path(
    ruta_out,
    "modelado_variantes.rds"
  )
)


espacial <- readRDS(
  file.path(
    ruta_out,
    "espacial.rds"
  )
)


robusto <- readRDS(
  file.path(
    ruta_out,
    "robusto.rds"
  )
)


# ------------------------------------------------------------------------------
# Datos municipales
# ------------------------------------------------------------------------------

datos <- variantes$datos

n_dom <- nrow(datos)

DIRECTO <- datos$pobreza_monetaria


# ==============================================================================
# 4. PONDERADOR POBLACIONAL
# ==============================================================================


# ------------------------------------------------------------------------------
# Comprobar existencia de la columna
# ------------------------------------------------------------------------------

if (!COL_PESO %in% names(matriz)) {
  
  cat(
    "\nNo está la columna",
    COL_PESO,
    "en la matriz de insumos.\n"
  )
  
  cat(
    "\nColumnas numéricas disponibles:\n"
  )
  
  print(
    names(matriz)[
      vapply(
        matriz,
        is.numeric,
        logical(1)
      )
    ]
  )
  
  stop(
    "Indicar la columna de población en COL_PESO."
  )
  
}


cat(
  "\nPonderador:",
  COL_PESO,
  if (PESO_EN_LOGARITMO)
    "(se exponencia)"
  else
    "",
  "\n"
)


# ------------------------------------------------------------------------------
# Extraer población
# ------------------------------------------------------------------------------

pesos <- matriz %>%
  
  as.data.frame() %>%
  
  select(
    cod_mun,
    peso = all_of(COL_PESO)
  ) %>%
  
  mutate(
    
    cod_mun =
      as.character(cod_mun),
    
    peso =
      if (PESO_EN_LOGARITMO)
        exp(peso)
    else
      peso
    
  )


# ------------------------------------------------------------------------------
# Validación de identificadores y población
# ------------------------------------------------------------------------------

stopifnot(
  
  "Hay cod_mun duplicados en la matriz de ponderadores" =
    !anyDuplicated(
      pesos$cod_mun
    ),
  
  "Hay poblaciones no finitas" =
    all(
      is.finite(
        pesos$peso
      )
    ),
  
  "Hay poblaciones no positivas" =
    all(
      pesos$peso > 0
    )
  
)


# ------------------------------------------------------------------------------
# Integración con los datos de modelación
# ------------------------------------------------------------------------------

datos <- datos %>%
  
  mutate(
    cod_mun =
      as.character(cod_mun)
  ) %>%
  
  left_join(
    pesos,
    by = "cod_mun"
  ) %>%
  
  mutate(
    cod_dpto =
      substr(
        cod_mun,
        1,
        2
      )
  )


# ------------------------------------------------------------------------------
# Validaciones posteriores
# ------------------------------------------------------------------------------

stopifnot(
  
  "Hay municipios sin ponderador" =
    all(
      is.finite(
        datos$peso
      )
    ),
  
  "Hay ponderadores no positivos" =
    all(
      datos$peso > 0
    ),
  
  "Hay municipios fuera de los departamentos del estudio" =
    all(
      datos$cod_dpto %in%
        DEPARTAMENTOS$cod_dpto
    ),
  
  "Hay identificadores municipales duplicados" =
    !anyDuplicated(
      datos$cod_mun
    )
  
)


# ==============================================================================
# 5. VERIFICACIÓN DEL PONDERADOR
# ==============================================================================


cat(
  "\n-- Cobertura municipal y población recuperada --\n\n"
)


cobertura <- datos %>%
  
  group_by(
    cod_dpto
  ) %>%
  
  summarise(
    
    municipios =
      n(),
    
    poblacion =
      round(
        sum(peso)
      ),
    
    .groups = "drop"
    
  ) %>%
  
  left_join(
    
    DEPARTAMENTOS %>%
      
      select(
        cod_dpto,
        departamento
      ),
    
    by = "cod_dpto"
    
  ) %>%
  
  select(
    cod_dpto,
    departamento,
    municipios,
    poblacion
  )


print(
  cobertura,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Verificar que ambos departamentos están representados
# ------------------------------------------------------------------------------

stopifnot(
  
  "Falta un departamento del estudio en los datos municipales" =
    setequal(
      cobertura$cod_dpto,
      DEPARTAMENTOS$cod_dpto
    )
  
)


# ------------------------------------------------------------------------------
# Verificar cobertura completa de cada departamento
# ------------------------------------------------------------------------------
# El contraste con la cifra departamental solo es válido si el agregado
# incluye todos los municipios del departamento.
# ------------------------------------------------------------------------------

cobertura_completa <- cobertura %>%
  left_join(
    DEPARTAMENTOS %>% select(cod_dpto, municipios_esperados),
    by = "cod_dpto"
  )

stopifnot(
  "Algún departamento no tiene todos sus municipios" =
    all(cobertura_completa$municipios == cobertura_completa$municipios_esperados)
)


# ------------------------------------------------------------------------------
# Municipios extremos por población
# ------------------------------------------------------------------------------

extremos <- datos %>%
  
  select(
    cod_mun,
    Municipio,
    peso
  ) %>%
  
  mutate(
    poblacion =
      round(peso)
  ) %>%
  
  arrange(
    desc(poblacion)
  ) %>%
  
  select(
    cod_mun,
    Municipio,
    poblacion
  )


cat(
  "\nMunicipios de mayor y menor población recuperada:\n\n"
)


print(
  
  rbind(
    head(extremos, 5),
    tail(extremos, 5)
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------------------------
# Control de valores extremos
# ------------------------------------------------------------------------------

implausibles <-
  sum(
    datos$peso < 100 |
      datos$peso > 5e6
  )


if (implausibles > 0) {
  
  cat(
    "\nAVISO:",
    implausibles,
    "municipios con población fuera de un rango amplio de control.\n"
  )
  
  cat(
    "Verificar si",
    COL_PESO,
    "está estandarizada además de transformada.\n"
  )
  
}


# ==============================================================================
# 6. ESPECIFICACIONES A AGREGAR
# ==============================================================================


AJUSTES <- list()


# ------------------------------------------------------------------------------
# Función para registrar modelos
# ------------------------------------------------------------------------------

agregar_ajuste <- function(
    etiqueta,
    objeto,
    familia,
    papel
) {
  
  if (is.null(objeto)) {
    
    cat(
      "AUSENTE:",
      etiqueta,
      "\n"
    )
    
    return(
      invisible(NULL)
    )
    
  }
  
  
  AJUSTES[[etiqueta]] <<-
    list(
      
      aj =
        objeto,
      
      familia =
        familia,
      
      papel =
        papel
      
    )
  
}


# ------------------------------------------------------------------------------
# Modelos FH base
# ------------------------------------------------------------------------------

agregar_ajuste(
  
  "M1 original",
  
  variantes$ajustes[["M1_reml_original"]],
  
  "base",
  
  "referencia"
  
)


agregar_ajuste(
  
  "M1 log",
  
  variantes$ajustes[["M1_reml_log"]],
  
  "base",
  
  "candidata"
  
)


agregar_ajuste(
  
  "M1 logit",
  
  variantes$ajustes[["M1_reml_logit"]],
  
  "base",
  
  "candidata"
  
)


# ------------------------------------------------------------------------------
# Modelos espaciales
# ------------------------------------------------------------------------------

for (nm in ESPACIALES) {
  
  agregar_ajuste(
    
    paste(
      "Espacial",
      nm
    ),
    
    espacial$ajustes_esp[[nm]],
    
    "espacial",
    
    "candidata"
    
  )
  
}


# ------------------------------------------------------------------------------
# Modelos robustos
# ------------------------------------------------------------------------------

for (nm in c("M1", "M2", "M3")) {
  
  for (met in c("reblup", "reblupbc")) {
    
    clave <-
      paste(
        nm,
        met,
        K_ROBUSTO,
        sep = "_"
      )
    
    
    agregar_ajuste(
      
      paste(
        met,
        nm
      ),
      
      robusto$ajustes_robustos[[clave]],
      
      "robusta",
      
      "candidata"
      
    )
    
  }
  
}


cat(
  "\nEspecificaciones a agregar:",
  length(AJUSTES),
  "\n"
)


# ==============================================================================
# 7. EXTRACCIÓN DE ESTIMACIONES Y MSE
# ==============================================================================


extraer_ajuste <- function(aj) {
  
  est <-
    tryCatch(
      as.numeric(
        aj$ind$FH
      ),
      error =
        function(e)
          NULL
    )
  
  
  mse <-
    tryCatch(
      as.numeric(
        aj$MSE$FH
      ),
      error =
        function(e)
          NULL
    )
  
  
  # ---------------------------------------------------------------------------
  # Validación de existencia
  # ---------------------------------------------------------------------------
  
  if (
    is.null(est) ||
    is.null(mse)
  ) {
    
    return(
      NULL
    )
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Validación de longitud
  # ---------------------------------------------------------------------------
  
  if (
    length(est) != n_dom ||
    length(mse) != n_dom
  ) {
    
    return(
      NULL
    )
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Identificadores de dominio
  # ---------------------------------------------------------------------------
  
  dom <-
    tryCatch(
      as.character(
        aj$ind$Domain
      ),
      error =
        function(e)
          NULL
    )
  
  
  if (
    is.null(dom) ||
    !identical(
      dom,
      datos$cod_mun
    )
  ) {
    
    stop(
      "Los dominios del ajuste no están alineados con 'datos'."
    )
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Validación numérica
  # ---------------------------------------------------------------------------
  
  if (
    any(!is.finite(est)) ||
    any(!is.finite(mse)) ||
    any(mse < 0)
  ) {
    
    return(
      NULL
    )
    
  }
  
  
  list(
    est = est,
    mse = mse
  )
  
}


# ==============================================================================
# 8. AGREGACIÓN
# ==============================================================================


filas <- list()


# ------------------------------------------------------------------------------
# 8.1 Estimador directo
# ------------------------------------------------------------------------------

ag <- agregar_por_dpto(
  datos = datos,
  est = DIRECTO
)


filas[[1]] <-
  ag %>%
  
  mutate(
    
    ajuste =
      "Estimador directo",
    
    familia =
      "directo",
    
    papel =
      "insumo"
    
  )


# ------------------------------------------------------------------------------
# 8.2 Especificaciones de modelos
# ------------------------------------------------------------------------------

for (nm in names(AJUSTES)) {
  
  ex <-
    extraer_ajuste(
      AJUSTES[[nm]]$aj
    )
  
  
  if (is.null(ex)) {
    
    cat(
      "NO EVALUABLE:",
      nm,
      "\n"
    )
    
    next
    
  }
  
  
  ag <-
    agregar_por_dpto(
      
      datos =
        datos,
      
      est =
        ex$est
      
    )
  
  
  filas[[length(filas) + 1]] <-
    ag %>%
    
    mutate(
      
      ajuste =
        nm,
      
      familia =
        AJUSTES[[nm]]$familia,
      
      papel =
        AJUSTES[[nm]]$papel
      
    )
  
}


comparacion <-
  do.call(
    rbind,
    filas
  )


# ==============================================================================
# 9. CONTRASTE CON LA REFERENCIA OFICIAL
# ==============================================================================


comparacion <-
  comparacion %>%
  
  left_join(
    
    referencia %>%
      
      select(
        
        cod_dpto,
        
        departamento,
        
        oficial =
          estimacion,
        
        oficial_ee =
          error_estandar,
        
        oficial_li =
          limite_inferior,
        
        oficial_ls =
          limite_superior
        
      ),
    
    by = "cod_dpto"
    
  ) %>%
  
  mutate(
    
    # --------------------------------------------------------------------------
    # Diferencia en escala de proporción
    # --------------------------------------------------------------------------
    
    brecha =
      agregado - oficial,
    
    
    # --------------------------------------------------------------------------
    # Diferencia en puntos porcentuales
    # --------------------------------------------------------------------------
    
    brecha_pp =
      round(
        (agregado - oficial) * 100,
        3
      ),
    
    
    # --------------------------------------------------------------------------
    # Diferencia relativa respecto a la referencia
    # --------------------------------------------------------------------------
    
    brecha_rel =
      round(
        (agregado - oficial) /
          oficial * 100,
        2
      ),
    
    
    # --------------------------------------------------------------------------
    # Distancia expresada en errores estándar de la referencia
    #
    # Medida descriptiva; no es una prueba formal.
    # --------------------------------------------------------------------------
    
    brecha_en_ee =
      round(
        (agregado - oficial) /
          oficial_ee,
        2
      ),
    
    
    # --------------------------------------------------------------------------
    # Compatibilidad con el intervalo reportado oficialmente
    #
    # No se interpreta como prueba de igualdad.
    # --------------------------------------------------------------------------
    
    dentro_ic_oficial =
      agregado >= oficial_li &
      agregado <= oficial_ls,
    
    
    # --------------------------------------------------------------------------
    # Factor informativo de benchmarking
    # --------------------------------------------------------------------------
    
    factor_benchmark =
      round(
        oficial / agregado,
        4
      )
    
  ) %>%
  
  select(
    
    departamento,
    
    ajuste,
    
    familia,
    
    papel,
    
    agregado,
    
    oficial,
    
    oficial_li,
    
    oficial_ls,
    
    brecha_pp,
    
    brecha_rel,
    
    brecha_en_ee,
    
    dentro_ic_oficial,
    
    factor_benchmark
    
  ) %>%
  
  arrange(
    
    departamento,
    
    papel,
    
    ajuste
    
  )


# ==============================================================================
# 10. RESULTADOS DEL CONTRASTE
# ==============================================================================


cat(
  "\n============================================================\n"
)

cat(
  "AGREGADO DEPARTAMENTAL FRENTE A LA CIFRA OFICIAL\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nbrecha_pp     : diferencia en puntos porcentuales\n"
)

cat(
  "brecha_rel    : diferencia relativa respecto a la referencia\n"
)

cat(
  "brecha_en_ee  : distancia respecto al EE oficial\n"
)

cat(
  "dentro_ic_oficial : el agregado cae dentro del intervalo oficial\n"
)

cat(
  "\nEstas medidas son descriptivas y no constituyen una prueba formal\n"
)

cat(
  "de igualdad entre el agregado del modelo y la referencia oficial.\n\n"
)


# ------------------------------------------------------------------------------
# Resultados por departamento
# ------------------------------------------------------------------------------

for (dp in referencia$departamento) {
  
  cat(
    "--",
    dp,
    "--\n\n"
  )
  
  
  print(
    
    comparacion %>%
      
      filter(
        departamento == dp
      ) %>%
      
      select(
        
        ajuste,
        
        familia,
        
        papel,
        
        agregado,
        
        oficial,
        
        brecha_pp,
        
        brecha_rel,
        
        brecha_en_ee,
        
        dentro_ic_oficial
        
      ),
    
    row.names = FALSE
    
  )
  
  
  cat(
    "\n"
  )
  
}


# ==============================================================================
# 11. COMPROBACIÓN DEL ESTIMADOR DIRECTO
# ==============================================================================
#
# Esta sección evalúa si el agregado de los estimadores directos municipales
# es consistente con la referencia departamental oficial.
#
# La interpretación correcta es "evidencia de consistencia", no una demostración
# de que ambas cifras hayan sido obtenidas exactamente mediante idéntico
# procedimiento.
# ==============================================================================


directo_cmp <-
  comparacion %>%
  
  filter(
    ajuste ==
      "Estimador directo"
  )


cat(
  "-- Agregado del estimador directo --\n"
)

cat(
  "Se verifica la consistencia entre la agregación de los estimadores\n"
)

cat(
  "municipales y la referencia departamental oficial.\n\n"
)


print(
  
  directo_cmp %>%
    
    select(
      
      departamento,
      
      agregado,
      
      oficial,
      
      brecha_pp,
      
      brecha_rel,
      
      brecha_en_ee,
      
      dentro_ic_oficial
      
    ),
  
  row.names = FALSE
  
)


# ==============================================================================
# 12. RESUMEN DE COMPATIBILIDAD
# ==============================================================================
#
# Se calcula únicamente sobre las especificaciones candidatas.
#
# No se incluyen:
#   - estimador directo;
#   - M1 original.
#
# Esto evita interpretar el insumo o la especificación de referencia como
# modelos competidores dentro del conjunto de candidatas.
# ==============================================================================


candidatas <-
  comparacion %>%
  
  filter(
    papel == "candidata"
  )


resumen_compatibilidad <-
  candidatas %>%
  
  group_by(
    departamento
  ) %>%
  
  summarise(
    
    especificaciones =
      n(),
    
    dentro_intervalo =
      sum(
        dentro_ic_oficial,
        na.rm = TRUE
      ),
    
    proporcion_dentro =
      round(
        mean(
          dentro_ic_oficial,
          na.rm = TRUE
        ),
        3
      ),
    
    .groups = "drop"
    
  )


cat(
  "\n============================================================\n"
)

cat(
  "RESUMEN DE COMPATIBILIDAD CON LA REFERENCIA OFICIAL\n"
)

cat(
  "Solo especificaciones candidatas. Excluye el estimador directo\n"
)

cat(
  "y M1 en escala original.\n"
)

cat(
  "============================================================\n\n"
)


print(
  resumen_compatibilidad,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Candidatas fuera del intervalo oficial
# ------------------------------------------------------------------------------

fuera <-
  candidatas %>%
  
  filter(
    !dentro_ic_oficial
  ) %>%
  
  select(
    
    departamento,
    
    ajuste,
    
    familia,
    
    agregado,
    
    oficial,
    
    brecha_pp,
    
    brecha_rel,
    
    brecha_en_ee
    
  )


if (nrow(fuera) > 0) {
  
  cat(
    "\n-- Candidatas fuera del intervalo oficial --\n\n"
  )
  
  print(
    fuera,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "\nTodas las candidatas caen dentro del intervalo oficial en ambos\n"
  )
  
  cat(
    "departamentos.\n"
  )
  
}


# ==============================================================================
# 13. EXPORTACIÓN
# ==============================================================================


write.csv(
  
  referencia,
  
  file.path(
    ruta_out,
    "agr_referencia_departamental.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


write.csv(
  
  comparacion,
  
  file.path(
    ruta_out,
    "agr_comparacion.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


write.csv(
  
  resumen_compatibilidad,
  
  file.path(
    ruta_out,
    "agr_resumen_compatibilidad.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


# ==============================================================================
# 14. FINAL
# ==============================================================================


cat(
  "\n============================================================\n"
)

cat(
  "CONTRASTE EXTERNO COMPLETADO\n"
)

cat(
  "============================================================\n"
)