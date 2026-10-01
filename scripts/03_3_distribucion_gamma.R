# ==============================================================================
# 03_3_distribucion_gamma.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA III.4 — Distribución del factor de encogimiento
#
# Objetivo
# --------
# Describir la distribución del factor de encogimiento en los 84 dominios
# para comparar la intensidad y heterogeneidad del uso de la información
# directa y sintética entre especificaciones.
#
# El promedio de gamma no es suficiente para caracterizar el comportamiento
# del encogimiento. Por ello se consideran medidas de posición, dispersión,
# extremos y recuentos de dominios bajo distintos niveles descriptivos de
# encogimiento.
#
# ----------------------------------------------------------------------
# INTERPRETACIÓN
# ----------------------------------------------------------------------
#
# En el FH estándar:
#
#   gamma_i = sigma_u^2 / (sigma_u^2 + psi_i)
#
# representa el peso del estimador directo del propio dominio.
#
# Por tanto:
#
#   1 - gamma_i
#
# representa el peso sintético.
#
# En los modelos robustos, el valor calculado para gamma se obtiene a partir
# de la estructura de varianza del modelo de trabajo y se utiliza solamente
# como medida aproximada de intensidad de encogimiento; no debe interpretarse
# como el peso exacto aplicado por el predictor robusto.
#
# En los modelos espaciales, la matriz de encogimiento no es diagonal.
# La diagonal representa únicamente el peso asociado al estimador directo
# del propio dominio. Para evaluar el peso total sobre los estimadores
# directos se reporta adicionalmente la suma de filas de la matriz.
#
# Por esa razón el peso sintético medio, 1 - media de gamma, NO se reporta en
# la tabla principal. Sobre la diagonal del modelo espacial esa cantidad no es
# el peso sintético: el predictor toma además el aporte de los municipios
# vecinos, de modo que el complemento de la diagonal lo sobrestima. El peso
# sintético se reporta únicamente en la tabla de suma de filas, donde el vector
# sí agota el peso sobre los estimadores directos.
#
# Los umbrales gamma < 0.50, gamma < 0.25 y gamma < 0.10 se emplean
# exclusivamente como referencias descriptivas de intensidad creciente
# del encogimiento. No constituyen criterios de aceptación o rechazo.
#
# ----------------------------------------------------------------------
# ENTRADAS
# ----------------------------------------------------------------------
# output/prec_cv_por_dominio.csv
# output/espacial.rds
#
# ----------------------------------------------------------------------
# SALIDAS
# ----------------------------------------------------------------------
# output/gam_distribucion.csv
# output/gam_distribucion.tex
# output/gam_espacial_filas.csv
# ==============================================================================


# ==============================================================================
# 0. PAQUETES Y PARÁMETROS
# ==============================================================================

library(dplyr)
library(here)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")


# ------------------------------------------------------------------------------
# Umbrales descriptivos
# ------------------------------------------------------------------------------

UMBRALES <- c(
  0.50,
  0.25,
  0.10
)


# ------------------------------------------------------------------------------
# Orden de presentación
# ------------------------------------------------------------------------------

ORDEN <- c(
  "M1 log",
  "M1 logit",
  "Espacial M1",
  "Espacial M2",
  "reblup M1",
  "reblup M2",
  "reblup M3",
  "reblupbc M1",
  "reblupbc M2",
  "reblupbc M3",
  "M1 original"
)


# ------------------------------------------------------------------------------
# Especificaciones espaciales
# ------------------------------------------------------------------------------

ESPACIALES <- c(
  "M1",
  "M2"
)


# ==============================================================================
# 1. CARGA DE LOS FACTORES DE ENCOGIMIENTO
# ==============================================================================

ruta_cv <- file.path(
  ruta_out,
  "prec_cv_por_dominio.csv"
)

stopifnot(
  "No se encuentra prec_cv_por_dominio.csv" =
    file.exists(ruta_cv)
)

cv <- read.csv(
  ruta_cv,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8",
  check.names = TRUE
)


# ------------------------------------------------------------------------------
# Identificación explícita de las columnas de gamma
# ------------------------------------------------------------------------------
# Se evita identificar gamma simplemente por descarte, para que futuras
# columnas auxiliares del archivo no sean interpretadas erróneamente como
# factores de encogimiento.
# ------------------------------------------------------------------------------

cols_gamma_esperadas <- make.names(
  ORDEN
)

cols_gamma <- intersect(
  cols_gamma_esperadas,
  names(cv)
)

stopifnot(
  "No se identificaron columnas de gamma" =
    length(cols_gamma) > 0
)

n_dom <- nrow(cv)


# ------------------------------------------------------------------------------
# Red de seguridad
# ------------------------------------------------------------------------------
# ORDEN cumple dos funciones a la vez: ordena y selecciona. Una especificación
# que se agregue al pipeline y no se añada a ORDEN desaparecería de la tabla sin
# producir ningún aviso. Se comparan por tanto las columnas seleccionadas contra
# las que quedan al descartar identificadores y coeficientes de variación.
# ------------------------------------------------------------------------------

cols_cv <- grep(
  "^CV_",
  names(cv),
  value = TRUE
)

omitidas <- setdiff(
  setdiff(
    names(cv),
    c("cod_mun", "Municipio", cols_cv)
  ),
  cols_gamma
)

if (length(omitidas) > 0) {
  
  cat(
    "\nAVISO: columnas de gamma en el archivo y ausentes de ORDEN:\n  ",
    paste(omitidas, collapse = ", "),
    "\n"
  )
  
  cat(
    "No se incluyen en la tabla. Agregarlas a ORDEN si corresponde.\n"
  )
  
}


cat("\n============================================================\n")
cat("FACTORES DE ENCOGIMIENTO\n")
cat("============================================================\n")

cat(
  "Dominios:",
  n_dom,
  "\n"
)

cat(
  "Especificaciones:",
  length(cols_gamma),
  "\n"
)


# ------------------------------------------------------------------------------
# Verificación básica
# ------------------------------------------------------------------------------

for (cl in cols_gamma) {
  
  if (all(!is.finite(cv[[cl]]))) {
    
    cat(
      "\nSin gamma disponible:",
      cl,
      "\n"
    )
    
  }
  
}


# ==============================================================================
# 2. FUNCIÓN DE RESUMEN
# ==============================================================================

resumir_gamma <- function(
    g,
    etiqueta,
    incluir_peso_sintetico = FALSE
) {
  
  g <- g[
    is.finite(g)
  ]
  
  if (length(g) == 0) {
    return(NULL)
  }
  
  q <- quantile(
    g,
    c(
      0.10,
      0.25,
      0.75,
      0.90
    ),
    names = FALSE
  )
  
  fila <- data.frame(
    
    modelo = etiqueta,
    
    n = length(g),
    
    media = mean(g),
    
    mediana = median(g),
    
    p10 = q[1],
    
    p25 = q[2],
    
    p75 = q[3],
    
    p90 = q[4],
    
    minimo = min(g),
    
    maximo = max(g),
    
    RIC = q[3] - q[2],
    
    stringsAsFactors = FALSE
  )
  
  
  # --------------------------------------------------------------------------
  # Peso sintético medio
  # --------------------------------------------------------------------------
  # Solo se calcula cuando el vector recibido agota el peso sobre los
  # estimadores directos. Eso ocurre en los modelos de matriz de encogimiento
  # diagonal y en la suma de filas del modelo espacial, pero no en la diagonal
  # del modelo espacial, donde 1 - gamma_i sobrestima el componente sintético.
  # --------------------------------------------------------------------------
  
  if (incluir_peso_sintetico) {
    
    fila$peso_sintetico_medio <- 1 - mean(g)
    
  }
  
  
  # --------------------------------------------------------------------------
  # Recuentos descriptivos
  # --------------------------------------------------------------------------
  
  for (u in UMBRALES) {
    
    nombre <- paste0(
      "n_gamma_menor_",
      format(
        u,
        trim = TRUE
      )
    )
    
    fila[[nombre]] <- sum(
      g < u
    )
    
  }
  
  
  fila
}


# ==============================================================================
# 3. TABLA PRINCIPAL
# ==============================================================================

tabla_gamma <- do.call(
  rbind,
  lapply(
    cols_gamma,
    function(cl) {
      
      resumir_gamma(
        cv[[cl]],
        gsub(
          "\\.",
          " ",
          cl
        )
      )
      
    }
  )
)


# ------------------------------------------------------------------------------
# Orden de presentación
# ------------------------------------------------------------------------------

tabla_gamma <- tabla_gamma %>%
  
  mutate(
    .orden = match(
      modelo,
      ORDEN
    )
  ) %>%
  
  arrange(
    is.na(.orden),
    .orden,
    modelo
  ) %>%
  
  select(
    -.orden
  )


cat(
  "\n-- Distribución de gamma por especificación --\n\n"
)

print(
  tabla_gamma %>%
    
    mutate(
      across(
        where(is.numeric) &
          !c(
            n,
            starts_with("n_gamma_menor_")
          ),
        ~ round(
          .x,
          3
        )
      )
    ),
  
  row.names = FALSE
)


# ==============================================================================
# 4. ESPECIFICACIONES CON GAMMA IDÉNTICO
# ==============================================================================
# reblup y reblupbc utilizan los mismos componentes de varianza en cada
# especificación, por lo que el gamma calculado coincide.
# ------------------------------------------------------------------------------

claves <- vapply(
  cols_gamma,
  
  function(cl) {
    
    paste(
      round(
        cv[[cl]],
        10
      ),
      collapse = "|"
    )
    
  },
  
  character(1)
)


grupos <- split(
  cols_gamma,
  claves
)

grupos <- grupos[
  vapply(
    grupos,
    length,
    integer(1)
  ) > 1
]


if (
  length(grupos) > 0
) {
  
  cat(
    "\n-- Especificaciones con gamma idéntico --\n"
  )
  
  for (g in grupos) {
    
    cat(
      "  ",
      paste(
        gsub(
          "\\.",
          " ",
          g
        ),
        collapse = " = "
      ),
      "\n",
      sep = ""
    )
    
  }
  
}


# ==============================================================================
# 5. COMPLEMENTO ESPACIAL: SUMA DE FILAS
# ==============================================================================
#
# En el FH espacial la matriz de encogimiento es no diagonal.
#
# Por ello se reportan:
#
#   a) diagonal: peso sobre el estimador directo del propio municipio;
#
#   b) suma de filas: peso total sobre los estimadores directos, incluyendo
#      la contribución procedente de municipios vecinos.
#
# Esta segunda cantidad es la referencia apropiada para describir el peso
# total de la información directa en el modelo espacial.
# ==============================================================================

ruta_esp <- file.path(
  ruta_out,
  "espacial.rds"
)


if (
  file.exists(ruta_esp)
) {
  
  espacial <- readRDS(
    ruta_esp
  )
  
  filas_esp <- list()
  
  
  for (
    nm in ESPACIALES
  ) {
    
    if (
      nm %in%
      names(
        espacial$gamma_esp_fila
      )
    ) {
      
      g <- espacial$gamma_esp_fila[[nm]]
      
      if (
        !is.null(g)
      ) {
        
        filas_esp[[nm]] <-
          resumir_gamma(
            as.numeric(g),
            paste0(
              "Espacial ",
              nm,
              " (suma de filas)"
            ),
            incluir_peso_sintetico = TRUE
          )
        
      }
      
    }
    
  }
  
  
  if (
    length(filas_esp) > 0
  ) {
    
    tabla_esp <- do.call(
      rbind,
      filas_esp
    )
    
    
    cat(
      "\n============================================================\n"
    )
    
    cat(
      "PESO TOTAL SOBRE LOS ESTIMADORES DIRECTOS — MODELOS ESPACIALES\n"
    )
    
    cat(
      "============================================================\n\n"
    )
    
    
    print(
      tabla_esp %>%
        
        mutate(
          across(
            where(is.numeric) &
              !c(
                n,
                starts_with("n_gamma_menor_")
              ),
            ~ round(
              .x,
              3
            )
          )
        ),
      
      row.names = FALSE
    )
    
    
    write.csv(
      tabla_esp,
      file.path(
        ruta_out,
        "gam_espacial_filas.csv"
      ),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    
    
  } else {
    
    cat(
      "\nespacial.rds no contiene gamma_esp_fila.\n"
    )
    
  }
  
} else {
  
  cat(
    "\nNo se encuentra espacial.rds; se omite el complemento espacial.\n"
  )
  
}


# ==============================================================================
# 6. EXPORTACIÓN PRINCIPAL
# ==============================================================================

write.csv(
  tabla_gamma,
  file.path(
    ruta_out,
    "gam_distribucion.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ==============================================================================
# 7. EXPORTACIÓN LATEX
# ==============================================================================
# Se genera una tabla compacta para el documento.
#
# Se priorizan:
#   - media
#   - mediana
#   - P25 y P75
#   - mínimo y máximo
#   - recuento bajo los tres umbrales
#
# El peso sintético medio no se incluye: es redundante con la media y, sobre la
# diagonal del modelo espacial, no corresponde al componente sintético.
#
# Los percentiles P10 y P90 quedan disponibles en CSV para análisis detallado.
# ==============================================================================

cols_tabla <- c(
  "media",
  "mediana",
  "p25",
  "p75",
  "minimo",
  "maximo"
)

cols_tabla <- cols_tabla[
  cols_tabla %in%
    names(tabla_gamma)
]


cols_cnt <- grep(
  "^n_gamma_menor_",
  names(tabla_gamma),
  value = TRUE
)


con <- file(
  file.path(
    ruta_out,
    "gam_distribucion.tex"
  ),
  open = "w",
  encoding = "UTF-8"
)


writeLines(
  paste0(
    "\\begin{tabular}{l",
    strrep(
      "r",
      length(cols_tabla) +
        length(cols_cnt)
    ),
    "}"
  ),
  con
)


writeLines(
  "\\hline",
  con
)


encabezados <- c(
  
  media =
    "Media",
  
  mediana =
    "Mediana",
  
  p25 =
    "P25",
  
  p75 =
    "P75",
  
  minimo =
    "Mín.",
  
  maximo =
    "Máx."
)


# Se generan a partir de cols_cnt y no se escriben a mano: si UMBRALES cambia,
# un encabezado fijo etiquetaría mal la tabla sin producir error.

encabezados_cnt <- sprintf(
  "n($\\gamma<%.2f$)",
  as.numeric(
    sub(
      "^n_gamma_menor_",
      "",
      cols_cnt
    )
  )
)


writeLines(
  
  paste0(
    
    "Modelo & ",
    
    paste(
      
      c(
        unname(
          encabezados[
            cols_tabla
          ]
        ),
        
        encabezados_cnt
      ),
      
      collapse = " & "
      
    ),
    
    " \\\\"
    
  ),
  
  con
  
)


writeLines(
  "\\hline",
  con
)


for (
  i in seq_len(
    nrow(tabla_gamma)
  )
) {
  
  valores <- c(
    
    sprintf(
      "%.3f",
      unlist(
        tabla_gamma[
          i,
          cols_tabla
        ]
      )
    ),
    
    sprintf(
      "%d",
      unlist(
        tabla_gamma[
          i,
          cols_cnt
        ]
      )
    )
    
  )
  
  
  writeLines(
    
    paste0(
      
      tabla_gamma$modelo[i],
      
      " & ",
      
      paste(
        valores,
        collapse = " & "
      ),
      
      " \\\\"
      
    ),
    
    con
  )
  
}


writeLines(
  "\\hline",
  con
)

writeLines(
  "\\end{tabular}",
  con
)

close(con)


# ==============================================================================
# 8. MENSAJE FINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "DISTRIBUCIÓN DE GAMMA COMPLETADA\n"
)

cat(
  "============================================================\n"
)

cat(
  "Salida principal:\n"
)

cat(
  "  - output/gam_distribucion.csv\n"
)

cat(
  "  - output/gam_distribucion.tex\n"
)

cat(
  "Salida espacial:\n"
)

cat(
  "  - output/gam_espacial_filas.csv\n"
)

cat(
  "============================================================\n"
)