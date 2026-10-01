# ==============================================================================
# 05_3_shapley_cv_m1_log.R
# ------------------------------------------------------------------------------
# Shapley-CV exacto sobre las variables de M1
#
# OBJETIVO
#   Cuantificar la contribución marginal promedio de cada variable de M1
#   a la reducción del CV promedio de los predictores SAE.
#
# FUNCIÓN CARACTERÍSTICA
#
#   v(S) = CV_promedio(modelo vacío) - CV_promedio(modelo S)
#
#   Por tanto:
#
#       v(emptyset) = 0
#
#   y un valor Shapley positivo indica una contribución promedio
#   a la reducción del CV.
#
# CONFIGURACIÓN
#   - FH
#   - REML
#   - Transformación logarítmica
#   - Backtransformation: bc_crude
#   - MSE analítico
#   - 6 variables de M1
#   - 2^6 = 64 coaliciones, incluida la vacía
#   - Shapley exacto
#
# INTERPRETACIÓN
#
#   Shapley_CV > 0:
#       la variable contribuye, en promedio, a reducir el CV.
#
#   Shapley_CV < 0:
#       la variable contribuye, en promedio, a aumentar el CV.
#
#   La unidad es puntos porcentuales de CV.
#
# IMPORTANTE
#   Este procedimiento NO es LOOCV ni validación fuera de muestra.
#   Es una medida de contribución marginal a la precisión interna
#   de los predictores SAE.
# ==============================================================================


# ==============================================================================
# 0. PAQUETES
# ==============================================================================

if (!requireNamespace("emdi", quietly = TRUE)) {
  stop("Instale el paquete 'emdi'.")
}

if (!requireNamespace("TUvalues", quietly = TRUE)) {
  stop("Instale el paquete 'TUvalues'.")
}

library(emdi)
library(TUvalues)
library(dplyr)
library(here)


# ==============================================================================
# 1. RUTAS Y DATOS
# ==============================================================================

ruta_out <- here("output")

matriz <- readRDS(
  file.path(
    ruta_out,
    "matriz_sae_transformada_v2.rds"
  )
)

datos <- as.data.frame(matriz)


# ==============================================================================
# 2. VARIABLES DE M1
# ==============================================================================

vars_m1 <- c(
  "alcantarillado_2018",
  "pct_rural",
  "partos_calificados_2020",
  "transferencias_pc",
  "saber11_lectura_2022",
  "tasa_violencia_intra_2019"
)

p <- length(vars_m1)

stopifnot(p == 6)


# ==============================================================================
# 3. VALIDACIONES
# ==============================================================================

vars_necesarias <- c(
  "cod_mun",
  "pobreza_monetaria",
  "varianza_pobreza",
  vars_m1
)

faltantes <- setdiff(
  vars_necesarias,
  names(datos)
)

if (length(faltantes) > 0) {
  
  stop(
    "Faltan variables en los datos: ",
    paste(faltantes, collapse = ", ")
  )
  
}


stopifnot(
  !anyNA(
    datos[
      ,
      vars_necesarias,
      drop = FALSE
    ]
  )
)

stopifnot(
  all(
    is.finite(
      datos$pobreza_monetaria
    )
  )
)

stopifnot(
  all(
    datos$pobreza_monetaria > 0
  )
)

stopifnot(
  all(
    is.finite(
      datos$varianza_pobreza
    )
  )
)

stopifnot(
  all(
    datos$varianza_pobreza > 0
  )
)

stopifnot(
  !anyDuplicated(
    datos$cod_mun
  )
)


# ==============================================================================
# 4. FÓRMULA
# ==============================================================================

formula_fh <- function(vars) {
  
  if (length(vars) == 0) {
    
    return(
      pobreza_monetaria ~ 1
    )
    
  }
  
  as.formula(
    paste(
      "pobreza_monetaria ~",
      paste(
        vars,
        collapse = " + "
      )
    )
  )
  
}


# ==============================================================================
# 5. EXTRAER PREDICCIÓN, MSE Y CV
# ==============================================================================

extraer_metricas <- function(
    ajuste,
    etiqueta = ""
) {
  
  # --------------------------------------------------------------------------
  # Verificar orden de municipios
  # --------------------------------------------------------------------------
  
  dominios_ajuste <- tryCatch(
    as.character(
      ajuste$ind$Domain
    ),
    error = function(e) NULL
  )
  
  dominios_datos <- as.character(
    datos$cod_mun
  )
  
  if (
    is.null(dominios_ajuste) ||
    !identical(
      dominios_ajuste,
      dominios_datos
    )
  ) {
    
    stop(
      "El orden de dominios no coincide en ",
      etiqueta
    )
    
  }
  
  
  # --------------------------------------------------------------------------
  # Extraer columna FH
  # --------------------------------------------------------------------------
  
  extraer_columna_fh <- function(
    objeto,
    patron = "^FH$"
  ) {
    
    if (is.null(objeto)) {
      return(NULL)
    }
    
    i <- grep(
      patron,
      names(objeto)
    )
    
    if (length(i) == 0) {
      return(NULL)
    }
    
    as.numeric(
      objeto[[i[1]]]
    )
    
  }
  
  
  # --------------------------------------------------------------------------
  # Predicción FH
  # --------------------------------------------------------------------------
  
  pred <- extraer_columna_fh(
    ajuste$ind
  )
  
  if (
    is.null(pred) ||
    length(pred) != nrow(datos)
  ) {
    
    stop(
      "No fue posible extraer la estimación FH de ",
      etiqueta
    )
    
  }
  
  
  # --------------------------------------------------------------------------
  # MSE
  # --------------------------------------------------------------------------
  
  mse <- extraer_columna_fh(
    ajuste$MSE
  )
  
  if (
    is.null(mse) ||
    length(mse) != nrow(datos)
  ) {
    
    stop(
      "No fue posible extraer el MSE de ",
      etiqueta
    )
    
  }
  
  
  # --------------------------------------------------------------------------
  # CV
  #
  # CV_i = sqrt(MSE_i) / theta_hat_i * 100
  # --------------------------------------------------------------------------
  
  admisible <- (
    is.finite(pred) &
      pred > 0 &
      pred <= 1 &
      is.finite(mse) &
      mse >= 0
  )
  
  cv <- rep(
    NA_real_,
    length(pred)
  )
  
  cv[admisible] <-
    sqrt(
      mse[admisible]
    ) /
    pred[admisible] *
    100
  
  
  list(
    pred = pred,
    mse = mse,
    cv = cv,
    admisible = admisible
  )
  
}


# ==============================================================================
# 6. AJUSTAR UN MODELO FH
# ==============================================================================

ajustar_fh <- function(
    vars
) {
  
  emdi::fh(
    
    fixed =
      formula_fh(vars),
    
    vardir =
      "varianza_pobreza",
    
    combined_data =
      datos,
    
    domains =
      "cod_mun",
    
    method =
      "reml",
    
    transformation =
      "log",
    
    backtransformation =
      "bc_crude",
    
    MSE =
      TRUE,
    
    mse_type =
      "analytical",
    
    B =
      c(500, 0),
    
    tol =
      1e-8,
    
    maxit =
      1000
  )
  
}


# ==============================================================================
# 7. GENERAR LAS 64 COALICIONES
# ------------------------------------------------------------------------------

# Incluye:
#
#   - coalición vacía
#   - 63 coaliciones no vacías
#
# ==============================================================================

coal <- TUvalues::coalitions(
  p
)

bin_coal <- as.matrix(
  coal$Binary
)

n_modelos <- nrow(
  bin_coal
)

stopifnot(
  n_modelos == 2^p
)


# Índices de las variables pertenecientes a cada coalición
all_comb <- lapply(
  seq_len(nrow(bin_coal)),
  function(r) {
    unname(
      which(
        bin_coal[r, ] == 1
      )
    )
  }
)


# Claves para identificar coaliciones
claves_tuv <- vapply(
  seq_len(nrow(bin_coal)),
  function(r) {
    
    idx <- which(
      bin_coal[r, ] == 1
    )
    
    if (length(idx) == 0) {
      ""
    } else {
      paste(
        idx,
        collapse = "|"
      )
    }
    
  },
  character(1)
)


cat(
  "\n============================================================\n"
)

cat(
  "SHAPLEY-CV - M1 LOG\n"
)

cat(
  "============================================================\n"
)

cat(
  "Variables:",
  p,
  "\n"
)

cat(
  "Coaliciones:",
  n_modelos,
  "\n"
)

cat(
  "Coalición vacía incluida: sí\n"
)


# ==============================================================================
# 8. AJUSTAR LAS 64 COALICIONES
# ==============================================================================

resultados <- vector(
  "list",
  n_modelos
)

cv_municipal <- matrix(
  NA_real_,
  nrow = nrow(datos),
  ncol = n_modelos
)

colnames(
  cv_municipal
) <- paste0(
  "modelo_",
  seq_len(n_modelos)
)


t0 <- Sys.time()


for (
  j in seq_along(all_comb)
) {
  
  idx <- all_comb[[j]]
  
  vars_j <- vars_m1[
    idx
  ]
  
  
  clave_j <- claves_tuv[j]
  
  
  cat(
    sprintf(
      "Modelo %2d / %2d | %d variables | %s\n",
      j,
      n_modelos,
      length(vars_j),
      if (
        length(vars_j) == 0
      ) {
        "(intercepto)"
      } else {
        paste(
          vars_j,
          collapse = " + "
        )
      }
    )
  )
  
  
  ajuste_j <- tryCatch(
    
    ajustar_fh(
      vars_j
    ),
    
    error = function(e) {
      
      stop(
        "Error en modelo ",
        j,
        " (",
        ifelse(
          clave_j == "",
          "vacío",
          clave_j
        ),
        "): ",
        conditionMessage(e)
      )
      
    }
    
  )
  
  
  met_j <- extraer_metricas(
    ajuste_j,
    paste(
      "modelo",
      j
    )
  )
  
  
  # --------------------------------------------------------------------------
  # MUY IMPORTANTE:
  #
  # Para que la comparación entre coaliciones sea limpia,
  # exigimos que los 84 municipios sean admisibles en TODOS los modelos.
  # --------------------------------------------------------------------------
  
  if (
    !all(
      met_j$admisible
    )
  ) {
    
    no_adm <- sum(
      !met_j$admisible
    )
    
    stop(
      "El modelo ",
      j,
      " tiene ",
      no_adm,
      " municipios con CV no admisible. ",
      "No se continúa para evitar comparar medias calculadas ",
      "sobre diferentes conjuntos de municipios."
    )
    
  }
  
  
  cv_municipal[
    ,
    j
  ] <-
    met_j$cv
  
  
  resultados[[j]] <- data.frame(
    
    modelo =
      j,
    
    clave =
      clave_j,
    
    variables =
      if (
        length(vars_j) == 0
      ) {
        "(intercepto)"
      } else {
        paste(
          vars_j,
          collapse = " + "
        )
      },
    
    n_variables =
      length(vars_j),
    
    media_CV =
      mean(
        met_j$cv
      ),
    
    mediana_CV =
      median(
        met_j$cv
      ),
    
    p90_CV =
      as.numeric(
        quantile(
          met_j$cv,
          0.90
        )
      ),
    
    max_CV =
      max(
        met_j$cv
      ),
    
    media_MSE =
      mean(
        met_j$mse
      ),
    
    mediana_MSE =
      median(
        met_j$mse
      ),
    
    n_admisibles =
      sum(
        met_j$admisible
      ),
    
    stringsAsFactors =
      FALSE
  )
  
}


cat(
  "\nTiempo total:",
  round(
    as.numeric(
      difftime(
        Sys.time(),
        t0,
        units = "mins"
      )
    ),
    2
  ),
  "minutos\n"
)


# ==============================================================================
# 9. TABLA DE LAS 64 COALICIONES
# ==============================================================================

tabla_modelos <- bind_rows(
  resultados
)


stopifnot(
  nrow(tabla_modelos) == 64
)


# Verificar que el orden sea exactamente el de TUvalues
stopifnot(
  identical(
    tabla_modelos$clave,
    claves_tuv
  )
)


# Verificar que todos tengan 84 municipios
stopifnot(
  all(
    tabla_modelos$n_admisibles ==
      nrow(datos)
  )
)


# ==============================================================================
# 10. FUNCIÓN CARACTERÍSTICA DEL JUEGO
# ------------------------------------------------------------------------------
#
# v(S) =
#
#   CV_promedio(modelo vacío)
#   -
#   CV_promedio(modelo S)
#
# Por construcción:
#
#   v(emptyset) = 0
#
# y valores positivos indican reducción del CV.
# ==============================================================================

clave_vacia <- ""

cv_vacio <- tabla_modelos$media_CV[
  tabla_modelos$clave == clave_vacia
]

if (
  length(cv_vacio) != 1
) {
  
  stop(
    "No se pudo identificar el modelo vacío."
  )
  
}


tabla_modelos <- tabla_modelos %>%
  
  mutate(
    valor_CV =
      cv_vacio -
      media_CV
  )


# ==============================================================================
# 11. SHAPLEY-CV EXACTO
# ------------------------------------------------------------------------------

# TUvalues::shapley() recibe las 63 coaliciones NO vacías.
#
# La coalición vacía tiene:
#
#   v(emptyset) = 0
#
# porque hemos centrado la función característica alrededor
# del CV del modelo vacío.
# ==============================================================================

valor_CV_no_vacio <- tabla_modelos$valor_CV[
  tabla_modelos$clave != ""
]


stopifnot(
  length(
    valor_CV_no_vacio
  ) == 2^p - 1
)


shapley_cv <- TUvalues::shapley(
  valor_CV_no_vacio,
  method = "exact"
)


stopifnot(
  length(
    shapley_cv
  ) == p
)


# ==============================================================================
# 12. RESULTADO SHAPLEY-CV
# ==============================================================================

resultado_shapley <- data.frame(
  
  variable =
    vars_m1,
  
  Shapley_CV_pp =
    as.numeric(
      shapley_cv
    ),
  
  stringsAsFactors =
    FALSE
  
) %>%
  
  arrange(
    desc(
      Shapley_CV_pp
    )
  ) %>%
  
  mutate(
    posicion =
      row_number()
  ) %>%
  
  select(
    posicion,
    variable,
    Shapley_CV_pp
  )


# ==============================================================================
# 13. PROPIEDAD DE EFICIENCIA
# ------------------------------------------------------------------------------

cv_full <- tabla_modelos$media_CV[
  tabla_modelos$clave ==
    paste(
      seq_len(p),
      collapse = "|"
    )
]


ganancia_total_CV <-
  cv_vacio -
  cv_full


suma_shapley <-
  sum(
    resultado_shapley$Shapley_CV_pp
  )


diferencia_eficiencia <-
  suma_shapley -
  ganancia_total_CV


cat(
  "\n============================================================\n"
)

cat(
  "EFICIENCIA DE SHAPLEY-CV\n"
)

cat(
  "============================================================\n"
)

cat(
  "CV modelo vacío: ",
  cv_vacio,
  "%\n",
  sep = ""
)

cat(
  "CV modelo completo: ",
  cv_full,
  "%\n",
  sep = ""
)

cat(
  "Reducción total del CV: ",
  ganancia_total_CV,
  " pp\n",
  sep = ""
)

cat(
  "Suma Shapley: ",
  suma_shapley,
  " pp\n",
  sep = ""
)

cat(
  "Diferencia: ",
  format(
    diferencia_eficiencia,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)


if (
  abs(
    diferencia_eficiencia
  ) > 1e-8
) {
  
  warning(
    "La propiedad de eficiencia no se cumple dentro de la tolerancia."
  )
  
}


# ==============================================================================
# 14. MOSTRAR RESULTADOS
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "RANKING SHAPLEY-CV\n"
)

cat(
  "============================================================\n"
)

print(
  resultado_shapley,
  row.names = FALSE
)


# ==============================================================================
# 15. GUARDAR RESULTADOS
# ==============================================================================

write.csv(
  
  tabla_modelos,
  
  file.path(
    ruta_out,
    "shapley_cv_m1_log_modelos.csv"
  ),
  
  row.names =
    FALSE,
  
  fileEncoding =
    "UTF-8"
  
)


write.csv(
  
  resultado_shapley,
  
  file.path(
    ruta_out,
    "shapley_cv_m1_log_ranking.csv"
  ),
  
  row.names =
    FALSE,
  
  fileEncoding =
    "UTF-8"
  
)


write.csv(
  
  data.frame(
    cod_mun =
      datos$cod_mun,
    cv_municipal,
    check.names =
      FALSE
  ),
  
  file.path(
    ruta_out,
    "shapley_cv_m1_log_cv_municipal.csv"
  ),
  
  row.names =
    FALSE,
  
  fileEncoding =
    "UTF-8"
  
)


saveRDS(
  
  list(
    
    variables_m1 =
      vars_m1,
    
    tabla_modelos =
      tabla_modelos,
    
    resultado_shapley =
      resultado_shapley,
    
    cv_municipal =
      cv_municipal,
    
    CV_vacio =
      cv_vacio,
    
    CV_full =
      cv_full,
    
    ganancia_total_CV =
      ganancia_total_CV,
    
    suma_shapley =
      suma_shapley,
    
    diferencia_eficiencia =
      diferencia_eficiencia
    
  ),
  
  file.path(
    ruta_out,
    "shapley_cv_m1_log.rds"
  )
  
)


cat(
  "\nProceso terminado correctamente.\n"
)