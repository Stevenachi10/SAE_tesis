# ==============================================================================
# 04_1_benchmarking_m1log.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA IV.1: Benchmarking de la especificación seleccionada (M1 log)
# ==============================================================================

# ------------------------------------------------------------------------------
# OBJETIVO
# ------------------------------------------------------------------------------
# Ajustar las estimaciones municipales de M1 log para que su agregado
# departamental, ponderado por población, coincida con la estimación oficial
# de Cauca y Valle del Cauca.
#
# Se compara:
#   1. estimación municipal antes y después del benchmarking;
#   2. magnitud del ajuste municipal;
#   3. agregado departamental antes y después del ajuste;
#   4. discrepancia municipal frente al estimador directo;
#   5. diagnóstico de consistencia de Brown, en la dirección directa y en la
#      formulación inversa de SAEval.
#
# No se recalcula el MSE ni el CV después del benchmarking. El MSE del EBLUP
# se conserva en la base solo como referencia de la estimación sin ajuste.
# ==============================================================================


# ==============================================================================
# 0. PAQUETES Y PARÁMETROS
# ==============================================================================

library(emdi)
library(dplyr)
library(here)
library(car)

stopifnot(
  "Falta el paquete sandwich (errores HC3)." =
    requireNamespace("sandwich", quietly = TRUE),
  "Falta el paquete SAEval." =
    requireNamespace("SAEval", quietly = TRUE)
)

select <- dplyr::select

ruta_out <- here("output")

AJUSTE_FINAL <- "M1_reml_log"
TIPO_BENCH   <- "ratio"
COL_PESO     <- "log_pob_total"
PESO_EN_LOGARITMO <- TRUE


# ==============================================================================
# 1. CARGA DE RESULTADOS
# ==============================================================================

variantes <- readRDS(
  file.path(ruta_out, "modelado_variantes.rds")
)

matriz <- readRDS(
  file.path(ruta_out, "matriz_sae_transformada_v2.rds")
)

ruta_ref <- file.path(
  ruta_out,
  "agr_referencia_departamental.csv"
)

stopifnot(
  "No se encuentra la referencia departamental." =
    file.exists(ruta_ref)
)

referencia <- read.csv(
  ruta_ref,
  colClasses = c(cod_dpto = "character"),
  fileEncoding = "UTF-8"
)

aj <- variantes$ajustes[[AJUSTE_FINAL]]

stopifnot(
  "No se encuentra el ajuste final." = !is.null(aj),
  "El ajuste final no tiene MSE." = !is.null(aj$MSE)
)

cat("\n============================================================\n")
cat("BENCHMARKING DE", AJUSTE_FINAL, "\n")
cat("============================================================\n")
# transformation y method son listas en los objetos de emdi
cat("Transformación:",
    paste(unlist(aj$transformation), collapse = " / "), "\n")
cat("MSE utilizado por el ajuste:",
    paste(unlist(aj$method$MSE_method), collapse = " / "), "\n")
cat("Tipo de benchmarking:", TIPO_BENCH, "\n")


# ==============================================================================
# 2. DATOS MUNICIPALES Y PONDERADORES
# ==============================================================================

pesos <- matriz %>%
  as.data.frame() %>%
  select(
    cod_mun,
    peso = all_of(COL_PESO)
  ) %>%
  mutate(
    cod_mun = as.character(cod_mun),
    peso = if (PESO_EN_LOGARITMO) exp(peso) else peso
  )

stopifnot(
  "Hay municipios con código duplicado en los ponderadores." =
    !anyDuplicated(pesos$cod_mun)
)

datos <- variantes$datos %>%
  mutate(
    cod_mun = as.character(cod_mun)
  ) %>%
  left_join(
    pesos,
    by = "cod_mun"
  ) %>%
  mutate(
    cod_dpto = substr(cod_mun, 1, 2)
  )

stopifnot(
  "Municipios sin ponderador." =
    all(is.finite(datos$peso) & datos$peso > 0),
  
  "Dominios del ajuste no alineados con los datos." =
    identical(
      as.character(aj$ind$Domain),
      datos$cod_mun
    ),
  
  "Existe un departamento sin cifra oficial." =
    all(
      unique(datos$cod_dpto) %in% referencia$cod_dpto
    )
)


# ==============================================================================
# 3. BASE MUNICIPAL
# ==============================================================================

base <- data.frame(
  cod_mun   = datos$cod_mun,
  Municipio = datos$Municipio,
  cod_dpto  = datos$cod_dpto,
  peso      = datos$peso,
  directo   = datos$pobreza_monetaria,
  psi       = datos$varianza_pobreza,
  fh        = as.numeric(aj$ind$FH),
  mse       = as.numeric(aj$MSE$FH),
  stringsAsFactors = FALSE
)

# Verificación del estimador directo
if (!isTRUE(
  all.equal(
    base$directo,
    as.numeric(aj$ind$Direct)
  )
)) {
  warning(
    "El estimador directo del ajuste no coincide con pobreza_monetaria."
  )
}

stopifnot(
  "Estimaciones FH o MSE no válidos." =
    all(
      is.finite(base$fh) &
        base$fh > 0 &
        is.finite(base$mse) &
        base$mse >= 0
    )
)


# ==============================================================================
# 4. FUNCIÓN DE BENCHMARKING
# ==============================================================================

benchmark_dpto <- function(aj, idx, valor_oficial, share) {
  
  obj <- aj
  
  obj$ind <- aj$ind[idx, , drop = FALSE]
  obj$MSE <- aj$MSE[idx, , drop = FALSE]
  
  res <- emdi::benchmark(
    obj,
    benchmark = valor_oficial,
    share     = share,
    type      = TIPO_BENCH
  )
  
  as.numeric(res$FH_Bench)
}


# ==============================================================================
# 5. BENCHMARKING POR DEPARTAMENTO
# ==============================================================================

base$fh_b <- NA_real_

for (dp in unique(base$cod_dpto)) {
  
  idx <- which(
    base$cod_dpto == dp
  )
  
  # Participación poblacional municipal dentro del departamento
  share <- base$peso[idx] /
    sum(base$peso[idx])
  
  # Cifra oficial departamental
  oficial <- referencia$estimacion[
    referencia$cod_dpto == dp
  ]
  
  # Benchmarking ratio
  base$fh_b[idx] <- benchmark_dpto(
    aj            = aj,
    idx           = idx,
    valor_oficial = oficial,
    share         = share
  )
}


# ==============================================================================
# 6. MEDIDAS ANTES Y DESPUÉS DEL BENCHMARKING
# ==============================================================================

base <- base %>%
  mutate(
    ajuste_pp   = (fh_b - fh) * 100,
    ajuste_rel  = (fh_b / fh - 1) * 100,
    
    discrepancia_antes =
      (fh - directo) * 100,
    
    discrepancia_despues =
      (fh_b - directo) * 100
  )


# ==============================================================================
# 7. RESUMEN DEPARTAMENTAL
# ==============================================================================

resumen_dpto <- lapply(
  unique(base$cod_dpto),
  function(dp) {
    
    idx <- which(
      base$cod_dpto == dp
    )
    
    share <- base$peso[idx] /
      sum(base$peso[idx])
    
    oficial <- referencia$estimacion[
      referencia$cod_dpto == dp
    ]
    
    agregado_antes <- sum(
      share * base$fh[idx]
    )
    
    agregado_despues <- sum(
      share * base$fh_b[idx]
    )
    
    data.frame(
      cod_dpto         = dp,
      departamento     = referencia$departamento[
        referencia$cod_dpto == dp
      ],
      oficial          = oficial,
      agregado_antes   = agregado_antes,
      agregado_despues = agregado_despues,
      
      diferencia_antes =
        (agregado_antes - oficial) * 100,
      
      diferencia_despues =
        (agregado_despues - oficial) * 100,
      
      ajuste_medio_pp =
        mean(base$ajuste_pp[idx]),
      
      ajuste_min_pp =
        min(base$ajuste_pp[idx]),
      
      ajuste_max_pp =
        max(base$ajuste_pp[idx]),
      
      stringsAsFactors = FALSE
    )
  }
)

resumen_dpto <- do.call(
  rbind,
  resumen_dpto
)


# ==============================================================================
# 8. VERIFICACIONES
# ==============================================================================

# El agregado benchmarkeado debe coincidir exactamente con la cifra oficial
stopifnot(
  "El agregado benchmarkeado no reproduce la cifra oficial." =
    all(
      abs(
        resumen_dpto$agregado_despues -
          resumen_dpto$oficial
      ) < 1e-10
    )
)


# ------------------------------------------------------------------------------
# Verificación manual del benchmarking ratio
# ------------------------------------------------------------------------------

fh_ratio_manual <- base$fh

for (dp in unique(base$cod_dpto)) {
  
  idx <- which(
    base$cod_dpto == dp
  )
  
  share <- base$peso[idx] /
    sum(base$peso[idx])
  
  oficial <- referencia$estimacion[
    referencia$cod_dpto == dp
  ]
  
  agregado <- sum(
    share * base$fh[idx]
  )
  
  fh_ratio_manual[idx] <-
    base$fh[idx] *
    oficial / agregado
}

dif_ratio <- max(
  abs(
    fh_ratio_manual -
      base$fh_b
  )
)

cat(
  "\nDiferencia máxima entre emdi::benchmark() y cálculo manual:",
  signif(dif_ratio, 6),
  "\n"
)

stopifnot(
  "El benchmarking ratio de emdi no coincide con el cálculo manual." =
    dif_ratio < 1e-10
)


# ------------------------------------------------------------------------------
# Verificación del rango de las estimaciones
# ------------------------------------------------------------------------------

fuera_rango <- sum(
  base$fh_b <= 0 |
    base$fh_b >= 1
)

cat(
  "Estimaciones benchmarkeadas fuera de (0,1):",
  fuera_rango,
  "\n"
)


# ==============================================================================
# 9. RESUMEN GLOBAL DEL AJUSTE
# ==============================================================================

resumen_global <- data.frame(
  
  ajuste_medio_abs_pp =
    mean(abs(base$ajuste_pp)),
  
  ajuste_mediano_abs_pp =
    median(abs(base$ajuste_pp)),
  
  ajuste_max_abs_pp =
    max(abs(base$ajuste_pp)),
  
  n_ajuste_positivo =
    sum(base$ajuste_pp > 0),
  
  n_ajuste_negativo =
    sum(base$ajuste_pp < 0),
  
  n_sin_ajuste =
    sum(abs(base$ajuste_pp) < 1e-12),
  
  disc_media_abs_antes_pp =
    mean(abs(base$discrepancia_antes)),
  
  disc_media_abs_despues_pp =
    mean(abs(base$discrepancia_despues)),
  
  disc_max_abs_antes_pp =
    max(abs(base$discrepancia_antes)),
  
  disc_max_abs_despues_pp =
    max(abs(base$discrepancia_despues))
)


# ==============================================================================
# 10. DIAGNÓSTICO DE CONSISTENCIA ANTES Y DESPUÉS DEL BENCHMARKING
#     Brown directo + formulación inversa de SAEval
# ==============================================================================

# ------------------------------------------------------------------------------
# 10.1 Brown en dirección directa
#
# Directo_i = beta0 + beta1 * SAE_i + error_i
#
# H0: beta0 = 0 y beta1 = 1
# Se reportan errores estándar convencionales y HC3.
# ------------------------------------------------------------------------------

brown_directo <- function(data, variable_sae, etiqueta) {
  
  d <- data.frame(
    directo = data$directo,
    sae = data[[variable_sae]]
  )
  
  d <- d[
    is.finite(d$directo) &
      is.finite(d$sae),
    ,
    drop = FALSE
  ]
  
  fit <- lm(
    directo ~ sae,
    data = d
  )
  
  H0 <- c(
    "(Intercept) = 0",
    "sae = 1"
  )
  
  # Contraste con errores estándar convencionales
  test_clasico <- car::linearHypothesis(
    fit,
    H0,
    test = "F"
  )
  
  # Contraste con errores estándar HC3
  test_hc3 <- car::linearHypothesis(
    fit,
    H0,
    vcov. = sandwich::vcovHC(
      fit,
      type = "HC3"
    ),
    test = "F"
  )
  
  data.frame(
    momento = etiqueta,
    direccion = "Directa (Brown)",
    intercepto = unname(coef(fit)[1]),
    pendiente = unname(coef(fit)[2]),
    R2 = summary(fit)$r.squared,
    p_clasico = test_clasico[2, "Pr(>F)"],
    p_HC3 = test_hc3[2, "Pr(>F)"],
    n = nrow(d),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# 10.2 Formulación inversa mediante SAEval
#
# SAE_i = beta0 + beta1 * Directo_i + error_i
#
# Esta es la dirección utilizada por SAEval::bias().
# ------------------------------------------------------------------------------

brown_sa_eval <- function(data, variable_sae, etiqueta) {
  
  d <- data.frame(
    directo = data$directo,
    sae = data[[variable_sae]]
  )
  
  d <- d[
    is.finite(d$directo) &
      is.finite(d$sae),
    ,
    drop = FALSE
  ]
  
  # SAEval: directo en X, SAE en Y
  resultado <- SAEval::bias(
    data = d,
    dir = ~directo,
    sae = ~sae,
    scatterplot = FALSE
  )
  
  # output1 corresponde a la regresión sin transformación
  out <- resultado$output1
  
  out <- out[
    out$methods == "sae",
    ,
    drop = FALSE
  ]
  
  if (nrow(out) == 0) {
    # Por seguridad, si la etiqueta interna de SAEval difiere
    out <- resultado$output1[1, , drop = FALSE]
  }
  
  # También calculamos manualmente la formulación inversa,
  # para tener el p-valor conjunto de forma explícita y comparable.
  fit <- lm(
    sae ~ directo,
    data = d
  )
  
  H0 <- c(
    "(Intercept) = 0",
    "directo = 1"
  )
  
  test_clasico <- car::linearHypothesis(
    fit,
    H0,
    test = "F"
  )
  
  test_hc3 <- car::linearHypothesis(
    fit,
    H0,
    vcov. = sandwich::vcovHC(
      fit,
      type = "HC3"
    ),
    test = "F"
  )
  
  data.frame(
    momento = etiqueta,
    direccion = "Inversa (SAEval)",
    intercepto = unname(coef(fit)[1]),
    pendiente = unname(coef(fit)[2]),
    R2 = summary(fit)$r.squared,
    
    # p-valores del contraste conjunto de la formulación inversa, con los
    # mismos nombres que la dirección directa para poder unir las tablas
    p_clasico = test_clasico[2, "Pr(>F)"],
    p_HC3 = test_hc3[2, "Pr(>F)"],
    
    # Decisión textual producida por SAEval, conservada para verificación
    SAEval_decision = if ("F" %in% names(out)) as.character(out$F) else NA_character_,
    
    n = nrow(d),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# 10.3 Aplicación antes y después del benchmarking
# ------------------------------------------------------------------------------

brown_directo_resultados <- rbind(
  
  brown_directo(
    data = base,
    variable_sae = "fh",
    etiqueta = "Antes"
  ),
  
  brown_directo(
    data = base,
    variable_sae = "fh_b",
    etiqueta = "Después"
  )
)


brown_sa_eval_resultados <- rbind(
  
  brown_sa_eval(
    data = base,
    variable_sae = "fh",
    etiqueta = "Antes"
  ),
  
  brown_sa_eval(
    data = base,
    variable_sae = "fh_b",
    etiqueta = "Después"
  )
)


# ------------------------------------------------------------------------------
# 10.4 Unificación de resultados
# ------------------------------------------------------------------------------

brown_benchmark <- bind_rows(
  brown_directo_resultados,
  brown_sa_eval_resultados
)


# ------------------------------------------------------------------------------
# 10.5 Mostrar resultados
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("BROWN: ANTES Y DESPUÉS DEL BENCHMARKING\n")
cat("============================================================\n\n")

print(
  brown_benchmark %>%
    mutate(
      across(
        where(is.numeric),
        ~ round(.x, 4)
      )
    ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 10.6 Verificación: la decisión de SAEval coincide con el cálculo manual
# ------------------------------------------------------------------------------

inv <- brown_sa_eval_resultados

decision_manual <- ifelse(inv$p_clasico < 0.05, "Reject", "Accept")
decision_saeval <- ifelse(grepl("^Reject", inv$SAEval_decision), "Reject",
                          ifelse(grepl("^Accept", inv$SAEval_decision), "Accept", NA))

if (any(is.na(decision_saeval)) || any(decision_manual != decision_saeval)) {
  warning("La decisión de SAEval no coincide con el cálculo manual en la formulación inversa.")
} else {
  cat("\nLa decisión de SAEval coincide con el cálculo manual.\n")
}



# ==============================================================================
# 11. RESULTADOS
# ==============================================================================

cat("\n-- Agregados departamentales --\n\n")

print(
  resumen_dpto %>%
    mutate(
      across(
        where(is.numeric),
        ~ round(.x, 4)
      )
    ),
  row.names = FALSE
)

cat("\n-- Resumen global del ajuste --\n\n")

print(
  resumen_global %>%
    mutate(
      across(
        where(is.numeric),
        ~ round(.x, 4)
      )
    ),
  row.names = FALSE
)

cat("\n-- Municipios con mayor ajuste absoluto --\n\n")

print(
  base %>%
    arrange(
      desc(abs(ajuste_pp))
    ) %>%
    select(
      Municipio,
      cod_dpto,
      directo,
      fh,
      fh_b,
      ajuste_pp,
      ajuste_rel
    ) %>%
    mutate(
      across(
        where(is.numeric),
        ~ round(.x, 3)
      )
    ) %>%
    head(10),
  row.names = FALSE
)


# ==============================================================================
# 12. EXPORTACIÓN
# ==============================================================================

write.csv(
  base,
  file.path(
    ruta_out,
    "bench_municipal.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  resumen_dpto,
  file.path(
    ruta_out,
    "bench_resumen_dpto.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  resumen_global,
  file.path(
    ruta_out,
    "bench_resumen_global.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  brown_benchmark,
  file.path(
    ruta_out,
    "bench_brown.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

saveRDS(
  list(
    ajuste          = AJUSTE_FINAL,
    tipo_bench      = TIPO_BENCH,
    municipal       = base,
    resumen_dpto    = resumen_dpto,
    resumen_global  = resumen_global,
    brown_benchmark = brown_benchmark
  ),
  file.path(
    ruta_out,
    "bench_m1log.rds"
  )
)


cat("\n============================================================\n")
cat("BENCHMARKING COMPLETADO\n")
cat("============================================================\n")