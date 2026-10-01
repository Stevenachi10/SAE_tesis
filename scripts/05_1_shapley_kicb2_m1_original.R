# ==============================================================================
# 05_1_shapley_kicb2_m1_original.R
# ------------------------------------------------------------------------------
# Shapley-KICb2 exacto sobre las variables de M1.
#
# IMPORTANTE:
#   - El juego se define mediante KICb2 de sae::eblupFH(), bajo ML.
#   - KICb2 corresponde al FH en escala original, que es la escala y el
#     criterio con los que se seleccionó M1.
#   - M1 permanece fijo: solo se estudian sus seis variables.
#   - No se vuelve a seleccionar entre las variables originales.
#
# DECISIONES
#   - Covariables estandarizadas (media 0, desviación 1). Con intercepto, la
#     matriz de diseño estandarizada genera el mismo espacio de columnas, de
#     modo que la verosimilitud, los valores ajustados, las réplicas bootstrap
#     y el KICb2 no cambian. Solo mejora el condicionamiento numérico: sin
#     estandarizar, solve() falla dentro del bootstrap en algunas coaliciones.
#   - Una sola semilla para todas las coaliciones (números aleatorios comunes).
#     El valor de Shapley se construye con diferencias de KICb2 entre
#     coaliciones; con la misma semilla, esas diferencias no mezclan ruido de
#     Monte Carlo de corridas distintas.
#   - Réplicas que fallan. sae::eblupFH() descarta y reemplaza las réplicas
#     bootstrap cuyo reajuste no converge, pero se detiene por completo si el
#     reajuste produce un sistema singular. Aquí se reproduce el mismo bucle
#     bootstrap de sae (versión 1.3) y ambas situaciones se tratan igual: la
#     réplica se descarta y se genera otra. El número de réplicas descartadas
#     se registra por coalición. Sin réplicas descartadas, el resultado es
#     idéntico al de sae::eblupFH() con la misma semilla (se verifica en la
#     sección 9b).
#   - B = 500, el mismo número de réplicas de la selección stepwise.
#
# SALIDAS (output/)
#   shapley_kicb2_m1_modelos.csv
#   shapley_kicb2_m1_ranking.csv
#   shapley_kicb2_m1.rds
# ==============================================================================


# ==============================================================================
# 0. PAQUETES
# ==============================================================================

if (!requireNamespace("sae", quietly = TRUE)) {
  stop("Instale el paquete 'sae'.")
}

if (!requireNamespace("TUvalues", quietly = TRUE)) {
  stop("Instale el paquete 'TUvalues'.")
}

library(sae)
library(TUvalues)
library(dplyr)
library(here)

select <- dplyr::select
filter <- dplyr::filter


# ==============================================================================
# 1. CONFIGURACIÓN
# ==============================================================================

ruta_out <- here("output")

B_BOOT    <- 500
MAXITER   <- 5000
PRECISION <- 1e-8
SEMILLA   <- 2906

# Límite de réplicas descartadas por coalición. Si se supera, el script se
# detiene: descartar demasiadas réplicas deja de ser un tratamiento marginal.
MAX_DESCARTES <- 50


# ==============================================================================
# 2. DATOS
# ==============================================================================

matriz <- readRDS(
  file.path(ruta_out, "matriz_sae_transformada_v2.rds")
)

datos <- as.data.frame(matriz)


# ==============================================================================
# 3. VARIABLES DE M1
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
# 4. VALIDACIONES
# ==============================================================================

vars_necesarias <- c("pobreza_monetaria", "varianza_pobreza", vars_m1)

faltantes <- setdiff(vars_necesarias, names(datos))

if (length(faltantes) > 0) {
  stop("Faltan variables en los datos: ", paste(faltantes, collapse = ", "))
}

stopifnot(
  "Existen valores faltantes." =
    !anyNA(datos[, vars_necesarias, drop = FALSE]),
  "pobreza_monetaria contiene valores no finitos." =
    all(is.finite(datos$pobreza_monetaria)),
  "La varianza de muestreo contiene valores no finitos." =
    all(is.finite(datos$varianza_pobreza)),
  "Existen varianzas de muestreo no positivas." =
    all(datos$varianza_pobreza > 0)
)


# ==============================================================================
# 5. ESTANDARIZACIÓN DE LAS COVARIABLES
# ------------------------------------------------------------------------------
# Reparametrización invariante: no cambia el KICb2 (ver encabezado).
# ==============================================================================

datos_std <- datos %>%
  mutate(across(all_of(vars_m1), ~ as.numeric(scale(.x))))

cat("\nNúmero de condición de X (modelo completo):\n")
cat("  sin estandarizar:",
    signif(kappa(model.matrix(
      as.formula(paste("pobreza_monetaria ~", paste(vars_m1, collapse = " + "))),
      datos), exact = TRUE), 4), "\n")
cat("  estandarizada:   ",
    signif(kappa(model.matrix(
      as.formula(paste("pobreza_monetaria ~", paste(vars_m1, collapse = " + "))),
      datos_std), exact = TRUE), 4), "\n")


# ==============================================================================
# 6. FÓRMULA
# ==============================================================================

formula_fh <- function(vars) {
  if (length(vars) == 0) return(pobreza_monetaria ~ 1)
  as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
}


# ==============================================================================
# 7. GENERAR LAS 63 COALICIONES
# ==============================================================================

coal     <- TUvalues::coalitions(p)
bin_coal <- as.matrix(coal$Binary)

all_comb <- lapply(
  2:nrow(bin_coal),
  function(r) unname(which(bin_coal[r, ] == 1))
)

n_modelos <- length(all_comb)

cat("\nNúmero de variables:", p, "\n")
cat("Modelos no vacíos:", n_modelos, "\n")

stopifnot(n_modelos == 2^p - 1)


# ==============================================================================
# 8. COMPROBAR RANGO DE LAS MATRICES DE DISEÑO
# ==============================================================================

for (j in seq_along(all_comb)) {
  
  vars_j <- vars_m1[all_comb[[j]]]
  X_j    <- model.matrix(formula_fh(vars_j), data = datos_std)
  
  if (qr(X_j)$rank < ncol(X_j)) {
    stop("La coalición ", j, " presenta una matriz de diseño singular: ",
         paste(vars_j, collapse = ", "))
  }
}

cat("Todas las matrices de diseño tienen rango completo.\n")


# ==============================================================================
# 9. AJUSTE Y KICb2
# ------------------------------------------------------------------------------
# sae::eblupFH() busca vardir por el NOMBRE escrito en la llamada
# (deparse(substitute(vardir))) dentro de 'data'. Por eso se escribe
# vardir = varianza_pobreza, sin 'datos$': con 'datos$varianza_pobreza' busca
# una columna con ese nombre literal, no la encuentra y el ajuste falla.
#
# El bucle bootstrap reproduce el de sae::eblupFH() (sae 1.3, R/eblupFH.R),
# con las mismas fórmulas y el mismo orden de generación de números
# aleatorios. La única diferencia es el tryCatch del reajuste.
# ==============================================================================

kicb2_boot <- function(vars, B = B_BOOT) {
  
  f       <- formula_fh(vars)
  X       <- model.matrix(f, data = datos_std)
  y       <- datos_std$pobreza_monetaria
  sigma2d <- datos_std$varianza_pobreza
  p_x     <- ncol(X)
  D       <- nrow(X)
  
  # Ajuste original, sin bootstrap
  aj <- sae::eblupFH(
    formula   = f,
    vardir    = varianza_pobreza,
    method    = "ML",
    MAXITER   = MAXITER,
    PRECISION = PRECISION,
    data      = datos_std
  )
  
  if (!isTRUE(aj$fit$convergence)) {
    stop("El ajuste original no convergió.")
  }
  
  loglike          <- unname(aj$fit$goodness["loglike"])
  min2loglike      <- -2 * loglike
  lambdahat        <- aj$fit$refvar
  betahat          <- matrix(aj$fit$estcoef[, "beta"], ncol = 1)
  Xbetahat         <- X %*% betahat
  lambdahatsigma2d <- lambdahat + sigma2d
  
  B1hatast  <- 0
  B3ast     <- 0
  B5ast     <- 0
  sum_logf_y    <- 0
  sum_logf_yast <- 0
  
  n_no_conv  <- 0
  n_singular <- 0
  
  b <- 1
  
  while (b <= B) {
    
    uastb <- sqrt(lambdahat) * matrix(rnorm(D, mean = 0, sd = 1), nrow = D, ncol = 1)
    eastb <- sqrt(sigma2d)   * matrix(rnorm(D, mean = 0, sd = 1), nrow = D, ncol = 1)
    yastb <- Xbetahat + uastb + eastb
    
    rb <- tryCatch(
      sae::eblupFH(yastb ~ X - 1, sigma2d, method = "ML",
                   MAXITER = MAXITER, PRECISION = PRECISION),
      error = function(e) NULL
    )
    
    if (is.null(rb)) {
      n_singular <- n_singular + 1
    } else if (!isTRUE(rb$fit$convergence)) {
      n_no_conv <- n_no_conv + 1
    }
    
    if (n_singular + n_no_conv > MAX_DESCARTES) {
      stop("Más de ", MAX_DESCARTES, " réplicas descartadas.")
    }
    
    # Réplica descartada: se genera otra sin avanzar el contador
    if (is.null(rb) || !isTRUE(rb$fit$convergence)) next
    
    betahatastb   <- matrix(rb$fit$estcoef[, "beta"], ncol = 1)
    lambdahatastb <- rb$fit$refvar
    
    Xbetahathatastb2     <- (X %*% (betahat - betahatastb))^2
    yastbXbetahatastb2   <- (yastb - X %*% betahatastb)^2
    lambdahatastbsigma2d <- lambdahatastb + sigma2d
    
    B1hatast <- B1hatast + sum((lambdahatsigma2d + Xbetahathatastb2 -
                                  yastbXbetahatastb2) / lambdahatastbsigma2d)
    
    logf <- (-0.5) * sum(log(2 * pi * lambdahatastbsigma2d) +
                           ((y - X %*% betahatastb)^2) / lambdahatastbsigma2d)
    
    sum_logf_y    <- sum_logf_y + logf
    sum_logf_yast <- sum_logf_yast + unname(rb$fit$goodness["loglike"])
    
    B3ast <- B3ast + sum((lambdahatastbsigma2d + Xbetahathatastb2) /
                           lambdahatsigma2d)
    B5ast <- B5ast + sum(log(lambdahatastbsigma2d) +
                           yastbXbetahatastb2 / lambdahatastbsigma2d)
    
    b <- b + 1
  }
  
  B2ast <- sum(log(lambdahatsigma2d)) + B3ast / B - B5ast / B
  AICc  <- min2loglike + B1hatast / B
  AICb2 <- as.vector(min2loglike - 4 / B * (sum_logf_y - loglike * B))
  
  list(
    loglike    = loglike,
    KIC        = min2loglike + 3 * (p_x + 1),
    KICc       = AICc + B2ast,
    KICb2      = AICb2 + B2ast,
    n_no_conv  = n_no_conv,
    n_singular = n_singular
  )
}


# ------------------------------------------------------------------------------
# 9b. Verificación: sin réplicas descartadas, la función reproduce a sae
# ------------------------------------------------------------------------------
# Se compara sobre el modelo completo. Si en esa corrida hubo réplicas
# descartadas, la comparación no aplica y solo se informa.

cat("\n============================================================\n")
cat("VERIFICACIÓN FRENTE A sae::eblupFH()\n")
cat("============================================================\n")

set.seed(SEMILLA)
ver_propia <- kicb2_boot(vars_m1)

if (ver_propia$n_singular + ver_propia$n_no_conv == 0) {
  
  set.seed(SEMILLA)
  ver_sae <- sae::eblupFH(
    formula   = formula_fh(vars_m1),
    vardir    = varianza_pobreza,
    method    = "ML",
    MAXITER   = MAXITER,
    PRECISION = PRECISION,
    B         = B_BOOT,
    data      = datos_std
  )
  
  dif_ver <- abs(ver_propia$KICb2 - unname(ver_sae$fit$goodness["KICb2"]))
  
  cat("KICb2 función propia:", ver_propia$KICb2, "\n")
  cat("KICb2 sae:           ", unname(ver_sae$fit$goodness["KICb2"]), "\n")
  cat("Diferencia:          ", format(dif_ver, scientific = TRUE), "\n")
  
  stopifnot("La función propia no reproduce el KICb2 de sae." = dif_ver < 1e-6)
  
} else {
  
  cat("El modelo completo tuvo réplicas descartadas;",
      "la comparación exacta con sae no aplica.\n")
}


# ------------------------------------------------------------------------------
# 9c. Las 63 coaliciones
# ------------------------------------------------------------------------------

resultados <- vector("list", n_modelos)

t0 <- Sys.time()

for (j in seq_along(all_comb)) {
  
  idx    <- all_comb[[j]]
  vars_j <- vars_m1[idx]
  
  cat(sprintf("Modelo %2d / %d | %s\n", j, n_modelos,
              paste(vars_j, collapse = " + ")))
  
  # Misma semilla en todas las coaliciones
  set.seed(SEMILLA)
  
  r <- tryCatch(kicb2_boot(vars_j),
                error = function(e) list(error = conditionMessage(e)))
  
  fila <- list(
    modelo       = j,
    clave        = paste(idx, collapse = "|"),
    variables    = paste(vars_j, collapse = " + "),
    n_variables  = length(vars_j),
    loglike      = NA_real_,
    KIC          = NA_real_,
    KICc         = NA_real_,
    KICb2        = NA_real_,
    n_singular   = NA_integer_,
    n_no_conv    = NA_integer_,
    mensaje      = NA_character_
  )
  
  if (!is.null(r$error)) {
    
    fila$mensaje <- r$error
    cat("  Error:", r$error, "\n")
    
  } else {
    
    fila$loglike    <- r$loglike
    fila$KIC        <- r$KIC
    fila$KICc       <- r$KICc
    fila$KICb2      <- r$KICb2
    fila$n_singular <- r$n_singular
    fila$n_no_conv  <- r$n_no_conv
    
    if (r$n_singular + r$n_no_conv > 0) {
      cat("  Réplicas descartadas: ", r$n_singular, " singulares, ",
          r$n_no_conv, " sin convergencia\n", sep = "")
    }
  }
  
  resultados[[j]] <- fila
}

cat("\nTiempo total:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    "minutos\n")


# ==============================================================================
# 10. TABLA DE RESULTADOS Y CONTROL
# ==============================================================================

tabla_modelos <- bind_rows(lapply(resultados, as.data.frame,
                                  stringsAsFactors = FALSE))

tabla_modelos$valido <- is.finite(tabla_modelos$KICb2)

cat("\n============================================================\n")
cat("CONTROL FINAL\n")
cat("============================================================\n")

print(tabla_modelos %>%
        select(modelo, n_variables, KIC, KICb2, n_singular, n_no_conv, valido),
      row.names = FALSE)

cat("\nRéplicas descartadas en total:",
    sum(tabla_modelos$n_singular, na.rm = TRUE), "singulares y",
    sum(tabla_modelos$n_no_conv, na.rm = TRUE), "sin convergencia, en",
    sum((tabla_modelos$n_singular + tabla_modelos$n_no_conv) > 0, na.rm = TRUE),
    "de", n_modelos, "coaliciones\n")

n_invalidos <- sum(!tabla_modelos$valido)

cat("\nModelos válidos:  ", sum(tabla_modelos$valido), "\n")
cat("Modelos inválidos:", n_invalidos, "\n")

if (n_invalidos > 0) {
  print(tabla_modelos %>% filter(!valido) %>% select(modelo, variables, mensaje),
        row.names = FALSE)
  stop("No se calcula Shapley: hay ", n_invalidos,
       " coaliciones sin KICb2 válido. No se eliminan en silencio.")
}


# ==============================================================================
# 11. COMPROBAR ORDEN CON TUvalues
# ==============================================================================

claves_tuv <- vapply(
  2:nrow(bin_coal),
  function(r) paste(which(bin_coal[r, ] == 1), collapse = "|"),
  character(1)
)

stopifnot(
  "El orden de las coaliciones no coincide con TUvalues." =
    identical(tabla_modelos$clave, claves_tuv)
)


# ==============================================================================
# 12. SHAPLEY EXACTO
# ==============================================================================

shapley <- TUvalues::shapley(as.numeric(tabla_modelos$KICb2),
                             method = "exact")

stopifnot(length(shapley) == p)

resultado_shapley <- data.frame(
  variable      = vars_m1,
  Shapley_KICb2 = as.numeric(shapley),
  stringsAsFactors = FALSE
) %>%
  arrange(Shapley_KICb2) %>%
  mutate(posicion = row_number()) %>%
  select(posicion, variable, Shapley_KICb2)

cat("\n============================================================\n")
cat("RANKING SHAPLEY-KICb2 - M1, ESCALA ORIGINAL\n")
cat("============================================================\n")

print(resultado_shapley, row.names = FALSE)


# ==============================================================================
# 13. PROPIEDAD DE EFICIENCIA
# ==============================================================================

KICb2_full <- tabla_modelos$KICb2[
  tabla_modelos$clave == paste(seq_len(p), collapse = "|")
]

suma_shapley <- sum(resultado_shapley$Shapley_KICb2)
diferencia   <- suma_shapley - KICb2_full

cat("\nKICb2 modelo completo:", KICb2_full, "\n")
cat("Suma Shapley:         ", suma_shapley, "\n")
cat("Diferencia:           ", format(diferencia, scientific = TRUE), "\n")

if (abs(diferencia) > 1e-8) {
  warning("La propiedad de eficiencia no se cumple.")
}


# ==============================================================================
# 14. GUARDAR
# ==============================================================================

write.csv(tabla_modelos,
          file.path(ruta_out, "shapley_kicb2_m1_modelos.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(resultado_shapley,
          file.path(ruta_out, "shapley_kicb2_m1_ranking.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(
  list(
    variables_m1          = vars_m1,
    tabla_modelos         = tabla_modelos,
    resultado_shapley     = resultado_shapley,
    B                     = B_BOOT,
    semilla               = SEMILLA,
    KICb2_full            = KICb2_full,
    suma_shapley          = suma_shapley,
    diferencia_eficiencia = diferencia,
    session               = sessionInfo()
  ),
  file.path(ruta_out, "shapley_kicb2_m1.rds")
)

cat("\nProceso terminado.\n")