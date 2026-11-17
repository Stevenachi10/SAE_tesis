# ==============================================================================
# 05_shapley_etapa4_completo.R
# ------------------------------------------------------------------------------
# ETAPA IV: contribución de las variables auxiliares (valores de Shapley).
# Un solo script que reemplaza a 05_1, 05_1b, 05_2 y 05_3:
#
#   A. Shapley-KICb2, escala ORIGINAL   (sae::eblupFH, ML, bootstrap B = 500)
#   B. Shapley-KICb2, escala LOGARÍTMICA (log(y), psi / y^2, igual que emdi)
#   C. Shapley-MSE y Shapley-CV          (emdi::fh: REML, log, bc_crude,
#                                         MSE analítico)
#   D. Resumen de los cuatro rankings
#   E. Auditoría de los resultados anteriores
#
# JUEGOS (M1 fijo; jugadores = sus 6 variables; 64 coaliciones con la vacía)
#   KICb2: v(S) = KICb2(S) - KICb2(vacía)          (negativo = favorece al criterio)
#   MSE:   v(S) = MSE_prom(vacía) - MSE_prom(S)    (positivo = reduce el MSE)
#   CV:    v(S) = CV_prom(vacía)  - CV_prom(S)     (positivo = reduce el CV)
#   En los tres, v(vacía) = 0 y la suma de los valores de Shapley es la
#   variación total entre el modelo con intercepto únicamente y M1.
#
# DECISIONES (KICb2)
#   - Covariables estandarizadas: no cambian la verosimilitud ni el KICb2, solo
#     el condicionamiento numérico (sin estandarizar, solve() falla en el
#     bootstrap de algunas coaliciones).
#   - Una sola semilla para todas las coaliciones (números aleatorios comunes):
#     las diferencias de KICb2 entre coaliciones no mezclan ruido Monte Carlo.
#   - Se reproduce el bucle bootstrap de sae::eblupFH (sae 1.3). Las réplicas
#     cuyo reajuste falla o no converge se descartan y se genera otra; el número
#     de descartes se registra por coalición. Sin descartes, el resultado es
#     idéntico al de sae::eblupFH con la misma semilla (se verifica en la
#     auditoría, prueba 6).
#   - Escala logarítmica: y* = log(y), psi* = psi / y^2 (método delta), la misma
#     transformación de emdi::fh(transformation = "log"). El criterio se
#     calcula bajo ML, mientras que la estimación final usa REML. El KICb2 en
#     escala log NO es comparable con el de la selección (escala original).
#
# DECISIONES (MSE y CV)
#   - Una sola corrida de emdi::fh por coalición alimenta MSE y CV.
#   - En emdi, tol y maxit solo afectan a los métodos espaciales y de
#     Ybarra-Lohr; se dejan por consistencia con 05_2 y 05_3.
#   - CV_i = sqrt(MSE_i) / theta_i * 100. Con log + bc_crude,
#     MSE_i = theta_i^2 * mse_log_i, de modo que CV_i = sqrt(mse_log_i) * 100 y
#     el MSE retransformado mezcla precisión con nivel (ver auditoría, prueba 4).
#   - Todos los municipios deben ser admisibles en todas las coaliciones.
#
# ENTRADAS (output/)
#   matriz_sae_transformada_v2.rds
#   final_estimaciones_municipales_m1_log.csv   (solo para la auditoría)
#
# SALIDAS (output/)
#   shapley_kicb2_m1_{modelos,ranking}.csv, shapley_kicb2_m1.rds          (original)
#   shapley_kicb2_m1_log_{modelos,ranking}.csv, shapley_kicb2_m1_log.rds  (log)
#   shapley_mse_m1_log_{modelos,ranking,mse_municipal}.csv, shapley_mse_m1_log.rds
#   shapley_cv_m1_log_{modelos,ranking,cv_municipal}.csv,   shapley_cv_m1_log.rds
#   shapley_etapa4_resumen.csv
#   auditoria_shapley_*.csv, auditoria_shapley_etapa4.rds
# ==============================================================================


# ==============================================================================
# 0. CONFIGURACIÓN
# ==============================================================================

ruta_out <- here::here("output")

B_BOOT        <- 500      # réplicas bootstrap del KICb2
MAXITER       <- 5000
PRECISION     <- 1e-8
SEMILLA       <- 2906
MAX_DESCARTES <- 50       # réplicas descartadas permitidas por coalición

# Auditoría
RUN_AUDITORIA        <- TRUE
AUD_N_SEMILLAS_MC    <- 20      # semillas para medir la variabilidad de KICb2(M1)
AUD_SEMILLAS_SHAPLEY <- 1:3     # semillas adicionales para el Shapley-KICb2
# (cada una = 64 coaliciones x B en cada escala;
#  usar integer(0) para omitir esta prueba)


# ==============================================================================
# 1. PAQUETES
# ==============================================================================

for (pk in c("sae", "emdi", "TUvalues", "here")) {
  if (!requireNamespace(pk, quietly = TRUE)) stop("Instale el paquete '", pk, "'.")
}

library(sae)
library(emdi)
library(TUvalues)
library(here)


# ==============================================================================
# 2. DATOS
# ==============================================================================

datos <- as.data.frame(
  readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))
)

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

vars_necesarias <- c("cod_mun", "pobreza_monetaria", "varianza_pobreza", vars_m1)
faltantes <- setdiff(vars_necesarias, names(datos))
if (length(faltantes) > 0) {
  stop("Faltan variables en los datos: ", paste(faltantes, collapse = ", "))
}

stopifnot(
  "Existen valores faltantes."                      = !anyNA(datos[, vars_necesarias]),
  "pobreza_monetaria contiene valores no finitos."  = all(is.finite(datos$pobreza_monetaria)),
  "pobreza_monetaria debe ser positiva (log)."      = all(datos$pobreza_monetaria > 0),
  "varianza_pobreza contiene valores no finitos."   = all(is.finite(datos$varianza_pobreza)),
  "Existen varianzas de muestreo no positivas."     = all(datos$varianza_pobreza > 0),
  "cod_mun duplicado."                              = !anyDuplicated(datos$cod_mun)
)

# Datos para el KICb2 en cada escala: covariables estandarizadas, respuesta y
# varianza de muestreo con nombres comunes (y_resp, psi_resp).
datos_std <- datos
datos_std[vars_m1] <- lapply(datos[vars_m1], function(x) as.numeric(scale(x)))

armar_escala <- function(escala) {
  d <- datos_std
  if (escala == "original") {
    d$y_resp   <- d$pobreza_monetaria
    d$psi_resp <- d$varianza_pobreza
  } else {
    d$y_resp   <- log(d$pobreza_monetaria)
    d$psi_resp <- d$varianza_pobreza / d$pobreza_monetaria^2
  }
  stopifnot(all(is.finite(d$y_resp)), all(is.finite(d$psi_resp)), all(d$psi_resp > 0))
  d
}

DAT <- list(original = armar_escala("original"), log = armar_escala("log"))


# ==============================================================================
# 3. COALICIONES
# ------------------------------------------------------------------------------
# 2^6 = 64 coaliciones, la vacía primero. TUvalues::shapley() recibe las 63 no
# vacías en este orden y asume v(vacía) = 0.
# ==============================================================================

coal     <- TUvalues::coalitions(p)
bin_coal <- as.matrix(coal$Binary)

n_modelos <- nrow(bin_coal)
stopifnot(n_modelos == 2^p)

all_comb <- lapply(seq_len(n_modelos), function(r) unname(which(bin_coal[r, ] == 1)))

claves_tuv <- vapply(
  all_comb,
  function(idx) if (length(idx) == 0) "" else paste(idx, collapse = "|"),
  character(1)
)

clave_full <- paste(seq_len(p), collapse = "|")
stopifnot(claves_tuv[1] == "", claves_tuv[n_modelos] == clave_full)

etiqueta <- function(vars) if (length(vars) == 0) "(intercepto)" else paste(vars, collapse = " + ")
etiquetas_modelos <- vapply(all_comb, function(idx) etiqueta(vars_m1[idx]), character(1))

formula_fh <- function(vars, resp = "pobreza_monetaria") {
  rhs <- if (length(vars) == 0) "1" else paste(vars, collapse = " + ")
  as.formula(paste(resp, "~", rhs))
}

# Shapley exacto de un juego con v(vacía) = 0. v64 está en el orden de all_comb.
shapley_exacto <- function(v64) {
  v64 <- unname(as.numeric(v64))
  stopifnot(length(v64) == n_modelos, abs(v64[1]) < 1e-12)
  phi <- as.numeric(TUvalues::shapley(v64[-1], method = "exact"))
  stopifnot(length(phi) == p)
  setNames(phi, vars_m1)
}


# ==============================================================================
# 4. KICb2 POR BOOTSTRAP (reproduce el bucle de sae::eblupFH, sae 1.3)
# ------------------------------------------------------------------------------
# sae::eblupFH() busca vardir por el NOMBRE escrito en la llamada
# (deparse(substitute(vardir))) dentro de 'data'. Por eso se escribe
# vardir = psi_resp, sin 'datos$'.
# ==============================================================================

kicb2_boot <- function(vars, escala, B = B_BOOT) {
  
  dat     <- DAT[[escala]]
  f       <- formula_fh(vars, "y_resp")
  X       <- model.matrix(f, data = dat)
  y       <- dat$y_resp
  sigma2d <- dat$psi_resp
  D       <- nrow(X)
  
  aj <- sae::eblupFH(
    formula   = f,
    vardir    = psi_resp,
    method    = "ML",
    MAXITER   = MAXITER,
    PRECISION = PRECISION,
    data      = dat
  )
  
  if (!isTRUE(aj$fit$convergence)) stop("El ajuste original no convergió.")
  
  loglike          <- unname(aj$fit$goodness["loglike"])
  min2loglike      <- -2 * loglike
  lambdahat        <- aj$fit$refvar
  betahat          <- matrix(aj$fit$estcoef[, "beta"], ncol = 1)
  Xbetahat         <- X %*% betahat
  lambdahatsigma2d <- lambdahat + sigma2d
  
  B3ast      <- 0
  B5ast      <- 0
  sum_logf_y <- 0
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
    
    logf <- (-0.5) * sum(log(2 * pi * lambdahatastbsigma2d) +
                           ((y - X %*% betahatastb)^2) / lambdahatastbsigma2d)
    
    sum_logf_y <- sum_logf_y + logf
    
    B3ast <- B3ast + sum((lambdahatastbsigma2d + Xbetahathatastb2) /
                           lambdahatsigma2d)
    B5ast <- B5ast + sum(log(lambdahatastbsigma2d) +
                           yastbXbetahatastb2 / lambdahatastbsigma2d)
    
    b <- b + 1
  }
  
  B2ast <- sum(log(lambdahatsigma2d)) + B3ast / B - B5ast / B
  AICb2 <- as.vector(min2loglike - 4 / B * (sum_logf_y - loglike * B))
  
  list(
    KICb2      = AICb2 + B2ast,
    n_singular = n_singular,
    n_no_conv  = n_no_conv
  )
}


# ==============================================================================
# 5. SHAPLEY-KICb2 PARA UNA ESCALA Y UNA SEMILLA
# ==============================================================================

shapley_kicb2 <- function(escala, semilla = SEMILLA, verbose = TRUE) {
  
  kic       <- rep(NA_real_, n_modelos)
  n_sing    <- integer(n_modelos)
  n_nconv   <- integer(n_modelos)
  
  for (j in seq_len(n_modelos)) {
    
    if (verbose) {
      cat(sprintf("[KICb2 %s] Modelo %2d / %d | %s\n",
                  escala, j, n_modelos, etiquetas_modelos[j]))
    }
    
    set.seed(semilla)   # misma semilla en todas las coaliciones
    
    r <- tryCatch(
      kicb2_boot(vars_m1[all_comb[[j]]], escala),
      error = function(e) {
        stop("Error en la coalición ", j, " (", etiquetas_modelos[j], ", escala ",
             escala, "): ", conditionMessage(e), call. = FALSE)
      }
    )
    
    kic[j]     <- r$KICb2
    n_sing[j]  <- r$n_singular
    n_nconv[j] <- r$n_no_conv
    
    if (verbose && r$n_singular + r$n_no_conv > 0) {
      cat("   Réplicas descartadas:", r$n_singular, "singulares,",
          r$n_no_conv, "sin convergencia\n")
    }
  }
  
  if (anyNA(kic) || any(!is.finite(kic))) {
    stop("Hay coaliciones sin KICb2 válido: no se calcula Shapley.")
  }
  
  v   <- kic - kic[1]                       # juego centrado en el modelo vacío
  phi <- shapley_exacto(v)
  
  dif_ef <- sum(phi) - v[n_modelos]
  if (abs(dif_ef) > 1e-8) warning("La propiedad de eficiencia no se cumple (KICb2, ", escala, ").")
  
  tabla <- data.frame(
    modelo      = seq_len(n_modelos),
    clave       = claves_tuv,
    variables   = etiquetas_modelos,
    n_variables = lengths(all_comb),
    KICb2       = kic,
    v_centrado  = v,
    n_singular  = n_sing,
    n_no_conv   = n_nconv,
    stringsAsFactors = FALSE
  )
  
  ranking <- data.frame(
    variable      = vars_m1,
    Shapley_KICb2 = as.numeric(phi),
    stringsAsFactors = FALSE
  )
  ranking <- ranking[order(ranking$Shapley_KICb2), ]    # más negativo = mayor aporte
  ranking$posicion      <- seq_len(p)
  ranking$participacion <- 100 * ranking$Shapley_KICb2 / sum(ranking$Shapley_KICb2)
  ranking <- ranking[, c("posicion", "variable", "Shapley_KICb2", "participacion")]
  rownames(ranking) <- NULL
  
  list(
    escala                = escala,
    semilla               = semilla,
    B                     = B_BOOT,
    variables_m1          = vars_m1,
    tabla_modelos         = tabla,
    resultado_shapley     = ranking,
    phi                   = phi,
    KICb2_vacia           = kic[1],
    KICb2_full            = kic[n_modelos],
    suma_shapley          = sum(phi),
    diferencia_eficiencia = dif_ef,
    descartes_total       = sum(n_sing + n_nconv),
    n_coaliciones_con_descartes = sum((n_sing + n_nconv) > 0)
  )
}

imprimir_kicb2 <- function(res) {
  cat("\n============================================================\n")
  cat("SHAPLEY-KICb2 - M1, ESCALA", toupper(res$escala), "\n")
  cat("============================================================\n")
  print(res$resultado_shapley, row.names = FALSE, digits = 5)
  cat("\nKICb2 vacío:", round(res$KICb2_vacia, 4),
      "| KICb2 M1:", round(res$KICb2_full, 4),
      "| variación:", round(res$KICb2_full - res$KICb2_vacia, 4), "\n")
  cat("Suma Shapley:", round(res$suma_shapley, 4),
      "| diferencia de eficiencia:", format(res$diferencia_eficiencia, scientific = TRUE), "\n")
  cat("Réplicas descartadas:", res$descartes_total, "de", n_modelos * B_BOOT,
      "en", res$n_coaliciones_con_descartes, "de", n_modelos, "coaliciones\n")
}

guardar_kicb2 <- function(res, sufijo) {
  write.csv(res$tabla_modelos,
            file.path(ruta_out, paste0("shapley_kicb2_m1", sufijo, "_modelos.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(res$resultado_shapley,
            file.path(ruta_out, paste0("shapley_kicb2_m1", sufijo, "_ranking.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
  saveRDS(c(res, list(session = sessionInfo())),
          file.path(ruta_out, paste0("shapley_kicb2_m1", sufijo, ".rds")))
}


# ==============================================================================
# A. SHAPLEY-KICb2, ESCALA ORIGINAL
# ==============================================================================

t0 <- Sys.time()
res_kic_orig <- shapley_kicb2("original")
imprimir_kicb2(res_kic_orig)
guardar_kicb2(res_kic_orig, "")


# ==============================================================================
# B. SHAPLEY-KICb2, ESCALA LOGARÍTMICA
# ==============================================================================

res_kic_log <- shapley_kicb2("log")
imprimir_kicb2(res_kic_log)
guardar_kicb2(res_kic_log, "_log")

cat("\nTiempo KICb2 (ambas escalas):",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutos\n")


# ==============================================================================
# C. SHAPLEY-MSE Y SHAPLEY-CV (emdi: REML, log, bc_crude, MSE analítico)
# ==============================================================================

ajustar_fh <- function(vars) {
  emdi::fh(
    fixed              = formula_fh(vars),
    vardir             = "varianza_pobreza",
    combined_data      = datos,
    domains            = "cod_mun",
    method             = "reml",
    transformation     = "log",
    backtransformation = "bc_crude",
    MSE                = TRUE,
    mse_type           = "analytical",
    tol                = 1e-8,
    maxit              = 1000
  )
}

extraer_metricas <- function(ajuste, etiqueta = "") {
  
  dominios_ajuste <- tryCatch(as.character(ajuste$ind$Domain), error = function(e) NULL)
  
  if (is.null(dominios_ajuste) ||
      !identical(dominios_ajuste, as.character(datos$cod_mun))) {
    stop("El orden de dominios no coincide en ", etiqueta)
  }
  
  columna_fh <- function(objeto) {
    if (is.null(objeto)) return(NULL)
    i <- grep("^FH$", names(objeto))
    if (length(i) == 0) return(NULL)
    as.numeric(objeto[[i[1]]])
  }
  
  pred <- columna_fh(ajuste$ind)
  mse  <- columna_fh(ajuste$MSE)
  
  if (is.null(pred) || length(pred) != nrow(datos)) {
    stop("No fue posible extraer la estimación FH de ", etiqueta)
  }
  if (is.null(mse) || length(mse) != nrow(datos)) {
    stop("No fue posible extraer el MSE de ", etiqueta)
  }
  
  admisible <- is.finite(pred) & pred > 0 & pred <= 1 & is.finite(mse) & mse >= 0
  
  cv <- rep(NA_real_, length(pred))
  cv[admisible] <- sqrt(mse[admisible]) / pred[admisible] * 100
  
  list(pred = pred, mse = mse, cv = cv, admisible = admisible)
}

cat("\n============================================================\n")
cat("SHAPLEY-MSE Y SHAPLEY-CV - M1 LOG (", n_modelos, " ajustes emdi)\n", sep = "")
cat("============================================================\n")

t1 <- Sys.time()

pred_mun <- matrix(NA_real_, nrow(datos), n_modelos)
mse_mun  <- matrix(NA_real_, nrow(datos), n_modelos)
cv_mun   <- matrix(NA_real_, nrow(datos), n_modelos)
colnames(pred_mun) <- colnames(mse_mun) <- colnames(cv_mun) <- paste0("modelo_", seq_len(n_modelos))

for (j in seq_len(n_modelos)) {
  
  cat(sprintf("[emdi] Modelo %2d / %d | %s\n", j, n_modelos, etiquetas_modelos[j]))
  
  ajuste_j <- tryCatch(
    ajustar_fh(vars_m1[all_comb[[j]]]),
    error = function(e) stop("Error en el modelo ", j, " (", etiquetas_modelos[j], "): ",
                             conditionMessage(e), call. = FALSE)
  )
  
  met_j <- extraer_metricas(ajuste_j, paste("modelo", j))
  
  if (!all(met_j$admisible)) {
    stop("El modelo ", j, " tiene ", sum(!met_j$admisible), " municipios no admisibles. ",
         "No se continúa para evitar comparar medias sobre conjuntos distintos.")
  }
  
  pred_mun[, j] <- met_j$pred
  mse_mun[, j]  <- met_j$mse
  cv_mun[, j]   <- met_j$cv
}

cat("Tiempo emdi:", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 2), "minutos\n")

tabla_pre <- data.frame(
  modelo       = seq_len(n_modelos),
  clave        = claves_tuv,
  variables    = etiquetas_modelos,
  n_variables  = lengths(all_comb),
  media_CV     = colMeans(cv_mun),
  mediana_CV   = apply(cv_mun, 2, median),
  p90_CV       = apply(cv_mun, 2, quantile, probs = 0.90),
  max_CV       = apply(cv_mun, 2, max),
  media_MSE    = colMeans(mse_mun),
  mediana_MSE  = apply(mse_mun, 2, median),
  p90_MSE      = apply(mse_mun, 2, quantile, probs = 0.90),
  max_MSE      = apply(mse_mun, 2, max),
  n_admisibles = nrow(datos),
  stringsAsFactors = FALSE,
  row.names = NULL
)

tabla_pre$reduccion_MSE <- tabla_pre$media_MSE[1] - tabla_pre$media_MSE
tabla_pre$valor_CV      <- tabla_pre$media_CV[1]  - tabla_pre$media_CV

# --- Shapley-MSE ---
phi_mse <- shapley_exacto(tabla_pre$reduccion_MSE)
red_mse <- tabla_pre$reduccion_MSE[n_modelos]

res_mse <- data.frame(variable = vars_m1, Shapley_MSE = as.numeric(phi_mse),
                      stringsAsFactors = FALSE)
res_mse <- res_mse[order(-res_mse$Shapley_MSE), ]
res_mse$posicion         <- seq_len(p)
res_mse$Shapley_MSE_pp2  <- res_mse$Shapley_MSE * 10000      # prop^2 -> pp^2
res_mse$participacion_pct <- 100 * res_mse$Shapley_MSE / red_mse
res_mse <- res_mse[, c("posicion", "variable", "Shapley_MSE", "Shapley_MSE_pp2", "participacion_pct")]
rownames(res_mse) <- NULL

# --- Shapley-CV ---
phi_cv <- shapley_exacto(tabla_pre$valor_CV)
red_cv <- tabla_pre$valor_CV[n_modelos]

res_cv <- data.frame(variable = vars_m1, Shapley_CV_pp = as.numeric(phi_cv),
                     stringsAsFactors = FALSE)
res_cv <- res_cv[order(-res_cv$Shapley_CV_pp), ]
res_cv$posicion          <- seq_len(p)
res_cv$participacion_pct <- 100 * res_cv$Shapley_CV_pp / red_cv
res_cv <- res_cv[, c("posicion", "variable", "Shapley_CV_pp", "participacion_pct")]
rownames(res_cv) <- NULL

dif_ef_mse <- sum(phi_mse) - red_mse
dif_ef_cv  <- sum(phi_cv)  - red_cv

if (abs(dif_ef_mse) > 1e-8) warning("La propiedad de eficiencia no se cumple (MSE).")
if (abs(dif_ef_cv)  > 1e-8) warning("La propiedad de eficiencia no se cumple (CV).")

cat("\n------------------------------------------------------------\n")
cat("SHAPLEY-MSE (pp^2)\n")
cat("------------------------------------------------------------\n")
cat("MSE vacío:", round(tabla_pre$media_MSE[1] * 1e4, 4), "pp^2 | MSE M1:",
    round(tabla_pre$media_MSE[n_modelos] * 1e4, 4), "pp^2 | reducción:",
    round(red_mse * 1e4, 4), "pp^2\n")
print(res_mse, row.names = FALSE, digits = 5)
cat("Diferencia de eficiencia:", format(dif_ef_mse, scientific = TRUE), "\n")

cat("\n------------------------------------------------------------\n")
cat("SHAPLEY-CV (pp)\n")
cat("------------------------------------------------------------\n")
cat("CV vacío:", round(tabla_pre$media_CV[1], 4), "% | CV M1:",
    round(tabla_pre$media_CV[n_modelos], 4), "% | reducción:", round(red_cv, 4), "pp\n")
print(res_cv, row.names = FALSE, digits = 5)
cat("Diferencia de eficiencia:", format(dif_ef_cv, scientific = TRUE), "\n")

# --- Guardar ---
tab_mse <- tabla_pre[, c("modelo", "clave", "variables", "n_variables", "media_MSE",
                         "mediana_MSE", "p90_MSE", "max_MSE", "n_admisibles", "reduccion_MSE")]
tab_cv  <- tabla_pre[, c("modelo", "clave", "variables", "n_variables", "media_CV",
                         "mediana_CV", "p90_CV", "max_CV", "media_MSE", "mediana_MSE",
                         "n_admisibles", "valor_CV")]

write.csv(tab_mse, file.path(ruta_out, "shapley_mse_m1_log_modelos.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(res_mse, file.path(ruta_out, "shapley_mse_m1_log_ranking.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(cod_mun = datos$cod_mun, mse_mun, check.names = FALSE),
          file.path(ruta_out, "shapley_mse_m1_log_mse_municipal.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(variables_m1 = vars_m1, tabla_modelos = tab_mse, resultado_shapley = res_mse,
             mse_municipal = mse_mun, MSE_vacio = tabla_pre$media_MSE[1],
             MSE_full = tabla_pre$media_MSE[n_modelos], reduccion_total_MSE = red_mse,
             suma_shapley = sum(phi_mse), diferencia_eficiencia = dif_ef_mse),
        file.path(ruta_out, "shapley_mse_m1_log.rds"))

write.csv(tab_cv, file.path(ruta_out, "shapley_cv_m1_log_modelos.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(res_cv, file.path(ruta_out, "shapley_cv_m1_log_ranking.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(cod_mun = datos$cod_mun, cv_mun, check.names = FALSE),
          file.path(ruta_out, "shapley_cv_m1_log_cv_municipal.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(variables_m1 = vars_m1, tabla_modelos = tab_cv, resultado_shapley = res_cv,
             cv_municipal = cv_mun, CV_vacio = tabla_pre$media_CV[1],
             CV_full = tabla_pre$media_CV[n_modelos], ganancia_total_CV = red_cv,
             suma_shapley = sum(phi_cv), diferencia_eficiencia = dif_ef_cv),
        file.path(ruta_out, "shapley_cv_m1_log.rds"))


# ==============================================================================
# D. RESUMEN DE LOS CUATRO RANKINGS
# ==============================================================================

f_rank <- function(res, col) {
  res[match(vars_m1, res$variable), c(col, "posicion")]
}

kio <- f_rank(res_kic_orig$resultado_shapley, "Shapley_KICb2")
kil <- f_rank(res_kic_log$resultado_shapley,  "Shapley_KICb2")
msr <- f_rank(res_mse, "Shapley_MSE_pp2")
cvr <- f_rank(res_cv,  "Shapley_CV_pp")

pct <- function(x) 100 * x / sum(x)

resumen <- data.frame(
  variable          = vars_m1,
  KICb2_orig        = kio$Shapley_KICb2, pct_KICb2_orig = pct(kio$Shapley_KICb2), pos_KICb2_orig = kio$posicion,
  KICb2_log         = kil$Shapley_KICb2, pct_KICb2_log  = pct(kil$Shapley_KICb2), pos_KICb2_log  = kil$posicion,
  MSE_pp2           = msr$Shapley_MSE_pp2, pct_MSE      = pct(msr$Shapley_MSE_pp2), pos_MSE      = msr$posicion,
  CV_pp             = cvr$Shapley_CV_pp,   pct_CV       = pct(cvr$Shapley_CV_pp),   pos_CV       = cvr$posicion,
  stringsAsFactors = FALSE
)
resumen <- resumen[order(resumen$pos_KICb2_orig), ]

cat("\n============================================================\n")
cat("RESUMEN DE LOS CUATRO SHAPLEY (orden: KICb2 en escala original)\n")
cat("============================================================\n")
print(resumen, row.names = FALSE, digits = 4)

write.csv(resumen, file.path(ruta_out, "shapley_etapa4_resumen.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# E. AUDITORÍA
# ------------------------------------------------------------------------------
# 1. El ajuste completo (coalición 64) coincide con la estimación final (06_0).
# 2. Réplica independiente de emdi::fh en R base (log, REML, bc_crude, MSE
#    analítico) frente a los 64 ajustes de emdi.
# 3. Ruido del optimizador en Shapley-CV: tolerancia por defecto vs estricta.
# 4. Efecto del nivel en el MSE: precisión a nivel fijo vs nivel de la predicción.
# 5. Variabilidad Monte Carlo de KICb2(M1) entre semillas, en ambas escalas.
# 6. La función bootstrap propia reproduce sae::eblupFH(B = 500).
# 7. Estabilidad del Shapley-KICb2 entre semillas, en ambas escalas.
# 8. Eficiencia de los cuatro juegos.
# ==============================================================================

if (RUN_AUDITORIA) {
  
  aud <- list()
  
  cat("\n\n############################################################\n")
  cat("AUDITORÍA\n")
  cat("############################################################\n")
  
  # ----------------------------------------------------------------------------
  # 1. Ajuste completo frente a la estimación final
  # ----------------------------------------------------------------------------
  cat("\n=== 1. M1 completo (emdi, coalición 64) frente a la estimación final ===\n")
  
  ruta_fin <- file.path(ruta_out, "final_estimaciones_municipales_m1_log.csv")
  
  if (file.exists(ruta_fin)) {
    
    fin <- read.csv(ruta_fin, colClasses = c(cod_mun = "character"))
    o   <- match(as.character(datos$cod_mun), fin$cod_mun)
    stopifnot(!anyNA(o))
    
    d1 <- c(estimacion_pp = max(abs(pred_mun[, n_modelos] * 100 - fin$fh_pct[o])),
            CV_pp         = max(abs(cv_mun[, n_modelos] - fin$cv_fh_pct[o])),
            MSE           = max(abs(mse_mun[, n_modelos] - fin$mse_fh[o])))
    
    cat("Diferencia máxima absoluta con el archivo final:\n")
    print(signif(d1, 3))
    cat("CV medio M1 (Shapley) / archivo final:",
        round(mean(cv_mun[, n_modelos]), 5), "/", round(mean(fin$cv_fh_pct), 5), "\n")
    aud$ajuste_final <- d1
    
  } else {
    cat("No se encontró final_estimaciones_municipales_m1_log.csv: prueba omitida.\n")
  }
  
  # ----------------------------------------------------------------------------
  # 2. Réplica independiente de emdi::fh en R base
  # ----------------------------------------------------------------------------
  cat("\n=== 2. Réplica en R base de emdi::fh frente a los 64 ajustes ===\n")
  
  X_de <- function(vars) {
    if (length(vars) == 0) matrix(1, nrow(datos), 1)
    else model.matrix(as.formula(paste("~", paste(vars, collapse = "+"))), datos)
  }
  
  # log(y), psi/y^2; REML por optimize() sobre c(0, var(log y)); theta = exp(mu +
  # 0.5 mse_log) (bc_crude); mse_log de Prasad-Rao (g1 + g2 + 2 g3); MSE =
  # theta^2 mse_log.
  fh_log_rep <- function(vars, tol = NULL) {
    y  <- datos$pobreza_monetaria; psi <- datos$varianza_pobreza
    yl <- log(y); pl <- psi / y^2
    X  <- X_de(vars); m <- length(yl)
    ll <- function(A) {
      Vi <- 1 / (A + pl); XtVi <- t(Vi * X); M <- XtVi %*% X; Q <- solve(M)
      Py <- Vi * yl - Vi * as.vector(X %*% (Q %*% (XtVi %*% yl)))
      -(m / 2) * log(2 * pi) + 0.5 * sum(log(Vi)) -
        0.5 * as.numeric(determinant(M)$modulus) - 0.5 * sum(yl * Py)
    }
    intervalo <- c(0, var(yl))
    A <- if (is.null(tol)) optimize(ll, intervalo, maximum = TRUE)$maximum
    else optimize(ll, intervalo, maximum = TRUE, tol = tol)$maximum
    Vi <- 1 / (A + pl); XtVi <- t(Vi * X); Q <- solve(XtVi %*% X)
    b   <- Q %*% (XtVi %*% yl); ajuste <- as.vector(X %*% b)
    mu  <- ajuste + A / (A + pl) * (yl - ajuste)
    Bd  <- pl / (A + pl)
    mse_log <- pl * (1 - Bd) + Bd^2 * rowSums((X %*% Q) * X) +
      2 * Bd^2 * (2 / sum(Vi^2)) / (A + pl)
    theta <- exp(mu + 0.5 * mse_log)
    list(A = A, theta = theta, mse = theta^2 * mse_log, cv = sqrt(mse_log) * 100)
  }
  
  rep_def <- lapply(all_comb, function(idx) fh_log_rep(vars_m1[idx]))
  
  d2 <- c(
    estimacion = max(sapply(seq_len(n_modelos), function(j) max(abs(rep_def[[j]]$theta - pred_mun[, j])))),
    CV_pp      = max(sapply(seq_len(n_modelos), function(j) max(abs(rep_def[[j]]$cv    - cv_mun[, j])))),
    MSE_rel    = max(sapply(seq_len(n_modelos), function(j) max(abs(rep_def[[j]]$mse / mse_mun[, j] - 1))))
  )
  cat("Diferencia máxima entre la réplica y emdi (64 coaliciones x 84 municipios):\n")
  print(signif(d2, 3))
  aud$replica_emdi <- d2
  
  # ----------------------------------------------------------------------------
  # 3. Ruido del optimizador (Shapley-CV)
  # ----------------------------------------------------------------------------
  cat("\n=== 3. Ruido del optimizador en Shapley-CV ===\n")
  
  rep_estricto <- lapply(all_comb, function(idx) fh_log_rep(vars_m1[idx], tol = 1e-12))
  
  cv_bar_def <- sapply(rep_def, function(r) mean(r$cv))
  cv_bar_str <- sapply(rep_estricto, function(r) mean(r$cv))
  
  dif_A  <- max(abs(sapply(seq_len(n_modelos), function(j) rep_def[[j]]$A - rep_estricto[[j]]$A)) /
                  sapply(rep_estricto, function(r) r$A))
  dif_cv <- max(abs(cv_bar_def - cv_bar_str))
  
  phi_cv_estricto <- shapley_exacto(cv_bar_str[1] - cv_bar_str)
  
  cat("Diferencia relativa máxima de sigma2_u:", signif(dif_A, 3),
      "| diferencia máxima del CV medio:", signif(dif_cv, 3), "pp\n")
  cat("Diferencia máxima en Shapley-CV (estricto vs emdi):",
      signif(max(abs(phi_cv_estricto - phi_cv)), 3), "pp\n")
  aud$optimizador <- c(dif_rel_sigma2u = dif_A, dif_CV_medio = dif_cv,
                       dif_shapley_CV = max(abs(phi_cv_estricto - phi_cv)))
  
  # ----------------------------------------------------------------------------
  # 4. Efecto del nivel en el MSE retransformado
  # ----------------------------------------------------------------------------
  cat("\n=== 4. Precisión y nivel en el Shapley-MSE ===\n")
  
  th_full <- pred_mun[, n_modelos]
  v_fijo  <- apply(cv_mun, 2, function(cv) mean(th_full^2 * (cv / 100)^2))  # MSE a nivel fijo
  
  phi_mse_fijo <- shapley_exacto(v_fijo[1] - v_fijo)
  
  tab4 <- rbind(
    `Shapley-CV (pp)`                = phi_cv,
    `Shapley-MSE (pp^2)`             = phi_mse * 1e4,
    `Shapley-MSE, nivel fijo (pp^2)` = phi_mse_fijo * 1e4
  )
  print(signif(tab4, 4))
  cat("Posiciones (1 = mayor aporte):\n")
  print(t(apply(-tab4, 1, rank)))
  
  tot  <- red_mse * 1e4
  prec <- (v_fijo[1] - v_fijo[n_modelos]) * 1e4
  cat("\nReducción total del MSE:", round(tot, 3), "pp^2 = precisión a nivel fijo",
      round(prec, 3), "+ efecto del nivel", round(tot - prec, 3), "\n")
  cat("DE entre municipios de la predicción: vacío", round(sd(pred_mun[, 1]), 4),
      "| M1", round(sd(th_full), 4), "\n")
  
  write.csv(tab4, file.path(ruta_out, "auditoria_shapley_cv_mse.csv"))
  aud$nivel <- list(tabla = tab4, total = tot, precision_nivel_fijo = prec)
  
  # ----------------------------------------------------------------------------
  # 5. Variabilidad Monte Carlo de KICb2(M1)
  # ----------------------------------------------------------------------------
  cat("\n=== 5. Variabilidad Monte Carlo de KICb2(M1), B =", B_BOOT, "===\n")
  
  mc <- sapply(c(original = "original", log = "log"), function(esc) {
    vapply(seq_len(AUD_N_SEMILLAS_MC), function(s) {
      set.seed(s)
      kicb2_boot(vars_m1, esc)$KICb2
    }, numeric(1))
  })
  
  res_mc <- t(apply(mc, 2, function(x) c(media = mean(x), DE = sd(x), minimo = min(x), maximo = max(x))))
  cat(AUD_N_SEMILLAS_MC, "semillas:\n"); print(round(res_mc, 3))
  
  kfull <- c(original = res_kic_orig$KICb2_full, log = res_kic_log$KICb2_full)
  cat("KICb2(M1) de la corrida principal (semilla", SEMILLA, "):", round(kfull, 3), "\n")
  cat("Distancia a la media, en desviaciones estándar:",
      round((kfull - res_mc[, "media"]) / res_mc[, "DE"], 2), "\n")
  
  write.csv(cbind(semilla = seq_len(AUD_N_SEMILLAS_MC), mc),
            file.path(ruta_out, "auditoria_shapley_kicb2_mc.csv"), row.names = FALSE)
  aud$mc_kicb2 <- list(valores = mc, resumen = res_mc)
  
  # ----------------------------------------------------------------------------
  # 6. La función propia reproduce sae::eblupFH(B = 500)
  # ----------------------------------------------------------------------------
  cat("\n=== 6. Función bootstrap propia frente a sae::eblupFH(B =", B_BOOT, ") ===\n")
  
  # sae::eblupFH(B = ) se detiene si una réplica produce un sistema singular, de
  # modo que la comparación solo es posible en una semilla sin réplicas
  # descartadas: se usa la semilla principal y, si tuvo descartes, la primera
  # semilla de 1 a 20 que no los tenga.
  verif_sae <- sapply(c(original = "original", log = "log"), function(esc) {
    
    for (s_ver in c(SEMILLA, 1:20)) {
      set.seed(s_ver)
      mio <- kicb2_boot(vars_m1, esc)
      if (mio$n_singular + mio$n_no_conv == 0) break
    }
    
    if (mio$n_singular + mio$n_no_conv > 0) {
      return(c(semilla = NA_real_, propia = mio$KICb2, sae = NA_real_, dif = NA_real_))
    }
    
    set.seed(s_ver)
    ref <- sae::eblupFH(
      formula   = formula_fh(vars_m1, "y_resp"),
      vardir    = psi_resp,
      method    = "ML",
      MAXITER   = MAXITER,
      PRECISION = PRECISION,
      B         = B_BOOT,
      data      = DAT[[esc]]
    )
    
    k_sae <- unname(ref$fit$goodness["KICb2"])
    c(semilla = s_ver, propia = mio$KICb2, sae = k_sae, dif = abs(mio$KICb2 - k_sae))
  })
  
  print(signif(t(verif_sae), 8))
  if (any(!is.na(verif_sae["dif", ]) & verif_sae["dif", ] > 1e-6, na.rm = TRUE)) {
    warning("La función propia no reproduce el KICb2 de sae::eblupFH.")
  }
  aud$verif_sae <- verif_sae
  
  # ----------------------------------------------------------------------------
  # 7. Estabilidad del Shapley-KICb2 entre semillas
  # ----------------------------------------------------------------------------
  if (length(AUD_SEMILLAS_SHAPLEY) > 0) {
    
    cat("\n=== 7. Estabilidad del Shapley-KICb2 entre semillas ===\n")
    
    rob <- list()
    
    for (esc in c("original", "log")) {
      
      principal <- if (esc == "original") res_kic_orig else res_kic_log
      
      filas <- lapply(AUD_SEMILLAS_SHAPLEY, function(s) {
        cat("  escala", esc, "| semilla", s, "\n")
        r <- shapley_kicb2(esc, s, verbose = FALSE)
        c(r$phi, kicb2_vacia = r$KICb2_vacia, kicb2_m1 = r$KICb2_full)
      })
      
      M <- rbind(c(principal$phi, kicb2_vacia = principal$KICb2_vacia, kicb2_m1 = principal$KICb2_full),
                 do.call(rbind, filas))
      rownames(M) <- c(SEMILLA, AUD_SEMILLAS_SHAPLEY)
      
      phi_m <- M[, vars_m1, drop = FALSE]
      medias <- colMeans(phi_m)
      
      cat("\nShapley-KICb2, escala", esc, ", por semilla:\n"); print(round(M, 2))
      cat("Resumen entre", nrow(M), "semillas:\n")
      print(round(rbind(media = medias, DE = apply(phi_m, 2, sd),
                        `% del total` = 100 * medias / sum(medias)), 2))
      
      rangos <- t(apply(phi_m, 1, rank))     # 1 = mayor aporte (más negativo)
      cat("Posiciones por semilla (1 = mayor aporte):\n"); print(rangos)
      cat("¿El orden se mantiene en todas las semillas? ",
          all(apply(rangos, 2, function(x) length(unique(x)) == 1)), "\n", sep = "")
      
      write.csv(cbind(semilla = rownames(M), as.data.frame(M)),
                file.path(ruta_out, paste0("auditoria_shapley_kicb2_semillas_", esc, ".csv")),
                row.names = FALSE)
      
      rob[[esc]] <- list(valores = M, rangos = rangos)
    }
    
    aud$robustez_semillas <- rob
  }
  
  # ----------------------------------------------------------------------------
  # 8. Eficiencia de los cuatro juegos
  # ----------------------------------------------------------------------------
  cat("\n=== 8. Propiedad de eficiencia ===\n")
  
  ef <- data.frame(
    juego = c("KICb2 original", "KICb2 log", "MSE (pp^2)", "CV (pp)"),
    variacion_total = c(res_kic_orig$KICb2_full - res_kic_orig$KICb2_vacia,
                        res_kic_log$KICb2_full  - res_kic_log$KICb2_vacia,
                        red_mse * 1e4, red_cv),
    suma_shapley    = c(res_kic_orig$suma_shapley, res_kic_log$suma_shapley,
                        sum(phi_mse) * 1e4, sum(phi_cv)),
    stringsAsFactors = FALSE
  )
  ef$diferencia <- ef$suma_shapley - ef$variacion_total
  print(ef, row.names = FALSE, digits = 6)
  write.csv(ef, file.path(ruta_out, "auditoria_shapley_eficiencia.csv"), row.names = FALSE)
  aud$eficiencia <- ef
  
  saveRDS(c(aud, list(session = sessionInfo())),
          file.path(ruta_out, "auditoria_shapley_etapa4.rds"))
}

cat("\nProceso terminado.\n")