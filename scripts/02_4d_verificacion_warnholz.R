# ==============================================================================
# 02_4d_verificacion_warnholz.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.4d. Verificación de la representación pseudolineal y del error
# cuadrático medio por pseudolinealización de REBLUP y REBLUP-BC
#
# Requiere output/robusto.rds, producido por 02_4_robusto.R.
#
# ------------------------------------------------------------------------------
# FUENTES
# ------------------------------------------------------------------------------
# Warnholz (2016) adapta al modelo de nivel de área la estimación del error
# cuadrático medio basada en la representación pseudolineal propuesta por
# Chambers, Chandra y Tzavidis (2011) para modelos de nivel de unidad. Esa
# adaptación es la que implementa saeRobust, del que emdi toma la estimación
# robusta. [Verificar en la tesis la numeración de las ecuaciones antes de
# citarlas.]
#
# ------------------------------------------------------------------------------
# QUE SE VERIFICA
# ------------------------------------------------------------------------------
# El predictor se escribe como W y, con W una matriz de pesos que depende de
# los datos. El error cuadrático medio por pseudolinealización trata esos pesos
# como fijos:
#
#     ECM = diag(A G A' + W Ve W') + (W X beta - X beta)^2,    A = W - I,
#     G = Z Vu Z'
#
# Pruebas:
#   0. El reajuste con saeRobust reproduce el de emdi.
#   1. W y reproduce REBLUP y Wbc y reproduce REBLUP-BC.
#   2. W está calibrada sobre X (W X = X). De la construcción de saeRobust,
#      W = XA + ZB(I - XA) con A X = I, de modo que W X = X y el término de
#      sesgo del REBLUP es nulo.
#   3. La fórmula matricial reproduce el ECM entregado por emdi.
#   4. saeRobust::mse() reproduce el ECM entregado por emdi.
#
# En REBLUP-BC las filas de los municipios recortados se reducen a un único
# peso sobre el directo, 1 -/+ c*ee_i/y_i, que no suma uno. La calibración no
# se exige para Wbc: ese desvío es el que genera el término de sesgo, y se
# reporta de forma descriptiva.
#
# ------------------------------------------------------------------------------
# SOBRE EL REAJUSTE
# ------------------------------------------------------------------------------
# El objeto que devuelve fh() no conserva el objeto interno de saeRobust. Se
# reajusta el modelo con saeRobust::rfh() usando la misma llamada que emdi
# (eblup_robust.R) y se comprueba que coeficientes, varianza del efecto de área
# y predicciones coincidan antes de continuar.
#
# Se usan tol y maxit guardados en robusto.rds, los mismos del ajuste en
# 02_4_robusto.R. Los ajustes verificados son los de pseudolinealización, que en
# 02_4_robusto.R se estiman como análisis de sensibilidad del ECM.
#
# saeRobust::mse() calcula los pesos de REBLUP-BC con c = 1 con independencia
# de mult_constant. La verificación solo es válida con mult_constant = 1.
#
# ------------------------------------------------------------------------------
# SALIDAS (output/)
# ------------------------------------------------------------------------------
#   warnholz_verificacion_resumen.csv
#   warnholz_verificacion_detalle.csv
# ==============================================================================

library(dplyr)
library(here)
library(Matrix)
library(saeRobust)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")

rob <- readRDS(file.path(ruta_out, "robusto.rds"))

datos <- rob$datos
K     <- rob$k_principal
C     <- rob$mult_constant

stopifnot("La verificación supone mult_constant = 1" = isTRUE(all.equal(C, 1)))

# Mismos valores del ajuste en 02_4_robusto.R
TOL_AJUSTE <- rob$tol
MAXIT      <- rob$maxit

y   <- as.numeric(datos$pobreza_monetaria)
psi <- as.numeric(datos$varianza_pobreza)
n   <- length(y)

# Tolerancias de la verificación
TOL_REL <- 1e-8    # diferencias relativas
TOL_ABS <- 1e-10   # diferencias absolutas en escala de la proporción

# W y reproduce el REBLUP solo hasta la precisión de la convergencia, y W X = X
# hasta el error numérico de invertir matrices. Con tol = 1e-8 ambas
# diferencias quedan del orden de 1e-7.
TOL_CONV <- 1e-6

formula_de <- function(vars) {
  as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
}

dif_rel <- function(a, b) max(abs(a / b - 1))


# ==============================================================================
# CICLO POR ESPECIFICACION
# ==============================================================================

resumen <- list()
detalle <- list()

for (nm in names(rob$especificaciones)) {
  
  aj_rb <- rob$ajustes_robustos_pseudo[[paste(nm, "reblup", format(K), sep = "_")]]
  aj_bc <- rob$ajustes_robustos_pseudo[[paste(nm, "reblupbc", format(K), sep = "_")]]
  
  for (aj in list(aj_rb, aj_bc)) {
    stopifnot("Orden de dominios distinto de 'datos'" =
                identical(as.character(aj$ind$Domain),
                          as.character(datos$cod_mun)))
  }
  
  reblup_emdi   <- as.numeric(aj_rb$ind$FH)
  reblupbc_emdi <- as.numeric(aj_bc$ind$FH)
  mse_rb_emdi   <- as.numeric(aj_rb$MSE$FH)
  mse_bc_emdi   <- as.numeric(aj_bc$MSE$FH)
  
  # ---------------------------------------------------------------------------
  # 0. Reajuste con saeRobust, misma llamada que emdi
  # ---------------------------------------------------------------------------
  fit <- saeRobust::rfh(
    formula_de(rob$especificaciones[[nm]]),
    data       = datos,
    samplingVar = "varianza_pobreza",
    k          = K,
    tol        = TOL_AJUSTE,
    maxIter    = MAXIT)
  
  p0_coef <- dif_rel(as.numeric(fit$coefficients),
                     as.numeric(aj_rb$model$coefficients$coefficients))
  p0_s2u  <- dif_rel(as.numeric(fit$variance)[1],
                     as.numeric(aj_rb$model$variance))
  p0_pred <- max(abs(as.numeric(fit$reblup) - reblup_emdi))
  
  stopifnot(
    "El reajuste no reproduce los coeficientes de emdi" = p0_coef < TOL_REL,
    "El reajuste no reproduce sigma2_u de emdi"         = p0_s2u  < TOL_REL,
    "El reajuste no reproduce el REBLUP de emdi"        = p0_pred < TOL_ABS)
  
  # ---------------------------------------------------------------------------
  # Matrices de varianza y pesos
  # ---------------------------------------------------------------------------
  matV <- saeRobust::variance(fit)
  
  Z  <- as.matrix(matV$Z())
  Vu <- as.matrix(matV$Vu())
  Ve <- as.matrix(matV$Ve())
  G  <- Z %*% Vu %*% t(Z)
  
  # Método S3 de saeRobust para el genérico weights() de stats
  wt  <- weights(fit, c = C)
  W   <- as.matrix(wt$W)
  Wbc <- as.matrix(wt$Wbc)
  
  X  <- as.matrix(fit$x)
  xb <- as.numeric(X %*% fit$coefficients)
  I_n <- diag(n)
  
  # Municipios en que REBLUP-BC recorta: su fila difiere de la del REBLUP
  recortado <- rowSums(abs(Wbc - W)) > 1e-12
  
  # ---------------------------------------------------------------------------
  # 1. Representación pseudolineal
  # ---------------------------------------------------------------------------
  p1_rb <- max(abs(as.numeric(W %*% y)   - reblup_emdi))
  p1_bc <- max(abs(as.numeric(Wbc %*% y) - reblupbc_emdi))
  
  # ---------------------------------------------------------------------------
  # 2. Calibración de W sobre X
  # ---------------------------------------------------------------------------
  p2_calibracion_W <- max(abs(W %*% X - X))
  suma_W   <- rowSums(W)
  suma_Wbc <- rowSums(Wbc)
  
  # ---------------------------------------------------------------------------
  # 3. Fórmula matricial del ECM
  # ---------------------------------------------------------------------------
  descomponer <- function(Wm) {
    Am <- Wm - I_n
    list(var_u = diag(Am %*% G %*% t(Am)),
         var_e = diag(Wm %*% Ve %*% t(Wm)),
         sesgo = as.numeric(Wm %*% xb - xb))
  }
  
  d_rb <- descomponer(W)
  d_bc <- descomponer(Wbc)
  
  mse_rb_formula <- d_rb$var_u + d_rb$var_e + d_rb$sesgo^2
  mse_bc_formula <- d_bc$var_u + d_bc$var_e + d_bc$sesgo^2
  
  p3_rb <- dif_rel(mse_rb_formula, mse_rb_emdi)
  p3_bc <- dif_rel(mse_bc_formula, mse_bc_emdi)
  
  # ---------------------------------------------------------------------------
  # 4. saeRobust::mse() frente a emdi
  # ---------------------------------------------------------------------------
  mse_rb_pkg <- saeRobust::mse(fit, type = "pseudo", predType = "reblup")$pseudo
  mse_bc_pkg <- saeRobust::mse(fit, type = "pseudo", predType = "reblupbc")$pseudobc
  
  p4_rb <- dif_rel(as.numeric(mse_rb_pkg), mse_rb_emdi)
  p4_bc <- dif_rel(as.numeric(mse_bc_pkg), mse_bc_emdi)
  
  # ---------------------------------------------------------------------------
  # Peso del término de sesgo en el ECM de REBLUP-BC
  # ---------------------------------------------------------------------------
  cuota_sesgo_bc <- d_bc$sesgo^2 / mse_bc_formula
  
  resumen[[nm]] <- data.frame(
    especificacion           = nm,
    p0_coef_dif_rel          = signif(p0_coef, 3),
    p0_sigma2u_dif_rel       = signif(p0_s2u, 3),
    p0_reblup_dif_abs        = signif(p0_pred, 3),
    p1_reblup_dif_abs        = signif(p1_rb, 3),
    p1_reblupbc_dif_abs      = signif(p1_bc, 3),
    p2_calibracion_W_dif     = signif(p2_calibracion_W, 3),
    p3_ecm_reblup_dif_rel    = signif(p3_rb, 3),
    p3_ecm_reblupbc_dif_rel  = signif(p3_bc, 3),
    p4_ecm_reblup_dif_rel    = signif(p4_rb, 3),
    p4_ecm_reblupbc_dif_rel  = signif(p4_bc, 3),
    n_recortados             = sum(recortado),
    sesgo_max_reblup         = signif(max(abs(d_rb$sesgo)), 3),
    cuota_sesgo_bc_recortados    = if (any(recortado))
      round(median(cuota_sesgo_bc[recortado]), 3)
    else NA_real_,
    cuota_sesgo_bc_no_recortados = round(median(cuota_sesgo_bc[!recortado]), 3),
    stringsAsFactors = FALSE)
  
  detalle[[nm]] <- data.frame(
    especificacion   = nm,
    cod_mun          = datos$cod_mun,
    Municipio        = datos$Municipio,
    directo          = y,
    reblup           = reblup_emdi,
    reblupbc         = reblupbc_emdi,
    recortado        = recortado,
    suma_pesos_W     = suma_W,
    suma_pesos_Wbc   = suma_Wbc,
    mse_reblup_emdi  = mse_rb_emdi,
    mse_reblup_form  = mse_rb_formula,
    var_u_reblup     = d_rb$var_u,
    var_e_reblup     = d_rb$var_e,
    sesgo2_reblup    = d_rb$sesgo^2,
    mse_reblupbc_emdi = mse_bc_emdi,
    mse_reblupbc_form = mse_bc_formula,
    var_u_reblupbc   = d_bc$var_u,
    var_e_reblupbc   = d_bc$var_e,
    sesgo2_reblupbc  = d_bc$sesgo^2,
    cuota_sesgo_bc   = cuota_sesgo_bc,
    stringsAsFactors = FALSE)
}

resumen <- do.call(rbind, resumen)
detalle <- do.call(rbind, detalle)


# ==============================================================================
# RESULTADOS
# ==============================================================================

marca <- function(ok) if (isTRUE(ok)) "OK" else "REVISAR"

cat("\n============================================================\n")
cat("VERIFICACION DE LA PSEUDOLINEALIZACION (WARNHOLZ, 2016)\n")
cat("============================================================\n\n")

cat("0. Reajuste con saeRobust reproduce emdi:          OK (verificado con stopifnot)\n\n")

cat("1. Representación pseudolineal\n")
cat("   REBLUP     W y   :", marca(all(resumen$p1_reblup_dif_abs   < TOL_CONV)), "\n")
cat("   REBLUP-BC  Wbc y :", marca(all(resumen$p1_reblupbc_dif_abs < TOL_ABS)), "\n\n")

cat("2. Calibración W X = X (sesgo del REBLUP nulo)\n")
cat("   REBLUP           :", marca(all(resumen$p2_calibracion_W_dif < TOL_CONV)), "\n\n")

cat("3. Fórmula matricial frente a emdi\n")
cat("   REBLUP           :", marca(all(resumen$p3_ecm_reblup_dif_rel   < TOL_REL)), "\n")
cat("   REBLUP-BC        :", marca(all(resumen$p3_ecm_reblupbc_dif_rel < TOL_REL)), "\n\n")

cat("4. saeRobust::mse() frente a emdi\n")
cat("   REBLUP           :", marca(all(resumen$p4_ecm_reblup_dif_rel   < TOL_REL)), "\n")
cat("   REBLUP-BC        :", marca(all(resumen$p4_ecm_reblupbc_dif_rel < TOL_REL)), "\n\n")

print(resumen, row.names = FALSE)

cat("\nPeso del sesgo al cuadrado en el ECM de REBLUP-BC (mediana):\n")
print(resumen %>%
        select(especificacion, n_recortados,
               cuota_sesgo_bc_recortados, cuota_sesgo_bc_no_recortados),
      row.names = FALSE)

cat("\nMunicipios recortados, descomposición del ECM de REBLUP-BC:\n\n")
print(detalle %>%
        filter(recortado) %>%
        transmute(especificacion, Municipio, cod_mun,
                  suma_pesos_Wbc = round(suma_pesos_Wbc, 4),
                  var_u  = signif(var_u_reblupbc, 3),
                  var_e  = signif(var_e_reblupbc, 3),
                  sesgo2 = signif(sesgo2_reblupbc, 3),
                  cuota_sesgo = round(cuota_sesgo_bc, 3)) %>%
        arrange(especificacion, desc(cuota_sesgo)),
      row.names = FALSE)


# ==============================================================================
# GUARDADO
# ==============================================================================

write.csv(resumen, file.path(ruta_out, "warnholz_verificacion_resumen.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(detalle, file.path(ruta_out, "warnholz_verificacion_detalle.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n============================================================\n")
cat("VERIFICACION COMPLETADA\n")
cat("============================================================\n")