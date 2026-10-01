# ==============================================================================
# 03_1_precision_diagnosticos.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA III.2 — Diagnósticos de consistencia frente al estimador directo
#
# Diagnósticos de Brown, Chambers, Heady y Heasman (2001), implementados
# principalmente mediante el paquete SAEval.
#
# Las estimaciones y errores cuadráticos medios provienen de los objetos
# de emdi ajustados en la Etapa II.
# ==============================================================================


# ==============================================================================
# 0. PAQUETES Y PARÁMETROS
# ==============================================================================

library(emdi)
library(SAEval)
library(dplyr)
library(here)
library(car)
library(lmtest)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")

HAY_SANDWICH <- requireNamespace("sandwich", quietly = TRUE)

if (!HAY_SANDWICH) {
  cat("AVISO: el paquete 'sandwich' no está instalado.\n")
  cat("El contraste HC3 no será calculado.\n")
}


# ------------------------------------------------------------------------------
# Especificaciones espaciales que se evalúan
# ------------------------------------------------------------------------------
# Sus factores de encogimiento se recuperan de espacial.rds en el bloque 1, una
# vez cargado el archivo.
# ------------------------------------------------------------------------------

ESPACIALES <- c("M1", "M2")


# ==============================================================================
# 1. CARGA Y DATOS DE REFERENCIA
# ==============================================================================

variantes <- readRDS(
  file.path(ruta_out, "modelado_variantes.rds")
)

espacial <- readRDS(
  file.path(ruta_out, "espacial.rds")
)

robusto <- readRDS(
  file.path(ruta_out, "robusto.rds")
)

datos <- variantes$datos


# ------------------------------------------------------------------------------
# Verificación de alineación
# ------------------------------------------------------------------------------

stopifnot(
  
  "cod_mun difiere entre fuentes" =
    identical(
      datos$cod_mun,
      robusto$datos$cod_mun
    ),
  
  "pobreza_monetaria difiere entre fuentes" =
    isTRUE(
      all.equal(
        datos$pobreza_monetaria,
        robusto$datos$pobreza_monetaria
      )
    ),
  
  "varianza_pobreza difiere entre fuentes" =
    isTRUE(
      all.equal(
        datos$varianza_pobreza,
        robusto$datos$varianza_pobreza
      )
    )
)


DIRECTO <- datos$pobreza_monetaria
PSI     <- datos$varianza_pobreza
n_dom   <- nrow(datos)


cat("\n============================================================\n")
cat("REFERENCIA\n")
cat("============================================================\n")

cat("Dominios:", n_dom, "\n")

cat(
  "Directo: media",
  round(mean(DIRECTO), 4),
  "| recorrido [",
  round(min(DIRECTO), 4),
  ",",
  round(max(DIRECTO), 4),
  "]\n"
)

cat(
  "CV directo (%): media",
  round(
    mean(
      sqrt(PSI) / DIRECTO * 100
    ),
    3
  ),
  "| máximo",
  round(
    max(
      sqrt(PSI) / DIRECTO * 100
    ),
    3
  ),
  "\n"
)


# ------------------------------------------------------------------------------
# Factor de encogimiento de las especificaciones espaciales
# ------------------------------------------------------------------------------
# emdi no lo devuelve en el caso espacial (ver la nota del bloque 3). La Etapa
# II lo calcula a partir de los componentes de varianza estimados y lo guarda en
# espacial.rds, por dominio, en dos formas:
#
#   gamma_esp        diagonal de la matriz de encogimiento, peso sobre el
#                    estimador directo del propio dominio
#   gamma_esp_fila   suma de cada fila, peso total sobre los estimadores
#                    directos, incluidos los de los dominios vecinos
#
# Si el archivo procede de una corrida anterior que no las guardaba, las listas
# quedan vacías y gamma se registra como no disponible para esa familia.
# ------------------------------------------------------------------------------

if (!is.null(espacial$cod_mun)) {
  
  stopifnot(
    "El orden de dominios de espacial.rds no coincide con 'datos'" =
      identical(
        as.character(espacial$cod_mun),
        as.character(datos$cod_mun)
      )
  )
}

GAMMA_ESPACIAL      <- list()
GAMMA_ESPACIAL_FILA <- list()

for (nm in ESPACIALES) {
  
  GAMMA_ESPACIAL[[paste("Espacial", nm)]] <-
    espacial$gamma_esp[[nm]]
  
  GAMMA_ESPACIAL_FILA[[paste("Espacial", nm)]] <-
    espacial$gamma_esp_fila[[nm]]
}

cat(
  "Encogimiento espacial recuperado de espacial.rds:",
  length(GAMMA_ESPACIAL),
  "de",
  length(ESPACIALES),
  "especificaciones\n"
)

if (length(GAMMA_ESPACIAL) < length(ESPACIALES)) {
  
  cat(
    "Volver a ejecutar 02_3_espacial.R para generarlo.\n"
  )
}


# ==============================================================================
# 2. AJUSTES A EVALUAR
# ------------------------------------------------------------------------------
# M1 original se conserva únicamente como referencia de contraste.
# ==============================================================================

AJUSTES <- list()


agregar <- function(
    etiqueta,
    objeto,
    familia,
    escala,
    mse_type,
    papel) {
  
  if (is.null(objeto)) {
    
    cat(
      "AUSENTE:",
      etiqueta,
      "\n"
    )
    
    return(invisible(NULL))
  }
  
  AJUSTES[[etiqueta]] <<- list(
    aj = objeto,
    familia = familia,
    escala = escala,
    mse_type = mse_type,
    papel = papel
  )
}


# ------------------------------------------------------------------------------
# FH estándar
# ------------------------------------------------------------------------------

agregar(
  "M1 original",
  variantes$ajustes[["M1_reml_original"]],
  "base",
  "original",
  "analitico",
  "referencia"
)

agregar(
  "M1 log",
  variantes$ajustes[["M1_reml_log"]],
  "base",
  "log",
  "analitico",
  "candidata"
)

agregar(
  "M1 logit",
  variantes$ajustes[["M1_reml_logit"]],
  "base",
  "logit",
  "bootstrap",
  "candidata"
)


# ------------------------------------------------------------------------------
# FH espacial
# ------------------------------------------------------------------------------

for (nm in c("M1", "M2")) {
  
  agregar(
    paste("Espacial", nm),
    espacial$ajustes_esp[[nm]],
    "espacial",
    "original",
    "analitico",
    "candidata"
  )
}


# ------------------------------------------------------------------------------
# FH robusto
# ------------------------------------------------------------------------------

for (nm in c("M1", "M2", "M3")) {
  
  for (met in c("reblup", "reblupbc")) {
    
    clave <- paste(
      nm,
      met,
      format(1.345),
      sep = "_"
    )
    
    aj_rob <- robusto$ajustes_robustos[[clave]]
    
    # El estimador del ECM se lee del ajuste. Desde la redefinicion de la
    # etapa II, ajustes_robustos contiene los ajustes con bootstrap
    # parametrico (B = 500); los de pseudolinealizacion quedan en
    # ajustes_robustos_pseudo y solo se usan como sensibilidad.
    mse_rob <- if (is.null(aj_rob)) NA_character_ else aj_rob$method$MSE_method
    
    if (!is.null(aj_rob) && !identical(mse_rob, "bootstrap")) {
      stop("El ajuste ", clave, " no usa ECM bootstrap: ", mse_rob,
           ". Volver a correr 02_4_robusto.R.")
    }
    
    agregar(
      paste(met, nm),
      aj_rob,
      "robusta",
      "original",
      "bootstrap",
      "candidata"
    )
  }
}


cat(
  "\nAjustes disponibles:",
  length(AJUSTES),
  "\n"
)


# ==============================================================================
# 3. EXTRACCIÓN
# ------------------------------------------------------------------------------
# Se extraen la estimación FH, su error cuadrático medio y el factor de
# encogimiento.
#
# DISPONIBILIDAD DE GAMMA.
# Verificado sobre el código fuente de emdi: el objeto devuelto por fh() incluye
# model$gamma únicamente cuando correlation = "no", tanto en la rama estándar
# como en la rama con transformación. La rama espacial llama a eblup_SFH(), que
# no devuelve gamma, y la rama robusta construye su lista model sin ese
# elemento. De las once especificaciones evaluadas, por tanto, solo las tres de
# la familia base lo entregan.
#
# Para las restantes se procede de forma distinta según la familia.
#
#   ESPACIAL. No se calcula aquí. Bajo el modelo SAR la varianza del efecto de
#   área es Omega = sigma_u^2 [(I - rho W)'(I - rho W)]^(-1), y el encogimiento
#   corresponde a la matriz Omega (Omega + Psi)^(-1), cuya diagonal no coincide
#   con la expresión escalar sigma_u^2 / (sigma_u^2 + psi). Trasladar esa
#   expresión produciría una cantidad ajena al modelo. El valor lo aporta la
#   Etapa II a través de GAMMA_ESPACIAL.
#
#   ROBUSTA. Se calcula mediante la expresión escalar y se marca como
#   aproximado. El estimador robusto conserva la estructura de covarianza del
#   modelo de trabajo, Var(u) = sigma_u^2 I y Var(e) = Psi, y modifica las
#   ecuaciones de estimación mediante una función de influencia acotada. La
#   expresión describe, en consecuencia, el encogimiento del modelo ajustado y
#   no el encogimiento que el predictor robusto aplica efectivamente.
#
# DOS FORMAS DEL ENCOGIMIENTO.
# Se reportan gamma_medio, resumen de la diagonal de la matriz de encogimiento,
# y gamma_fila_medio, resumen de la suma de sus filas. En el modelo no espacial
# esa matriz es diagonal y las dos cantidades coinciden, de modo que la segunda
# columna solo aporta información en la familia espacial: allí la diagonal es el
# peso sobre el estimador directo del propio dominio y la suma de la fila es el
# peso total sobre los estimadores directos, que es la cantidad comparable con
# el gamma del modelo estándar.
#
# La columna gamma_fuente registra la procedencia en cada especificación.
# ==============================================================================

extraer <- function(aj) {
  
  est <- tryCatch(
    as.numeric(aj$ind$FH),
    error = function(e) NULL
  )
  
  mse <- tryCatch(
    as.numeric(aj$MSE$FH),
    error = function(e) NULL
  )
  
  
  if (
    is.null(est) ||
    is.null(mse)
  ) {
    return(NULL)
  }
  
  
  if (
    length(est) != n_dom ||
    length(mse) != n_dom
  ) {
    return(NULL)
  }
  
  
  # --------------------------------------------------------------------------
  # Verificación de dominios
  # --------------------------------------------------------------------------
  
  dom <- tryCatch(
    as.character(aj$ind$Domain),
    error = function(e) NULL
  )
  
  if (
    is.null(dom) ||
    !identical(
      dom,
      as.character(datos$cod_mun)
    )
  ) {
    
    stop(
      "Los dominios del ajuste no están alineados con 'datos'."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Gamma entregado por emdi
  # --------------------------------------------------------------------------
  
  gam <- tryCatch({
    
    g <- aj$model$gamma
    
    if (
      !is.null(g) &&
      "Gamma" %in% names(g)
    ) {
      
      as.numeric(g$Gamma)
      
    } else {
      
      NULL
    }
    
  }, error = function(e) NULL)
  
  
  if (
    !is.null(gam) &&
    length(gam) == n_dom
  ) {
    
    return(
      list(
        est = est,
        mse = mse,
        gamma = gam,
        gamma_fila = gam,
        gamma_fuente = "emdi"
      )
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Gamma no entregado por emdi
  #
  # Se distingue la familia a partir del propio objeto, mediante
  # model$correlation, y no a partir de la etiqueta asignada en el bloque 2.
  # --------------------------------------------------------------------------
  
  es_espacial <- tryCatch(
    identical(
      as.character(aj$model$correlation),
      "spatial"
    ),
    error = function(e) FALSE
  )
  
  
  if (es_espacial) {
    
    return(
      list(
        est = est,
        mse = mse,
        gamma = rep(NA_real_, n_dom),
        gamma_fila = rep(NA_real_, n_dom),
        gamma_fuente = "no disponible"
      )
    )
  }
  
  
  # Familia robusta: expresión escalar del modelo de trabajo
  s2u <- tryCatch({
    
    v <- aj$model$variance
    
    if (
      !is.null(names(v)) &&
      "variance" %in% names(v)
    ) {
      as.numeric(v[["variance"]])
    } else {
      as.numeric(v)[1]
    }
    
  }, error = function(e) NA_real_)
  
  
  if (
    is.finite(s2u) &&
    s2u > 0
  ) {
    
    # La matriz de encogimiento es diagonal, de modo que la suma de cada fila
    # coincide con el elemento diagonal correspondiente.
    gam_ap <- s2u / (s2u + PSI)
    
    list(
      est = est,
      mse = mse,
      gamma = gam_ap,
      gamma_fila = gam_ap,
      gamma_fuente = "aproximado"
    )
    
  } else {
    
    list(
      est = est,
      mse = mse,
      gamma = rep(NA_real_, n_dom),
      gamma_fila = rep(NA_real_, n_dom),
      gamma_fuente = "no disponible"
    )
  }
}



# ==============================================================================
# 4. BASE COMÚN PARA SAEval
# ==============================================================================

eval_data <- data.frame(
  cod_mun = datos$cod_mun,
  Municipio = datos$Municipio,
  DIRECTO = DIRECTO,
  PSI = PSI,
  stringsAsFactors = FALSE
)

metadata <- list()
gammas   <- list()


for (nm in names(AJUSTES)) {
  
  ex <- extraer(
    AJUSTES[[nm]]$aj
  )
  
  
  if (is.null(ex)) {
    
    cat(
      "NO EVALUABLE:",
      nm,
      "(estimación o MSE ausente)\n"
    )
    
    next
  }
  
  
  # --------------------------------------------------------------------------
  # Gamma espacial suministrado por la Etapa II
  # --------------------------------------------------------------------------
  
  if (nm %in% names(GAMMA_ESPACIAL)) {
    
    g_ext <- as.numeric(
      GAMMA_ESPACIAL[[nm]]
    )
    
    if (length(g_ext) == n_dom) {
      
      ex$gamma <- g_ext
      ex$gamma_fuente <- "etapa II"
      
      # La suma de filas solo se adopta si viene completa. No se sustituye por
      # la diagonal: en el modelo espacial son cantidades distintas.
      f_ext <- as.numeric(
        GAMMA_ESPACIAL_FILA[[nm]]
      )
      
      if (length(f_ext) == n_dom) {
        ex$gamma_fila <- f_ext
      }
      
    } else {
      
      cat(
        "AVISO: GAMMA_ESPACIAL de",
        nm,
        "no tiene longitud",
        n_dom,
        "- se ignora\n"
      )
    }
  }
  
  
  col <- make.names(nm)
  
  
  eval_data[[col]] <-
    ex$est
  
  eval_data[[paste0("MSE_", col)]] <-
    ex$mse
  
  
  gammas[[col]] <-
    ex$gamma
  
  
  # --------------------------------------------------------------------------
  # Valores fuera de la escala de una proporción
  # --------------------------------------------------------------------------
  
  fuera <- sum(
    ex$est <= 0 |
      ex$est >= 1,
    na.rm = TRUE
  )
  
  
  # --------------------------------------------------------------------------
  # Brown W manual
  # --------------------------------------------------------------------------
  
  den <- PSI + ex$mse
  
  ok <- is.finite(ex$est) &
    is.finite(ex$mse) &
    is.finite(den) &
    den > 0
  
  
  W_manual <- sum(
    (DIRECTO[ok] - ex$est[ok])^2 /
      den[ok]
  )
  
  
  # --------------------------------------------------------------------------
  # Grados de libertad del estadístico W
  # --------------------------------------------------------------------------
  # Se usa la convención de SAEval::gof(), m - 1, porque es la implementación
  # con la que se calcula W en este script. Brown et al. (2001) contrastan W
  # con m grados de libertad; esa convención se reporta como comparación en el
  # bloque 6.
  # --------------------------------------------------------------------------
  
  m_dom <- sum(ok)
  
  gl_w <- as.integer(m_dom - 1L)
  
  chi2_95 <- qchisq(0.95, df = gl_w)
  
  p_w <- pchisq(
    W_manual,
    df = gl_w,
    lower.tail = FALSE
  )
  
  
  
  # --------------------------------------------------------------------------
  # Resumen de gamma
  #
  # Si gamma no está disponible, queda NA.
  # --------------------------------------------------------------------------
  
  gamma_disponible <-
    any(
      is.finite(ex$gamma)
    )
  
  
  if (gamma_disponible) {
    
    gamma_medio <-
      mean(
        ex$gamma,
        na.rm = TRUE
      )
    
    gamma_min <-
      min(
        ex$gamma,
        na.rm = TRUE
      )
    
    gamma_max <-
      max(
        ex$gamma,
        na.rm = TRUE
      )
    
  } else {
    
    gamma_medio <- NA_real_
    gamma_min   <- NA_real_
    gamma_max   <- NA_real_
  }
  
  
  gamma_fila_medio <-
    if (any(is.finite(ex$gamma_fila))) {
      mean(
        ex$gamma_fila,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }
  
  
  # --------------------------------------------------------------------------
  # Metadata
  # --------------------------------------------------------------------------
  
  metadata[[nm]] <- data.frame(
    
    ajuste = nm,
    
    variable = col,
    
    papel =
      AJUSTES[[nm]]$papel,
    
    familia =
      AJUSTES[[nm]]$familia,
    
    escala =
      AJUSTES[[nm]]$escala,
    
    mse_type =
      AJUSTES[[nm]]$mse_type,
    
    gamma_fuente =
      ex$gamma_fuente,
    
    gamma_medio =
      round(
        gamma_medio,
        4
      ),
    
    gamma_min =
      round(
        gamma_min,
        4
      ),
    
    gamma_max =
      round(
        gamma_max,
        4
      ),
    
    gamma_fila_medio =
      round(
        gamma_fila_medio,
        4
      ),
    
    W_manual =
      round(
        W_manual,
        2
      ),
    
    m_dominios =
      m_dom,
    
    W_sobre_m =
      round(
        W_manual / m_dom,
        4
      ),
    
    gl_w =
      gl_w,
    
    chi2_95_w =
      round(
        chi2_95,
        2
      ),
    
    p_w =
      round(
        p_w,
        4
      ),
    
    n_fuera_rango =
      fuera,
    
    stringsAsFactors = FALSE
  )
}


metadata <- do.call(
  rbind,
  metadata
)


modelo_cols <-
  metadata$variable

mse_cols <-
  paste0(
    "MSE_",
    modelo_cols
  )


# ------------------------------------------------------------------------------
# Verificaciones
# ------------------------------------------------------------------------------

stopifnot(
  
  "Faltan columnas de estimación" =
    all(
      modelo_cols %in%
        names(eval_data)
    ),
  
  "Faltan columnas de MSE" =
    all(
      mse_cols %in%
        names(eval_data)
    ),
  
  "Hay valores no finitos en el directo" =
    all(
      is.finite(
        eval_data$DIRECTO
      )
    ),
  
  "Hay valores no finitos en la varianza" =
    all(
      is.finite(
        eval_data$PSI
      )
    )
)


if (
  any(
    metadata$n_fuera_rango > 0
  )
) {
  
  cat(
    "\n-- Predicciones fuera del intervalo unitario --\n"
  )
  
  print(
    metadata %>%
      filter(
        n_fuera_rango > 0
      ) %>%
      select(
        ajuste,
        n_fuera_rango
      ),
    row.names = FALSE
  )
}


# ------------------------------------------------------------------------------
# Procedencia de gamma
# ------------------------------------------------------------------------------

cat(
  "\n-- Procedencia del factor de encogimiento --\n"
)

print(
  table(
    metadata$familia,
    metadata$gamma_fuente
  )
)


formula_sae <- as.formula(
  paste(
    "~",
    paste(
      modelo_cols,
      collapse = " + "
    )
  )
)


formula_mse <- as.formula(
  paste(
    "~",
    paste(
      mse_cols,
      collapse = " + "
    )
  )
)


# ==============================================================================
# 5. DIAGNÓSTICO DE SESGO — SAEval
# ------------------------------------------------------------------------------
# Se conserva como referencia de la implementación del paquete.
#
# IMPORTANTE:
# SAEval plantea la regresión en una dirección que no coincide con la utilizada
# en Brown et al. (2001). Por ello, esta salida NO es la utilizada como
# diagnóstico principal de la tesis.
#
# scatterplot = FALSE evita readline() durante la ejecución del script.
# ==============================================================================

res_bias <- SAEval::bias(
  data = eval_data,
  dir = ~DIRECTO,
  sae = formula_sae,
  scatterplot = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "DIAGNÓSTICO DE SESGO — SAEval\n"
)

cat(
  "============================================================\n\n"
)

print(
  res_bias$output1,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# output2
# ------------------------------------------------------------------------------

if (
  !is.null(
    res_bias$output2
  )
) {
  
  cat(
    "\n------------------------------------------------------------\n"
  )
  
  cat(
    "SAEval — transformación raíz cuadrada\n"
  )
  
  cat(
    "------------------------------------------------------------\n\n"
  )
  
  print(
    res_bias$output2,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "\noutput2 no se construyó: la regla interna de SAEval no activó la\n"
  )
  
  cat(
    "transformación. Esto no constituye evidencia de homocedasticidad.\n"
  )
  
  cat(
    "El contraste complementario ordenado por PSI se reporta en 5b.\n"
  )
}


# ==============================================================================
# 5b. DIAGNÓSTICO DE SESGO — DIRECCIÓN DE BROWN
# ------------------------------------------------------------------------------
# Regresión:
#
#     DIRECTO = beta0 + beta1 * MODELO + error
#
# Hipótesis conjunta:
#
#     H0: beta0 = 0 y beta1 = 1
#
# Se reporta:
#   - contraste clásico
#   - contraste HC3 como sensibilidad
#   - GQ ordenado por PSI como diagnóstico complementario
#
# El GQ NO se utiliza como prueba aislada de validez/invalidez del FH.
# ==============================================================================

brown_bias <- do.call(
  
  rbind,
  
  lapply(
    
    modelo_cols,
    
    function(col) {
      
      
      d <- data.frame(
        
        y = eval_data$DIRECTO,
        
        x = eval_data[[col]],
        
        psi = eval_data$PSI
      )
      
      
      d <- d[
        is.finite(d$y) &
          is.finite(d$x) &
          is.finite(d$psi),
        ,
        drop = FALSE
      ]
      
      
      # ------------------------------------------------------------------------
      # Modelo auxiliar
      # ------------------------------------------------------------------------
      
      fit <- lm(
        y ~ x,
        data = d
      )
      
      smr <- summary(fit)
      
      
      # ------------------------------------------------------------------------
      # Contraste clásico
      # ------------------------------------------------------------------------
      
      lh <- car::linearHypothesis(
        fit,
        diag(2),
        c(0, 1)
      )
      
      p_ols <-
        lh[["Pr(>F)"]][2]
      
      
      # ------------------------------------------------------------------------
      # Contraste HC3
      # ------------------------------------------------------------------------
      
      if (HAY_SANDWICH) {
        
        lh3 <-
          car::linearHypothesis(
            fit,
            diag(2),
            c(0, 1),
            vcov. =
              sandwich::vcovHC(
                fit,
                type = "HC3"
              )
          )
        
        p_hc3 <-
          lh3[["Pr(>F)"]][2]
        
      } else {
        
        p_hc3 <-
          NA_real_
      }
      
      
      # ------------------------------------------------------------------------
      # Goldfeld-Quandt ordenado por PSI
      #
      # Interpretación:
      # diagnóstico complementario de la heterogeneidad de la dispersión
      # asociada al ordenamiento por la varianza de muestreo.
      # ------------------------------------------------------------------------
      
      gq <- tryCatch(
        
        lmtest::gqtest(
          fit,
          order.by = d$psi,
          alternative = "two.sided"
        )$p.value,
        
        error = function(e)
          NA_real_
      )
      
      
      # ------------------------------------------------------------------------
      # Resultado
      # ------------------------------------------------------------------------
      
      data.frame(
        
        variable = col,
        
        b0 =
          round(
            unname(
              coef(fit)[1]
            ),
            4
          ),
        
        b1 =
          round(
            unname(
              coef(fit)[2]
            ),
            4
          ),
        
        ee_b1 =
          round(
            smr$coefficients[2, 2],
            4
          ),
        
        R2 =
          round(
            smr$r.squared,
            4
          ),
        
        p_conjunto =
          round(
            p_ols,
            4
          ),
        
        p_conjunto_hc3 =
          round(
            as.numeric(
              p_hc3
            ),
            4
          ),
        
        gq_p_psi =
          round(
            as.numeric(
              gq
            ),
            4
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)


brown_bias <-
  metadata %>%
  select(
    ajuste,
    variable,
    papel,
    familia
  ) %>%
  left_join(
    brown_bias,
    by = "variable"
  )


cat(
  "\n============================================================\n"
)

cat(
  "DIAGNÓSTICO DE SESGO — DIRECCIÓN DE BROWN\n"
)

cat(
  "H0 conjunta: beta0 = 0 y beta1 = 1\n"
)

cat(
  "============================================================\n\n"
)

print(
  brown_bias %>%
    select(
      -variable
    ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Modelos sensibles a HC3
# ------------------------------------------------------------------------------

if (HAY_SANDWICH) {
  
  discrepa <-
    with(
      brown_bias,
      
      is.finite(p_conjunto) &
        is.finite(p_conjunto_hc3) &
        
        (
          (p_conjunto < 0.05) !=
            (p_conjunto_hc3 < 0.05)
        )
    )
  
  
  if (
    any(discrepa)
  ) {
    
    cat(
      "\n-- La conclusión cambia al utilizar HC3 --\n"
    )
    
    print(
      brown_bias[
        discrepa,
        c(
          "ajuste",
          "p_conjunto",
          "p_conjunto_hc3"
        )
      ],
      row.names = FALSE
    )
    
  } else {
    
    cat(
      "\nLa conclusión del contraste no cambia al emplear HC3.\n"
    )
  }
}


# ==============================================================================
# 5c. DISPERSIÓN FRENTE A IDENTIDAD
# ==============================================================================

pdf(
  file.path(
    ruta_out,
    "prec_bias_scatter.pdf"
  ),
  width = 6,
  height = 6
)


for (
  i in seq_len(
    nrow(brown_bias)
  )
) {
  
  col <-
    brown_bias$variable[i]
  
  x <-
    eval_data[[col]]
  
  y <-
    eval_data$DIRECTO
  
  
  lim <-
    range(
      c(x, y),
      finite = TRUE
    )
  
  
  plot(
    x,
    y,
    xlim = lim,
    ylim = lim,
    pch = 19,
    cex = 0.7,
    col = "grey30",
    xlab = "Estimador del modelo",
    ylab = "Estimador directo",
    main = brown_bias$ajuste[i]
  )
  
  
  abline(
    0,
    1,
    col = "red",
    lwd = 1.5
  )
  
  
  abline(
    lm(y ~ x),
    col = "blue",
    lwd = 1.5,
    lty = 2
  )
  
  
  legend(
    "topleft",
    bty = "n",
    cex = 0.8,
    legend = c(
      "Identidad",
      "Ajuste"
    ),
    col = c(
      "red",
      "blue"
    ),
    lty = c(
      1,
      2
    ),
    lwd = 1.5
  )
}


dev.off()


# ==============================================================================
# 6. BONDAD DEL AJUSTE — BROWN
# ==============================================================================

res_gof <- SAEval::gof(
  data = eval_data,
  dir = ~DIRECTO,
  sae = formula_sae,
  v.dir = ~PSI,
  mse.sae = formula_mse,
  alfa = 0.05
)


cat(
  "\n============================================================\n"
)

cat(
  "BONDAD DEL AJUSTE — BROWN\n"
)

cat(
  "============================================================\n\n"
)

print(
  res_gof,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Información complementaria
#
# gamma_fuente indica si el factor de encogimiento procede de emdi, de la Etapa
# II, de una aproximación escalar, o no está disponible.
#
# gamma_fila_medio difiere de gamma_medio únicamente en la familia espacial.
# Ver la nota del bloque 3.
# ------------------------------------------------------------------------------

cat(
  "\n-- Información complementaria --\n"
)

print(
  metadata %>%
    select(
      ajuste,
      mse_type,
      W_manual,
      W_sobre_m,
      gl_w,
      chi2_95_w,
      p_w
    ),
  row.names = FALSE
)


cat(
  "\n-- Factor de encogimiento --\n"
)

print(
  metadata %>%
    select(
      ajuste,
      gamma_fuente,
      gamma_medio,
      gamma_min,
      gamma_max,
      gamma_fila_medio
    ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Verificación de W
# ------------------------------------------------------------------------------

W_pkg <-
  suppressWarnings(
    as.numeric(
      res_gof$W
    )
  )


if (
  length(W_pkg) ==
  nrow(metadata)
) {
  
  dif <-
    abs(
      W_pkg -
        metadata$W_manual
    )
  
  
  if (
    any(
      is.finite(dif) &
      dif > 0.05
    )
  ) {
    
    cat(
      "\n-- AVISO: discrepancia entre W manual y SAEval --\n"
    )
    
    print(
      data.frame(
        
        ajuste =
          metadata$ajuste,
        
        W_saeval =
          round(
            W_pkg,
            2
          ),
        
        W_manual =
          metadata$W_manual,
        
        diferencia =
          round(
            dif,
            3
          )
        
      ),
      row.names = FALSE
    )
    
  } else {
    
    cat(
      "\nEl W de SAEval coincide con el cálculo directo.\n"
    )
  }
}


# ------------------------------------------------------------------------------
# Convención de grados de libertad
# ------------------------------------------------------------------------------
# Se contrasta el valor observado frente a los dos valores críticos: el de
# Brown et al. (2001), con m grados de libertad, y el que emplea SAEval, con
# m - 1. Ver la nota del bloque 4.
# ------------------------------------------------------------------------------

gl_brown <- n_dom

c2_brown <- qchisq(0.95, df = gl_brown)


comparacion_gl <- data.frame(
  
  ajuste = metadata$ajuste,
  
  W = metadata$W_manual,
  
  gl_saeval = metadata$gl_w,
  
  c2_saeval = metadata$chi2_95_w,
  
  gl_brown = gl_brown,
  
  c2_brown = round(c2_brown, 2),
  
  rechaza_saeval = metadata$W_manual > metadata$chi2_95_w,
  
  rechaza_brown = metadata$W_manual > c2_brown,
  
  stringsAsFactors = FALSE
)


cat(
  "\n-- Grados de libertad del estadístico W --\n"
)

print(
  comparacion_gl,
  row.names = FALSE
)


if (
  all(
    comparacion_gl$rechaza_saeval ==
    comparacion_gl$rechaza_brown,
    na.rm = TRUE
  )
) {
  
  cat(
    "\nLa conclusión no depende de la convención de grados de libertad.\n"
  )
  
} else {
  
  cat(
    "\nAVISO: la conclusión cambia según la convención de grados de libertad.\n"
  )
}



# ==============================================================================
# 7. COBERTURA
# ==============================================================================

res_coverage <-
  SAEval::coverage(
    data = eval_data,
    dir = ~DIRECTO,
    sae = formula_sae,
    v.dir = ~PSI,
    mse.sae = formula_mse,
    alfa = 0.05
  )


cat(
  "\n============================================================\n"
)

cat(
  "DIAGNÓSTICO DE COBERTURA\n"
)

cat(
  "============================================================\n\n"
)

print(
  res_coverage,
  row.names = FALSE
)


# ==============================================================================
# 8. INTERVALOS
# ==============================================================================

res_cinterval <-
  SAEval::cinterval(
    data = eval_data,
    dir = ~DIRECTO,
    sae = formula_sae,
    v.dir = ~PSI,
    mse.sae = formula_mse,
    level = 0.95,
    plot = FALSE
  )


cat(
  "\n============================================================\n"
)

cat(
  "ANÁLISIS DE INTERVALOS\n"
)

cat(
  "============================================================\n\n"
)

print(
  res_cinterval,
  row.names = FALSE
)


# ==============================================================================
# 8b. VERIFICACIÓN MANUAL DE cinterval() Y coverage()
# ------------------------------------------------------------------------------
# La condición de solapamiento que emplea SAEval::cinterval() es
#
#     (yUP > xLW | xLW > yLW) & (yUP > xUP | xUP > yLW)
#
# y se cumple siempre. Si xLW <= yLW, entonces yUP > yLW >= xLW, de modo que la
# primera mitad es verdadera; si xUP <= yLW, entonces yUP > yLW >= xUP, y la
# segunda también lo es. La columna overlap devuelve por tanto el número total
# de dominios, con independencia de los intervalos.
#
# La condición correcta es xLW < yUP & yLW < xUP.
#
# Este bloque:
#   a) reproduce el error con la propia librería en un ejemplo mínimo;
#   b) recalcula sobre los datos included, el solapamiento correcto al 95 % y
#      el diagnóstico de cobertura de Brown et al. (2001), y los compara con lo
#      que devolvió el paquete.
# ==============================================================================


# ------------------------------------------------------------------------------
# a) Ejemplo mínimo
# ------------------------------------------------------------------------------
# Tres dominios con errores estándar de 0,01. Los intervalos al 95 % tienen una
# semiamplitud cercana a 0,02 y el estimador del modelo se ubica a más de 0,25
# del directo, de modo que ningún par de intervalos se solapa.
# ------------------------------------------------------------------------------

ejemplo <- data.frame(
  dir  = c(0.10, 0.30, 0.50),
  vdir = c(1e-4, 1e-4, 1e-4),
  mod  = c(0.60, 0.80, 0.05),
  mse  = c(1e-4, 1e-4, 1e-4)
)

ej_cint <- SAEval::cinterval(
  data = ejemplo,
  dir = ~dir,
  sae = ~mod,
  v.dir = ~vdir,
  mse.sae = ~mse,
  level = 0.95,
  plot = FALSE
)

ej_cov <- SAEval::coverage(
  data = ejemplo,
  dir = ~dir,
  sae = ~mod,
  v.dir = ~vdir,
  mse.sae = ~mse,
  alfa = 0.05
)

cat("\n============================================================\n")
cat("VERIFICACIÓN DE cinterval() Y coverage()\n")
cat("============================================================\n\n")

cat("Ejemplo mínimo: 3 dominios, ningún par de intervalos se solapa\n\n")

cat(
  "  cinterval(): overlap =", ej_cint$overlap,
  "de", nrow(ejemplo), "\n"
)

cat(
  "  coverage():  dominios sin solapamiento =", ej_cov$non_coverage,
  "de", nrow(ejemplo), "\n"
)

cat(
  "  solapamientos reales =",
  sum(
    (ejemplo$mod - qnorm(0.975) * sqrt(ejemplo$mse)) <
      (ejemplo$dir + qnorm(0.975) * sqrt(ejemplo$vdir)) &
      (ejemplo$dir - qnorm(0.975) * sqrt(ejemplo$vdir)) <
      (ejemplo$mod + qnorm(0.975) * sqrt(ejemplo$mse))
  ),
  "\n"
)


# ------------------------------------------------------------------------------
# b) Recálculo sobre los datos
# ------------------------------------------------------------------------------

z95 <- qnorm(0.975)

verif_int <- do.call(
  rbind,
  lapply(seq_along(modelo_cols), function(i) {
    
    col  <- modelo_cols[i]
    mcol <- mse_cols[i]
    
    y  <- eval_data$DIRECTO
    sy <- sqrt(eval_data$PSI)
    x  <- eval_data[[col]]
    sx <- sqrt(eval_data[[mcol]])
    
    ok <- is.finite(y) & is.finite(sy) & sy > 0 &
      is.finite(x) & is.finite(sx)
    
    y  <- y[ok]
    sy <- sy[ok]
    x  <- x[ok]
    sx <- sx[ok]
    
    # Intervalos al 95 % con el multiplicador estándar, como en cinterval()
    yUP <- y + z95 * sy
    yLW <- y - z95 * sy
    xUP <- x + z95 * sx
    xLW <- x - z95 * sx
    
    cond_saeval   <- (yUP > xLW | xLW > yLW) & (yUP > xUP | xUP > yLW)
    cond_correcta <- (xLW < yUP) & (yLW < xUP)
    incluido      <- (y > xLW) & (y < xUP)
    
    # Diagnóstico de cobertura de Brown et al. (2001): multiplicador ajustado
    # para que, bajo independencia, la probabilidad de no solapamiento sea 5 %
    zb <- z95 * sqrt(sx^2 + sy^2) / (sx + sy)
    
    solapa_brown <- (x - zb * sx < y + zb * sy) &
      (y - zb * sy < x + zb * sx)
    
    bt <- binom.test(
      sum(solapa_brown),
      length(solapa_brown),
      p = 0.95
    )
    
    data.frame(
      variable                 = col,
      dominios                 = length(y),
      included_manual          = sum(incluido),
      solap_saeval_reproducido = sum(cond_saeval),
      solap_correcto_95        = sum(cond_correcta),
      sin_solap_brown_manual   = sum(!solapa_brown),
      p_brown_manual           = round(bt$p.value, 4),
      stringsAsFactors = FALSE
    )
  })
)


verif_int <- verif_int %>%
  left_join(
    data.frame(
      variable        = as.character(res_cinterval$methods),
      included_saeval = res_cinterval$included,
      overlap_saeval  = res_cinterval$overlap,
      stringsAsFactors = FALSE
    ),
    by = "variable"
  ) %>%
  left_join(
    data.frame(
      variable         = as.character(res_coverage$methods),
      sin_solap_saeval = res_coverage$non_coverage,
      p_saeval         = round(res_coverage$p_value, 4),
      stringsAsFactors = FALSE
    ),
    by = "variable"
  )


verif_int <- metadata %>%
  select(ajuste, variable) %>%
  left_join(verif_int, by = "variable")


cat("\nRecálculo sobre los datos\n\n")

print(
  verif_int %>%
    select(
      ajuste,
      dominios,
      included_manual,
      included_saeval,
      solap_correcto_95,
      solap_saeval_reproducido,
      overlap_saeval,
      sin_solap_brown_manual,
      sin_solap_saeval,
      p_brown_manual,
      p_saeval
    ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Comprobaciones
# ------------------------------------------------------------------------------

chequeo <- function(condicion, si, no) {
  cat(if (isTRUE(condicion)) si else no, "\n")
}

cat("\n")

chequeo(
  all(verif_int$included_manual == verif_int$included_saeval),
  "OK  included: el recálculo coincide con SAEval.",
  "AVISO  included: el recálculo NO coincide con SAEval."
)

chequeo(
  all(verif_int$solap_saeval_reproducido == verif_int$dominios),
  "OK  la condición de SAEval se cumple en todos los dominios (error confirmado).",
  "AVISO  la condición de SAEval no se cumple en todos los dominios."
)

chequeo(
  all(verif_int$overlap_saeval == verif_int$solap_saeval_reproducido),
  "OK  overlap de SAEval coincide con la condición reproducida.",
  "AVISO  overlap de SAEval no coincide con la condición reproducida."
)

chequeo(
  all(verif_int$sin_solap_brown_manual == verif_int$sin_solap_saeval),
  "OK  cobertura de Brown: el recálculo coincide con SAEval.",
  "AVISO  cobertura de Brown: el recálculo NO coincide con SAEval."
)

cat(
  "\nSolapamiento correcto al 95 %: entre",
  min(verif_int$solap_correcto_95),
  "y",
  max(verif_int$solap_correcto_95),
  "de",
  max(verif_int$dominios),
  "dominios.\n"
)


write.csv(
  verif_int,
  file.path(ruta_out, "prec_intervalos_verificacion.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)



# ==============================================================================
# 9. COEFICIENTE DE VARIACIÓN
# ==============================================================================

cv_de <- function(
    mse,
    est
) {
  
  ifelse(
    
    is.finite(est) &
      est > 0 &
      is.finite(mse),
    
    sqrt(
      pmax(
        mse,
        0
      )
    ) /
      est,
    
    NA_real_
  )
}


# ------------------------------------------------------------------------------
# CV directo
# ------------------------------------------------------------------------------

cv_data <-
  data.frame(
    
    CV_directo =
      cv_de(
        eval_data$PSI,
        eval_data$DIRECTO
      )
  )


# ------------------------------------------------------------------------------
# CV por modelo
# ------------------------------------------------------------------------------

for (
  col in modelo_cols
) {
  
  cv_data[[paste0("CV_", col)]] <-
    cv_de(
      eval_data[[paste0("MSE_", col)]],
      eval_data[[col]]
    )
}


formula_cv <-
  as.formula(
    paste(
      "~",
      paste(
        names(cv_data),
        collapse = " + "
      )
    )
  )


res_cv <-
  SAEval::cv_table(
    data = cv_data,
    cv = formula_cv,
    boxplot = FALSE
  )


cat(
  "\n============================================================\n"
)

cat(
  "DISTRIBUCIÓN DE COEFICIENTES DE VARIACIÓN\n"
)

cat(
  "============================================================\n\n"
)

print(
  res_cv,
  row.names = FALSE
)


# ==============================================================================
# 10. RESUMEN DE PRECISIÓN
# ==============================================================================

precision <-
  do.call(
    
    rbind,
    
    lapply(
      
      modelo_cols,
      
      function(col) {
        
        
        cvm <-
          cv_data[[paste0("CV_", col)]]
        
        cvd <-
          cv_data$CV_directo
        
        
        ok <-
          is.finite(cvm) &
          is.finite(cvd)
        
        
        data.frame(
          
          variable =
            col,
          
          cv_promedio =
            round(
              mean(cvm[ok]) * 100,
              3
            ),
          
          cv_mediana =
            round(
              median(cvm[ok]) * 100,
              3
            ),
          
          cv_p90 =
            round(
              quantile(
                cvm[ok],
                0.90,
                names = FALSE
              ) * 100,
              3
            ),
          
          cv_max =
            round(
              max(cvm[ok]) * 100,
              3
            ),
          
          cv_directo =
            round(
              mean(cvd[ok]) * 100,
              3
            ),
          
          reduccion_cv =
            round(
              (
                mean(cvd[ok]) -
                  mean(cvm[ok])
              ) * 100,
              3
            ),
          
          n_menor_directo =
            sum(
              cvm[ok] <
                cvd[ok]
            ),
          
          n_cv_na =
            sum(
              !ok
            ),
          
          n_dominios =
            sum(
              ok
            ),
          
          stringsAsFactors = FALSE
        )
      }
    )
  )


precision <-
  metadata %>%
  select(
    ajuste,
    variable,
    papel,
    familia,
    mse_type
  ) %>%
  left_join(
    precision,
    by = "variable"
  )


cat(
  "\n============================================================\n"
)

cat(
  "RESUMEN DE PRECISIÓN\n"
)

cat(
  "============================================================\n\n"
)

print(
  precision %>%
    select(
      -variable
    ),
  row.names = FALSE
)


# ==============================================================================
# 11. EXPORTACIÓN
# ==============================================================================

cv_por_dominio <-
  cbind(
    
    data.frame(
      
      cod_mun =
        datos$cod_mun,
      
      Municipio =
        datos$Municipio,
      
      stringsAsFactors = FALSE
    ),
    
    cv_data,
    
    as.data.frame(
      gammas
    )
  )


# ------------------------------------------------------------------------------
# Metadata
# ------------------------------------------------------------------------------

write.csv(
  metadata,
  file.path(
    ruta_out,
    "prec_modelos_metadata.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------------------------
# SAEval bias
# ------------------------------------------------------------------------------

write.csv(
  res_bias$output1,
  file.path(
    ruta_out,
    "prec_bias_saeval.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------------------------
# Brown bias
# ------------------------------------------------------------------------------

write.csv(
  brown_bias,
  file.path(
    ruta_out,
    "prec_bias_brown.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------------------------
# Brown GOF
# ------------------------------------------------------------------------------

write.csv(
  res_gof,
  file.path(
    ruta_out,
    "prec_brown_gof.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------------------------
# Coverage
# ------------------------------------------------------------------------------

write.csv(
  res_coverage,
  file.path(
    ruta_out,
    "prec_coverage.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------------------------
# Intervalos
# ------------------------------------------------------------------------------

write.csv(
  res_cinterval,
  file.path(
    ruta_out,
    "prec_intervalos.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------------------------
# CV resumen
# ------------------------------------------------------------------------------

write.csv(
  precision,
  file.path(
    ruta_out,
    "prec_cv_resumen.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------------------------
# CV por dominio
# ------------------------------------------------------------------------------

write.csv(
  cv_por_dominio,
  file.path(
    ruta_out,
    "prec_cv_por_dominio.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ==============================================================================
# 11b. FORMULACIÓN INVERSA DEL DIAGNÓSTICO DE SESGO (ANEXO)
# ------------------------------------------------------------------------------
# SAEval::bias() (bloque 5) plantea la regresión en la dirección inversa a la
# de Brown et al. (2001):
#
#     MODELO = b0 + b1 * DIRECTO
#
# y solo entrega la decisión al 5 % ("Reject"/"Accept"), sin valores p. Aquí se
# repite esa formulación con los mismos contrastes del bloque 5b: prueba
# conjunta H0: b0 = 0 y b1 = 1 con la matriz de covarianzas convencional y HC3,
# y Goldfeld-Quandt ordenado por PSI. La pendiente debe coincidir con la de
# SAEval.
#
# Salida: output/prec_bias_inversa.csv
# ==============================================================================

inversa <- do.call(rbind, lapply(modelo_cols, function(col) {

  d <- data.frame(
    y   = eval_data[[col]],        # estimación modelada como respuesta
    x   = eval_data$DIRECTO,       # estimador directo como regresor
    psi = eval_data$PSI
  )
  d <- d[is.finite(d$y) & is.finite(d$x) & is.finite(d$psi), , drop = FALSE]

  fit <- lm(y ~ x, data = d)
  smr <- summary(fit)

  p_ols <- car::linearHypothesis(fit, diag(2), c(0, 1))[["Pr(>F)"]][2]

  p_hc3 <- if (HAY_SANDWICH) {
    car::linearHypothesis(
      fit, diag(2), c(0, 1),
      vcov. = sandwich::vcovHC(fit, type = "HC3")
    )[["Pr(>F)"]][2]
  } else {
    NA_real_
  }

  gq <- tryCatch(
    lmtest::gqtest(fit, order.by = d$psi, alternative = "two.sided")$p.value,
    error = function(e) NA_real_
  )

  data.frame(
    variable       = col,
    b0             = round(unname(coef(fit)[1]), 4),
    b1             = round(unname(coef(fit)[2]), 4),
    ee_b1          = round(smr$coefficients[2, 2], 4),
    R2             = round(smr$r.squared, 4),
    p_conjunto     = round(p_ols, 4),
    p_conjunto_hc3 = round(as.numeric(p_hc3), 4),
    gq_p_psi       = round(as.numeric(gq), 4),
    stringsAsFactors = FALSE
  )
}))


# Control: la pendiente debe coincidir con la de SAEval
control_inv <- merge(
  inversa[, c("variable", "b1")],
  data.frame(variable  = res_bias$output1$methods,
             b1_saeval = round(res_bias$output1$b1, 4)),
  by = "variable", all.x = TRUE
)

if (any(abs(control_inv$b1 - control_inv$b1_saeval) > 1e-4, na.rm = TRUE)) {
  warning("La pendiente de la formulación inversa no coincide con la de SAEval.")
}


cat("\n============================================================\n")
cat("FORMULACIÓN INVERSA: MODELO = b0 + b1 * DIRECTO\n")
cat("============================================================\n\n")

print(inversa, row.names = FALSE)

cat("\nRango de pendientes:", min(inversa$b1), "a", max(inversa$b1), "\n")

write.csv(
  inversa,
  file.path(ruta_out, "prec_bias_inversa.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ==============================================================================
# 12. FINAL
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "EVALUACIÓN COMPLETADA\n"
)

cat(
  "Ajustes evaluados:",
  nrow(metadata),
  "\n"
)

cat(
  "============================================================\n"
)


