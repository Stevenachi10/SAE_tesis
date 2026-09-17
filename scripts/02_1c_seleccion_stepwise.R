# ==============================================================================
# 02_1c_seleccion_stepwise.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.1c — Selección secuencial bajo criterios específicos del
#               modelo Fay-Herriot
#
# ------------------------------------------------------------------------------
# OBJETIVO
# ------------------------------------------------------------------------------
# Seleccionar covariables directamente sobre el modelo Fay-Herriot mediante
# búsqueda secuencial, empleando los criterios KIC, KICc y KICb2. Los tres
# pertenecen a la misma familia y difieren en la corrección aplicada a la
# penalización, de menor a mayor exigencia.
#
# La búsqueda se realiza bajo máxima verosimilitud ordinaria, requisito para
# que los criterios de información resulten comparables entre modelos con
# distinto número de covariables. La especificación seleccionada se reestima
# bajo verosimilitud restringida para obtener la componente de varianza, los
# factores de contracción y las medidas de precisión.
#
# ------------------------------------------------------------------------------
# DIRECCION
# ------------------------------------------------------------------------------
# "forward" : cada paso únicamente incorpora covariables.
# "both"    : cada paso evalúa la incorporación de una covariable no presente
#             y la eliminación de una ya seleccionada, y ejecuta el movimiento
#             que produzca la mayor reducción del criterio.
#
# ------------------------------------------------------------------------------
# REPLICAS BOOTSTRAP
# ------------------------------------------------------------------------------
# emdi calcula los criterios corregidos únicamente cuando el bootstrap está
# activo. Esto incluye a KICc, pese a tratarse de una corrección de muestra
# finita: con B = c(0, 0) el elemento correspondiente no se genera.
#
# El valor de KICc depende del número de réplicas de manera apreciable, de
# modo que B_BOOT no es un parámetro de cómputo sino una decisión que puede
# incidir en el punto de parada de la búsqueda. Se declara como tal y se
# registra junto a cada especificación resultante.
#
# ------------------------------------------------------------------------------
# SEMILLA Y DETERMINISMO
# ------------------------------------------------------------------------------
# fh() fija internamente el generador de números aleatorios mediante su
# argumento seed, cuyo valor por defecto es 123. Por ese motivo la semilla se
# pasa como argumento y no mediante set.seed(), que no tendría efecto sobre el
# remuestreo interno.
#
# Verificado empíricamente: dos ajustes de la misma especificación, con el
# mismo B y la misma semilla, devuelven valores idénticos de KICc y KICb2. El
# criterio es por tanto una función determinista de la especificación, dados B
# y la semilla.
#
# ------------------------------------------------------------------------------
# TERMINACION
# ------------------------------------------------------------------------------
# Un movimiento se acepta solo si reduce el criterio en más de TOL_MEJORA.
# Como el criterio decrece estrictamente en cada movimiento aceptado y el
# número de subconjuntos de covariables es finito, ninguna especificación
# puede visitarse dos veces y la búsqueda no puede entrar en ciclo. La
# propiedad se apoya en el determinismo verificado arriba. TOL_MEJORA opera
# entonces como protección frente al error numérico del optimizador.
#
# ------------------------------------------------------------------------------
# RESTRICCIONES
# ------------------------------------------------------------------------------
# VIF_MAX opera como condición de admisibilidad durante la búsqueda: un
# candidato que llevaría el diseño por encima del umbral se descarta antes de
# ajustarse. Se aplica solo a las incorporaciones, porque retirar una columna
# no puede elevar el factor de inflación de las restantes.
#
# COR_MAX descarta covariables prácticamente duplicadas antes de iniciar la
# búsqueda.
#
# Ninguna de las dos constituye un supuesto del modelo Fay-Herriot.
#
# ------------------------------------------------------------------------------
# ALGORITMO PROPIO
# ------------------------------------------------------------------------------
# Se implementa la búsqueda en lugar de emplear step.fh() porque, con un
# universo amplio de covariables, algunas especificaciones candidatas producen
# matrices de diseño mal condicionadas. step.fh() detiene toda la búsqueda ante
# el primero de esos errores. Aquí un candidato no estimable se descarta y la
# búsqueda continúa; los descartes quedan registrados.
#
# ------------------------------------------------------------------------------
# ENTRADAS
# ------------------------------------------------------------------------------
#   output/matriz_sae_transformada_v2.rds
#
# ------------------------------------------------------------------------------
# SALIDAS (output/)
# ------------------------------------------------------------------------------
#   seleccion_stepwise<sufijo>.csv
#   seleccion_stepwise_historial<sufijo>.csv
#   seleccion_stepwise_singulares<sufijo>.csv
#   seleccion_stepwise_descartes<sufijo>.csv
#   seleccion_stepwise<sufijo>.rds
# ==============================================================================


# ==============================================================================
# PAQUETES
# ==============================================================================

library(emdi)
library(car)
library(dplyr)
library(here)

select <- dplyr::select
filter <- dplyr::filter


# ==============================================================================
# RUTAS
# ==============================================================================

ruta_out <- here("output")
matriz   <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))


# ==============================================================================
# CONFIGURACION
# ==============================================================================

# ---- Universo de candidatas -------------------------------------------------
# "reducido" : las 14 covariables definidas explicitamente
# "completo" : todas las numericas admisibles tras los filtros
UNIVERSO <- "completo"

# ---- Direccion de la busqueda -----------------------------------------------
DIRECCION <- "both"

# ---- Metodos de estimacion --------------------------------------------------
METODO_SELECCION <- "ml"      # comparacion de especificaciones
METODO_AJUSTE    <- "reml"    # ajuste final del modelo seleccionado

# ---- Criterios --------------------------------------------------------------
# Orden de ejecucion. Conviene situar primero los que no requieren bootstrap:
# cualquier error de configuracion se manifiesta en segundos y no al cabo de
# varios minutos.
CRITERIOS <- c("KIC", "KICc", "KICb2")

# ---- Replicas bootstrap -----------------------------------------------------
B_BOOT <- 500

# ---- Semilla ----------------------------------------------------------------
SEMILLA <- 2906

# ---- Tolerancia de mejora ---------------------------------------------------
TOL_MEJORA <- 1e-6

# ---- Multicolinealidad ------------------------------------------------------
VIF_MAX     <- 10
FILTRAR_VIF <- TRUE

# ---- Correlacion extrema ----------------------------------------------------
COR_MAX <- 0.99

# ---- Limite de movimientos --------------------------------------------------
# NULL lo deriva del tamaño del universo. No limita el numero de covariables:
# la busqueda se detiene por si sola cuando ningun movimiento mejora.
MOVS_MAX <- NULL


# ==============================================================================
# CRITERIOS DISPONIBLES Y REPLICAS REQUERIDAS
# ==============================================================================

CRITERIOS_DISPONIBLES <- c("AIC", "AICc", "AICb1", "AICb2",
                           "BIC",
                           "KIC", "KICc", "KICb1", "KICb2")

# Criterios que emdi solo calcula con bootstrap activo.
CRITERIOS_CON_BOOT <- c("AICc", "AICb1", "AICb2",
                        "KICc", "KICb1", "KICb2")

B_de <- function(criterio) {
  if (criterio %in% CRITERIOS_CON_BOOT) c(0, B_BOOT) else c(0, 0)
}


# ==============================================================================
# CANDIDATAS DEL UNIVERSO REDUCIDO
# ==============================================================================

CANDIDATAS_REDUCIDO <- c(
  "extension",
  "saber11_lectura_2022",
  "partos_calificados_2020",
  "dist_popayan_km",
  "dist_cali_km",
  "pct_contributivo",
  "rural_disperso",
  "tiene_coca",
  "desplaz_tasa_1000",
  "pct_rural",
  "categoria_617",
  "alcantarillado_2018",
  "mort_infantil_1a_2020",
  "pct_terciario_2023"
)


# ==============================================================================
# VARIABLES QUE NO PUEDEN SER COVARIABLES
# ------------------------------------------------------------------------------
# Identificadores, la variable de interes y toda magnitud derivada de ella o de
# su varianza de muestreo. Su incorporacion constituiria una fuga de
# informacion: el criterio premiaria el uso de la propia respuesta.
# ==============================================================================

NO_COVARIABLES <- c(
  "cod_mun", "Municipio", "cod_dpto", "Departamento",
  "pobreza_monetaria", "varianza_pobreza",
  "ee_pobreza", "cv_directo", "n_eff", "n_muestra",
  "error_estandar", "cvlog"
)


# ==============================================================================
# CONSTRUCCION DEL UNIVERSO DE CANDIDATAS
# ==============================================================================

if (UNIVERSO == "reducido") {
  
  CANDIDATAS <- CANDIDATAS_REDUCIDO
  
  faltan <- setdiff(CANDIDATAS, names(matriz))
  if (length(faltan) > 0) {
    stop("Covariables ausentes en la matriz: ",
         paste(faltan, collapse = ", "))
  }
  
} else if (UNIVERSO == "completo") {
  
  # ---- Solo columnas numericas ----------------------------------------------
  # Cualquier variable almacenada como factor o caracter queda fuera del
  # universo sin figurar en los descartes. Se reporta para dejar constancia.
  numericas    <- names(matriz)[sapply(matriz, is.numeric)]
  no_numericas <- setdiff(names(matriz), numericas)
  
  if (length(no_numericas) > 0) {
    cat("Columnas no numericas, fuera del universo: ",
        paste(no_numericas, collapse = ", "), "\n", sep = "")
  }
  
  CANDIDATAS <- setdiff(numericas, NO_COVARIABLES)
  
  # ---- Constantes -----------------------------------------------------------
  constantes <- CANDIDATAS[sapply(CANDIDATAS, function(v) {
    x <- matriz[[v]]
    length(unique(x[!is.na(x)])) < 2
  })]
  
  if (length(constantes) > 0) {
    cat("Descartadas por ser constantes: ",
        paste(constantes, collapse = ", "), "\n", sep = "")
    CANDIDATAS <- setdiff(CANDIDATAS, constantes)
  }
  
  # ---- Valores faltantes ----------------------------------------------------
  # Se exige informacion completa para mantener constante el conjunto de
  # dominios entre especificaciones.
  con_na <- CANDIDATAS[sapply(CANDIDATAS, function(v) {
    any(is.na(matriz[[v]]))
  })]
  
  if (length(con_na) > 0) {
    cat("Descartadas por valores faltantes: ",
        paste(con_na, collapse = ", "), "\n", sep = "")
    CANDIDATAS <- setdiff(CANDIDATAS, con_na)
  }
  
  # ---- Correlacion extrema --------------------------------------------------
  # Ante un par por encima del umbral se conserva la de menor indice, es decir
  # la que aparece antes en la matriz. Los pares se imprimen para que la
  # decision pueda revisarse.
  if (length(CANDIDATAS) > 1) {
    
    M  <- as.matrix(matriz[, CANDIDATAS, drop = FALSE])
    cm <- abs(cor(M, use = "pairwise.complete.obs"))
    cm[!upper.tri(cm)] <- 0
    
    redundantes <- character(0)
    
    for (j in seq_along(CANDIDATAS)) {
      
      if (CANDIDATAS[j] %in% redundantes) next
      
      pares <- CANDIDATAS[which(cm[j, ] > COR_MAX)]
      pares <- setdiff(pares, redundantes)
      
      if (length(pares) > 0) {
        cat("Redundante con ", CANDIDATAS[j], ": ",
            paste(pares, collapse = ", "), "\n", sep = "")
        redundantes <- c(redundantes, pares)
      }
    }
    
    CANDIDATAS <- setdiff(CANDIDATAS, redundantes)
  }
  
  # ---- Dependencia lineal exacta --------------------------------------------
  # Un conjunto completo de indicadoras de una misma categorica suma la unidad
  # y resulta colineal con el intercepto. La descomposicion QR identifica las
  # columnas que no aportan rango.
  if (length(CANDIDATAS) > 0) {
    
    X <- cbind(1, as.matrix(matriz[, CANDIDATAS, drop = FALSE]))
    colnames(X)[1] <- "(Intercept)"
    
    qrX <- qr(X)
    
    if (qrX$rank < ncol(X)) {
      
      conservar <- colnames(X)[qrX$pivot[seq_len(qrX$rank)]]
      aliadas   <- setdiff(CANDIDATAS, conservar)
      
      if (length(aliadas) > 0) {
        cat("Descartadas por dependencia lineal exacta: ",
            paste(aliadas, collapse = ", "), "\n", sep = "")
      }
      
      CANDIDATAS <- intersect(CANDIDATAS, conservar)
    }
  }
  
} else {
  stop("UNIVERSO debe ser 'reducido' o 'completo'")
}


# ==============================================================================
# VALIDACIONES
# ==============================================================================

if (!DIRECCION %in% c("forward", "both")) {
  stop("DIRECCION debe ser 'forward' o 'both'")
}

no_validos <- setdiff(CRITERIOS, CRITERIOS_DISPONIBLES)

if (length(no_validos) > 0) {
  stop("Criterios no reconocidos: ", paste(no_validos, collapse = ", "),
       ". Disponibles: ", paste(CRITERIOS_DISPONIBLES, collapse = ", "))
}

if (length(CANDIDATAS) == 0) {
  stop("El universo de candidatas quedo vacio tras los filtros.")
}

if (is.null(MOVS_MAX)) {
  MOVS_MAX <- if (DIRECCION == "both") {
    3 * length(CANDIDATAS)
  } else {
    length(CANDIDATAS)
  }
}


# ==============================================================================
# INFORMACION INICIAL
# ==============================================================================

cat("\n============================================================\n")
cat("SELECCION SECUENCIAL — FAY-HERRIOT\n")
cat("============================================================\n")
cat("Universo:            ", UNIVERSO, "\n", sep = "")
cat("Direccion:           ", DIRECCION, "\n", sep = "")
cat("Candidatas:          ", length(CANDIDATAS), "\n", sep = "")
cat("Dominios:            ", nrow(matriz), "\n", sep = "")
cat("Metodo seleccion:    ", METODO_SELECCION, "\n", sep = "")
cat("Metodo ajuste final: ", METODO_AJUSTE, "\n", sep = "")
cat("Criterios:           ", paste(CRITERIOS, collapse = ", "), "\n", sep = "")
cat("Con bootstrap:       ",
    paste(intersect(CRITERIOS, CRITERIOS_CON_BOOT), collapse = ", "),
    "\n", sep = "")
cat("Replicas B:          ", B_BOOT, "\n", sep = "")
cat("Semilla:             ", SEMILLA, "\n", sep = "")
cat("Tolerancia mejora:   ", TOL_MEJORA, "\n", sep = "")
cat("Maximo movimientos:  ", MOVS_MAX, "\n", sep = "")
cat("Filtro VIF < ", VIF_MAX, ":     ",
    if (FILTRAR_VIF) "activo" else "inactivo", "\n", sep = "")
cat("============================================================\n\n")


# ==============================================================================
# DATOS
# ==============================================================================

datos <- matriz %>%
  select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza,
         all_of(CANDIDATAS)) %>%
  as.data.frame()


# ==============================================================================
# FUNCIONES AUXILIARES
# ==============================================================================

formula_de <- function(vars) {
  if (length(vars) == 0) {
    return(as.formula("pobreza_monetaria ~ 1"))
  }
  as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
}

escalar <- function(x) {
  v <- tryCatch(as.numeric(x), error = function(e) NA_real_)
  if (length(v) != 1 || !is.finite(v)) return(NA_real_)
  v
}

# Identificador de una especificacion, invariante al orden de las covariables.
clave_de <- function(vars) {
  if (length(vars) == 0) return("<intercepto>")
  paste(sort(vars), collapse = "|")
}

# ------------------------------------------------------------------------------
# VIF del diseño
# ------------------------------------------------------------------------------
# Todas las covariables son numericas o indicadoras, de modo que car::vif()
# devuelve un vector con nombres. La rama matricial contempla el caso GVIF por
# seguridad y no se activa en esta aplicacion.
vif_de <- function(vars) {
  
  if (length(vars) < 2) return(1)
  
  tryCatch({
    
    v <- car::vif(lm(formula_de(vars), data = datos))
    
    if (is.matrix(v)) {
      
      col <- if ("GVIF^(1/(2*Df))" %in% colnames(v)) {
        "GVIF^(1/(2*Df))"
      } else {
        colnames(v)[1]
      }
      
      escalar(max(v[, col], na.rm = TRUE))
      
    } else {
      
      escalar(max(v, na.rm = TRUE))
    }
    
  }, error = function(e) NA_real_)
}


# ==============================================================================
# AJUSTE DURANTE LA BUSQUEDA
# ------------------------------------------------------------------------------
# El bootstrap se activa unicamente cuando el criterio lo exige. Un candidato
# no estimable devuelve NULL y la busqueda continua.
# ==============================================================================

ajustar_candidato <- function(vars, criterio) {
  
  argumentos <- list(
    formula_de(vars),
    vardir        = "varianza_pobreza",
    combined_data = datos,
    domains       = "cod_mun",
    method        = METODO_SELECCION,
    MSE           = FALSE,
    B             = B_de(criterio),
    seed          = SEMILLA
  )
  
  tryCatch(
    withCallingHandlers(
      do.call(fh, argumentos),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
}


# ==============================================================================
# EXTRACCION DEL CRITERIO
# ==============================================================================

extraer_criterio <- function(ajuste, criterio) {
  
  if (is.null(ajuste)) return(NA_real_)
  
  tryCatch({
    x <- ajuste$model$model_select[[criterio]]
    x <- as.numeric(x)[1]
    if (!is.finite(x)) return(NA_real_)
    x
  }, error = function(e) NA_real_)
}


# ==============================================================================
# AJUSTE FINAL BAJO MAXIMA VEROSIMILITUD
# ------------------------------------------------------------------------------
# Sirve unicamente para recuperar los criterios de informacion de la
# especificacion seleccionada. El bootstrap se activa para que todos queden
# disponibles en el mismo objeto. No se solicita MSE.
# ==============================================================================

ajustar_final_ml <- function(vars) {
  
  limite <- FALSE
  
  ajuste <- withCallingHandlers(
    tryCatch(
      fh(formula_de(vars),
         vardir        = "varianza_pobreza",
         combined_data = datos,
         domains       = "cod_mun",
         method        = METODO_SELECCION,
         MSE           = FALSE,
         B             = c(0, B_BOOT),
         seed          = SEMILLA),
      error = function(e) NULL
    ),
    warning = function(w) {
      if (grepl("interval limit|singular", conditionMessage(w),
                ignore.case = TRUE)) {
        limite <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  
  list(ajuste = ajuste, limite = limite)
}


# ==============================================================================
# AJUSTE FINAL BAJO VEROSIMILITUD RESTRINGIDA
# ------------------------------------------------------------------------------
# De aqui salen la componente de varianza, los factores de contraccion, el
# EBLUP y las medidas de precision. Como el MSE es analitico, no se activa
# bootstrap.
# ==============================================================================

ajustar_final_reml <- function(vars) {
  
  limite <- FALSE
  
  ajuste <- withCallingHandlers(
    tryCatch(
      fh(formula_de(vars),
         vardir        = "varianza_pobreza",
         combined_data = datos,
         domains       = "cod_mun",
         method        = METODO_AJUSTE,
         MSE           = TRUE,
         mse_type      = "analytical",
         B             = c(0, 0),
         seed          = SEMILLA),
      error = function(e) NULL
    ),
    warning = function(w) {
      if (grepl("interval limit|singular", conditionMessage(w),
                ignore.case = TRUE)) {
        limite <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  
  list(ajuste = ajuste, limite = limite)
}


# ==============================================================================
# RESUMEN DEL MODELO SELECCIONADO
# ------------------------------------------------------------------------------
# Los criterios de informacion provienen del ajuste bajo maxima verosimilitud;
# los indicadores de precision, del ajuste bajo verosimilitud restringida. Si
# la reestimacion bajo REML falla, no se sustituye por ML: la fila no se
# produce y se emite una advertencia.
# ==============================================================================

resumir_modelo <- function(vars, modelo, criterio, etapa) {
  
  res_ml <- ajustar_final_ml(vars)
  r_ml   <- res_ml$ajuste
  
  if (is.null(r_ml)) {
    warning("El modelo seleccionado por ", criterio,
            " no pudo ajustarse bajo ", METODO_SELECCION, ".")
    return(NULL)
  }
  
  res_re <- ajustar_final_reml(vars)
  r_re   <- res_re$ajuste
  
  if (is.null(r_re)) {
    warning("El modelo seleccionado por ", criterio,
            " no pudo reestimarse bajo ", METODO_AJUSTE, ".")
    return(NULL)
  }
  
  toma <- function(nm) {
    tryCatch(as.numeric(r_ml$model$model_select[[nm]])[1],
             error = function(e) NA_real_)
  }
  
  s2u <- escalar(r_re$model$variance)
  gam <- s2u / (s2u + datos$varianza_pobreza)
  cvs <- (sqrt(r_re$MSE$FH) / r_re$ind$FH) * 100
  
  sh <- tryCatch(
    summary(r_re)$normality["Random_effects", "Shapiro_p"],
    error = function(e) NA_real_
  )
  
  data.frame(
    modelo       = modelo,
    criterio     = criterio,
    direccion    = DIRECCION,
    etapa        = etapa,
    n_vars       = length(vars),
    kic          = round(toma("KIC"),   4),
    kicc         = round(toma("KICc"),  4),
    kicb2        = round(toma("KICb2"), 4),
    aic          = round(toma("AIC"),   4),
    bic          = round(toma("BIC"),   4),
    fh_r2        = round(toma("FH_R2"), 4),
    adj_r2       = round(toma("AdjR2"), 4),
    metodo_sel   = METODO_SELECCION,
    metodo_ind   = METODO_AJUSTE,
    b_boot       = B_BOOT,
    semilla      = SEMILLA,
    vif_max      = round(vif_de(vars), 2),
    sigma2_u     = signif(s2u, 6),
    gamma_medio  = round(mean(gam, na.rm = TRUE), 4),
    gamma_min    = round(min(gam,  na.rm = TRUE), 4),
    cv_medio     = round(mean(cvs, na.rm = TRUE), 3),
    cv_max       = round(max(cvs,  na.rm = TRUE), 3),
    n_cv_20      = sum(cvs < 20, na.rm = TRUE),
    shapiro_ef_p = signif(sh, 4),
    var_limite   = res_ml$limite | res_re$limite,
    covariables  = paste(vars, collapse = " + "),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# ALGORITMO SECUENCIAL
# ==============================================================================

stepwise_fh <- function(criterio, candidatas, movs_max,
                        direccion = DIRECCION) {
  
  usar_boot <- criterio %in% CRITERIOS_CON_BOOT
  
  seleccionadas    <- character(0)
  historial        <- list()
  modelos_singular <- list()
  descartes        <- list()
  
  # ---- Memoria de especificaciones ya evaluadas -----------------------------
  # Con movimientos en ambas direcciones una misma especificacion puede
  # proponerse varias veces. El valor del criterio se conserva para no repetir
  # el ajuste, lo que importa sobre todo con los criterios que usan bootstrap.
  memoria   <- new.env(parent = emptyenv())
  n_ajustes <- 0L
  
  evaluar <- function(vars) {
    
    k <- clave_de(vars)
    
    if (exists(k, envir = memoria, inherits = FALSE)) {
      return(get(k, envir = memoria, inherits = FALSE))
    }
    
    ajuste <- ajustar_candidato(vars, criterio)
    n_ajustes <<- n_ajustes + 1L
    
    res <- if (is.null(ajuste)) {
      list(valor = NA_real_, estado = "singular")
    } else {
      v <- extraer_criterio(ajuste, criterio)
      if (is.na(v)) {
        list(valor = NA_real_, estado = "sin_criterio")
      } else {
        list(valor = v, estado = "ok")
      }
    }
    
    assign(k, res, envir = memoria)
    res
  }
  
  cat("\n============================================================\n")
  cat("INICIO BUSQUEDA — ", criterio, " (", direccion, ")",
      if (usar_boot) paste0(" bootstrap B = ", B_BOOT) else "",
      "\n", sep = "")
  cat("============================================================\n")
  
  res0 <- evaluar(character(0))
  
  if (res0$estado == "singular") {
    stop("El modelo nulo no pudo estimarse bajo ", METODO_SELECCION,
         ". Verificar los argumentos de fh().")
  }
  
  if (is.na(res0$valor)) {
    stop("El criterio ", criterio, " no esta disponible en el objeto fh. ",
         "Verificar si requiere replicas bootstrap (argumento B).")
  }
  
  criterio_actual <- res0$valor
  
  cat("Paso 0: ", criterio, " = ", round(criterio_actual, 4), "\n", sep = "")
  
  for (paso in seq_len(movs_max)) {
    
    candidatos_paso <- list()
    
    cat("\n------------------------------------------------------------\n")
    cat("Paso ", paso, " | Modelo actual: ",
        if (length(seleccionadas) == 0) "intercepto"
        else paste(seleccionadas, collapse = " + "),
        "\n", sep = "")
    
    tiempo_paso <- Sys.time()
    n_vif <- 0L
    
    # ------------------------------------------------------------------------
    # Movimientos de incorporacion
    # ------------------------------------------------------------------------
    
    restantes <- setdiff(candidatas, seleccionadas)
    
    hay_elim <- direccion == "both" && length(seleccionadas) >= 2
    
    cat("Evaluando ", length(restantes), " incorporaciones",
        if (hay_elim) {
          paste0(" y ", length(seleccionadas), " eliminaciones")
        } else "",
        "...\n", sep = "")
    
    for (variable in restantes) {
      
      vars_candidato <- c(seleccionadas, variable)
      
      # -- Admisibilidad: multicolinealidad ----------------------------------
      # Se evalua antes del ajuste. El calculo del VIF sobre el diseño es
      # varios ordenes de magnitud mas rapido que la estimacion del modelo.
      if (FILTRAR_VIF) {
        
        vm <- vif_de(vars_candidato)
        
        if (is.na(vm) || vm >= VIF_MAX) {
          
          descartes[[length(descartes) + 1]] <- data.frame(
            criterio = criterio,
            paso     = paso,
            variable = variable,
            motivo   = if (is.na(vm)) "VIF_no_calculable" else "VIF",
            valor    = round(vm, 2),
            stringsAsFactors = FALSE
          )
          
          n_vif <- n_vif + 1L
          next
        }
      }
      
      res <- evaluar(vars_candidato)
      
      if (res$estado == "singular") {
        
        modelos_singular[[length(modelos_singular) + 1]] <- data.frame(
          criterio    = criterio,
          paso        = paso,
          movimiento  = "+",
          variable    = variable,
          covariables = paste(vars_candidato, collapse = " + "),
          stringsAsFactors = FALSE
        )
        
        cat("  + ", variable, " -> SINGULAR / NO ESTIMABLE\n", sep = "")
        next
      }
      
      if (res$estado == "sin_criterio") {
        cat("  + ", variable, " -> criterio no disponible\n", sep = "")
        next
      }
      
      candidatos_paso[[paste0("+", variable)]] <- list(
        vars       = vars_candidato,
        criterio   = res$valor,
        movimiento = "+",
        variable   = variable
      )
      
      cat(sprintf("  + %-32s %.4f\n", variable, res$valor))
    }
    
    # ------------------------------------------------------------------------
    # Movimientos de eliminacion
    # ------------------------------------------------------------------------
    # Se evaluan a partir de dos covariables seleccionadas. Con una sola, la
    # unica eliminacion posible devuelve al modelo nulo, cuyo criterio ya
    # resulto superior al de la especificacion vigente.
    #
    # No se aplica el filtro de VIF: retirar una columna del diseño no puede
    # elevar el factor de inflacion de las restantes.
    
    if (hay_elim) {
      
      for (variable in seleccionadas) {
        
        vars_candidato <- setdiff(seleccionadas, variable)
        
        res <- evaluar(vars_candidato)
        
        if (res$estado == "singular") {
          
          modelos_singular[[length(modelos_singular) + 1]] <- data.frame(
            criterio    = criterio,
            paso        = paso,
            movimiento  = "-",
            variable    = variable,
            covariables = paste(vars_candidato, collapse = " + "),
            stringsAsFactors = FALSE
          )
          
          cat("  - ", variable, " -> SINGULAR / NO ESTIMABLE\n", sep = "")
          next
        }
        
        if (res$estado == "sin_criterio") {
          cat("  - ", variable, " -> criterio no disponible\n", sep = "")
          next
        }
        
        candidatos_paso[[paste0("-", variable)]] <- list(
          vars       = vars_candidato,
          criterio   = res$valor,
          movimiento = "-",
          variable   = variable
        )
        
        cat(sprintf("  - %-32s %.4f\n", variable, res$valor))
      }
    }
    
    if (FILTRAR_VIF) {
      cat("\nDescartados por VIF en este paso: ", n_vif,
          " | movimientos evaluados: ", length(candidatos_paso), "\n",
          sep = "")
    }
    
    if (length(candidatos_paso) == 0) {
      cat("\nNingun movimiento satisface las condiciones de admisibilidad.\n")
      cat("La busqueda termina con ", length(seleccionadas),
          " covariables.\n", sep = "")
      break
    }
    
    valores      <- sapply(candidatos_paso, function(x) x$criterio)
    mejor_nombre <- names(valores)[which.min(valores)]
    mejor        <- candidatos_paso[[mejor_nombre]]
    
    nuevo_criterio <- mejor$criterio
    mejora         <- (criterio_actual - nuevo_criterio) > TOL_MEJORA
    
    cat("\nMejor movimiento: ", mejor$movimiento, " ", mejor$variable,
        "\n", sep = "")
    cat(criterio, " actual:    ", round(criterio_actual, 4), "\n", sep = "")
    cat(criterio, " candidato: ", round(nuevo_criterio, 4), "\n", sep = "")
    
    if (!mejora) {
      cat("\nNingun movimiento reduce el criterio por encima de la ",
          "tolerancia.\n", sep = "")
      cat("La busqueda termina con ", length(seleccionadas),
          " covariables.\n", sep = "")
      break
    }
    
    criterio_previo <- criterio_actual
    
    seleccionadas   <- mejor$vars
    criterio_actual <- nuevo_criterio
    
    tiempo_paso_min <- as.numeric(
      difftime(Sys.time(), tiempo_paso, units = "mins")
    )
    
    historial[[length(historial) + 1]] <- data.frame(
      criterio       = criterio,
      paso           = paso,
      movimiento     = mejor$movimiento,
      variable       = mejor$variable,
      n_vars         = length(seleccionadas),
      valor_criterio = round(nuevo_criterio, 4),
      mejora         = round(criterio_previo - nuevo_criterio, 4),
      vif_max        = round(vif_de(seleccionadas), 2),
      n_descartes    = n_vif,
      tiempo_min     = round(tiempo_paso_min, 3),
      covariables    = paste(seleccionadas, collapse = " + "),
      stringsAsFactors = FALSE
    )
    
    cat("\n>>> MOVIMIENTO ACEPTADO: ", mejor$movimiento, " ",
        mejor$variable, "\n", sep = "")
    cat(">>> ", criterio, " = ", round(nuevo_criterio, 4),
        " (mejora ", round(criterio_previo - nuevo_criterio, 4), ")\n",
        sep = "")
    cat(">>> Covariables: ", length(seleccionadas), "\n", sep = "")
    cat(">>> Tiempo del paso: ", round(tiempo_paso_min, 2),
        " minutos\n", sep = "")
  }
  
  cat("\nAjustes efectivos (sin repeticiones): ", n_ajustes, "\n", sep = "")
  
  historial_df <- if (length(historial) > 0) {
    do.call(rbind, historial)
  } else data.frame()
  
  singular_df <- if (length(modelos_singular) > 0) {
    do.call(rbind, modelos_singular)
  } else data.frame()
  
  descartes_df <- if (length(descartes) > 0) {
    do.call(rbind, descartes)
  } else data.frame()
  
  list(
    variables      = seleccionadas,
    criterio       = criterio,
    valor_criterio = criterio_actual,
    n_ajustes      = n_ajustes,
    historial      = historial_df,
    singular       = singular_df,
    descartes      = descartes_df
  )
}


# ==============================================================================
# EJECUCION
# ==============================================================================

resultados  <- list()
historiales <- list()
singulares  <- list()
descartados <- list()

ajustes_por_criterio <- integer(0)

tiempo_total <- Sys.time()

for (criterio in CRITERIOS) {
  
  cat("\n\n############################################################\n")
  cat("# CRITERIO: ", criterio, "\n", sep = "")
  cat("############################################################\n")
  
  tiempo_criterio <- Sys.time()
  
  resultado <- tryCatch(
    stepwise_fh(criterio   = criterio,
                candidatas = CANDIDATAS,
                movs_max   = MOVS_MAX,
                direccion  = DIRECCION),
    error = function(e) {
      cat("\nERROR EN ", criterio, ": ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  
  if (is.null(resultado)) next
  
  vars_finales <- resultado$variables
  
  if (length(vars_finales) == 0) {
    cat("\nEl criterio no selecciono covariables.\n")
    next
  }
  
  etiqueta <- paste0("F_", criterio, "_", DIRECCION)
  
  resumen <- resumir_modelo(vars     = vars_finales,
                            modelo   = etiqueta,
                            criterio = criterio,
                            etapa    = "modelo_final")
  
  if (!is.null(resumen)) resultados[[criterio]] <- resumen
  
  historiales[[criterio]] <- resultado$historial
  singulares[[criterio]]  <- resultado$singular
  descartados[[criterio]] <- resultado$descartes
  
  ajustes_por_criterio[criterio] <- resultado$n_ajustes
  
  cat("\nTiempo total ", criterio, ": ",
      round(as.numeric(difftime(Sys.time(), tiempo_criterio,
                                units = "mins")), 2),
      " minutos\n", sep = "")
}

if (length(resultados) == 0) {
  stop("Ninguna variante produjo un modelo final valido bajo ",
       METODO_AJUSTE, ".")
}

tabla_sel <- do.call(rbind, resultados)


# ==============================================================================
# RESULTADOS
# ==============================================================================

cat("\n\n============================================================\n")
cat("MODELOS FINALES SELECCIONADOS\n")
cat("============================================================\n")
cat("Criterios bajo ", METODO_SELECCION,
    " | indicadores bajo ", METODO_AJUSTE,
    " | B = ", B_BOOT, " | semilla = ", SEMILLA, "\n\n", sep = "")

print(
  tabla_sel %>%
    select(modelo, criterio, n_vars, kic, kicc, kicb2, bic,
           fh_r2, vif_max, sigma2_u, gamma_medio, gamma_min,
           cv_medio, cv_max, n_cv_20, shapiro_ef_p, var_limite),
  row.names = FALSE
)

cat("\n\n============================================================\n")
cat("COVARIABLES SELECCIONADAS\n")
cat("============================================================\n")

for (i in seq_len(nrow(tabla_sel))) {
  cat("\n", tabla_sel$modelo[i], " (", tabla_sel$n_vars[i], "):\n", sep = "")
  cat("  ", gsub(" \\+ ", ", ", tabla_sel$covariables[i]), "\n", sep = "")
}


# ==============================================================================
# DIAGNOSTICO DE MULTICOLINEALIDAD
# ==============================================================================

cat("\n\n============================================================\n")
cat("DIAGNOSTICO DE MULTICOLINEALIDAD\n")
cat("============================================================\n")

if (FILTRAR_VIF) {
  cat("El VIF opera como condicion de admisibilidad durante la busqueda: ",
      "ninguna\nespecificacion con VIF mayor o igual a ", VIF_MAX,
      " llega a ajustarse.\n\n", sep = "")
} else {
  cat("El VIF no participa en la seleccion; se reporta como diagnostico.\n\n")
}

excede <- tabla_sel %>% filter(!is.na(vif_max), vif_max >= VIF_MAX)

if (nrow(excede) > 0) {
  print(excede %>% select(modelo, n_vars, vif_max), row.names = FALSE)
} else {
  cat("Ningun modelo seleccionado presenta VIF mayor o igual a ",
      VIF_MAX, ".\n", sep = "")
}


# ==============================================================================
# HISTORIAL
# ==============================================================================

historial_total <- if (length(historiales) > 0) {
  do.call(rbind, historiales)
} else data.frame()

cat("\n\n============================================================\n")
cat("HISTORIAL DE SELECCION\n")
cat("============================================================\n\n")

if (nrow(historial_total) > 0) {
  
  print(historial_total %>%
          select(criterio, paso, movimiento, variable, n_vars,
                 valor_criterio, mejora, vif_max, n_descartes, tiempo_min),
        row.names = FALSE)
  
  retiradas <- historial_total %>% filter(movimiento == "-")
  
  if (nrow(retiradas) > 0) {
    cat("\n-- Covariables retiradas por la direccion hacia atras --\n")
    print(as.data.frame(
      retiradas %>% select(criterio, paso, variable, mejora)),
      row.names = FALSE)
  } else {
    cat("\nNinguna covariable fue retirada: la busqueda hacia atras no ",
        "encontro\nespecificaciones mejores que las de la ruta hacia ",
        "adelante.\n", sep = "")
  }
}


# ==============================================================================
# CANDIDATOS NO ESTIMABLES
# ==============================================================================

singulares_total <- if (length(singulares) > 0) {
  x <- singulares[sapply(singulares, nrow) > 0]
  if (length(x) > 0) do.call(rbind, x) else data.frame()
} else data.frame()

cat("\n\n============================================================\n")
cat("CANDIDATOS NO ESTIMABLES\n")
cat("============================================================\n\n")

if (nrow(singulares_total) > 0) {
  
  cat("Especificaciones descartadas: ", nrow(singulares_total), "\n\n",
      sep = "")
  
  cat("-- Frecuencia por variable --\n")
  print(as.data.frame(
    singulares_total %>%
      group_by(variable) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n))),
    row.names = FALSE)
  
} else {
  cat("No se detectaron candidatos no estimables.\n")
}


# ==============================================================================
# DESCARTES POR ADMISIBILIDAD
# ==============================================================================

descartes_total <- if (length(descartados) > 0) {
  x <- descartados[sapply(descartados, nrow) > 0]
  if (length(x) > 0) do.call(rbind, x) else data.frame()
} else data.frame()

if (nrow(descartes_total) > 0) {
  
  cat("\n\n============================================================\n")
  cat("DESCARTES POR CONDICIONES DE ADMISIBILIDAD\n")
  cat("============================================================\n\n")
  
  print(as.data.frame(
    descartes_total %>%
      group_by(criterio, motivo) %>%
      summarise(descartes = n(), .groups = "drop")),
    row.names = FALSE)
  
  cat("\n-- Variables mas descartadas --\n")
  print(as.data.frame(
    descartes_total %>%
      group_by(variable, motivo) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n)) %>% head(15)),
    row.names = FALSE)
}


# ==============================================================================
# COMPARACION ENTRE CRITERIOS
# ==============================================================================

if (nrow(tabla_sel) > 1) {
  
  cat("\n\n============================================================\n")
  cat("COMPARACION ENTRE CRITERIOS\n")
  cat("============================================================\n\n")
  
  conjuntos <- strsplit(tabla_sel$covariables, " + ", fixed = TRUE)
  names(conjuntos) <- tabla_sel$criterio
  
  for (i in seq_along(conjuntos)) {
    cat(names(conjuntos)[i], " (", length(conjuntos[[i]]), "): ",
        paste(conjuntos[[i]], collapse = ", "), "\n", sep = "")
  }
  
  comunes <- Reduce(intersect, conjuntos)
  
  cat("\nSeleccionadas por todos los criterios (", length(comunes), "): ",
      if (length(comunes) == 0) "ninguna"
      else paste(comunes, collapse = ", "), "\n", sep = "")
  
  for (i in seq_along(conjuntos)) {
    propias <- setdiff(conjuntos[[i]], comunes)
    cat("Solo en ", names(conjuntos)[i], ": ",
        if (length(propias) == 0) "ninguna"
        else paste(propias, collapse = ", "), "\n", sep = "")
  }
  
  # ---- Anidamiento ----------------------------------------------------------
  # Si las especificaciones forman una cadena, la diferencia entre criterios se
  # reduce al punto de parada y no a la composicion del modelo.
  ordenados <- conjuntos[order(sapply(conjuntos, length))]
  
  cadena <- if (length(ordenados) < 2) TRUE else {
    all(sapply(seq_len(length(ordenados) - 1), function(i) {
      all(ordenados[[i]] %in% ordenados[[i + 1]])
    }))
  }
  
  if (cadena) {
    cat("\nLas especificaciones forman una cadena anidada: la diferencia\n")
    cat("entre criterios se reduce al punto de parada.\n")
  } else {
    cat("\nLas especificaciones NO son anidadas: difieren en composicion y\n")
    cat("no solo en el numero de covariables.\n")
  }
}


# ==============================================================================
# GUARDADO
# ==============================================================================

sufijo <- paste0("_", UNIVERSO, "_", DIRECCION)

write.csv(tabla_sel,
          file.path(ruta_out,
                    paste0("seleccion_stepwise", sufijo, ".csv")),
          row.names = FALSE, fileEncoding = "UTF-8")

if (nrow(historial_total) > 0) {
  write.csv(historial_total,
            file.path(ruta_out,
                      paste0("seleccion_stepwise_historial", sufijo, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}

if (nrow(singulares_total) > 0) {
  write.csv(singulares_total,
            file.path(ruta_out,
                      paste0("seleccion_stepwise_singulares", sufijo, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}

if (nrow(descartes_total) > 0) {
  write.csv(descartes_total,
            file.path(ruta_out,
                      paste0("seleccion_stepwise_descartes", sufijo, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}

saveRDS(
  list(
    tabla                = tabla_sel,
    historial            = historial_total,
    singulares           = singulares_total,
    descartes            = descartes_total,
    direccion            = DIRECCION,
    filtrar_vif          = FILTRAR_VIF,
    universo             = UNIVERSO,
    candidatas           = CANDIDATAS,
    criterios            = CRITERIOS,
    criterios_con_boot   = intersect(CRITERIOS, CRITERIOS_CON_BOOT),
    metodo_seleccion     = METODO_SELECCION,
    metodo_ajuste        = METODO_AJUSTE,
    B_boot               = B_BOOT,
    semilla              = SEMILLA,
    tol_mejora           = TOL_MEJORA,
    movs_max             = MOVS_MAX,
    cor_max              = COR_MAX,
    vif_max              = VIF_MAX,
    ajustes_por_criterio = ajustes_por_criterio,
    session              = sessionInfo()
  ),
  file.path(ruta_out, paste0("seleccion_stepwise", sufijo, ".rds"))
)


# ==============================================================================
# FINAL
# ==============================================================================

cat("\n\n============================================================\n")
cat("SELECCION SECUENCIAL COMPLETADA\n")
cat("============================================================\n")
cat("Universo:              ", UNIVERSO, "\n", sep = "")
cat("Direccion:             ", DIRECCION, "\n", sep = "")
cat("Candidatas evaluadas:  ", length(CANDIDATAS), "\n", sep = "")
cat("Criterios ejecutados:  ", paste(names(resultados), collapse = ", "),
    "\n", sep = "")
cat("Ajustes por criterio:  ",
    paste(names(ajustes_por_criterio), ajustes_por_criterio,
          sep = "=", collapse = ", "), "\n", sep = "")
cat("Tiempo total:          ",
    round(as.numeric(difftime(Sys.time(), tiempo_total, units = "mins")), 2),
    " minutos\n", sep = "")
cat("============================================================\n")


f6 <- pobreza_monetaria ~ alcantarillado_2018 + pct_rural +
  partos_calificados_2020 + transferencias_pc +
  saber11_lectura_2022 + tasa_violencia_intra_2019

f7 <- update(f6, . ~ . + pct_secundario_2023)

arg <- list(vardir = "varianza_pobreza", combined_data = datos,
            domains = "cod_mun", method = "ml", MSE = FALSE,
            B = c(0, 500))

sapply(c(123, 2906, 7, 42, 2024), function(s) {
  a6 <- do.call(fh, c(list(f6), arg, list(seed = s)))
  a7 <- do.call(fh, c(list(f7), arg, list(seed = s)))
  a6$model$model_select$KICc - a7$model$model_select$KICc
})

args(emdi::fh)

?fh
browseVignettes("emdi")
browseVignettes("emdi")

file.path(find.package("emdi"), "doc", "vignette_fh.pdf")
