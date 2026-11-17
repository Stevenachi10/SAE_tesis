# ==============================================================================
# 02_1e_estabilidad_seleccion_kicb2.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.1e: Estabilidad de la selección de M1 (KICb2). Revisión, punto 20.
#
# ------------------------------------------------------------------------------
# OBJETIVO
# ------------------------------------------------------------------------------
# Documentar la estabilidad de la selección con KICb2, sin reemplazar el
# procedimiento utilizado. Se repite la búsqueda de 02_1c (KICb2, ambas
# direcciones, VIF < 10, mismo universo de 76 candidatas, ML, emdi::fh) en
# tres sensibilidades:
#
#   1. Semillas:    629, 2960, 296 y 6290, con B = 500 y los 84 municipios.
#   2. Réplicas:    B = 750 y B = 1000, con semilla 2906 y los 84 municipios.
#   3. Submuestreo: 50 submuestras de 68 municipios (34 por departamento),
#                   sin reposición, con semilla 2906 y B = 500. Se reporta la
#                   frecuencia con que se selecciona cada covariable.
#
# Además se corre la configuración original (semilla 2906, B = 500, 84
# municipios) para verificar que este script reproduce M1.
#
# ------------------------------------------------------------------------------
# SUBMUESTREO SIN REPOSICIÓN
# ------------------------------------------------------------------------------
# En Fay-Herriot cada municipio es un dominio con su propio efecto aleatorio.
# Con reposición, un municipio repetido contaría como dos dominios distintos
# con datos idénticos, lo que distorsiona la estimación de sigma2_u. Por eso se
# submuestrea sin reposición (De Bin et al., 2016), estratificando por
# departamento para conservar el balance entre Cauca y Valle del Cauca. Las 50
# submuestras se sortean una sola vez con semilla 2906 y se guardan.
#
# ------------------------------------------------------------------------------
# INTERRUPCIÓN Y REANUDACIÓN
# ------------------------------------------------------------------------------
# - Cada corrida terminada se guarda en output/estabilidad_kicb2/corridas/.
# - Durante cada corrida, cada especificación evaluada se guarda en
#   output/estabilidad_kicb2/cache/. Como emdi::fh fija la semilla en cada
#   ajuste, el criterio es una función determinista de la especificación: al
#   reanudar, la búsqueda recorre el mismo camino con los valores guardados y
#   continúa desde donde iba.
#
# PARA DETENER: crear un archivo vacío llamado DETENER (cualquier extensión,
# por ejemplo DETENER.txt) dentro de output/estabilidad_kicb2/. Todos los
# procesos se detienen en menos de un minuto, sin perder lo avanzado.
#
# PARA REANUDAR: volver a ejecutar el script completo. Al iniciar, borra el
# archivo DETENER si existe.
#
# PARA VER EL AVANCE: abrir los archivos de output/estabilidad_kicb2/log/.
# Cada corrida escribe ahí cada paso aceptado.
#
# ------------------------------------------------------------------------------
# PARALELIZACIÓN
# ------------------------------------------------------------------------------
# Las corridas se reparten entre N_NUCLEOS procesos de R. Cada proceso toma una
# corrida completa; cuando termina, toma la siguiente pendiente.
#
# ------------------------------------------------------------------------------
# ENTRADAS
# ------------------------------------------------------------------------------
#   output/matriz_sae_transformada_v2.rds
#   output/seleccion_stepwise_completo_both.rds   (universo y M1)
#
# ------------------------------------------------------------------------------
# SALIDAS (output/estabilidad_kicb2/)
# ------------------------------------------------------------------------------
#   submuestras.rds
#   corridas/<id>.rds
#   estabilidad_semillas_B.csv
#   estabilidad_frecuencias.csv
#   estabilidad_modelos.csv
#   estabilidad_resumen.rds
# ==============================================================================


# ==============================================================================
# PAQUETES
# ==============================================================================

library(emdi)
library(car)
library(dplyr)
library(here)
library(parallel)

select <- dplyr::select
filter <- dplyr::filter


# ==============================================================================
# CONFIGURACIÓN
# ==============================================================================

SEMILLA      <- 2906                       # semilla principal
SEMILLAS_ALT <- c(629, 2960, 296, 6290)    # 0629 y 0296 equivalen a 629 y 296
B_BASE       <- 500
B_ALT        <- c(750, 1000)
N_SUB        <- 50                         # número de submuestras
N_POR_DPTO   <- 34                         # municipios por departamento

# Por defecto deja 2 núcleos libres (portátil). En una máquina dedicada se
# puede indicar otro valor al lanzar el script: N_NUCLEOS=32 Rscript ...
N_NUCLEOS <- as.integer(Sys.getenv("N_NUCLEOS",
                                   max(1, parallel::detectCores() - 2)))

# TRUE: no ejecuta corridas, solo resume las que ya terminaron.
SOLO_RESUMEN <- FALSE

# Idénticos a 02_1c
CRITERIO         <- "KICb2"
METODO_SELECCION <- "ml"
TOL_MEJORA       <- 1e-6
VIF_MAX          <- 10
COR_MAX          <- 0.99

NO_COVARIABLES <- c(
  "cod_mun", "Municipio", "cod_dpto", "Departamento",
  "pobreza_monetaria", "varianza_pobreza",
  "ee_pobreza", "cv_directo", "n_eff", "n_muestra",
  "error_estandar", "cvlog"
)


# ==============================================================================
# RUTAS
# ==============================================================================

ruta_out <- here("output")

DIR_BASE     <- file.path(ruta_out, "estabilidad_kicb2")
DIR_CORRIDAS <- file.path(DIR_BASE, "corridas")
DIR_CACHE    <- file.path(DIR_BASE, "cache")
DIR_LOG      <- file.path(DIR_BASE, "log")

for (d in c(DIR_BASE, DIR_CORRIDAS, DIR_CACHE, DIR_LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Si quedó una orden de detención de una sesión anterior, se elimina.
unlink(list.files(DIR_BASE, pattern = "^DETENER", full.names = TRUE))

matriz  <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))
corrida <- readRDS(file.path(ruta_out, "seleccion_stepwise_completo_both.rds"))


# ==============================================================================
# UNIVERSO DE CANDIDATAS
# ------------------------------------------------------------------------------
# Mismos filtros de 02_1c. Debe coincidir con el universo de la corrida
# original; si no coincide, el script se detiene.
# ==============================================================================

numericas  <- names(matriz)[sapply(matriz, is.numeric)]
CANDIDATAS <- setdiff(numericas, NO_COVARIABLES)

constantes <- CANDIDATAS[sapply(CANDIDATAS, function(v) {
  x <- matriz[[v]]
  length(unique(x[!is.na(x)])) < 2
})]
CANDIDATAS <- setdiff(CANDIDATAS, constantes)

con_na <- CANDIDATAS[sapply(CANDIDATAS, function(v) any(is.na(matriz[[v]])))]
CANDIDATAS <- setdiff(CANDIDATAS, con_na)

if (length(CANDIDATAS) > 1) {
  M  <- as.matrix(matriz[, CANDIDATAS, drop = FALSE])
  cm <- abs(cor(M, use = "pairwise.complete.obs"))
  cm[!upper.tri(cm)] <- 0
  redundantes <- character(0)
  for (j in seq_along(CANDIDATAS)) {
    if (CANDIDATAS[j] %in% redundantes) next
    pares <- setdiff(CANDIDATAS[which(cm[j, ] > COR_MAX)], redundantes)
    if (length(pares) > 0) redundantes <- c(redundantes, pares)
  }
  CANDIDATAS <- setdiff(CANDIDATAS, redundantes)
}

X   <- cbind(1, as.matrix(matriz[, CANDIDATAS, drop = FALSE]))
colnames(X)[1] <- "(Intercept)"
qrX <- qr(X)
if (qrX$rank < ncol(X)) {
  CANDIDATAS <- intersect(CANDIDATAS, colnames(X)[qrX$pivot[seq_len(qrX$rank)]])
}

if (!setequal(CANDIDATAS, corrida$candidatas)) {
  stop("El universo reconstruido no coincide con el de la corrida original.")
}

# Mismo orden que la corrida original, para que los empates se resuelvan igual
CANDIDATAS <- corrida$candidatas

VARS_M1 <- strsplit(
  corrida$tabla$covariables[corrida$tabla$criterio == "KICb2"],
  " + ", fixed = TRUE)[[1]]

KICB2_M1_ORIGINAL <- corrida$tabla$kicb2[corrida$tabla$criterio == "KICb2"]


# ==============================================================================
# DATOS Y DEPARTAMENTOS
# ==============================================================================

DATOS_FULL <- matriz %>%
  select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza,
         all_of(CANDIDATAS)) %>%
  as.data.frame()

# Departamento a partir del código DIVIPOLA (19 = Cauca, 76 = Valle del Cauca)
dpto <- substr(sprintf("%05d", as.integer(as.character(DATOS_FULL$cod_mun))),
               1, 2)

tabla_dpto <- table(dpto)
print(tabla_dpto)

stopifnot(
  "Se esperaban dos departamentos." = length(tabla_dpto) == 2,
  "Algún departamento tiene menos municipios que N_POR_DPTO." =
    all(tabla_dpto >= N_POR_DPTO)
)


# ==============================================================================
# SUBMUESTRAS (se sortean una sola vez y se guardan)
# ==============================================================================

archivo_sub <- file.path(DIR_BASE, "submuestras.rds")

if (file.exists(archivo_sub)) {
  
  SUBMUESTRAS <- readRDS(archivo_sub)
  stopifnot("submuestras.rds no tiene N_SUB submuestras." =
              length(SUBMUESTRAS) == N_SUB)
  
} else {
  
  set.seed(SEMILLA)
  
  SUBMUESTRAS <- lapply(seq_len(N_SUB), function(k) {
    idx <- unlist(lapply(split(seq_len(nrow(DATOS_FULL)), dpto),
                         function(i) sample(i, N_POR_DPTO)))
    DATOS_FULL$cod_mun[sort(idx)]
  })
  
  saveRDS(SUBMUESTRAS, archivo_sub)
}


# ==============================================================================
# LISTA DE CORRIDAS
# ------------------------------------------------------------------------------
# Las más largas primero (B mayor), para repartir mejor la carga.
# ==============================================================================

CORRIDAS <- rbind(
  data.frame(id = paste0("B", rev(B_ALT), "_s", SEMILLA), tipo = "B",
             semilla = SEMILLA, B = rev(B_ALT), sub = NA_integer_),
  data.frame(id = paste0("base_s", SEMILLA, "_B", B_BASE), tipo = "base",
             semilla = SEMILLA, B = B_BASE, sub = NA_integer_),
  data.frame(id = paste0("semilla_", SEMILLAS_ALT, "_B", B_BASE),
             tipo = "semilla", semilla = SEMILLAS_ALT, B = B_BASE,
             sub = NA_integer_),
  data.frame(id = sprintf("sub_%02d", seq_len(N_SUB)), tipo = "submuestra",
             semilla = SEMILLA, B = B_BASE, sub = seq_len(N_SUB)),
  stringsAsFactors = FALSE
)


# ==============================================================================
# FUNCIONES
# ==============================================================================

formula_de <- function(vars) {
  if (length(vars) == 0) return(as.formula("pobreza_monetaria ~ 1"))
  as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
}

clave_de <- function(vars) {
  if (length(vars) == 0) return("<intercepto>")
  paste(sort(vars), collapse = "|")
}

escalar <- function(x) {
  v <- tryCatch(as.numeric(x), error = function(e) NA_real_)
  if (length(v) != 1 || !is.finite(v)) return(NA_real_)
  v
}

# Igual que 02_1c, pero sobre los datos de la corrida (submuestra o completos)
vif_de <- function(vars, datos) {
  if (length(vars) < 2) return(1)
  tryCatch({
    v <- car::vif(lm(formula_de(vars), data = datos))
    if (is.matrix(v)) {
      col <- if ("GVIF^(1/(2*Df))" %in% colnames(v)) "GVIF^(1/(2*Df))" else
        colnames(v)[1]
      escalar(max(v[, col], na.rm = TRUE))
    } else {
      escalar(max(v, na.rm = TRUE))
    }
  }, error = function(e) NA_real_)
}

# Escritura atómica: si el proceso se corta, el archivo anterior no se daña.
guardar_atomico <- function(objeto, archivo) {
  tmp <- paste0(archivo, ".", Sys.getpid(), ".tmp")
  saveRDS(objeto, tmp)
  if (!file.rename(tmp, archivo)) {
    file.copy(tmp, archivo, overwrite = TRUE)
    unlink(tmp)
  }
}

detener_solicitado <- function() {
  length(list.files(DIR_BASE, pattern = "^DETENER")) > 0
}


# ------------------------------------------------------------------------------
# Una corrida completa de la búsqueda (misma lógica que 02_1c, solo KICb2)
# ------------------------------------------------------------------------------

ejecutar_corrida <- function(cfg) {
  
  archivo_res   <- file.path(DIR_CORRIDAS, paste0(cfg$id, ".rds"))
  archivo_cache <- file.path(DIR_CACHE,    paste0(cfg$id, ".rds"))
  archivo_log   <- file.path(DIR_LOG,      paste0(cfg$id, ".log"))
  
  if (file.exists(archivo_res)) return(paste("ya terminada:", cfg$id))
  
  log_ <- function(...) {
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n",
        sep = "", file = archivo_log, append = TRUE)
  }
  
  datos <- if (is.na(cfg$sub)) {
    DATOS_FULL
  } else {
    DATOS_FULL[DATOS_FULL$cod_mun %in% SUBMUESTRAS[[cfg$sub]], , drop = FALSE]
  }
  rownames(datos) <- NULL
  
  # ---- Memoria persistente ---------------------------------------------------
  memoria <- new.env(parent = emptyenv())
  
  if (file.exists(archivo_cache)) {
    previa <- readRDS(archivo_cache)
    for (k in names(previa)) assign(k, previa[[k]], envir = memoria)
  }
  
  n_cache_inicial <- length(ls(memoria))
  n_nuevos        <- 0L
  t0              <- Sys.time()
  
  evaluar <- function(vars) {
    
    k <- clave_de(vars)
    
    if (exists(k, envir = memoria, inherits = FALSE)) {
      return(get(k, envir = memoria, inherits = FALSE))
    }
    
    if (detener_solicitado()) stop("detenida por archivo DETENER")
    
    aj <- tryCatch(
      withCallingHandlers(
        fh(formula_de(vars),
           vardir        = "varianza_pobreza",
           combined_data = datos,
           domains       = "cod_mun",
           method        = METODO_SELECCION,
           MSE           = FALSE,
           B             = c(0, cfg$B),
           seed          = cfg$semilla),
        warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) NULL)
    
    v <- if (is.null(aj)) NA_real_ else
      tryCatch(escalar(aj$model$model_select[[CRITERIO]]),
               error = function(e) NA_real_)
    
    assign(k, v, envir = memoria)
    n_nuevos <<- n_nuevos + 1L
    guardar_atomico(as.list(memoria), archivo_cache)
    v
  }
  
  log_("Inicio | semilla ", cfg$semilla, " | B ", cfg$B,
       " | municipios ", nrow(datos),
       " | especificaciones ya guardadas: ", n_cache_inicial)
  
  # ---- Búsqueda --------------------------------------------------------------
  seleccionadas <- character(0)
  actual        <- evaluar(character(0))
  
  if (is.na(actual)) stop("No se pudo calcular el KICb2 del modelo nulo.")
  
  log_("Paso 0 | KICb2 = ", round(actual, 4))
  
  historial <- list()
  movs_max  <- 3 * length(CANDIDATAS)
  
  for (paso in seq_len(movs_max)) {
    
    candidatos <- list()
    
    # Incorporaciones, con filtro de VIF previo al ajuste
    for (variable in setdiff(CANDIDATAS, seleccionadas)) {
      vars_c <- c(seleccionadas, variable)
      vm <- vif_de(vars_c, datos)
      if (is.na(vm) || vm >= VIF_MAX) next
      v <- evaluar(vars_c)
      if (is.na(v)) next
      candidatos[[paste0("+", variable)]] <- list(
        vars = vars_c, valor = v, mov = "+", variable = variable)
    }
    
    # Eliminaciones, desde dos covariables seleccionadas
    if (length(seleccionadas) >= 2) {
      for (variable in seleccionadas) {
        vars_c <- setdiff(seleccionadas, variable)
        v <- evaluar(vars_c)
        if (is.na(v)) next
        candidatos[[paste0("-", variable)]] <- list(
          vars = vars_c, valor = v, mov = "-", variable = variable)
      }
    }
    
    if (length(candidatos) == 0) break
    
    valores <- vapply(candidatos, function(x) x$valor, numeric(1))
    mejor   <- candidatos[[which.min(valores)]]
    
    if ((actual - mejor$valor) <= TOL_MEJORA) break
    
    historial[[length(historial) + 1]] <- data.frame(
      paso       = paso,
      movimiento = mejor$mov,
      variable   = mejor$variable,
      n_vars     = length(mejor$vars),
      kicb2      = mejor$valor,
      mejora     = actual - mejor$valor,
      stringsAsFactors = FALSE)
    
    seleccionadas <- mejor$vars
    actual        <- mejor$valor
    
    log_("Paso ", paso, " | ", mejor$mov, " ", mejor$variable,
         " | KICb2 = ", round(actual, 4),
         " | ajustes nuevos en esta sesión: ", n_nuevos)
  }
  
  resultado <- list(
    id             = cfg$id,
    tipo           = cfg$tipo,
    semilla        = cfg$semilla,
    B              = cfg$B,
    sub            = cfg$sub,
    n_municipios   = nrow(datos),
    municipios     = datos$cod_mun,
    variables      = seleccionadas,
    kicb2          = actual,
    historial      = if (length(historial) > 0) do.call(rbind, historial) else
      data.frame(),
    n_ajustes      = length(ls(memoria)),
    minutos_ultima_sesion = round(as.numeric(
      difftime(Sys.time(), t0, units = "mins")), 1)
  )
  
  guardar_atomico(resultado, archivo_res)
  
  log_("Fin | ", length(seleccionadas), " covariables: ",
       paste(seleccionadas, collapse = ", "))
  
  paste("terminada:", cfg$id)
}

correr_seguro <- function(cfg) {
  tryCatch(
    ejecutar_corrida(cfg),
    error = function(e) {
      cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | DETENIDA O ERROR: ",
          conditionMessage(e), "\n", sep = "",
          file = file.path(DIR_LOG, paste0(cfg$id, ".log")), append = TRUE)
      paste0("sin terminar: ", cfg$id, " (", conditionMessage(e), ")")
    })
}


# ==============================================================================
# EJECUCIÓN
# ==============================================================================

terminadas <- file.exists(file.path(DIR_CORRIDAS, paste0(CORRIDAS$id, ".rds")))
pendientes <- CORRIDAS[!terminadas, , drop = FALSE]

cat("\n============================================================\n")
cat("ESTABILIDAD DE LA SELECCIÓN KICb2\n")
cat("============================================================\n")
cat("Corridas totales:    ", nrow(CORRIDAS), "\n")
cat("Ya terminadas:       ", sum(terminadas), "\n")
cat("Pendientes:          ", nrow(pendientes), "\n")
cat("Núcleos a usar:      ", min(N_NUCLEOS, max(1, nrow(pendientes))), "\n")
cat("Avance en:           ", DIR_LOG, "\n")
cat("Para detener, crear: ", file.path(DIR_BASE, "DETENER.txt"), "\n")
cat("============================================================\n\n")

if (!SOLO_RESUMEN && nrow(pendientes) > 0) {
  
  lista <- lapply(seq_len(nrow(pendientes)),
                  function(i) as.list(pendientes[i, , drop = FALSE]))
  
  n_nuc <- min(N_NUCLEOS, length(lista))
  t_ini <- Sys.time()
  
  if (n_nuc > 1) {
    
    cl <- parallel::makeCluster(n_nuc)
    
    salida <- tryCatch({
      
      parallel::clusterExport(cl, c(
        "DATOS_FULL", "SUBMUESTRAS", "CANDIDATAS", "CRITERIO",
        "METODO_SELECCION", "TOL_MEJORA", "VIF_MAX",
        "DIR_BASE", "DIR_CORRIDAS", "DIR_CACHE", "DIR_LOG",
        "formula_de", "clave_de", "escalar", "vif_de", "guardar_atomico",
        "detener_solicitado", "ejecutar_corrida", "correr_seguro"))
      
      parallel::clusterEvalQ(cl, {
        suppressPackageStartupMessages({
          library(emdi)
          library(car)
        })
        NULL
      })
      
      parallel::parLapplyLB(cl, lista, correr_seguro, chunk.size = 1)
      
    }, finally = parallel::stopCluster(cl))
    
  } else {
    
    salida <- lapply(lista, correr_seguro)
  }
  
  cat("Sesión terminada en",
      round(as.numeric(difftime(Sys.time(), t_ini, units = "hours")), 2),
      "horas\n\n")
  cat(unlist(salida), sep = "\n")
}


# ==============================================================================
# RESUMEN (con las corridas terminadas hasta ahora)
# ==============================================================================

archivos <- file.path(DIR_CORRIDAS, paste0(CORRIDAS$id, ".rds"))
archivos <- archivos[file.exists(archivos)]

if (length(archivos) == 0) {
  stop("Todavía no hay corridas terminadas para resumir.")
}

res <- lapply(archivos, readRDS)

tabla <- bind_rows(lapply(res, function(r) {
  data.frame(
    id           = r$id,
    tipo         = r$tipo,
    semilla      = r$semilla,
    B            = r$B,
    sub          = r$sub,
    n_municipios = r$n_municipios,
    n_vars       = length(r$variables),
    kicb2        = round(r$kicb2, 4),
    igual_M1     = setequal(r$variables, VARS_M1),
    comunes_M1   = length(intersect(r$variables, VARS_M1)),
    faltan_de_M1 = paste(setdiff(VARS_M1, r$variables), collapse = ", "),
    nuevas       = paste(setdiff(r$variables, VARS_M1), collapse = ", "),
    variables    = paste(r$variables, collapse = " + "),
    stringsAsFactors = FALSE)
}))

# ---- 1. Verificación y sensibilidades con los 84 municipios ------------------

cat("\n============================================================\n")
cat("1. CONFIGURACIÓN ORIGINAL, SEMILLAS Y B (84 municipios)\n")
cat("============================================================\n")
cat("M1: ", paste(VARS_M1, collapse = ", "), "\n", sep = "")
cat("KICb2 de M1 en 02_1c: ", KICB2_M1_ORIGINAL, "\n\n", sep = "")

tabla_84 <- tabla %>% filter(tipo != "submuestra")

print(tabla_84 %>%
        select(id, semilla, B, n_vars, kicb2, igual_M1, comunes_M1,
               faltan_de_M1, nuevas),
      row.names = FALSE)

base <- tabla_84 %>% filter(tipo == "base")

if (nrow(base) == 1) {
  if (base$igual_M1 && abs(base$kicb2 - KICB2_M1_ORIGINAL) < 1e-3) {
    cat("\nVerificación: la configuración original reproduce M1 y su KICb2.\n")
  } else {
    warning("La configuración original NO reproduce M1. Revisar antes de ",
            "interpretar las demás corridas.")
  }
}

# ---- 2. Submuestreo -----------------------------------------------------------

sub_res <- Filter(function(r) r$tipo == "submuestra", res)
n_sub   <- length(sub_res)

frecuencias <- NULL
modelos     <- NULL

if (n_sub > 0) {
  
  cat("\n============================================================\n")
  cat("2. SUBMUESTREO: ", n_sub, " de ", N_SUB, " submuestras terminadas\n",
      sep = "")
  cat("============================================================\n")
  
  if (n_sub < N_SUB) {
    cat("ATENCIÓN: resultados parciales; reportar solo con las ", N_SUB,
        " submuestras completas.\n", sep = "")
  }
  
  todas <- unlist(lapply(sub_res, `[[`, "variables"))
  
  frecuencias <- data.frame(variable = union(VARS_M1, unique(todas)),
                            stringsAsFactors = FALSE) %>%
    mutate(
      n_selecciones = vapply(variable, function(v) sum(todas == v), numeric(1)),
      porcentaje    = round(100 * n_selecciones / n_sub, 1),
      en_M1         = variable %in% VARS_M1
    ) %>%
    arrange(desc(n_selecciones), variable)
  
  cat("\nFrecuencia de selección por covariable:\n\n")
  print(frecuencias, row.names = FALSE)
  
  claves <- vapply(sub_res, function(r) clave_de(r$variables), character(1))
  
  modelos <- as.data.frame(table(modelo = claves), stringsAsFactors = FALSE) %>%
    rename(n = Freq) %>%
    mutate(porcentaje = round(100 * n / n_sub, 1),
           es_M1      = modelo == clave_de(VARS_M1)) %>%
    arrange(desc(n))
  
  cat("\nModelos distintos seleccionados: ", nrow(modelos), "\n", sep = "")
  cat("Veces que se seleccionó M1 exacto: ",
      sum(claves == clave_de(VARS_M1)), " de ", n_sub, "\n", sep = "")
  
  tam <- vapply(sub_res, function(r) length(r$variables), numeric(1))
  cat("Número de covariables seleccionadas: mediana ", median(tam),
      ", mínimo ", min(tam), ", máximo ", max(tam), "\n", sep = "")
  
  comunes <- vapply(sub_res, function(r) length(intersect(r$variables, VARS_M1)),
                    numeric(1))
  cat("Covariables de M1 presentes por submuestra: mediana ", median(comunes),
      " de 6\n", sep = "")
  
  cat("\nModelos más frecuentes:\n\n")
  print(head(modelos, 10), row.names = FALSE)
}


# ==============================================================================
# GUARDAR RESUMEN
# ==============================================================================

write.csv(tabla_84, file.path(DIR_BASE, "estabilidad_semillas_B.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

if (!is.null(frecuencias)) {
  write.csv(frecuencias, file.path(DIR_BASE, "estabilidad_frecuencias.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(modelos, file.path(DIR_BASE, "estabilidad_modelos.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}

saveRDS(
  list(
    tabla           = tabla,
    frecuencias     = frecuencias,
    modelos         = modelos,
    n_submuestras   = n_sub,
    vars_m1         = VARS_M1,
    configuracion   = list(semilla = SEMILLA, semillas_alt = SEMILLAS_ALT,
                           B_base = B_BASE, B_alt = B_ALT, n_sub = N_SUB,
                           n_por_dpto = N_POR_DPTO, vif_max = VIF_MAX,
                           tol = TOL_MEJORA, candidatas = CANDIDATAS),
    session         = sessionInfo()
  ),
  file.path(DIR_BASE, "estabilidad_resumen.rds")
)

cat("\nResumen guardado en ", DIR_BASE, "\n", sep = "")