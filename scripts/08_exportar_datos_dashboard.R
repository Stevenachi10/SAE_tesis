# ==============================================================================
# 08_exportar_datos_dashboard.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# EXPORTACIÓN DE DATOS PARA EL TABLERO INTERACTIVO
#
# Reúne en UN SOLO archivo (output/dashboard_datos.rds) lo que muestra el
# tablero, para todos los modelos estimados en la Etapa II:
#
#   Fay-Herriot, REML           M1, M2, M3 en escala original, log y logit   9
#   Fay-Herriot espacial, REML  M1, M2, M3 (queen, orden 1)                  3
#   Fay-Herriot robusto         REBLUP y REBLUPbc, k = 1,345, M1, M2, M3     6
#                                                                     Total 18
#
# M3s (M3 sin la variable de distancia) no se incluye: en la tesis es solo un
# contraste auxiliar, no una especificación candidata.
#
# ------------------------------------------------------------------------------
# PRINCIPIOS
# ------------------------------------------------------------------------------
# 1. El script solo LEE los objetos de la Etapa II y los CSV de la Etapa III.
#    No reajusta ningún modelo ni recalcula ninguna prueba.
#
# 2. Solo se exporta lo que el tablero muestra. El archivo de datos queda
#    incrustado completo dentro del HTML, de modo que cualquier columna de más
#    queda al alcance de quien abra el código de la página.
#
# 3. No se exporta la referencia municipal DANE ni nada que permita
#    reconstruirla municipio por municipio. Por eso:
#      - no salen el factor de encogimiento (gamma), el MSE ni el error
#        estándar por municipio; solo la estimación, el CV y el intervalo;
#      - los efectos aleatorios y los residuos de los QQ-plots salen
#        estandarizados (centrados y divididos por su desviación estándar),
#        que es como los grafica emdi. Con las series sin estandarizar y la
#        varianza del efecto de área se podían despejar la referencia y su
#        varianza de muestreo.
#
# 4. Cada objeto de ajuste se contrasta con lo que dice su rótulo (método,
#    transformación, estructura espacial y constante k). Si algo no coincide
#    se avisa en la consola: un modelo mal rotulado es el peor error posible
#    en un tablero de comparación.
#
# ------------------------------------------------------------------------------
# CONTENIDO DEL ARCHIVO (lista; todas las tablas son data.frame de tipos base)
# ------------------------------------------------------------------------------
#   meta              Fecha, semilla, regla del intervalo, réplicas bootstrap,
#                     constante k y versiones de R, emdi y saeRobust.
#   municipios        cod_mun, municipio, departamento (84 filas).
#   modelos           Una fila por modelo: identificador, rótulo, familia,
#                     especificación, escala, método, k, estimador del MSE y sus
#                     réplicas, sigma2_u, rho, convergencia, estado en la Etapa
#                     III y resumen del CV.
#   estimaciones      Formato largo: modelo x municipio. Estimación, CV (%) e
#                     intervalo al 95 %.
#   coeficientes      Formato largo: modelo x término. Coeficiente y error
#                     estándar. No salen z ni el valor p: el tablero no los
#                     muestra.
#   especificaciones  Covariables de M1, M2 y M3 y criterio de información con
#                     que se seleccionó cada una.
#   qq                Formato largo: residuos y efectos aleatorios
#                     estandarizados por municipio, solo de los nueve
#                     Fay-Herriot (el tablero no hace este diagnóstico en los
#                     espaciales ni en los robustos).
#   normalidad        Asimetría, curtosis y Shapiro-Wilk de los nueve
#                     Fay-Herriot (de output/sup_normalidad.csv), si existe.
#   pruebas           Tablas a nivel de modelo de las pruebas de concordancia
#                     con la referencia municipal DANE (CSV de 03_1), si
#                     existen: sesgo, bondad de ajuste y solapamiento de
#                     intervalos. Solo las columnas que el tablero explica.
#
# ------------------------------------------------------------------------------
# ENTRADAS (output/)
# ------------------------------------------------------------------------------
#   modelado_variantes.rds, espacial.rds, robusto.rds      (obligatorias)
#   sup_normalidad.csv, prec_bias_brown.csv, prec_brown_gof.csv,
#   prec_coverage.csv, prec_modelos_metadata.csv           (opcionales)
#
# SALIDA (output/)
#   dashboard_datos.rds
#
# El script solo usa R base.
# ==============================================================================


# ==============================================================================
# 0. CONFIGURACIÓN
# ==============================================================================

ruta_out <- if (requireNamespace("here", quietly = TRUE)) {
  here::here("output")
} else {
  file.path(getwd(), "output")
}

# Nivel del intervalo. Es el mismo del CSV final de M1 log: estimación más o
# menos 1,96 errores estándar, con el error estándar tomado del MSE que
# devuelve cada ajuste y recortado al intervalo [0, 1].
Z_IC <- 1.96

# Constante de sintonización del modelo robusto que se muestra. Si robusto.rds
# registra otra constante principal, se avisa.
K_ROBUSTO <- 1.345

# Modelo que se marca como especificación final.
MODELO_FINAL <- "FH_M1_log"

# Especificaciones que entran al tablero (M3s queda fuera).
ESPECIFICACIONES <- c("M1", "M2", "M3")

DEPARTAMENTOS <- c("19" = "Cauca", "76" = "Valle del Cauca")

# Decimales con que se guardan las series. El tablero muestra un decimal de
# porcentaje; guardar más de lo necesario solo engorda el HTML.
DEC_PROPORCION <- 6   # estimación e intervalo, como proporción (0,0001 pp)
DEC_CV         <- 4   # CV en porcentaje
DEC_QQ         <- 4   # series estandarizadas de los QQ-plots

# Archivos CSV de pruebas a nivel de modelo, con las columnas que se exportan.
# Solo salen estadísticos por modelo, nunca series por municipio.
CSV_PRUEBAS <- list(
  bias_brown = list(archivo = "prec_bias_brown.csv",
                    columnas = c("b0", "b1", "ee_b1", "R2", "p_conjunto",
                                 "p_conjunto_hc3", "gq_p_psi")),
  gof_brown  = list(archivo = "prec_brown_gof.csv",
                    columnas = c("W", "c2", "p_value", "results")),
  cobertura  = list(archivo = "prec_coverage.csv",
                    columnas = c("non_coverage", "domains", "non_overlap",
                                 "p_value", "results"))
)

# Metadatos de la Etapa III: de aquí salen el estado de cada modelo, los grados
# de libertad de la prueba de bondad de ajuste y, si el objeto de ajuste no lo
# trae, el tipo de MSE. La tabla en sí no se exporta.
CSV_METADATA <- "prec_modelos_metadata.csv"

# Columnas de la tabla de normalidad que usa el tablero.
COLUMNAS_NORMALIDAD <- c("ef_asimetria", "ef_curtosis", "ef_shapiro_W",
                         "ef_shapiro_p", "res_asimetria", "res_curtosis",
                         "res_shapiro_W", "res_shapiro_p")

# Ninguna tabla exportada puede traer columnas con estos nombres.
PATRON_REFERENCIA <- "direct|dane|referencia|gamma|psi|vardir"


# ==============================================================================
# 1. LECTURA
# ==============================================================================

leer_rds <- function(nombre) {
  f <- file.path(ruta_out, nombre)
  if (!file.exists(f)) {
    stop("No existe ", f, ". Ejecute antes el script que lo genera.")
  }
  readRDS(f)
}

variantes <- leer_rds("modelado_variantes.rds")
espacial  <- leer_rds("espacial.rds")
robusto   <- leer_rds("robusto.rds")

datos <- variantes$datos
cod   <- as.character(datos$cod_mun)
n_dom <- length(cod)

stopifnot(
  "Hay codigos municipales duplicados" = !any(duplicated(cod)),
  "El orden de dominios de espacial.rds no coincide con modelado_variantes.rds" =
    identical(as.character(espacial$cod_mun), cod),
  "El orden de dominios de robusto.rds no coincide con modelado_variantes.rds" =
    identical(as.character(robusto$datos$cod_mun), cod)
)

cat("\nDominios:", n_dom, "\n")

# Avisos de la exportación: se muestran al final, todos juntos.
avisos <- character(0)

avisar <- function(...) {
  msg <- paste0(...)
  avisos <<- c(avisos, msg)
  cat("  AVISO: ", msg, "\n", sep = "")
}

# Constante k y réplicas bootstrap registradas en la Etapa II.
k_rds <- tryCatch(as.numeric(robusto$k_principal)[1], error = function(e) NA_real_)

if (length(k_rds) == 1 && is.finite(k_rds) && !isTRUE(all.equal(k_rds, K_ROBUSTO))) {
  avisar("robusto.rds registra k principal = ", k_rds, " y este script usa ",
         K_ROBUSTO, ". Revise K_ROBUSTO.")
}

numero_o_na <- function(x) {
  v <- tryCatch(as.numeric(x)[1], error = function(e) NA_real_)
  if (length(v) != 1 || !is.finite(v)) NA_real_ else v
}

B_LOGIT   <- numero_o_na(variantes$B_mse)
B_ROBUSTO <- numero_o_na(robusto$b_boot)


# ==============================================================================
# 2. FUNCIONES
# ==============================================================================

escalar <- function(x) {
  v <- tryCatch(as.numeric(x), error = function(e) NA_real_)
  if (length(v) < 1 || !is.finite(v[1])) NA_real_ else v[1]
}

# Columna de un data.frame o NA si no existe.
col_o_na <- function(df, nombre, n) {
  if (is.data.frame(df) && nombre %in% names(df)) as.numeric(df[[nombre]])
  else rep(NA_real_, n)
}

# Texto de un elemento del objeto de ajuste, o NA si no existe.
texto_o_na <- function(x) {
  v <- tryCatch(as.character(x)[1], error = function(e) NA_character_)
  if (length(v) != 1 || is.na(v) || !nzchar(v)) NA_character_ else v
}

# Parámetros de varianza. En el modelo espacial model$variance trae los
# componentes por nombre; en los demás es un escalar.
componentes <- function(aj) {
  v  <- tryCatch(aj$model$variance, error = function(e) NULL)
  nm <- names(v)

  if (!is.null(nm) && all(c("correlation", "variance") %in% nm)) {
    return(list(s2u = as.numeric(v[["variance"]])[1],
                rho = as.numeric(v[["correlation"]])[1],
                convergencia =
                  if ("convergence" %in% nm) as.logical(v[["convergence"]])[1]
                  else NA))
  }

  list(s2u = escalar(v), rho = NA_real_, convergencia = NA)
}

# Nombres de los términos del objeto de coeficientes. Si el objeto no trae
# nombres utilizables se usa el orden intercepto + covariables.
terminos_de <- function(cf, vars, id) {
  rn <- rownames(cf)

  if (!is.null(rn) && !all(rn == as.character(seq_len(nrow(cf))))) return(rn)

  avisar(id, ": coeficientes sin nombres; se asignan por posicion.")
  c("(Intercept)", vars)[seq_len(nrow(cf))]
}

# Resumen de una serie de CV.
resumen_cv <- function(cv) {
  ok <- is.finite(cv)
  c(cv_medio   = if (any(ok)) mean(cv[ok])   else NA_real_,
    cv_mediana = if (any(ok)) median(cv[ok]) else NA_real_,
    cv_max     = if (any(ok)) max(cv[ok])    else NA_real_)
}

# Serie centrada y escalada con su media y su desviación estándar muestral.
# Es la transformación que aplica plot() de emdi antes de dibujar los QQ-plots.
estandarizar <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# ------------------------------------------------------------------------------
# Contraste del objeto de ajuste con su rótulo
# ------------------------------------------------------------------------------
# El objeto fh de emdi guarda el método, la transformación, la estructura de
# correlación y, en los robustos, la constante k. Si el objeto no trae alguno de
# esos elementos, ese contraste se omite.
verificar_ajuste <- function(aj, id, familia, escala, metodo, k) {

  met <- texto_o_na(aj$method$method)
  if (!is.na(met) && !identical(tolower(met), metodo)) {
    avisar(id, ": el objeto registra el metodo '", met, "' y el rotulo dice '",
           metodo, "'.")
  }

  tr       <- texto_o_na(aj$transformation$transformation)
  esperada <- c(original = "no", log = "log", logit = "logit")[[escala]]
  if (!is.na(tr) && !identical(tr, esperada)) {
    avisar(id, ": el objeto registra la transformacion '", tr,
           "' y el rotulo dice '", escala, "'.")
  }

  co       <- texto_o_na(aj$model$correlation)
  esperada <- if (familia == "SFH") "spatial" else "no"
  if (!is.na(co) && !identical(co, esperada)) {
    avisar(id, ": el objeto registra correlacion '", co, "' y se esperaba '",
           esperada, "'.")
  }

  if (familia == "RFH") {
    kk <- numero_o_na(aj$model$k)
    if (is.finite(kk) && !isTRUE(all.equal(kk, k))) {
      avisar(id, ": el objeto registra k = ", kk, " y el rotulo dice ", k, ".")
    }
  }

  invisible(NULL)
}

# ------------------------------------------------------------------------------
# Procesa un ajuste y devuelve sus tablas
# ------------------------------------------------------------------------------
procesar <- function(aj, id, rotulo, grupo, familia, espec, vars, escala,
                     metodo, k = NA_real_, con_qq = FALSE) {

  verificar_ajuste(aj, id, familia, escala, metodo, k)

  # ---- alineación por dominio: el orden se verifica, no se supone -------------
  dom <- as.character(aj$ind$Domain)
  ord <- match(cod, dom)

  if (anyNA(ord)) {
    stop(id, ": los dominios del ajuste no coinciden con los de 'datos'.")
  }

  est <- as.numeric(aj$ind$FH)[ord]

  mse <- if (is.data.frame(aj$MSE) && "FH" %in% names(aj$MSE)) {
    as.numeric(aj$MSE$FH)[match(cod, as.character(aj$MSE$Domain))]
  } else {
    avisar(id, ": el objeto no trae MSE; quedan sin CV ni intervalo.")
    rep(NA_real_, n_dom)
  }

  se <- sqrt(pmax(mse, 0))

  cv <- ifelse(is.finite(est) & est > 0 & is.finite(se), se / est * 100,
               NA_real_)

  inf_bruto <- est - Z_IC * se
  sup_bruto <- est + Z_IC * se
  inf <- pmax(inf_bruto, 0)
  sup <- pmin(sup_bruto, 1)

  comp <- componentes(aj)

  tipo_mse <- texto_o_na(aj$method$MSE_method)

  # Réplicas del bootstrap, cuando el MSE se obtuvo por ese procedimiento.
  mse_b <- if (!is.na(tipo_mse) && grepl("boot", tipo_mse, ignore.case = TRUE)) {
    if (familia == "RFH") B_ROBUSTO else B_LOGIT
  } else NA_real_

  estimaciones <- data.frame(
    modelo     = id,
    cod_mun    = cod,
    estimacion = round(est, DEC_PROPORCION),
    cv_pct     = round(cv, DEC_CV),
    ic_inf     = round(inf, DEC_PROPORCION),
    ic_sup     = round(sup, DEC_PROPORCION),
    stringsAsFactors = FALSE)

  # ---- coeficientes -----------------------------------------------------------
  cf <- tryCatch(as.data.frame(aj$model$coefficients), error = function(e) NULL)

  coeficientes <- if (is.null(cf) || nrow(cf) == 0) NULL else {
    n <- nrow(cf)
    data.frame(
      modelo  = id,
      termino = terminos_de(cf, vars, id),
      coef    = col_o_na(cf, "coefficients", n),
      ee      = col_o_na(cf, "std.error", n),
      stringsAsFactors = FALSE)
  }

  # Las covariables del ajuste deben ser las de su especificación.
  if (!is.null(coeficientes)) {
    en_ajuste <- setdiff(coeficientes$termino, "(Intercept)")
    if (!setequal(en_ajuste, vars)) {
      avisar(id, ": los terminos del ajuste no coinciden con las covariables ",
             "de ", espec, ".")
    }
  }

  # ---- series estandarizadas de los QQ-plots ----------------------------------
  # random_effects y std_real_residuals siguen el orden de los dominios del
  # ajuste, el mismo de aj$ind; se reordenan igual que las estimaciones.
  qq <- NULL

  if (con_qq) {
    ef <- tryCatch(as.numeric(aj$model$random_effects), error = function(e) NULL)
    rs <- tryCatch(as.numeric(aj$model$std_real_residuals),
                   error = function(e) NULL)

    if (length(ef) == n_dom && length(rs) == n_dom &&
        length(dom) == n_dom) {
      qq <- data.frame(
        modelo  = id,
        cod_mun = cod,
        res_z   = round(estandarizar(rs[ord]), DEC_QQ),
        ef_z    = round(estandarizar(ef[ord]), DEC_QQ),
        stringsAsFactors = FALSE)
    } else {
      avisar(id, ": sin efectos aleatorios o residuos por municipio; ",
             "no se dibujara su QQ-plot.")
    }
  }

  # ---- fila de metadatos -------------------------------------------------------
  rcv <- resumen_cv(cv)

  modelo <- data.frame(
    modelo          = id,
    rotulo          = rotulo,
    grupo           = grupo,
    familia         = familia,
    especificacion  = espec,
    n_covariables   = length(vars),
    escala          = escala,
    metodo          = metodo,
    k               = k,
    mse_tipo        = tipo_mse,
    mse_b           = mse_b,
    sigma2_u        = comp$s2u,
    rho             = comp$rho,
    convergencia    = comp$convergencia,
    es_final        = identical(id, MODELO_FINAL),
    cv_medio        = unname(rcv["cv_medio"]),
    cv_mediana      = unname(rcv["cv_mediana"]),
    cv_max          = unname(rcv["cv_max"]),
    n_ic_recortados = sum(inf_bruto < 0 | sup_bruto > 1, na.rm = TRUE),
    n_fuera_rango   = sum(est < 0 | est > 1, na.rm = TRUE),
    estimado        = TRUE,
    stringsAsFactors = FALSE)

  list(estimaciones = estimaciones, coeficientes = coeficientes,
       qq = qq, modelo = modelo)
}

# Fila de un modelo que no se pudo exportar, para que quede registrado.
fila_ausente <- function(id, rotulo, grupo, familia, espec, vars, escala,
                         metodo, k = NA_real_) {
  data.frame(modelo = id, rotulo = rotulo, grupo = grupo, familia = familia,
             especificacion = espec, n_covariables = length(vars),
             escala = escala, metodo = metodo, k = k, mse_tipo = NA_character_,
             mse_b = NA_real_, sigma2_u = NA_real_, rho = NA_real_,
             convergencia = NA, es_final = FALSE, cv_medio = NA_real_,
             cv_mediana = NA_real_, cv_max = NA_real_,
             n_ic_recortados = NA_integer_, n_fuera_rango = NA_integer_,
             estimado = FALSE, stringsAsFactors = FALSE)
}


# ==============================================================================
# 3. ESPECIFICACIONES
# ==============================================================================

# modelado_variantes.rds guarda cada especificación como lista con etiqueta,
# criterio y covariables; espacial.rds y robusto.rds guardan el vector de
# covariables. Se normaliza a un vector por especificación.
vars_de <- function(x) if (is.list(x) && !is.null(x$vars)) x$vars else x

criterio_de <- function(x) {
  if (is.list(x) && !is.null(x$criterio)) texto_o_na(x$criterio)
  else NA_character_
}

ESP_BASE <- lapply(variantes$especificaciones, vars_de)
ESP_ESP  <- lapply(espacial$especificaciones, vars_de)
ESP_ROB  <- lapply(robusto$especificaciones, vars_de)

faltan_esp <- setdiff(ESPECIFICACIONES, names(ESP_BASE))

if (length(faltan_esp) > 0) {
  stop("modelado_variantes.rds no trae las especificaciones: ",
       paste(faltan_esp, collapse = ", "))
}

# Los tres archivos de la Etapa II deben describir las mismas covariables. Si
# uno quedó de una corrida anterior de la selección, se avisa.
for (nm in ESPECIFICACIONES) {
  for (fuente in c("espacial.rds", "robusto.rds")) {
    otra <- if (fuente == "espacial.rds") ESP_ESP[[nm]] else ESP_ROB[[nm]]
    if (!is.null(otra) && !setequal(otra, ESP_BASE[[nm]])) {
      avisar("Las covariables de ", nm, " en ", fuente,
             " no coinciden con las de modelado_variantes.rds.")
    }
  }
}

especificaciones <- do.call(rbind, lapply(ESPECIFICACIONES, function(nm) {
  v <- ESP_BASE[[nm]]
  data.frame(especificacion = nm,
             criterio = criterio_de(variantes$especificaciones[[nm]]),
             orden = seq_along(v), variable = v, stringsAsFactors = FALSE)
}))


# ==============================================================================
# 4. RECORRIDO DE LOS MODELOS
# ==============================================================================

resultados <- list()

registrar_modelo <- function(r) {
  resultados[[r$modelo$modelo]] <<- r
}

registrar_ausente <- function(...) {
  fila <- fila_ausente(...)
  resultados[[fila$modelo]] <<- list(modelo = fila, estimaciones = NULL,
                                     coeficientes = NULL, qq = NULL)
  avisar("Modelo no disponible: ", fila$modelo)
}

# ---- 4.1 Fay-Herriot, REML, tres escalas ------------------------------------
cat("\n== FAY-HERRIOT ==\n")

rotulo_escala <- c(original = "escala original",
                   log      = "escala logar\u00edtmica",
                   logit    = "escala logit")

for (esp in ESPECIFICACIONES) {
  for (escala in c("original", "log", "logit")) {

    clave  <- paste(esp, "reml", escala, sep = "_")
    id     <- paste("FH", esp, escala, sep = "_")
    rotulo <- paste0("Fay-Herriot (", esp, ", ", rotulo_escala[[escala]], ")")
    vars   <- ESP_BASE[[esp]]
    aj     <- variantes$ajustes[[clave]]

    if (is.null(aj)) {
      registrar_ausente(id, rotulo, "Fay-Herriot", "FH", esp, vars, escala,
                        "reml")
      next
    }

    registrar_modelo(procesar(aj, id, rotulo, "Fay-Herriot", "FH", esp, vars,
                              escala, "reml", con_qq = TRUE))
    cat("  ", id, " -> ok\n", sep = "")
  }
}

# ---- 4.2 Fay-Herriot espacial -------------------------------------------------
cat("\n== FAY-HERRIOT ESPACIAL ==\n")

for (esp in ESPECIFICACIONES) {

  id     <- paste("SFH", esp, sep = "_")
  rotulo <- paste0("Fay-Herriot espacial (", esp, ")")
  vars   <- if (!is.null(ESP_ESP[[esp]])) ESP_ESP[[esp]] else ESP_BASE[[esp]]
  aj     <- espacial$ajustes_esp[[esp]]

  if (is.null(aj)) {
    registrar_ausente(id, rotulo, "Fay-Herriot espacial", "SFH", esp, vars,
                      "original", "reml")
    next
  }

  registrar_modelo(procesar(aj, id, rotulo, "Fay-Herriot espacial", "SFH", esp,
                            vars, "original", "reml"))
  cat("  ", id, " -> ok\n", sep = "")
}

# ---- 4.3 Fay-Herriot robusto --------------------------------------------------
cat("\n== FAY-HERRIOT ROBUSTO (k = ", K_ROBUSTO, ") ==\n", sep = "")

nombre_metodo <- c(reblup = "REBLUP", reblupbc = "REBLUPbc")

for (esp in ESPECIFICACIONES) {
  for (metodo in c("reblup", "reblupbc")) {

    clave  <- paste(esp, metodo, format(K_ROBUSTO), sep = "_")
    id     <- paste("RFH", esp, metodo, sep = "_")
    # La constante k es la misma en los seis robustos: se informa fuera del
    # nombre, en la nota del tablero.
    rotulo <- paste0("Fay-Herriot robusto ", nombre_metodo[[metodo]],
                     " (", esp, ")")
    vars   <- if (!is.null(ESP_ROB[[esp]])) ESP_ROB[[esp]] else ESP_BASE[[esp]]
    aj     <- robusto$ajustes_robustos[[clave]]

    if (is.null(aj)) {
      registrar_ausente(id, rotulo, "Fay-Herriot robusto", "RFH", esp, vars,
                        "original", metodo, K_ROBUSTO)
      next
    }

    registrar_modelo(procesar(aj, id, rotulo, "Fay-Herriot robusto", "RFH",
                              esp, vars, "original", metodo, K_ROBUSTO))
    cat("  ", id, " -> ok\n", sep = "")
  }
}


# ==============================================================================
# 5. ENSAMBLE
# ==============================================================================

unir <- function(campo) {
  piezas <- lapply(resultados, function(r) r[[campo]])
  piezas <- piezas[!vapply(piezas, is.null, logical(1))]
  if (length(piezas) == 0) return(NULL)
  out <- do.call(rbind, piezas)
  rownames(out) <- NULL
  out
}

modelos      <- unir("modelo")
estimaciones <- unir("estimaciones")
coeficientes <- unir("coeficientes")
qq           <- unir("qq")

# Orden de presentación: familia, especificación, escala (original, log, logit)
# y método. El modelo final no se adelanta: se reconoce por su marca.
orden_familia <- c(FH = 1, SFH = 2, RFH = 3)
orden_escala  <- c(original = 1, log = 2, logit = 3)
orden_metodo  <- c(reml = 1, reblup = 2, reblupbc = 3)

modelos <- modelos[order(orden_familia[modelos$familia],
                         modelos$especificacion,
                         orden_escala[modelos$escala],
                         orden_metodo[modelos$metodo]), ]
rownames(modelos) <- NULL
modelos$orden <- seq_len(nrow(modelos))

if (!any(modelos$es_final & modelos$estimado)) {
  avisar("El modelo final (", MODELO_FINAL, ") no esta entre los exportados; ",
         "el tablero no marcara ninguno.")
}

municipios <- data.frame(
  cod_mun      = cod,
  municipio    = as.character(datos$Municipio),
  departamento = unname(DEPARTAMENTOS[substr(cod, 1, 2)]),
  stringsAsFactors = FALSE)

stopifnot("Hay municipios fuera de Cauca y Valle del Cauca" =
            !anyNA(municipios$departamento))


# ==============================================================================
# 6. TABLAS DE PRUEBAS YA CALCULADAS
# ------------------------------------------------------------------------------
# Se leen los CSV de la Etapa III y se conservan solo las columnas que el
# tablero explica.
# ==============================================================================

leer_csv <- function(archivo) {
  f <- file.path(ruta_out, archivo)

  if (!file.exists(f)) {
    cat("  (no existe ", archivo, ")\n", sep = "")
    return(NULL)
  }

  tab <- read.csv(f, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  cat("  ", archivo, ": ", nrow(tab), " filas, ", ncol(tab), " columnas\n",
      sep = "")
  tab
}

# Clave de 03_1 ("M1.log", "Espacial.M2", "reblupbc.M3") a clave del tablero
# ("FH_M1_log", "SFH_M2", "RFH_M3_reblupbc"). Lo que no coincide queda NA.
id_tablero <- function(v) {
  v   <- gsub(" ", ".", trimws(as.character(v)), fixed = TRUE)
  out <- rep(NA_character_, length(v))

  i <- grepl("^M[123]\\.(original|log|logit)$", v)
  out[i] <- sub("^(M[123])\\.(original|log|logit)$", "FH_\\1_\\2", v[i])

  i <- grepl("^Espacial\\.M[123]$", v, ignore.case = TRUE)
  out[i] <- sub("^Espacial\\.(M[123])$", "SFH_\\1", v[i], ignore.case = TRUE)

  i <- grepl("^(reblup|reblupbc)\\.M[123]$", v)
  out[i] <- sub("^(reblup|reblupbc)\\.(M[123])$", "RFH_\\2_\\1", v[i])

  out
}

# Se queda con las columnas que el tablero explica y añade la clave 'modelo'.
curar_prueba <- function(tab, columnas, nombre) {
  if (is.null(tab)) return(NULL)

  clave <- intersect(c("variable", "methods"), names(tab))[1]
  if (is.na(clave)) {
    avisar(nombre, ": no tiene columna 'variable' ni 'methods'; se omite.")
    return(NULL)
  }

  id   <- id_tablero(tab[[clave]])
  cols <- intersect(columnas, names(tab))

  faltan <- setdiff(columnas, names(tab))
  if (length(faltan) > 0) {
    avisar(nombre, ": faltan las columnas ", paste(faltan, collapse = ", "), ".")
  }

  if (anyNA(id)) {
    avisar(nombre, ": filas sin modelo reconocido: ",
           paste(tab[[clave]][is.na(id)], collapse = ", "), ".")
  }

  out <- data.frame(modelo = id, tab[, cols, drop = FALSE],
                    stringsAsFactors = FALSE)
  out <- out[!is.na(out$modelo), , drop = FALSE]

  if (any(duplicated(out$modelo))) {
    avisar(nombre, ": hay modelos repetidos; se conserva la primera fila.")
    out <- out[!duplicated(out$modelo), , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

cat("\n== PRUEBAS ==\n")

pruebas <- lapply(names(CSV_PRUEBAS), function(nm) {
  curar_prueba(leer_csv(CSV_PRUEBAS[[nm]]$archivo),
               CSV_PRUEBAS[[nm]]$columnas, nm)
})
names(pruebas) <- names(CSV_PRUEBAS)
pruebas <- pruebas[!vapply(pruebas, is.null, logical(1))]


# ---- Normalidad: solo los nueve Fay-Herriot -----------------------------------
normalidad_csv <- leer_csv("sup_normalidad.csv")
normalidad     <- NULL

if (!is.null(normalidad_csv) &&
    all(c("familia", "especificacion", "escala") %in% names(normalidad_csv))) {

  nb <- normalidad_csv[tolower(normalidad_csv$familia) == "base" &
                         normalidad_csv$especificacion %in% ESPECIFICACIONES, ,
                       drop = FALSE]

  if ("disponible" %in% names(nb)) {
    nb <- nb[!is.na(nb$disponible) & as.logical(nb$disponible), , drop = FALSE]
  }

  if (nrow(nb) > 0) {
    cols <- intersect(COLUMNAS_NORMALIDAD, names(nb))
    normalidad <- data.frame(
      modelo = paste("FH", nb$especificacion, nb$escala, sep = "_"),
      nb[, cols, drop = FALSE], stringsAsFactors = FALSE)
    normalidad <- normalidad[normalidad$modelo %in% modelos$modelo &
                               !duplicated(normalidad$modelo), , drop = FALSE]
    rownames(normalidad) <- NULL
  }
}


# ---- Estado en la Etapa III, grados de libertad y respaldo del tipo de MSE ---
# El estado se toma de los modelos que 03_1 evaluó:
#   final         especificación final seleccionada
#   evaluado      incluido en la evaluación de concordancia y precisión
#   contraste     incluido en esa evaluación solo como contraste
#   no_evaluado   no incluido en esa evaluación
metadata_e3 <- leer_csv(CSV_METADATA)

modelos$etapa3 <- NA_character_

if (!is.null(metadata_e3) && "variable" %in% names(metadata_e3)) {

  ids_e3   <- id_tablero(metadata_e3$variable)
  papel_e3 <- if ("papel" %in% names(metadata_e3)) {
    as.character(metadata_e3$papel)
  } else rep("candidata", length(ids_e3))

  modelos$etapa3 <- ifelse(
    modelos$es_final, "final",
    ifelse(modelos$modelo %in% ids_e3[papel_e3 == "referencia"], "contraste",
           ifelse(modelos$modelo %in% ids_e3, "evaluado", "no_evaluado")))

  if (!MODELO_FINAL %in% ids_e3) {
    avisar("El modelo final no aparece en ", CSV_METADATA, ".")
  }

  # Grados de libertad con que SAEval contrasta W (m - 1).
  if (!is.null(pruebas$gof_brown) && "gl_w" %in% names(metadata_e3)) {
    pruebas$gof_brown$gl <-
      as.integer(metadata_e3$gl_w[match(pruebas$gof_brown$modelo, ids_e3)])
  }

  # Si el objeto de ajuste no trae el tipo de MSE, se toma el que registró
  # 03_1.
  if ("mse_type" %in% names(metadata_e3)) {
    mse_e3 <- setNames(as.character(metadata_e3$mse_type), ids_e3)
    sin    <- is.na(modelos$mse_tipo) & modelos$modelo %in% names(mse_e3)
    modelos$mse_tipo[sin] <- unname(mse_e3[modelos$modelo[sin]])
    if (any(sin)) {
      cat("  Tipo de MSE tomado de ", CSV_METADATA, " para: ",
          paste(modelos$modelo[sin], collapse = ", "), "\n", sep = "")
    }
  }

} else {
  avisar("Sin ", CSV_METADATA, "; no se marca el estado en la Etapa III.")
}

if (any(is.na(modelos$mse_tipo) & modelos$estimado)) {
  avisar("Modelos sin tipo de MSE: ",
         paste(modelos$modelo[is.na(modelos$mse_tipo) & modelos$estimado],
               collapse = ", "), ".")
}

# Las pruebas solo pueden referirse a modelos exportados.
for (nm in names(pruebas)) {
  sobran <- setdiff(pruebas[[nm]]$modelo, modelos$modelo[modelos$estimado])
  if (length(sobran) > 0) {
    avisar(nm, ": se descartan filas de modelos no exportados: ",
           paste(sobran, collapse = ", "), ".")
    pruebas[[nm]] <- pruebas[[nm]][!pruebas[[nm]]$modelo %in% sobran, ,
                                   drop = FALSE]
  }
}


# ==============================================================================
# 7. METADATOS Y VERIFICACIONES
# ==============================================================================

# Versiones con que se ajustaron los modelos: se leen de la sesión que guardó
# cada script de la Etapa II. Si no está, se usa la sesión actual.
version_de_sesion <- function(ses, paquete) {
  v <- tryCatch({
    x <- ses$otherPkgs[[paquete]]
    if (is.null(x)) x <- ses$loadedOnly[[paquete]]
    as.character(x$Version)[1]
  }, error = function(e) NA_character_)
  if (length(v) != 1 || is.na(v) || !nzchar(v)) NA_character_ else v
}

version_instalada <- function(paquete) {
  tryCatch(as.character(utils::packageVersion(paquete)),
           error = function(e) NA_character_)
}

r_de_sesion <- function(ses) {
  v <- tryCatch(paste(ses$R.version$major, ses$R.version$minor, sep = "."),
                error = function(e) NA_character_)
  if (length(v) != 1 || is.na(v) || v == ".") NA_character_ else v
}

ses_fh  <- variantes$session
ses_rob <- robusto$session

versiones <- c(R         = r_de_sesion(ses_fh),
               emdi      = version_de_sesion(ses_fh, "emdi"),
               saeRobust = version_de_sesion(ses_rob, "saeRobust"))
versiones_origen <- "sesion de ajuste"

if (anyNA(versiones)) {
  actual <- c(R         = paste(R.version$major, R.version$minor, sep = "."),
              emdi      = version_instalada("emdi"),
              saeRobust = version_instalada("saeRobust"))
  versiones[is.na(versiones)] <- actual[is.na(versiones)]
  versiones_origen <- "sesion de ajuste o de exportacion"
}

meta <- list(
  fecha             = format(Sys.time(), "%Y-%m-%d %H:%M"),
  n_municipios      = n_dom,
  n_modelos         = sum(modelos$estimado),
  modelos_ausentes  = modelos$rotulo[!modelos$estimado],
  modelo_final      = MODELO_FINAL,
  z_ic              = Z_IC,
  k_robusto         = K_ROBUSTO,
  b_mse_logit       = B_LOGIT,
  b_boot_robusto    = B_ROBUSTO,
  semilla           = numero_o_na(variantes$semilla),
  versiones         = versiones,
  versiones_origen  = versiones_origen
)

cat("\n============================================================\n")
cat("VERIFICACIONES\n")
cat("============================================================\n")

n_ok <- sum(modelos$estimado)

cat("Modelos exportados:", n_ok, "de", nrow(modelos), "\n")
cat("Filas de estimaciones:", nrow(estimaciones), "(esperadas",
    n_ok * n_dom, ")\n")

n_qq   <- if (is.null(qq)) 0 else length(unique(qq$modelo))
n_coef <- if (is.null(coeficientes)) 0 else length(unique(coeficientes$modelo))
n_fh   <- sum(modelos$familia == "FH" & modelos$estimado)

cat("Modelos con coeficientes:", n_coef, "\n")
cat("Modelos con series para los QQ-plots:", n_qq, "(esperados", n_fh,
    ": los Fay-Herriot)\n")
cat("Modelos con estadisticos de normalidad:",
    if (is.null(normalidad)) 0 else nrow(normalidad), "\n")

if (n_coef < n_ok) avisar("Hay modelos sin tabla de coeficientes.")

stopifnot(
  "Filas de estimaciones inesperadas" = nrow(estimaciones) == n_ok * n_dom,
  "Hay estimaciones duplicadas" =
    !any(duplicated(estimaciones[, c("modelo", "cod_mun")])),
  "Hay identificadores de modelo duplicados" = !any(duplicated(modelos$modelo))
)

# Ninguna tabla debe traer columnas de la referencia DANE ni cantidades que
# permitan reconstruirla.
todas <- c(list(modelos = modelos, estimaciones = estimaciones,
                coeficientes = coeficientes, qq = qq,
                municipios = municipios, normalidad = normalidad,
                especificaciones = especificaciones), pruebas)

for (nm in names(todas)) {
  if (is.null(todas[[nm]])) next
  malas <- grep(PATRON_REFERENCIA, names(todas[[nm]]), ignore.case = TRUE,
                value = TRUE)
  malas <- setdiff(malas, "gq_p_psi")   # valor p de una prueba, no una serie
  if (length(malas) > 0) {
    stop("La tabla ", nm, " trae columnas que no deben exportarse: ",
         paste(malas, collapse = ", "))
  }
}

cat("Sin columnas de la referencia DANE en ninguna tabla.\n")

cat("\nModelos:\n")
print(modelos[, c("orden", "modelo", "escala", "mse_tipo", "mse_b", "etapa3",
                  "sigma2_u", "rho", "cv_medio", "cv_max", "estimado")],
      row.names = FALSE, digits = 4)

if (any(modelos$n_fuera_rango > 0, na.rm = TRUE)) {
  cat("\nModelos con estimaciones fuera de [0, 1]:\n")
  print(modelos[which(modelos$n_fuera_rango > 0),
                c("modelo", "n_fuera_rango")], row.names = FALSE)
}

if (any(modelos$n_ic_recortados > 0, na.rm = TRUE)) {
  cat("\nModelos con intervalos recortados a [0, 1]:\n")
  print(modelos[which(modelos$n_ic_recortados > 0),
                c("modelo", "n_ic_recortados")], row.names = FALSE)
}


# ==============================================================================
# 8. GUARDADO
# ==============================================================================

salida <- list(
  meta             = meta,
  municipios       = municipios,
  modelos          = modelos,
  estimaciones     = estimaciones,
  coeficientes     = coeficientes,
  especificaciones = especificaciones,
  qq               = qq,
  normalidad       = normalidad,
  pruebas          = pruebas
)

archivo <- file.path(ruta_out, "dashboard_datos.rds")
saveRDS(salida, archivo)

cat("\n============================================================\n")
cat("EXPORTACION COMPLETADA\n")
cat("Archivo: ", archivo, "\n", sep = "")
cat("Tamano:  ", round(file.size(archivo) / 1024, 1), " KB\n", sep = "")

if (length(avisos) > 0) {
  cat("\nAVISOS (", length(avisos), "):\n", sep = "")
  cat(paste0("  - ", avisos, collapse = "\n"), "\n", sep = "")
} else {
  cat("Sin avisos.\n")
}

cat("============================================================\n")
