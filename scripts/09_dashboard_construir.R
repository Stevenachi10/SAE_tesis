# ==============================================================================
# 09_dashboard_construir.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# CONSTRUCCIÓN DEL TABLERO INTERACTIVO
#
# Genera UN SOLO archivo HTML autocontenido (output/dashboard_pobreza.html) que
# se abre con doble clic en cualquier navegador, sin servidor, sin conexión a
# internet y sin paquetes de R. Todo el código del tablero (HTML, estilos y
# JavaScript) y todos los datos van dentro del mismo archivo.
#
# ------------------------------------------------------------------------------
# QUÉ INCLUYE EL TABLERO
# ------------------------------------------------------------------------------
#   Mapa y tabla   Mapa coroplético de los 84 municipios, por clases. Se puede
#                  mostrar la estimación, el CV, la amplitud del intervalo, la
#                  diferencia con otro modelo o el rango de la estimación entre
#                  todos los modelos. Tabla buscable y ordenable, sin descarga.
#   Modelo         Ficha del modelo, covariables de M1, M2 y M3, coeficientes y
#                  tabla comparativa de los 18 modelos, agrupada por familia o
#                  por estimador del MSE.
#   Municipio      Estimación e intervalo del municipio elegido en todos los
#                  modelos, y municipios más sensibles a la especificación.
#   Supuestos      QQ-plots y estadísticos de normalidad, solo de los nueve
#                  modelos Fay-Herriot (sin espaciales ni robustos).
#   Concordancia   Pruebas de sesgo, bondad de ajuste y solapamiento de
#                  intervalos frente a la referencia municipal DANE, con la
#                  explicación de cada prueba y de cada columna.
#
# El modelo final se marca con una estrella y los que no entraron a la
# evaluación de concordancia y precisión de la Etapa III con una cruz. No
# incluye la referencia municipal DANE; solo los resultados de las pruebas por
# modelo.
#
# ------------------------------------------------------------------------------
# QUÉ QUEDA DENTRO DEL HTML
# ------------------------------------------------------------------------------
# El tablero no tiene botón de descarga, pero los datos van incrustados en el
# archivo y quien abra el código de la página puede leerlos. Por eso este script
# publica solo las columnas de la lista PUBLICAR, que son las que el tablero
# muestra. Cualquier otra columna que traiga dashboard_datos.rds se descarta y
# se informa en la consola.
#
# ------------------------------------------------------------------------------
# ENTRADAS
# ------------------------------------------------------------------------------
#   output/dashboard_datos.rds           Lo genera 08_exportar_datos_dashboard.R
#   scripts/dashboard/plantilla.html     Código del tablero (HTML, CSS y JS)
#   scripts/dashboard/geometria_municipios.json
#                                        Contornos de los 84 municipios, límite
#                                        de los dos departamentos y ubicación de
#                                        las capitales (ver la nota siguiente).
#
# La geometría se generó una sola vez a partir de la cartografía del proyecto
# (data/pivoteadas/MGN2025_MPIO_GRAFICO, DANE, MGN 2025) con la herramienta
# mapshaper: simplificación de Visvalingam ponderada al 4 %, conservando la
# topología entre municipios vecinos, coordenadas con cuatro decimales y
# disolución por departamento para obtener su límite. Las capitales usan las
# mismas coordenadas aproximadas de 06_1_mapas_m1_log.R. Si el archivo no trae
# el límite departamental ni las capitales, el mapa se dibuja sin ellos.
#
# SALIDA
#   output/dashboard_pobreza.html
#
# El script solo usa R base: no requiere jsonlite ni sf. Si jsonlite está
# instalado, se usa únicamente para comprobar que los datos son un JSON válido.
# ==============================================================================


# ==============================================================================
# 0. CONFIGURACIÓN
# ==============================================================================

raiz <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()

RUTA_DATOS      <- file.path(raiz, "output", "dashboard_datos.rds")
RUTA_PLANTILLA  <- file.path(raiz, "scripts", "dashboard", "plantilla.html")
RUTA_GEOMETRIA  <- file.path(raiz, "scripts", "dashboard", "geometria_municipios.json")
RUTA_SALIDA     <- file.path(raiz, "output", "dashboard_pobreza.html")

# Nombres legibles de las covariables. Para añadir una covariable nueva basta con
# agregar una línea. Si falta alguna, el tablero muestra el nombre del objeto.
# Las tildes van como escapes Unicode (\u00f3 = ó) para que el texto salga bien
# aunque la sesión de R no esté en UTF-8.
ETIQUETAS_COVARIABLES <- c(
  alcantarillado_2018       = "Cobertura de alcantarillado (2018)",
  pct_rural                 = "Porcentaje de poblaci\u00f3n rural (2024)",
  partos_calificados_2020   = "Partos con atenci\u00f3n calificada (2020)",
  transferencias_pc         = "Transferencias corrientes per c\u00e1pita (2024)",
  saber11_lectura_2022      = "Puntaje Saber 11, lectura cr\u00edtica (2022)",
  tasa_violencia_intra_2019 = "Tasa de violencia intrafamiliar (2019)",
  pct_secundario_2023       = "Participaci\u00f3n del sector secundario en el valor agregado (2023)",
  tiene_empresas_formales   = "Presencia de empresas formales (2016)",
  dist_cali_km              = "Distancia geod\u00e9sica a Cali (2024)"
)

# Columnas que se publican dentro del HTML. Son las que el tablero muestra.
PUBLICAR <- list(
  municipios       = c("cod_mun", "municipio", "departamento"),
  modelos          = c("modelo", "rotulo", "grupo", "familia", "especificacion",
                       "n_covariables", "escala", "metodo", "k", "mse_tipo",
                       "mse_b", "sigma2_u", "rho", "convergencia", "es_final",
                       "cv_medio", "cv_mediana", "cv_max", "n_ic_recortados",
                       "n_fuera_rango", "etapa3", "orden"),
  estimaciones     = c("modelo", "cod_mun", "estimacion", "cv_pct", "ic_inf",
                       "ic_sup"),
  coeficientes     = c("modelo", "termino", "coef", "ee"),
  especificaciones = c("especificacion", "criterio", "orden", "variable"),
  qq               = c("modelo", "cod_mun", "res_z", "ef_z"),
  normalidad       = c("modelo", "ef_asimetria", "ef_curtosis", "ef_shapiro_W",
                       "ef_shapiro_p", "res_asimetria", "res_curtosis",
                       "res_shapiro_W", "res_shapiro_p")
)

PUBLICAR_PRUEBAS <- list(
  bias_brown = c("modelo", "b0", "b1", "ee_b1", "R2", "p_conjunto",
                 "p_conjunto_hc3", "gq_p_psi"),
  gof_brown  = c("modelo", "W", "c2", "gl", "p_value", "results"),
  cobertura  = c("modelo", "non_coverage", "domains", "non_overlap", "p_value",
                 "results")
)

PUBLICAR_META <- c("fecha", "n_municipios", "n_modelos", "modelos_ausentes",
                   "modelo_final", "z_ic", "k_robusto", "b_mse_logit",
                   "b_boot_robusto", "semilla", "versiones", "versiones_origen")

MARCADOR_DATOS <- "__DATOS_JSON__"
MARCADOR_GEO   <- "__GEO_JSON__"


# ==============================================================================
# 1. SERIALIZADOR JSON EN R BASE
# ------------------------------------------------------------------------------
# Convenciones:
#   data.frame            objeto con un arreglo por columna. Ocupa mucho menos
#                         que un objeto por fila; la plantilla lo pasa a filas.
#   lista con nombres     objeto
#   lista sin nombres     arreglo
#   vector con nombres    objeto
#   vector sin nombres    escalar si tiene longitud 1; arreglo en otro caso
#   NA, NaN, Inf          null
#
# El texto que sale es ASCII puro: las tildes y cualquier otro carácter fuera de
# ASCII se escriben como \uXXXX, igual que los caracteres de control y los
# signos <, > y &. Así el resultado no depende de la codificación de la sesión
# de R y ningún texto de los datos puede cerrar ni alterar la etiqueta <script>
# donde queda incrustado.
# ==============================================================================

# Escapa una cadena UTF-8 y la devuelve entre comillas.
.escapar_cadena <- function(s) {
  cp <- utf8ToInt(s)

  if (anyNA(cp)) {
    # Bytes que no forman UTF-8 válido: se sustituyen para no cortar el JSON.
    s  <- iconv(s, from = "UTF-8", to = "UTF-8", sub = "?")
    cp <- utf8ToInt(s)
  }

  if (length(cp) == 0) return("\"\"")

  simple <- cp >= 32L & cp < 127L & !(cp %in% c(34L, 92L, 60L, 62L, 38L))
  ch     <- intToUtf8(cp, multiple = TRUE)

  if (!all(simple)) {
    ch[!simple] <- vapply(cp[!simple], function(u) {
      if (u == 34L) return("\\\"")
      if (u == 92L) return("\\\\")
      if (u == 10L) return("\\n")
      if (u == 13L) return("\\r")
      if (u == 9L)  return("\\t")
      if (u < 65536L) return(sprintf("\\u%04x", u))
      v <- u - 65536L
      sprintf("\\u%04x\\u%04x", 55296L + v %/% 1024L, 56320L + v %% 1024L)
    }, "")
  }

  paste0("\"", paste(ch, collapse = ""), "\"")
}

json_cadena <- function(x) {
  x   <- enc2utf8(as.character(x))
  out <- rep("null", length(x))
  ok  <- !is.na(x)

  if (!any(ok)) return(out)

  # Camino rápido: cadenas ASCII sin caracteres que escapar.
  directa <- ok & !grepl("[^ -~]|[\"\\\\<>&]", x, useBytes = TRUE)
  out[directa] <- paste0("\"", x[directa], "\"")

  resto <- which(ok & !directa)
  if (length(resto) > 0) {
    out[resto] <- vapply(x[resto], .escapar_cadena, "", USE.NAMES = FALSE)
  }

  out
}

# Un elemento JSON por cada valor del vector atómico.
json_vector <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (inherits(x, c("Date", "POSIXt"))) x <- as.character(x)

  if (is.logical(x)) {
    return(ifelse(is.na(x), "null", ifelse(x, "true", "false")))
  }

  if (is.integer(x)) {
    return(ifelse(is.na(x), "null", as.character(x)))
  }

  if (is.numeric(x)) {
    out <- rep("null", length(x))
    ok  <- is.finite(x)
    if (any(ok)) out[ok] <- trimws(formatC(x[ok], digits = 10, format = "g"))
    return(out)
  }

  json_cadena(x)
}

a_json <- function(x) {

  if (is.null(x)) return("null")

  if (is.data.frame(x)) {
    if (ncol(x) == 0) return("{}")
    cols <- vapply(x, function(v) {
      paste0("[", paste0(json_vector(v), collapse = ","), "]")
    }, "")
    return(paste0("{", paste0(json_cadena(names(x)), ":", cols, collapse = ","),
                  "}"))
  }

  if (is.list(x)) {
    if (length(x) == 0) return("[]")
    if (!is.null(names(x)) && all(nzchar(names(x)))) {
      return(paste0("{", paste0(json_cadena(names(x)), ":",
                                vapply(x, a_json, ""), collapse = ","), "}"))
    }
    return(paste0("[", paste0(vapply(x, a_json, ""), collapse = ","), "]"))
  }

  if (is.atomic(x)) {
    if (length(x) == 0) return("[]")
    el <- json_vector(x)
    if (!is.null(names(x)) && all(nzchar(names(x)))) {
      return(paste0("{", paste0(json_cadena(names(x)), ":", el, collapse = ","),
                    "}"))
    }
    if (length(x) == 1) return(el)
    return(paste0("[", paste0(el, collapse = ","), "]"))
  }

  stop("Tipo no serializable: ", class(x)[1])
}


# ==============================================================================
# 2. LECTURA
# ------------------------------------------------------------------------------
# La plantilla y la geometría se leen como bytes y el HTML se arma sobre bytes,
# de modo que la codificación del sistema (UTF-8 o no) no altera los acentos.
# ==============================================================================

leer_bytes <- function(ruta) {
  if (!file.exists(ruta)) stop("No existe ", ruta)
  b <- readBin(ruta, "raw", file.size(ruta))

  # Marca de orden de bytes (BOM) que añaden algunos editores de Windows: dentro
  # del HTML rompería el JSON incrustado.
  if (length(b) >= 3 && identical(as.integer(b[1:3]), c(239L, 187L, 191L))) {
    b <- b[-(1:3)]
  }

  b
}

# Número de veces que aparece un marcador en un texto (por bytes).
contar <- function(texto, marcador) {
  m <- gregexpr(marcador, texto, fixed = TRUE, useBytes = TRUE)[[1]]
  if (m[1] < 0) 0L else length(m)
}

if (!file.exists(RUTA_DATOS)) {
  stop("No existe ", RUTA_DATOS, ". Ejecute antes 08_exportar_datos_dashboard.R.")
}

datos         <- readRDS(RUTA_DATOS)
plantilla_raw <- leer_bytes(RUTA_PLANTILLA)
geometria_raw <- leer_bytes(RUTA_GEOMETRIA)

# Solo para comprobaciones y para localizar marcadores (texto ASCII).
plantilla_txt <- rawToChar(plantilla_raw)
geometria_txt <- trimws(rawToChar(geometria_raw))

componentes_esperados <- c("meta", "municipios", "modelos", "estimaciones",
                           "coeficientes", "especificaciones")
faltan <- setdiff(componentes_esperados, names(datos))

if (length(faltan) > 0) {
  stop("dashboard_datos.rds no trae: ", paste(faltan, collapse = ", "))
}

if ("diagnosticos_dominio" %in% names(datos)) {
  stop("dashboard_datos.rds viene de una version anterior de ",
       "08_exportar_datos_dashboard.R. Vuelva a ejecutar ese script.")
}

if (contar(plantilla_txt, MARCADOR_DATOS) != 1L ||
    contar(plantilla_txt, MARCADOR_GEO) != 1L) {
  stop("La plantilla debe contener exactamente una vez cada marcador (",
       MARCADOR_DATOS, " y ", MARCADOR_GEO, "). Revise que sea la plantilla ",
       "que acompana a este script.")
}

if (!grepl("id=\"datos-json\"", plantilla_txt, fixed = TRUE, useBytes = TRUE)) {
  stop("La plantilla es de una version anterior del tablero. Use la que ",
       "acompana a este script.")
}

stopifnot(
  "La geometria no parece un JSON" = startsWith(geometria_txt, "{"),
  "La geometria no trae la lista de municipios" =
    grepl("\"municipios\"", geometria_txt, fixed = TRUE, useBytes = TRUE)
)

if (!grepl("\"departamentos\"", geometria_txt, fixed = TRUE, useBytes = TRUE)) {
  cat("AVISO: la geometria es de una version anterior y no trae el limite ",
      "departamental ni las capitales. El mapa se dibuja sin ellos. Use el ",
      "archivo geometria_municipios.json que acompana a este script.\n", sep = "")
}

# Los códigos de la geometría deben coincidir con los de los datos.
cod_geo <- regmatches(geometria_txt,
                      gregexpr("\"cod_mun\"[[:space:]]*:[[:space:]]*\"[0-9]+\"",
                               geometria_txt, useBytes = TRUE))[[1]]
cod_geo <- gsub("[^0-9]", "", cod_geo)

if (any(duplicated(cod_geo))) {
  stop("La geometria trae municipios repetidos: ",
       paste(unique(cod_geo[duplicated(cod_geo)]), collapse = ", "))
}

if (!setequal(cod_geo, as.character(datos$municipios$cod_mun))) {
  solo_geo   <- setdiff(cod_geo, as.character(datos$municipios$cod_mun))
  solo_datos <- setdiff(as.character(datos$municipios$cod_mun), cod_geo)
  stop("Los codigos municipales de la geometria no coinciden con los de los ",
       "datos. Solo en la geometria: ",
       if (length(solo_geo)) paste(solo_geo, collapse = ", ") else "ninguno",
       ". Solo en los datos: ",
       if (length(solo_datos)) paste(solo_datos, collapse = ", ") else "ninguno",
       ".")
}


# ==============================================================================
# 3. SELECCIÓN DE LO QUE SE PUBLICA
# ==============================================================================

descartadas <- character(0)

# Se queda con las columnas permitidas de una tabla y anota las demás.
solo_publicables <- function(tab, permitidas, nombre) {
  if (is.null(tab)) return(NULL)
  sobran <- setdiff(names(tab), permitidas)
  if (length(sobran) > 0) {
    descartadas <<- c(descartadas, paste0(nombre, ": ",
                                          paste(sobran, collapse = ", ")))
  }
  tab[, intersect(permitidas, names(tab)), drop = FALSE]
}

# Solo los modelos que se exportaron completos.
estimado <- if ("estimado" %in% names(datos$modelos)) {
  !is.na(datos$modelos$estimado) & datos$modelos$estimado
} else rep(TRUE, nrow(datos$modelos))

ausentes <- datos$modelos$rotulo[!estimado]
modelos_ok <- datos$modelos$modelo[estimado]

publicar <- list()

publicar$meta <- datos$meta[intersect(PUBLICAR_META, names(datos$meta))]
publicar$meta$modelos_ausentes <- as.character(ausentes)
publicar$meta$construido <- format(Sys.time(), "%Y-%m-%d %H:%M")

datos$modelos$estimado <- NULL
datos$modelos <- datos$modelos[estimado, , drop = FALSE]

for (nm in names(PUBLICAR)) {
  tab <- solo_publicables(datos[[nm]], PUBLICAR[[nm]], nm)
  if (!is.null(tab) && "modelo" %in% names(tab) && nm != "modelos") {
    tab <- tab[tab$modelo %in% modelos_ok, , drop = FALSE]
  }
  if (!is.null(tab)) rownames(tab) <- NULL
  publicar[[nm]] <- tab
}

publicar$pruebas <- list()

for (nm in names(PUBLICAR_PRUEBAS)) {
  tab <- solo_publicables(datos$pruebas[[nm]], PUBLICAR_PRUEBAS[[nm]],
                          paste0("pruebas$", nm))
  if (!is.null(tab) && nrow(tab) > 0) {
    tab <- tab[tab$modelo %in% modelos_ok, , drop = FALSE]
    rownames(tab) <- NULL
    publicar$pruebas[[nm]] <- tab
  }
}

otras <- setdiff(names(datos), c(names(PUBLICAR), "meta", "pruebas"))
if (length(otras) > 0) {
  descartadas <- c(descartadas, paste0("tablas completas: ",
                                       paste(otras, collapse = ", ")))
}

publicar$etiquetas <- as.list(ETIQUETAS_COVARIABLES)

cat("Modelos:", length(modelos_ok), "| municipios:",
    nrow(publicar$municipios), "| poligonos:", length(cod_geo), "\n")

if (length(ausentes) > 0) {
  cat("AVISO: modelos que no se exportaron y no apareceran en el tablero:\n  ",
      paste(ausentes, collapse = "\n  "), "\n", sep = "")
}

if (length(descartadas) > 0) {
  cat("Columnas que no se publican (el tablero no las usa):\n  ",
      paste(descartadas, collapse = "\n  "), "\n", sep = "")
}

# Covariables sin etiqueta: se avisan para completarlas.
usadas <- unique(c(publicar$especificaciones$variable,
                   setdiff(publicar$coeficientes$termino, "(Intercept)")))
sin_etiqueta <- setdiff(usadas, names(ETIQUETAS_COVARIABLES))

if (length(sin_etiqueta) > 0) {
  cat("AVISO: covariables sin etiqueta (se mostrara el nombre del objeto):\n  ",
      paste(sin_etiqueta, collapse = ", "), "\n", sep = "")
}


# ==============================================================================
# 4. ENSAMBLE
# ==============================================================================

json_datos <- a_json(publicar)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  if (!isTRUE(jsonlite::validate(json_datos))) {
    stop("Los datos serializados no forman un JSON valido.")
  }
}

json_datos_raw <- charToRaw(json_datos)

# En la geometría, un signo < solo puede estar dentro de un texto: se escribe
# como \u003c para que no interfiera con la etiqueta <script>.
if (grepl("<", geometria_txt, fixed = TRUE, useBytes = TRUE)) {
  geometria_txt <- gsub("<", "\\u003c", geometria_txt, fixed = TRUE,
                        useBytes = TRUE)
}
geometria_raw <- charToRaw(geometria_txt)

# Se parte la plantilla en los marcadores en lugar de usar sub(): asi ningun
# caracter de los datos se interpreta como instruccion de reemplazo.
insertar <- function(bytes, marcador, contenido) {
  texto <- rawToChar(bytes)
  pos   <- as.integer(regexpr(marcador, texto, fixed = TRUE, useBytes = TRUE))
  if (pos < 0) stop("Marcador ausente: ", marcador)
  fin   <- pos + nchar(marcador, type = "bytes")
  c(bytes[seq_len(pos - 1)], contenido,
    if (fin <= length(bytes)) bytes[fin:length(bytes)] else raw(0))
}

html_raw <- insertar(plantilla_raw, MARCADOR_DATOS, json_datos_raw)
html_raw <- insertar(html_raw,      MARCADOR_GEO,   geometria_raw)


# ==============================================================================
# 5. GUARDADO
# ==============================================================================

dir.create(dirname(RUTA_SALIDA), showWarnings = FALSE, recursive = TRUE)

con <- file(RUTA_SALIDA, "wb")
writeBin(html_raw, con)
close(con)

cat("\n============================================================\n")
cat("TABLERO CONSTRUIDO\n")
cat("Archivo:   ", RUTA_SALIDA, "\n", sep = "")
cat("Tamano:    ", round(file.size(RUTA_SALIDA) / 1024), " KB (datos ",
    round(length(json_datos_raw) / 1024), " KB, geometria ",
    round(length(geometria_raw) / 1024), " KB)\n", sep = "")
cat("Abrir con doble clic en el navegador.\n")
cat("============================================================\n")
