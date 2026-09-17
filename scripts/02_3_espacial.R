# ==============================================================================
# 02_3_espacial.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.3 — Dependencia espacial y modelo Fay-Herriot espacial
#
# Construye la matriz de proximidad por contigüidad de primer orden, contrasta
# la autocorrelación del estimador directo y de las covariables, y ajusta el
# modelo espacial sobre las especificaciones seleccionadas.
#
# ------------------------------------------------------------------------------
# DEFINICION DE VECINDAD
# ------------------------------------------------------------------------------
# Se emplea una única definición: contigüidad tipo queen de primer orden,
# estandarizada por filas. Es la definición convencional en la literatura y no
# deja dominios aislados. No se evalúan definiciones alternativas.
#
# ------------------------------------------------------------------------------
# CRITERIO DE DECISION
# ------------------------------------------------------------------------------
# La pertinencia del modelo espacial se establece a partir de dos elementos:
#
#   1. DEPENDENCIA ESPACIAL EN LOS DATOS. Contrastes de Moran y Geary sobre el
#      estimador directo y sobre las covariables, siguiendo el procedimiento
#      del estudio de caso de Pratesi y Salvati (2009). Ambas son variables
#      observadas, de modo que la distribución de referencia de los contrastes
#      es la que corresponde.
#
#   2. MAGNITUD DEL PARAMETRO DE AUTOCORRELACION. La literatura sugiere
#      considerar el modelo espacial cuando rho supera 0,5 o cuando no hay
#      covariables geográficas que recojan la estructura territorial
#      (Bertarelli et al., 2021, citado en Harmening et al., 2023), en línea
#      con Pratesi y Salvati (2009) y Pfeffermann (2002). El script no
#      automatiza ese umbral: reporta el parámetro y la valoración se hace en
#      el documento. La especificación M3s, que retira de M3 la covariable de
#      distancia, permite establecer si el valor de rho depende de ella.
#
# El coeficiente de variación frente al modelo no espacial se reporta con
# carácter descriptivo; la comparación de precisión corresponde a la Etapa III.
#
# ------------------------------------------------------------------------------
# SOBRE LOS EFECTOS ALEATORIOS
# ------------------------------------------------------------------------------
# No se contrasta la autocorrelación de los efectos aleatorios estimados. El
# efecto aleatorio es una esperanza condicional, u_i = E(u_i | y), no una
# realización: su varianza es gamma_i * sigma2_u, que depende del factor de
# contracción y no es homogénea entre dominios. La distribución de referencia
# implementada en spdep corresponde a una variable observada, y la corrección
# de Cliff y Ord disponible en lm.morantest está derivada para residuos de
# mínimos cuadrados ordinarios, sin cubrir el caso de varianzas de muestreo
# conocidas y heterogéneas. No se identificó una derivación establecida para
# este caso, y la limitación se declara de forma explícita en el documento.
#
# ------------------------------------------------------------------------------
# SOBRE LA RAZON DE VEROSIMILITUDES
# ------------------------------------------------------------------------------
# No se emplea. La rutina que estima los componentes de varianza en el caso
# espacial (SREML) resuelve las ecuaciones de score de la verosimilitud
# RESTRINGIDA, que incluye el término -0,5 log|X'V^{-1}X|, mientras que la
# cantidad reportada como loglike corresponde a la verosimilitud marginal, que
# omite dicho término. La diferencia entre las verosimilitudes de ambos modelos
# no constituye, por tanto, un estadístico de razón de verosimilitudes.
#
# Adicionalmente, Pratesi y Salvati (2009) señalan que la función de
# verosimilitud presenta máximos locales además del máximo global, y por esa
# razón encadenan Nelder-Mead con el algoritmo de scoring. SREML inicia la
# búsqueda en rho = 0,5 con un punto de partida fijo.
#
# ------------------------------------------------------------------------------
# ENTRADAS
# ------------------------------------------------------------------------------
#   data/pivoteadas/MGN2025_MPIO_GRAFICO/MGN_ADM_MPIO_GRAFICO.shp
#   output/matriz_sae_transformada_v2.rds
#   output/seleccion_stepwise_completo_both.rds
#
# ------------------------------------------------------------------------------
# SALIDAS (output/)
# ------------------------------------------------------------------------------
#   matriz_proximidad.rds          Matriz W estandarizada por filas.
#   esp_autocorrelacion_queen.csv  Moran y Geary del directo y las covariables.
#   esp_ajuste_queen.csv           Ajuste espacial por especificación.
#   espacial.rds                   Objetos del ajuste para la Etapa III.
#   fig_esp_vecindad.png           Grafo de vecindad sobre el mapa.
# ==============================================================================

library(sf)
library(spdep)
library(emdi)
library(car)
library(dplyr)
library(ggplot2)
library(here)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")
ruta_car <- here("data", "pivoteadas", "MGN2025_MPIO_GRAFICO")

matriz  <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))
corrida <- readRDS(file.path(ruta_out, "seleccion_stepwise_completo_both.rds"))


# ==============================================================================
# CONFIGURACION
# ==============================================================================

SEMILLA <- 2906

# Replicas de las pruebas por permutaciones.
B_PERM <- 9999

# Diagnostico numerico: proximidad de rho al limite del espacio de parametros.
# No constituye criterio de admisibilidad.
RHO_LIMITE <- 0.99

# Tolerancia y maximo de iteraciones del algoritmo de estimacion de los
# componentes de varianza. Los valores por defecto de fh() son 1e-4 y 100.
TOL   <- 1e-8
MAXIT <- 1000

# Origen Nacional. La cartografia viene en coordenadas geograficas; los
# centroides requieren un sistema proyectado.
CRS_PROYECTADO <- 9377

tema_tesis <- theme_minimal(base_size = 13) +
  theme(plot.title = element_blank(), panel.grid.minor = element_blank())

COL_PUNTO <- "#b5482e"; COL_LINEA <- "#2c6e6b"; COL_BANDA <- "#c9c2b3"

escalar <- function(x) {
  v <- tryCatch(as.numeric(x), error = function(e) NA_real_)
  if (length(v) != 1 || !is.finite(v)) NA_real_ else v
}

formula_de <- function(vars) {
  as.formula(paste("pobreza_monetaria ~", paste(vars, collapse = " + ")))
}


# ==============================================================================
# ESPECIFICACIONES
# ------------------------------------------------------------------------------
# Las tres especificaciones seleccionadas, mas una cuarta que retira la unica
# covariable de distancia de M3. Esta ultima no constituye una especificacion
# candidata: permite establecer si el valor del parametro de autocorrelacion
# obedece a que dicha covariable recoge la estructura territorial.
# ==============================================================================

tabla_sel <- corrida$tabla
etiquetas <- c(KICb2 = "M1", KICc = "M2", KIC = "M3")

ESPECIFICACIONES <- list()

for (i in seq_len(nrow(tabla_sel))) {
  et <- unname(etiquetas[tabla_sel$criterio[i]])
  ESPECIFICACIONES[[et]] <- strsplit(tabla_sel$covariables[i],
                                     " + ", fixed = TRUE)[[1]]
}

ESPECIFICACIONES <- ESPECIFICACIONES[order(names(ESPECIFICACIONES))]

DISTANCIAS <- unique(grep("^dist_", unlist(ESPECIFICACIONES), value = TRUE))

if (length(DISTANCIAS) > 0 && "M3" %in% names(ESPECIFICACIONES)) {
  ESPECIFICACIONES[["M3s"]] <- setdiff(ESPECIFICACIONES[["M3"]], DISTANCIAS)
}

cat("\n== ESPECIFICACIONES ==\n")
for (nm in names(ESPECIFICACIONES)) {
  cat(nm, " (", length(ESPECIFICACIONES[[nm]]), "): ",
      paste(ESPECIFICACIONES[[nm]], collapse = ", "), "\n", sep = "")
}

vars_todas <- unique(unlist(ESPECIFICACIONES))

datos <- matriz %>%
  select(cod_mun, Municipio, pobreza_monetaria, varianza_pobreza,
         all_of(vars_todas)) %>%
  as.data.frame()

n_dom <- nrow(datos)
PSI   <- datos$varianza_pobreza

stopifnot(
  "Hay varianzas de muestreo no positivas" = all(PSI > 0),
  "Hay valores faltantes en las covariables" =
    !any(is.na(datos[, vars_todas, drop = FALSE])))


# ==============================================================================
# PARTE 1 — CARTOGRAFIA Y MATRIZ DE PROXIMIDAD
# ==============================================================================

poligonos <- st_read(file.path(ruta_car, "MGN_ADM_MPIO_GRAFICO.shp"),
                     quiet = TRUE)

cat("\n== CARTOGRAFIA MUNICIPAL ==\n")
cat("Poligonos en el archivo:", nrow(poligonos), "\n")

poligonos <- poligonos %>%
  filter(dpto_ccdgo %in% c("19", "76")) %>%
  mutate(cod_mun = as.character(mpio_cdpmp)) %>%
  select(cod_mun, nombre_mgn = mpio_cnmbr, geometry)

faltan <- setdiff(datos$cod_mun, poligonos$cod_mun)

if (length(faltan) > 0) {
  cat("\nMunicipios ausentes en la cartografia:\n")
  print(datos %>% filter(cod_mun %in% faltan) %>% select(cod_mun, Municipio),
        row.names = FALSE)
}

stopifnot("Faltan municipios del estudio en la cartografia" =
            length(faltan) == 0)

poligonos <- poligonos %>%
  filter(cod_mun %in% datos$cod_mun) %>%
  arrange(match(cod_mun, datos$cod_mun))

# La contiguidad es invariante a la proyeccion, pero los centroides no lo son.
cat("CRS de origen:", st_crs(poligonos)$input, "\n")
poligonos <- st_transform(poligonos, CRS_PROYECTADO)

# El orden de las filas de la matriz debe coincidir con el de los dominios: una
# discordancia produce ajustes sin sentido y sin ningun aviso.
stopifnot("El orden de la cartografia no coincide con el de los datos" =
            identical(poligonos$cod_mun, datos$cod_mun))

cat("Cartografia alineada con los", nrow(poligonos), "dominios.\n")

nb_queen <- poly2nb(poligonos, queen = TRUE, row.names = poligonos$cod_mun)
W        <- nb2mat(nb_queen, style = "W", zero.policy = FALSE)
lw_queen <- nb2listw(nb_queen, style = "W", zero.policy = FALSE)

cat("\n== ESTRUCTURA DE VECINDAD, CONTIGUIDAD DE PRIMER ORDEN ==\n")
cat("Vecinos por municipio: minimo", min(card(nb_queen)),
    "| media", round(mean(card(nb_queen)), 2),
    "| maximo", max(card(nb_queen)), "\n")
cat("Enlaces:", sum(card(nb_queen)) / 2,
    "| Dominios sin vecinos:", sum(card(nb_queen) == 0), "\n")

stopifnot(
  "Hay dominios sin vecinos" = all(card(nb_queen) > 0),
  "Dimensiones incorrectas de la matriz" =
    all(dim(W) == c(n_dom, n_dom)),
  "Las filas de la matriz no suman la unidad" =
    all(abs(rowSums(W) - 1) < 1e-10))

resumen_vec <- data.frame(
  cod_mun = poligonos$cod_mun, Municipio = datos$Municipio,
  n_vecinos = card(nb_queen), stringsAsFactors = FALSE) %>%
  arrange(desc(n_vecinos))

cat("\nMunicipios con mayor y menor numero de vecinos:\n")
print(rbind(head(resumen_vec, 5), tail(resumen_vec, 5)), row.names = FALSE)

saveRDS(W, file.path(ruta_out, "matriz_proximidad.rds"))

coords  <- st_coordinates(st_centroid(st_geometry(poligonos)))
enlaces <- nb2lines(nb_queen, coords = coords, as_sf = TRUE)
st_crs(enlaces) <- st_crs(poligonos)

ggsave(file.path(ruta_out, "fig_esp_vecindad.png"),
       ggplot() +
         geom_sf(data = poligonos, fill = "white", color = COL_BANDA,
                 linewidth = 0.3) +
         geom_sf(data = enlaces, color = COL_LINEA, linewidth = 0.3,
                 alpha = 0.6) +
         geom_point(data = as.data.frame(coords), aes(X, Y),
                    color = COL_PUNTO, size = 0.9) +
         labs(x = NULL, y = NULL) + tema_tesis +
         theme(axis.text = element_blank()),
       width = 8, height = 9, dpi = 400)


# ==============================================================================
# PARTE 2 — FUNCIONES
# ==============================================================================

# ------------------------------------------------------------------------------
# Contrastes de autocorrelacion
# ------------------------------------------------------------------------------
# Los dos estadisticos que emplea spatialcor.tests() del paquete: la I de Moran
# y la C de Geary. La I toma valores proximos a cero bajo ausencia de
# autocorrelacion, con esperanza -1/(m-1); la C toma valores proximos a la
# unidad, y valores inferiores indican autocorrelacion positiva. spdep invierte
# el signo del estadistico de Geary, de modo que alternative = "greater"
# corresponde a autocorrelacion positiva en ambos.
#
# Para la I se acompañan las versiones por permutaciones en cola superior y
# bilateral, que no descansan en la aproximacion normal.
pruebas_de <- function(x, lw) {
  
  vacio <- data.frame(
    I = NA_real_, I_esperanza = NA_real_, I_p_asint = NA_real_,
    I_p_perm = NA_real_, I_p_perm_bil = NA_real_,
    C = NA_real_, C_p_asint = NA_real_, stringsAsFactors = FALSE)
  
  if (is.null(x) || all(is.na(x))) return(vacio)
  
  mt <- tryCatch(moran.test(x, lw, randomisation = TRUE,
                            alternative = "greater", zero.policy = FALSE),
                 error = function(e) NULL)
  
  set.seed(SEMILLA)
  mc <- tryCatch(moran.mc(x, lw, nsim = B_PERM, alternative = "greater",
                          zero.policy = FALSE),
                 error = function(e) NULL)
  
  set.seed(SEMILLA)
  mb <- tryCatch(moran.mc(x, lw, nsim = B_PERM, alternative = "two.sided",
                          zero.policy = FALSE),
                 error = function(e) NULL)
  
  gt <- tryCatch(geary.test(x, lw, randomisation = TRUE,
                            alternative = "greater", zero.policy = FALSE),
                 error = function(e) NULL)
  
  data.frame(
    I            = if (is.null(mt)) NA_real_ else escalar(mt$estimate[1]),
    I_esperanza  = if (is.null(mt)) NA_real_ else escalar(mt$estimate[2]),
    I_p_asint    = if (is.null(mt)) NA_real_ else escalar(mt$p.value),
    I_p_perm     = if (is.null(mc)) NA_real_ else escalar(mc$p.value),
    I_p_perm_bil = if (is.null(mb)) NA_real_ else escalar(mb$p.value),
    C            = if (is.null(gt)) NA_real_ else escalar(gt$estimate[1]),
    C_p_asint    = if (is.null(gt)) NA_real_ else escalar(gt$p.value),
    stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------------------
# Ajuste no espacial
# ------------------------------------------------------------------------------
# Los avisos se registran en lugar de descartarse: con tol = 1e-8 y maxit = 1000
# una no convergencia por tope de iteraciones no produce error y afectaria a
# todas las cantidades derivadas.
ajustar_base <- function(vars) {
  
  avisos <- character(0)
  
  aj <- withCallingHandlers(
    tryCatch(fh(formula_de(vars), vardir = "varianza_pobreza",
                combined_data = datos, domains = "cod_mun",
                method = "reml", MSE = TRUE, mse_type = "analytical",
                tol = TOL, maxit = MAXIT, B = c(0, 0), seed = SEMILLA),
             error = function(e) NULL),
    warning = function(w) {
      avisos <<- c(avisos, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  
  list(ajuste = aj, avisos = unique(avisos))
}

# ------------------------------------------------------------------------------
# Alineacion de dominios
# ------------------------------------------------------------------------------
# W, gamma y las estimaciones se combinan por posicion. Una reordenacion
# silenciosa de la salida de fh() invalidaria el ajuste sin producir error.
alineado <- function(aj) {
  if (is.null(aj)) return(NA)
  dom <- tryCatch(as.character(aj$model$gamma$Domain), error = function(e) NULL)
  !is.null(dom) && identical(dom, as.character(datos$cod_mun))
}

# ------------------------------------------------------------------------------
# Componentes de varianza del ajuste espacial
# ------------------------------------------------------------------------------
# Se indexa por nombre. Si los componentes no aparecen la funcion falla, en
# lugar de suponer un orden: un intercambio de rho y sigma2_u produciria una
# matriz de encogimiento sin sentido y un gamma medio plausible.
extraer_componentes <- function(r) {
  
  v  <- r$model$variance
  nm <- names(v)
  
  if (is.null(nm) || !all(c("correlation", "variance") %in% nm)) {
    stop("model$variance no trae 'correlation' y 'variance'. ",
         "Revisar la version de emdi antes de interpretar el ajuste espacial.")
  }
  
  s2u  <- as.numeric(v[["variance"]])
  rho  <- as.numeric(v[["correlation"]])
  conv <- if ("convergence" %in% nm) as.logical(v[["convergence"]]) else NA
  
  list(s2u = s2u, rho = rho, convergencia = conv,
       coherente = is.finite(s2u) && is.finite(rho) && s2u > 0 && abs(rho) < 1)
}

# ------------------------------------------------------------------------------
# Medidas de precision
# ------------------------------------------------------------------------------
# Se conservan por separado el error cuadratico medio, el error estandar y el
# coeficiente de variacion, de modo que una anomalia pueda atribuirse al
# componente que la origina.
#
# El coeficiente de variacion no toma el valor absoluto del denominador: una
# estimacion negativa constituiria una prediccion fuera del rango admisible de
# una proporcion y debe detectarse, no encubrirse.
precision_de <- function(r) {
  
  est <- tryCatch(as.numeric(r$ind$FH), error = function(e) rep(NA_real_, n_dom))
  mse <- tryCatch(as.numeric(r$MSE$FH), error = function(e) rep(NA_real_, n_dom))
  
  se <- sqrt(pmax(mse, 0))
  
  cv <- ifelse(is.finite(est) & est > 0 & is.finite(se),
               se / est * 100, NA_real_)
  
  list(estimacion = est, mse = mse, se = se, cv = cv,
       n_fuera_rango = sum(est < 0 | est > 1, na.rm = TRUE),
       n_cv_na = sum(!is.finite(cv)))
}

# ------------------------------------------------------------------------------
# Comparacion de precision sobre dominios comunes
# ------------------------------------------------------------------------------
# El CV es NA donde la estimacion no es positiva. Si ese conjunto difiere entre
# el modelo espacial y el base, promediar por separado compararia medias sobre
# poblaciones distintas.
comparar_cv <- function(cv_esp, cv_base) {
  
  ok <- is.finite(cv_esp) & is.finite(cv_base)
  
  list(cv        = if (any(ok)) mean(cv_esp[ok])  else NA_real_,
       cv_base   = if (any(ok)) mean(cv_base[ok]) else NA_real_,
       dif       = if (any(ok)) mean(cv_esp[ok]) - mean(cv_base[ok]) else NA_real_,
       cv_max    = if (any(ok)) max(cv_esp[ok])   else NA_real_,
       n_comunes = sum(ok))
}

# ------------------------------------------------------------------------------
# Ajuste espacial
# ------------------------------------------------------------------------------
ajustar_esp <- function(vars, W_k) {
  
  limite <- FALSE
  
  r <- withCallingHandlers(
    tryCatch(fh(formula_de(vars), vardir = "varianza_pobreza",
                combined_data = datos, domains = "cod_mun",
                method = "reml", correlation = "spatial", corMatrix = W_k,
                MSE = TRUE, mse_type = "analytical",
                tol = TOL, maxit = MAXIT, B = c(0, 0), seed = SEMILLA),
             error = function(e) NULL),
    warning = function(w) {
      if (grepl("interval limit|singular", conditionMessage(w),
                ignore.case = TRUE)) limite <<- TRUE
      invokeRestart("muffleWarning")
    })
  
  if (is.null(r)) return(NULL)
  
  comp <- extraer_componentes(r)
  prec <- precision_de(r)
  
  # Diagonal de la matriz de encogimiento bajo estructura autorregresiva.
  # La covarianza de los efectos aleatorios es sigma2_u (A'A)^{-1} con
  # A = I - rho W. Se resuelve el sistema una sola vez sobre A'A en lugar de
  # invertir A y su traspuesta por separado.
  gam <- tryCatch({
    if (!comp$coherente) stop("componentes no coherentes")
    A     <- diag(n_dom) - comp$rho * W_k
    Omega <- comp$s2u * solve(t(A) %*% A)
    diag(Omega %*% solve(Omega + diag(PSI)))
  }, error = function(e) rep(NA_real_, n_dom))
  
  list(ajuste = r, limite = limite, alineado = alineado(r),
       s2u = comp$s2u, rho = comp$rho, convergencia = comp$convergencia,
       coherente = comp$coherente, gamma = gam,
       estimacion = prec$estimacion, mse = prec$mse, se = prec$se,
       cv = prec$cv, n_fuera_rango = prec$n_fuera_rango,
       n_cv_na = prec$n_cv_na)
}

vif_de <- function(vars) {
  if (length(vars) < 2) return(1)
  tryCatch(escalar(max(car::vif(lm(formula_de(vars), data = datos)))),
           error = function(e) NA_real_)
}

# ------------------------------------------------------------------------------
# Fila de resultados del ajuste espacial
# ------------------------------------------------------------------------------
fila_ajuste <- function(k, nm, res, cv_base_vec, segundos = NA_real_) {
  
  base <- data.frame(criterio_vecindad = k, especificacion = nm,
                     n_vars = length(ESPECIFICACIONES[[nm]]),
                     stringsAsFactors = FALSE)
  
  if (is.null(res)) {
    return(cbind(base, data.frame(
      converge = FALSE, rho = NA_real_, sigma2_u = NA_real_,
      gamma_medio = NA_real_, cv = NA_real_, cv_max = NA_real_,
      cv_base = NA_real_, dif_cv = NA_real_, n_comunes = NA_integer_,
      n_cv_na = NA_integer_, n_fuera_rango = NA_integer_,
      vif_max = vif_de(ESPECIFICACIONES[[nm]]), var_limite = NA,
      convergencia = NA, coherente = NA, alineado = NA, rho_al_limite = NA,
      segundos = round(segundos, 1), stringsAsFactors = FALSE)))
  }
  
  cc <- comparar_cv(res$cv, cv_base_vec)
  
  cbind(base, data.frame(
    converge       = TRUE,
    rho            = round(res$rho, 4),
    sigma2_u       = signif(res$s2u, 6),
    gamma_medio    = round(mean(res$gamma, na.rm = TRUE), 4),
    cv             = round(cc$cv, 3),
    cv_max         = round(cc$cv_max, 3),
    cv_base        = round(cc$cv_base, 3),
    dif_cv         = round(cc$dif, 3),
    n_comunes      = cc$n_comunes,
    n_cv_na        = res$n_cv_na,
    n_fuera_rango  = res$n_fuera_rango,
    vif_max        = round(vif_de(ESPECIFICACIONES[[nm]]), 2),
    var_limite     = res$limite,
    convergencia   = res$convergencia,
    coherente      = res$coherente,
    alineado       = res$alineado,
    rho_al_limite  = !is.na(res$rho) && abs(res$rho) > RHO_LIMITE,
    segundos       = round(segundos, 1),
    stringsAsFactors = FALSE))
}


# ==============================================================================
# PARTE 3 — DEPENDENCIA ESPACIAL EN LOS DATOS
# ------------------------------------------------------------------------------
# Contrastes sobre el estimador directo y sobre las covariables, siguiendo el
# estudio de caso de Pratesi y Salvati (2009). Ambas son variables observadas.
# ==============================================================================

cat("\n== DEPENDENCIA ESPACIAL EN LOS DATOS ==\n\n")

filas <- list()

pr <- pruebas_de(datos$pobreza_monetaria, lw_queen)

filas[[1]] <- cbind(
  data.frame(criterio_vecindad = "Queen, orden 1",
             serie = "pobreza_monetaria", tipo = "estimador directo",
             stringsAsFactors = FALSE), pr)

cat(sprintf("  %-28s %-18s I = %7.4f (p = %.4f)   C = %6.4f (p = %.4f)\n",
            "pobreza_monetaria", "estimador directo",
            pr$I, pr$I_p_perm, pr$C, pr$C_p_asint))
cat("\n")

for (v in vars_todas) {
  
  pr <- pruebas_de(datos[[v]], lw_queen)
  
  filas[[length(filas) + 1]] <- cbind(
    data.frame(criterio_vecindad = "Queen, orden 1", serie = v,
               tipo = "covariable", stringsAsFactors = FALSE), pr)
  
  cat(sprintf("  %-28s %-18s I = %7.4f (p = %.4f)   C = %6.4f (p = %.4f)\n",
              v, "covariable", pr$I, pr$I_p_perm, pr$C, pr$C_p_asint))
}

autocor_queen <- do.call(rbind, filas)

write.csv(autocor_queen,
          file.path(ruta_out, "esp_autocorrelacion_queen.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# PARTE 4 — AJUSTES BASE Y MODELO ESPACIAL
# ==============================================================================

cat("\n== AJUSTES NO ESPACIALES ==\n")

base_ajustes <- list(); base_avisos <- list(); base_cv <- list()

for (nm in names(ESPECIFICACIONES)) {
  
  ab <- ajustar_base(ESPECIFICACIONES[[nm]])
  
  base_ajustes[[nm]] <- ab$ajuste
  base_avisos[[nm]]  <- ab$avisos
  base_cv[[nm]]      <- if (is.null(ab$ajuste)) rep(NA_real_, n_dom)
  else precision_de(ab$ajuste)$cv
  
  cat(sprintf("  %-4s  alineado: %-5s  %s\n", nm,
              as.character(alineado(ab$ajuste)),
              if (is.null(ab$ajuste)) "ERROR en el ajuste"
              else if (length(ab$avisos) == 0) "sin avisos"
              else paste(length(ab$avisos), "aviso(s):",
                         paste(ab$avisos, collapse = " | "))))
}

cat("\n== MODELO ESPACIAL ==\n")

filas <- list(); esp_ajustes <- list(); res_queen <- list()

for (nm in names(ESPECIFICACIONES)) {
  
  t0  <- Sys.time()
  res <- ajustar_esp(ESPECIFICACIONES[[nm]], W)
  tt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  res_queen[[nm]] <- res
  
  if (is.null(res)) {
    filas[[length(filas) + 1]] <- fila_ajuste("Queen, orden 1", nm, NULL,
                                              base_cv[[nm]], tt)
    cat("  ", nm, " -> NO CONVERGE\n", sep = "")
    next
  }
  
  esp_ajustes[[nm]] <- res$ajuste
  
  fa <- fila_ajuste("Queen, orden 1", nm, res, base_cv[[nm]], tt)
  filas[[length(filas) + 1]] <- fa
  
  cat(sprintf("  %-4s  rho = %7.4f   CV %.3f (base %.3f, dif %+.3f)\n",
              nm, res$rho, fa$cv, fa$cv_base, fa$dif_cv))
}

ajuste_queen <- do.call(rbind, filas)

cat("\n")
print(ajuste_queen %>%
        select(especificacion, n_vars, rho, sigma2_u, gamma_medio,
               cv, cv_base, dif_cv, n_comunes, n_fuera_rango,
               var_limite, convergencia, alineado, rho_al_limite),
      row.names = FALSE)

write.csv(ajuste_queen, file.path(ruta_out, "esp_ajuste_queen.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# GUARDADO
# ==============================================================================

saveRDS(
  list(
    W                = W,
    nb_queen         = nb_queen,
    poligonos        = poligonos,
    crs_proyectado   = CRS_PROYECTADO,
    resumen_vec      = resumen_vec,
    especificaciones = ESPECIFICACIONES,
    autocor_queen    = autocor_queen,
    ajuste_queen     = ajuste_queen,
    ajustes_base     = base_ajustes,
    avisos_base      = base_avisos,
    cv_base          = base_cv,
    ajustes_esp      = esp_ajustes,
    semilla          = SEMILLA,
    B_perm           = B_PERM,
    tol              = TOL,
    maxit            = MAXIT,
    session          = sessionInfo()
  ),
  file.path(ruta_out, "espacial.rds"))

cat("\n============================================================\n")
cat("EVALUACION ESPACIAL COMPLETADA\n")
cat("============================================================\n")