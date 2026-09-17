# ==============================================================================
# 03_1_normalidad.R
# DIAGNÓSTICOS DE COMPATIBILIDAD CON NORMALIDAD — MODELOS FH
# ==============================================================================
#
# Objetivo
# --------
# Evaluar descriptivamente la compatibilidad de las cantidades derivadas de
# modelos Fay-Herriot (FH) con supuestos de normalidad, utilizando:
#
#   1. QQ-plots
#   2. Asimetría (skewness)
#   3. Curtosis
#   4. Shapiro-Wilk
#
# IMPORTANTE:
# Los diagnósticos se realizan sobre cantidades estimadas/predichas por el
# modelo, no sobre las realizaciones latentes u_i y e_i, que no son observables.
#
# Por tanto, los resultados deben interpretarse como evidencia diagnóstica
# sobre la compatibilidad de las cantidades obtenidas con el comportamiento
# esperado bajo el modelo, y no como una observación directa de la normalidad
# de los efectos aleatorios o errores latentes.
#
# Todas las cifras provienen de summary() y todas las figuras de plot(), ambas
# del paquete emdi. El script no recalcula los estadísticos ni construye
# gráficos propios.
#
# ------------------------------------------------------------------------------
# AJUSTES EVALUADOS
# ------------------------------------------------------------------------------
# FH base, REML        M1, M2, M3 en escala original, log y logit       9
# FH espacial, orden 1 M1, M2                                           2
#
# M3 y M3s espaciales no se evalúan: la primera presenta un parámetro de
# autocorrelación prácticamente nulo y la segunda es un contraste auxiliar
# sobre la covariable de distancia.
#
# FH BASE
# -------
# Se examinan random_effects y std_real_residuals, con QQ-plots, asimetría,
# curtosis y Shapiro-Wilk leídos de forma conjunta.
#
# FH ESPACIAL
# -----------
# Se examinan las mismas cantidades, con especial atención a la estructura de
# dependencia espacial. Los random_effects son efectos de área espacialmente
# correlacionados; por tanto, un Shapiro-Wilk aplicado a esta serie NO
# constituye un contraste estándar de normalidad de innovaciones
# independientes. El valor-p se reporta como parte de la salida automática del
# paquete, pero NO se utiliza como criterio formal de aceptación o rechazo.
#
# FH ROBUSTO
# ----------
# No se incluyen contrastes de normalidad estándar. La estimación robusta
# modifica la influencia que las observaciones extremas tienen sobre la
# estimación, por lo que los diagnósticos de normalidad del FH clásico no se
# trasladan automáticamente a este modelo. Su comparación con el modelo
# convencional se presenta en 02_4_robusto_diagnostico.R.
#
# ==============================================================================
# INTERPRETACIÓN GENERAL
# ==============================================================================
#
# 1. random_effects
#
#    En el FH clásico:
#
#       y_i = x_i' beta + u_i + e_i,   u_i ~ N(0, sigma_u^2),  e_i ~ N(0, psi_i)
#
#    emdi proporciona predicciones de los efectos de área, uhat_i = E(u_i | y).
#    Por tanto, normalidad de u_i != normalidad de uhat_i. Los random_effects
#    presentan contracción hacia cero, y sus propiedades dependen de la varianza
#    de muestreo de cada dominio.
#
# 2. std_real_residuals
#
#    emdi utiliza real_residuals = y_i - x_i' beta_hat - uhat_i y construye
#    std_real_residuals = real_residuals / sqrt(psi_i). Estas cantidades no son
#    observaciones directas de e_i: normalidad de e_i != normalidad de
#    std_real_residuals.
#
# 3. Shapiro-Wilk
#
#    Un valor-p superior al nivel de significación indica que NO se encontró
#    evidencia estadística suficiente para rechazar la hipótesis nula. No
#    significa que se haya demostrado que la serie sea normal. La lectura debe
#    ser especialmente cautelosa cuando las cantidades analizadas son
#    predicciones, tienen varianzas que pueden diferir entre dominios o están
#    espacialmente correlacionadas.
#
# 4. QQ-plot
#
#    Se utiliza para identificar asimetría, colas pesadas o ligeras,
#    observaciones extremas y desviaciones sistemáticas. Tampoco demuestra por
#    sí solo que el supuesto sea verdadero.
#
#    emdi parametriza la distribución teórica con la media y la desviación
#    muestrales de la serie y usa la identidad como recta de referencia, por lo
#    que su QQ-plot no coincide con el de qqnorm() y qqline().
#
# 5. Curtosis
#
#    emdi la calcula con el paquete moments, sin restar 3: su valor de
#    referencia bajo normalidad es 3.
#
# ------------------------------------------------------------------------------
# ENTRADAS (output/)
# ------------------------------------------------------------------------------
#   modelado_variantes.rds
#   espacial.rds
#
# SALIDAS
# ------------------------------------------------------------------------------
#   output/sup_normalidad.csv
#   output/sup_normalidad.rds
#   output/figuras/qq_<clave>.png
#   output/figuras/dens_res_<clave>.png
#   output/figuras/dens_ran_<clave>.png
#
# Autor: Steven Venachi Pizo
# Proyecto: SAE — Pobreza Monetaria Municipal
# ==============================================================================

library(emdi)
library(ggplot2)
library(dplyr)
library(here)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")
ruta_fig <- file.path(ruta_out, "figuras")

if (!dir.exists(ruta_fig)) dir.create(ruta_fig, recursive = TRUE)

leer <- function(archivo) {
  ruta <- file.path(ruta_out, archivo)
  if (!file.exists(ruta)) {
    warning("No se encontro ", archivo, ". Esa familia queda fuera.")
    return(NULL)
  }
  readRDS(ruta)
}

mv  <- leer("modelado_variantes.rds")
esp <- leer("espacial.rds")


# ==============================================================================
# CONFIGURACION
# ==============================================================================

# Familia base. Claves de mv$ajustes: <espec>_<metodo>_<escala>.
BASE_ESPEC  <- c("M1", "M2", "M3")
BASE_METODO <- "reml"
BASE_ESCALA <- c("original", "log", "logit")

# Familia espacial. Claves de esp$ajustes_esp.
ESP_ESPEC <- c("M1", "M2")

# TRUE guarda las figuras de todos los ajustes; FALSE, solo las de
# GRAFICAR_CLAVES.
GRAFICAR        <- TRUE
GRAFICAR_TODO   <- TRUE
GRAFICAR_CLAVES <- c("M1_reml_original", "esp_M1")

PNG_ANCHO <- 2000
PNG_ALTO  <- 1100
PNG_RES   <- 200


# ==============================================================================
# FUNCIONES
# ==============================================================================

escalar <- function(x) {
  v <- suppressWarnings(tryCatch(as.numeric(x), error = function(e) NA_real_))
  if (length(v) < 1 || !is.finite(v[1])) NA_real_ else v[1]
}

# ------------------------------------------------------------------------------
# Extraccion de la matriz de normalidad de summary()
# ------------------------------------------------------------------------------
# summary()$normality es una tabla de estadisticos ya calculados, una fila por
# componente (residuos y efectos aleatorios) y una columna por estadistico. No
# contiene las series.
#
# El nombre de la columna del valor p se exige de manera explicita: un patron
# que solo busque "shapiro" devuelve primero el estadistico, cuyo rango
# habitual entre 0,9 y 1 puede confundirse con un valor p no significativo.
via_summary <- function(aj) {
  
  m <- tryCatch(summary(aj)$normality, error = function(e) NULL)
  
  if (is.null(m) || is.null(rownames(m)) || is.null(colnames(m))) return(NULL)
  
  col_de <- function(patron) {
    i <- grep(patron, colnames(m), ignore.case = TRUE)
    if (length(i) == 0) NA_character_ else colnames(m)[i[1]]
  }
  
  c_skew <- col_de("skew")
  c_kurt <- col_de("kurt")
  c_w    <- col_de("shapiro.*w$|^w$")
  c_p    <- if ("Shapiro_p" %in% colnames(m)) "Shapiro_p" else
    col_de("_p$|p.?value")
  
  f_ef  <- grep("random", rownames(m), ignore.case = TRUE)
  f_res <- grep("resid",  rownames(m), ignore.case = TRUE)
  
  saca <- function(fila, col) {
    if (length(fila) == 0 || is.na(col)) NA_real_ else escalar(m[fila[1], col])
  }
  
  list(ef_asim  = saca(f_ef,  c_skew), ef_curt  = saca(f_ef,  c_kurt),
       ef_W     = saca(f_ef,  c_w),    ef_p     = saca(f_ef,  c_p),
       res_asim = saca(f_res, c_skew), res_curt = saca(f_res, c_kurt),
       res_W    = saca(f_res, c_w),    res_p    = saca(f_res, c_p),
       n_filas  = nrow(m))
}

# ------------------------------------------------------------------------------
# Papel del contraste de Shapiro-Wilk segun la familia
# ------------------------------------------------------------------------------
INTERPRETACION_SHAPIRO <- function(familia) {
  if (familia == "base") "evaluacion conjunta" else "descriptivo"
}

# ------------------------------------------------------------------------------
# Filas de la tabla
# ------------------------------------------------------------------------------
fila_vacia <- function(familia, espec, escala, clave, nota) {
  data.frame(
    familia = familia, especificacion = espec, escala = escala,
    clave = clave, disponible = FALSE,
    interpretacion_shapiro = INTERPRETACION_SHAPIRO(familia), nota = nota,
    ef_asimetria = NA_real_, ef_curtosis = NA_real_,
    ef_shapiro_W = NA_real_, ef_shapiro_p = NA_real_,
    res_asimetria = NA_real_, res_curtosis = NA_real_,
    res_shapiro_W = NA_real_, res_shapiro_p = NA_real_,
    stringsAsFactors = FALSE)
}

fila_de <- function(aj, familia, espec, escala, clave) {
  
  if (is.null(aj) || !inherits(aj, "fh")) {
    cat(sprintf("  %-9s %-4s %-9s  AJUSTE NO DISPONIBLE\n",
                familia, espec, escala))
    return(fila_vacia(familia, espec, escala, clave, "ajuste ausente"))
  }
  
  r <- via_summary(aj)
  
  if (is.null(r) || all(is.na(unlist(r[names(r) != "n_filas"])))) {
    cat(sprintf("  %-9s %-4s %-9s  summary() SIN NORMALIDAD\n",
                familia, espec, escala))
    return(fila_vacia(familia, espec, escala, clave, "summary sin normalidad"))
  }
  
  aviso <- if (!is.na(r$n_filas) && r$n_filas > 2)
    paste0("summary devuelve ", r$n_filas, " filas") else NA_character_
  
  cat(sprintf("  %-9s %-4s %-9s  ok%s\n", familia, espec, escala,
              if (is.na(aviso)) "" else paste0("  [", aviso, "]")))
  
  data.frame(
    familia = familia, especificacion = espec, escala = escala,
    clave = clave, disponible = TRUE,
    interpretacion_shapiro = INTERPRETACION_SHAPIRO(familia), nota = aviso,
    ef_asimetria  = signif(r$ef_asim,  4),
    ef_curtosis   = signif(r$ef_curt,  4),
    ef_shapiro_W  = signif(r$ef_W,     4),
    ef_shapiro_p  = signif(r$ef_p,     4),
    res_asimetria = signif(r$res_asim, 4),
    res_curtosis  = signif(r$res_curt, 4),
    res_shapiro_W = signif(r$res_W,    4),
    res_shapiro_p = signif(r$res_p,    4),
    stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------------------
# Salida cruda de summary()$normality
# ------------------------------------------------------------------------------
cruda <- function(aj, etiqueta) {
  cat("\n-- ", etiqueta, " --\n", sep = "")
  if (is.null(aj)) { cat("ajuste ausente\n"); return(invisible(NULL)) }
  m <- tryCatch(summary(aj)$normality, error = function(e) conditionMessage(e))
  if (is.character(m)) {
    cat("summary() fallo: ", m, "\n", sep = ""); return(invisible(NULL))
  }
  if (is.null(m)) {
    cat("summary() no devolvio la matriz.\n"); return(invisible(NULL))
  }
  print(m)
}

# ------------------------------------------------------------------------------
# Figuras del paquete
# ------------------------------------------------------------------------------
# plot() sobre un objeto fh devuelve una lista: qq_plots es un gtable con los
# dos QQ combinados; density_res y density_ran son objetos ggplot; los demas
# elementos corresponden a modelos de nivel de unidad y valen FALSE.
#
# En sesion interactiva el paquete pide pulsar enter entre figuras; ejecutar
# con Rscript evita esas pausas. La version instalada de ggplot2 puede emitir
# avisos de fortify y descartar los colores previstos por el paquete: revisar
# las figuras antes de darlas por buenas.
graficar <- function(aj, clave) {
  
  if (!GRAFICAR) return(invisible(FALSE))
  
  if (is.null(aj) || !inherits(aj, "fh")) {
    cat("    figuras: ajuste no valido\n"); return(invisible(FALSE))
  }
  
  lst <- tryCatch(
    suppressWarnings(plot(aj, label = "no_title")),
    error = function(e) {
      cat("    figuras: plot() fallo -> ", conditionMessage(e), "\n", sep = "")
      NULL
    })
  
  if (is.null(lst)) return(invisible(FALSE))
  
  guardados <- character(0)
  
  # --- QQ-plots combinados (gtable) ---
  qq <- lst[["qq_plots"]]
  
  if (is.null(qq)) {
    cat("    qq_plots: ausente\n")
  } else if (!inherits(qq, "grob")) {
    cat("    qq_plots: clase inesperada -> ",
        paste(class(qq), collapse = "/"), "\n", sep = "")
  } else {
    f <- file.path(ruta_fig, paste0("qq_", clave, ".png"))
    ok <- tryCatch({
      grDevices::png(f, width = PNG_ANCHO, height = PNG_ALTO, res = PNG_RES)
      grid::grid.newpage()
      grid::grid.draw(qq)
      grDevices::dev.off()
      file.exists(f)
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      cat("    qq_plots: fallo al guardar -> ", conditionMessage(e), "\n",
          sep = "")
      FALSE
    })
    if (isTRUE(ok)) guardados <- c(guardados, basename(f))
  }
  
  # --- Densidades (ggplot) ---
  for (nm in c("density_res", "density_ran")) {
    
    g <- lst[[nm]]
    
    if (is.null(g)) { cat("    ", nm, ": ausente\n", sep = ""); next }
    
    if (!inherits(g, "gg")) {
      cat("    ", nm, ": clase inesperada -> ",
          paste(class(g), collapse = "/"), "\n", sep = ""); next
    }
    
    sufijo <- if (nm == "density_res") "dens_res_" else "dens_ran_"
    f <- file.path(ruta_fig, paste0(sufijo, clave, ".png"))
    
    ok <- tryCatch({
      suppressWarnings(
        ggplot2::ggsave(f, plot = g, width = PNG_ANCHO / PNG_RES,
                        height = PNG_ALTO / PNG_RES, dpi = PNG_RES))
      file.exists(f)
    }, error = function(e) {
      cat("    ", nm, ": fallo al guardar -> ", conditionMessage(e), "\n",
          sep = ""); FALSE
    })
    
    if (isTRUE(ok)) guardados <- c(guardados, basename(f))
  }
  
  if (length(guardados) > 0) {
    cat("    figuras guardadas: ", paste(guardados, collapse = ", "), "\n",
        sep = "")
  } else {
    cat("    figuras: ningun archivo generado\n")
  }
  
  invisible(length(guardados) > 0)
}

toca_graficar <- function(clave) {
  GRAFICAR && (GRAFICAR_TODO || clave %in% GRAFICAR_CLAVES)
}


# ==============================================================================
# BLOQUE A — ESTRUCTURA DE LA SALIDA DEL PAQUETE
# ==============================================================================

cat("\n============================================================\n")
cat("BLOQUE A — summary()$normality, salida cruda\n")
cat("============================================================\n")

cruda(mv$ajustes[["M1_reml_original"]], "base M1 original")

for (e in ESP_ESPEC) {
  cruda(if (e %in% names(esp$ajustes_esp)) esp$ajustes_esp[[e]] else NULL,
        paste("espacial", e))
}


# ==============================================================================
# BLOQUE B — FH BASE
# ==============================================================================

cat("\n\n============================================================\n")
cat("BLOQUE B — FH BASE, REML, tres escalas\n")
cat("============================================================\n")

filas <- list()

for (e_base in BASE_ESPEC) {
  for (s_base in BASE_ESCALA) {
    
    clave <- paste(e_base, BASE_METODO, s_base, sep = "_")
    aj    <- if (!is.null(mv) && clave %in% names(mv$ajustes))
      mv$ajustes[[clave]] else NULL
    
    filas[[length(filas) + 1]] <- fila_de(aj, "base", e_base, s_base, clave)
    
    if (toca_graficar(clave)) graficar(aj, clave)
  }
}


# ==============================================================================
# BLOQUE C — FH ESPACIAL
# ------------------------------------------------------------------------------
# La clave lleva prefijo para no colisionar con la del ajuste base de la misma
# especificacion en escala original.
# ==============================================================================

cat("\n============================================================\n")
cat("BLOQUE C — FH ESPACIAL, contiguidad de primer orden\n")
cat("============================================================\n")

for (e_esp in ESP_ESPEC) {
  
  clave <- paste0("esp_", e_esp)
  
  aj <- if (!is.null(esp) && e_esp %in% names(esp$ajustes_esp))
    esp$ajustes_esp[[e_esp]] else NULL
  
  filas[[length(filas) + 1]] <- fila_de(aj, "espacial", e_esp, "original",
                                        clave)
  
  if (toca_graficar(clave)) graficar(aj, clave)
}


# ==============================================================================
# TABLAS
# ==============================================================================

normalidad <- do.call(rbind, filas)

cat("\n\n============================================================\n")
cat("EFECTOS ALEATORIOS PREDICHOS\n")
cat("============================================================\n\n")

print(normalidad %>%
        select(familia, especificacion, escala, disponible,
               interpretacion_shapiro,
               ef_asimetria, ef_curtosis, ef_shapiro_W, ef_shapiro_p) %>%
        as.data.frame(),
      row.names = FALSE)

cat("\n\n============================================================\n")
cat("RESIDUOS REALIZADOS ESTANDARIZADOS\n")
cat("============================================================\n\n")

print(normalidad %>%
        select(familia, especificacion, escala, disponible,
               interpretacion_shapiro,
               res_asimetria, res_curtosis, res_shapiro_W, res_shapiro_p) %>%
        as.data.frame(),
      row.names = FALSE)

faltan <- normalidad %>% filter(!disponible)

if (nrow(faltan) > 0) {
  cat("\nAVISO: no se obtuvo normalidad para ", nrow(faltan),
      " ajuste(s).\n", sep = "")
  print(faltan %>% select(familia, especificacion, escala, nota),
        row.names = FALSE)
}

notas <- normalidad %>% filter(disponible, !is.na(nota))

if (nrow(notas) > 0) {
  cat("\nAVISO: summary() devuelve mas de dos filas en algunos ajustes.\n")
  cat("Verificar en el Bloque A cual corresponde a los residuos.\n")
  print(notas %>% select(familia, especificacion, escala, nota),
        row.names = FALSE)
}


# ==============================================================================
# GUARDADO
# ==============================================================================

write.csv(normalidad, file.path(ruta_out, "sup_normalidad.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(list(normalidad  = normalidad,
             base_espec  = BASE_ESPEC,
             base_metodo = BASE_METODO,
             base_escala = BASE_ESCALA,
             esp_espec   = ESP_ESPEC,
             session     = sessionInfo()),
        file.path(ruta_out, "sup_normalidad.rds"))

cat("\n============================================================\n")
cat("03_1_normalidad.R — FINALIZADO: ", sum(normalidad$disponible), " de ",
    nrow(normalidad), " ajustes con diagnostico disponible\n", sep = "")
cat("Figuras en ", ruta_fig, ": ",
    length(list.files(ruta_fig, pattern = "\\.png$")), " archivo(s)\n",
    sep = "")
cat("============================================================\n")