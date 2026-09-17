# ==============================================================================
# 02_4_robusto_diagnostico.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.4 — Extensión robusta del modelo Fay-Herriot (DIAGNOSTICO)
#
# Requiere output/robusto.rds, producido por 02_4_robusto.R.
#
# ------------------------------------------------------------------------------
# CRITERIOS
# ------------------------------------------------------------------------------
# COEFICIENTE DE VARIACION. Se calcula únicamente sobre dominios con estimación
# admisible: predicción finita en el intervalo (0, 1] y error cuadrático medio
# finito y no negativo. Los dominios excluidos se cuentan por causa y se
# reportan; no se promedian ni se descartan en silencio.
#
# La comparación del coeficiente de variación entre el modelo convencional y el
# robusto enfrenta estimadores distintos del error cuadrático medio: la
# aproximación analítica de Prasad-Rao en el primero y la linealización en el
# segundo. La comparación se reporta con esa salvedad y no sustenta por sí sola
# ninguna conclusión.
#
# Los umbrales del coeficiente de variación se reportan como indicadores
# descriptivos de la distribución. No constituyen criterios formales de
# publicación mientras no se adopte y se cite una referencia que los establezca.
#
# DESPLAZAMIENTO. La diferencia entre la estimación robusta y la convencional
# se denomina sensibilidad a la robustificación. Un desplazamiento apreciable
# indica que la estimación del dominio depende del tratamiento que reciben los
# residuos de magnitud elevada; no constituye por sí mismo evidencia de que el
# dominio sea influyente, afirmación que requeriría una medida de influencia
# específica.
#
# COEFICIENTES. Se reporta el valor bajo cada ajuste, la diferencia absoluta,
# esa diferencia expresada en unidades del error estándar del ajuste
# convencional y el cambio de signo. La cantidad expresada en unidades del
# error estándar es una medida descriptiva de magnitud relativa: no es el error
# estándar de la diferencia ni sustenta un contraste de significancia. No se
# emplea la diferencia relativa, que carece de sentido cuando el coeficiente de
# referencia se aproxima a cero.
#
# El emparejamiento de coeficientes entre ajustes se hace por nombre de término.
# Si alguno de los objetos no proporciona nombres, se recurre al emparejamiento
# por posición y se deja constancia en la columna correspondiente.
#
# PESOS DE LA FUNCION DE HUBER. No se reconstruyen a partir de la serie que el
# paquete devuelve como residuo estandarizado. Establecer sobre qué residuo
# escalado opera la robustificación exigiría verificar la implementación, y sin
# esa verificación los pesos reconstruidos podrían no ser los que el algoritmo
# aplicó.
#
# ------------------------------------------------------------------------------
# SALIDAS (output/)
# ------------------------------------------------------------------------------
#   rob_cv.csv             Distribución del coeficiente de variación.
#   rob_desplazamiento.csv Desplazamiento por dominio, constante principal.
#   rob_coeficientes.csv   Coeficientes bajo ambos ajustes.
#   rob_sensibilidad_k.csv Sensibilidad a la constante de sintonización.
#   robusto_diagnostico.rds
#   fig_rob_desplazamiento.png
#   fig_rob_cv.png
# ==============================================================================

library(emdi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")

rob    <- readRDS(file.path(ruta_out, "robusto.rds"))
matriz <- readRDS(file.path(ruta_out, "matriz_sae_transformada_v2.rds"))

datos            <- rob$datos
ESPECIFICACIONES <- rob$especificaciones
K                <- rob$k_principal
METODOS          <- rob$metodos
n_dom            <- nrow(datos)

# ---- Correspondencia de dominios --------------------------------------------
# Todas las comparaciones posteriores operan por posicion.
stopifnot(
  "El orden de dominios guardado no coincide con la matriz de entrada" =
    identical(as.character(datos$cod_mun), as.character(matriz$cod_mun)),
  "Hay codigos municipales duplicados" = !any(duplicated(datos$cod_mun)))

# Indicadores descriptivos de la distribucion del coeficiente de variacion.
CV_UMBRALES <- c(10, 15, 20)

tema_tesis <- theme_minimal(base_size = 13) +
  theme(plot.title = element_blank(), panel.grid.minor = element_blank())

COL_PUNTO <- "#b5482e"; COL_LINEA <- "#2c6e6b"; COL_BANDA <- "#c9c2b3"


# ==============================================================================
# FUNCIONES
# ==============================================================================

# ------------------------------------------------------------------------------
# Extraccion de predicciones y errores cuadraticos medios
# ------------------------------------------------------------------------------
# Se localiza la columna por nombre, se verifica la longitud y se comprueba que
# el orden de dominios del ajuste coincida con el de 'datos'. Cualquier
# discrepancia produce un fallo visible en lugar de una comparacion corrida.
extraer <- function(aj, etiqueta = "") {
  
  col_de <- function(df, patron) {
    if (is.null(df)) return(NULL)
    i <- grep(patron, names(df))
    if (length(i) == 0) return(NULL)
    as.numeric(df[[i[1]]])
  }
  
  dom <- tryCatch(as.character(aj$ind$Domain), error = function(e) NULL)
  
  if (is.null(dom) || !identical(dom, as.character(datos$cod_mun))) {
    stop("El orden de dominios del ajuste ", etiqueta,
         " no coincide con 'datos'.")
  }
  
  pred <- col_de(aj$ind, "^FH$")
  
  if (is.null(pred) || length(pred) != n_dom) {
    stop("No se pudo extraer la prediccion del ajuste ", etiqueta,
         " con la longitud esperada.")
  }
  
  mse <- col_de(aj$MSE, "^FH$")
  if (is.null(mse) || length(mse) != n_dom) mse <- rep(NA_real_, n_dom)
  
  list(pred = pred, mse = mse)
}

# ------------------------------------------------------------------------------
# Coeficiente de variacion con exclusion explicita
# ------------------------------------------------------------------------------
cv_de <- function(pred, mse) {
  
  admisible <- is.finite(pred) & pred > 0 & pred <= 1 &
    is.finite(mse) & mse >= 0
  
  cv <- rep(NA_real_, length(pred))
  cv[admisible] <- sqrt(mse[admisible]) / pred[admisible] * 100
  
  list(cv = cv,
       n_admisibles   = sum(admisible),
       n_pred_na      = sum(!is.finite(pred)),
       n_no_positiva  = sum(is.finite(pred) & pred <= 0),
       n_fuera_uno    = sum(is.finite(pred) & pred > 1),
       n_mse_invalido = sum(!is.finite(mse) | (is.finite(mse) & mse < 0)))
}

# ------------------------------------------------------------------------------
# Resumen de la distribucion del coeficiente de variacion
# ------------------------------------------------------------------------------
resumen_cv <- function(cv, umbrales = CV_UMBRALES) {
  
  v <- cv[is.finite(cv)]
  
  if (length(v) == 0) {
    base <- data.frame(n = 0L, media = NA_real_, mediana = NA_real_,
                       p25 = NA_real_, p75 = NA_real_, p90 = NA_real_,
                       minimo = NA_real_, maximo = NA_real_)
  } else {
    q <- unname(quantile(v, c(0.25, 0.75, 0.90)))
    base <- data.frame(n = length(v),
                       media   = round(mean(v), 3),
                       mediana = round(median(v), 3),
                       p25 = round(q[1], 3), p75 = round(q[2], 3),
                       p90 = round(q[3], 3),
                       minimo = round(min(v), 3), maximo = round(max(v), 3))
  }
  
  conteos <- as.data.frame(setNames(
    lapply(umbrales, function(u) sum(v > u)),
    paste0("n_cv_mayor_", umbrales)))
  
  cbind(base, conteos)
}

# ------------------------------------------------------------------------------
# Nombres de los terminos de un objeto de coeficientes
# ------------------------------------------------------------------------------
# Devuelve NULL si el objeto no trae nombres utilizables, de modo que el
# emparejamiento por posicion quede registrado en lugar de aplicarse en
# silencio.
nombres_coef <- function(cf) {
  if (is.null(cf)) return(NULL)
  rn <- rownames(cf)
  if (is.null(rn)) return(NULL)
  if (all(rn == as.character(seq_len(nrow(cf))))) return(NULL)
  rn
}

etiqueta_modelo <- function(nm, metodo, k) {
  paste(nm, metodo, format(k), sep = "_")
}


# ==============================================================================
# PARTE 1 — DISTRIBUCION DEL COEFICIENTE DE VARIACION
# ==============================================================================

cat("\n============================================================\n")
cat("DISTRIBUCION DEL COEFICIENTE DE VARIACION\n")
cat("============================================================\n")

filas <- list(); cv_largo <- list()

registrar_cv <- function(nm, metodo, mse_type, k, cc) {
  cbind(
    data.frame(especificacion = nm, metodo = metodo, mse_type = mse_type,
               k = k, stringsAsFactors = FALSE),
    resumen_cv(cc$cv),
    data.frame(n_excluidos    = n_dom - cc$n_admisibles,
               n_pred_na      = cc$n_pred_na,
               n_no_positiva  = cc$n_no_positiva,
               n_fuera_uno    = cc$n_fuera_uno,
               n_mse_invalido = cc$n_mse_invalido))
}

for (nm in names(ESPECIFICACIONES)) {
  
  aj <- rob$ajustes_base[[nm]]
  
  if (!is.null(aj)) {
    
    ex <- extraer(aj, paste(nm, rob$metodo_base))
    cc <- cv_de(ex$pred, ex$mse)
    
    filas[[length(filas) + 1]] <-
      registrar_cv(nm, rob$metodo_base, rob$mse_base, NA_real_, cc)
    
    cv_largo[[length(cv_largo) + 1]] <- data.frame(
      especificacion = nm, metodo = rob$metodo_base, cv = cc$cv,
      stringsAsFactors = FALSE)
  }
  
  for (metodo in METODOS) {
    
    mod <- etiqueta_modelo(nm, metodo, K)
    aj  <- rob$ajustes_robustos[[mod]]
    if (is.null(aj)) next
    
    ex <- extraer(aj, mod)
    cc <- cv_de(ex$pred, ex$mse)
    
    filas[[length(filas) + 1]] <-
      registrar_cv(nm, metodo, rob$mse_robusto, K, cc)
    
    cv_largo[[length(cv_largo) + 1]] <- data.frame(
      especificacion = nm, metodo = metodo, cv = cc$cv,
      stringsAsFactors = FALSE)
  }
}

tabla_cv <- do.call(rbind, filas)
cv_largo <- do.call(rbind, cv_largo)

cat("\n")
print(tabla_cv %>%
        select(especificacion, metodo, mse_type, n, media, mediana,
               p75, p90, maximo, n_excluidos),
      row.names = FALSE)

cat("\nDominios por encima de cada umbral (indicadores descriptivos):\n\n")
print(tabla_cv %>%
        select(especificacion, metodo, starts_with("n_cv_mayor_")),
      row.names = FALSE)

if (any(tabla_cv$n_excluidos > 0)) {
  cat("\n-- Dominios excluidos del calculo, por causa --\n")
  print(tabla_cv %>%
          filter(n_excluidos > 0) %>%
          select(especificacion, metodo, n_excluidos, n_pred_na,
                 n_no_positiva, n_fuera_uno, n_mse_invalido),
        row.names = FALSE)
}

cat("\nNota: el coeficiente de variacion del ajuste convencional procede de la\n")
cat("aproximacion analitica y el de los ajustes robustos de la linealizacion.\n")
cat("No son estimadores equivalentes del error cuadratico medio.\n")

write.csv(tabla_cv, file.path(ruta_out, "rob_cv.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# PARTE 2 — SENSIBILIDAD A LA ROBUSTIFICACION, POR DOMINIO
# ==============================================================================

cat("\n\n============================================================\n")
cat("SENSIBILIDAD A LA ROBUSTIFICACION\n")
cat("============================================================\n")

filas <- list()

for (nm in names(ESPECIFICACIONES)) {
  
  aj_base <- rob$ajustes_base[[nm]]
  if (is.null(aj_base)) next
  
  pred_base <- extraer(aj_base, paste(nm, rob$metodo_base))$pred
  
  for (metodo in METODOS) {
    
    mod <- etiqueta_modelo(nm, metodo, K)
    aj  <- rob$ajustes_robustos[[mod]]
    if (is.null(aj)) next
    
    pred_rob <- extraer(aj, mod)$pred
    
    filas[[length(filas) + 1]] <- data.frame(
      modelo = mod, especificacion = nm, metodo = metodo, k = K,
      cod_mun = datos$cod_mun, Municipio = datos$Municipio,
      directo = datos$pobreza_monetaria,
      varianza = datos$varianza_pobreza,
      base = pred_base, robusto = pred_rob,
      dif_pp = round((pred_rob - pred_base) * 100, 4),
      stringsAsFactors = FALSE)
  }
}

desplazamiento <- if (length(filas) > 0) do.call(rbind, filas) else data.frame()

if (nrow(desplazamiento) > 0) {
  
  cat("\nResumen del desplazamiento, en puntos porcentuales:\n\n")
  
  print(desplazamiento %>%
          group_by(especificacion, metodo) %>%
          summarise(media_abs   = round(mean(abs(dif_pp)), 4),
                    mediana_abs = round(median(abs(dif_pp)), 4),
                    maximo_abs  = round(max(abs(dif_pp)), 4),
                    n_mayor_1pp = sum(abs(dif_pp) > 1),
                    n_mayor_2pp = sum(abs(dif_pp) > 2),
                    .groups = "drop") %>%
          as.data.frame(),
        row.names = FALSE)
  
  cat("\nDominios con mayor desplazamiento:\n\n")
  
  print(desplazamiento %>%
          group_by(cod_mun, Municipio) %>%
          summarise(max_abs   = round(max(abs(dif_pp)), 4),
                    n_modelos = dplyr::n(),
                    varianza  = round(dplyr::first(varianza), 5),
                    .groups = "drop") %>%
          arrange(desc(max_abs)) %>%
          head(12) %>%
          as.data.frame(),
        row.names = FALSE)
  
  cat("\nEl desplazamiento mide sensibilidad al tratamiento de los residuos de\n")
  cat("magnitud elevada. No constituye por si mismo una medida de influencia.\n")
  
  write.csv(desplazamiento, file.path(ruta_out, "rob_desplazamiento.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}


# ==============================================================================
# PARTE 3 — COEFICIENTES
# ==============================================================================

cat("\n\n============================================================\n")
cat("COEFICIENTES\n")
cat("============================================================\n")

filas <- list()

for (nm in names(ESPECIFICACIONES)) {
  
  aj_base <- rob$ajustes_base[[nm]]
  if (is.null(aj_base)) next
  
  cf_base <- tryCatch(aj_base$model$coefficients, error = function(e) NULL)
  if (is.null(cf_base)) next
  
  nom_base <- nombres_coef(cf_base)
  
  for (metodo in METODOS) {
    
    mod <- etiqueta_modelo(nm, metodo, K)
    aj  <- rob$ajustes_robustos[[mod]]
    if (is.null(aj)) next
    
    cf_rob  <- tryCatch(aj$model$coefficients, error = function(e) NULL)
    if (is.null(cf_rob)) next
    
    nom_rob <- nombres_coef(cf_rob)
    
    # ---- Emparejamiento ------------------------------------------------------
    if (!is.null(nom_base) && !is.null(nom_rob) &&
        setequal(nom_base, nom_rob)) {
      
      orden   <- match(nom_base, nom_rob)
      terminos <- nom_base
      metodo_emparejamiento <- "nombre"
      
    } else {
      
      if (nrow(cf_base) != nrow(cf_rob)) {
        warning("Los objetos de coeficientes de ", mod,
                " no son emparejables; se omite la comparacion.")
        next
      }
      
      orden    <- seq_len(nrow(cf_rob))
      terminos <- if (!is.null(nom_base)) nom_base
      else c("(Intercept)", ESPECIFICACIONES[[nm]])
      metodo_emparejamiento <- "posicion"
      
      warning("Emparejamiento por posicion en ", mod,
              ": los objetos de coeficientes no traen nombres utilizables.")
    }
    
    b_base  <- as.numeric(cf_base$coefficients)
    se_base <- as.numeric(cf_base$std.error)
    b_rob   <- as.numeric(cf_rob$coefficients)[orden]
    
    filas[[length(filas) + 1]] <- data.frame(
      especificacion = nm, metodo = metodo, k = K,
      covariable = terminos,
      beta_base = signif(b_base, 4),
      beta_robusto = signif(b_rob, 4),
      dif_absoluta = signif(b_rob - b_base, 4),
      desplazamiento_ee_base = round((b_rob - b_base) / se_base, 3),
      cambia_signo = sign(b_base) != sign(b_rob),
      emparejamiento = metodo_emparejamiento,
      stringsAsFactors = FALSE)
  }
}

coeficientes <- if (length(filas) > 0) do.call(rbind, filas) else data.frame()

if (nrow(coeficientes) > 0) {
  
  cat("\n")
  print(coeficientes, row.names = FALSE)
  
  cat("\nCoeficientes que cambian de signo:",
      sum(coeficientes$cambia_signo), "de", nrow(coeficientes), "\n")
  cat("Desplazamiento maximo en unidades del error estandar convencional:",
      round(max(abs(coeficientes$desplazamiento_ee_base), na.rm = TRUE), 3),
      "\n")
  
  if (any(coeficientes$emparejamiento == "posicion")) {
    cat("\nAviso: hay comparaciones emparejadas por posicion.\n")
  }
  
  cat("\nEl desplazamiento en unidades del error estandar es una medida\n")
  cat("descriptiva de magnitud relativa; no es el error estandar de la\n")
  cat("diferencia ni sustenta un contraste de significancia.\n")
  
  write.csv(coeficientes, file.path(ruta_out, "rob_coeficientes.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}


# ==============================================================================
# PARTE 4 — SENSIBILIDAD A LA CONSTANTE DE SINTONIZACION
# ==============================================================================

cat("\n\n============================================================\n")
cat("SENSIBILIDAD A LA CONSTANTE DE SINTONIZACION\n")
cat("============================================================\n")

k_todas <- sort(unique(c(rob$k_principal, rob$k_sensibilidad)))

modelos_presentes <- names(rob$ajustes_robustos)

k_disponibles <- k_todas[vapply(k_todas, function(k) {
  any(grepl(paste0("_", format(k), "$"), modelos_presentes))
}, logical(1))]

if (length(k_disponibles) <= 1) {
  
  cat("\nNo hay ajustes con constantes adicionales en robusto.rds.\n")
  cat("Constantes presentes:", paste(k_disponibles, collapse = ", "), "\n")
  sensibilidad_k <- data.frame()
  
} else {
  
  filas <- list()
  
  for (nm in names(ESPECIFICACIONES)) {
    
    aj_base <- rob$ajustes_base[[nm]]
    if (is.null(aj_base)) next
    
    pred_base <- extraer(aj_base, paste(nm, rob$metodo_base))$pred
    
    for (metodo in METODOS) {
      for (k in k_disponibles) {
        
        mod <- etiqueta_modelo(nm, metodo, k)
        aj  <- rob$ajustes_robustos[[mod]]
        if (is.null(aj)) next
        
        ex <- extraer(aj, mod)
        cc <- cv_de(ex$pred, ex$mse)
        d  <- (ex$pred - pred_base) * 100
        
        filas[[length(filas) + 1]] <- data.frame(
          especificacion = nm, metodo = metodo, k = k,
          principal = isTRUE(all.equal(k, rob$k_principal)),
          cv_medio = round(mean(cc$cv, na.rm = TRUE), 3),
          cv_mediana = round(median(cc$cv, na.rm = TRUE), 3),
          dif_media_abs = round(mean(abs(d)), 4),
          dif_max_abs = round(max(abs(d)), 4),
          n_mayor_1pp = sum(abs(d) > 1),
          n_excluidos = n_dom - cc$n_admisibles,
          stringsAsFactors = FALSE)
      }
    }
  }
  
  sensibilidad_k <- do.call(rbind, filas)
  
  cat("\n")
  print(sensibilidad_k %>% arrange(especificacion, metodo, k),
        row.names = FALSE)
  
  cat("\nLa constante principal es k =", rob$k_principal,
      "; los restantes valores constituyen analisis secundario.\n")
  
  write.csv(sensibilidad_k, file.path(ruta_out, "rob_sensibilidad_k.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}


# ==============================================================================
# FIGURAS
# ==============================================================================

if (nrow(desplazamiento) > 0) {
  
  ggsave(file.path(ruta_out, "fig_rob_desplazamiento.png"),
         desplazamiento %>%
           ggplot(aes(x = base, y = robusto)) +
           geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                       color = COL_BANDA) +
           geom_point(color = COL_LINEA, alpha = 0.7, size = 1.8) +
           facet_grid(metodo ~ especificacion) +
           labs(x = "Estimación bajo el modelo convencional",
                y = "Estimación bajo el modelo robusto") +
           tema_tesis,
         width = 10, height = 7, dpi = 400)
}

if (!is.null(cv_largo) && nrow(cv_largo) > 0) {
  
  ggsave(file.path(ruta_out, "fig_rob_cv.png"),
         cv_largo %>% filter(is.finite(cv)) %>%
           ggplot(aes(x = metodo, y = cv)) +
           geom_boxplot(fill = "white", color = COL_LINEA,
                        outlier.color = COL_PUNTO, outlier.size = 1.4) +
           facet_wrap(~ especificacion) +
           labs(x = NULL, y = "Coeficiente de variación (%)") +
           tema_tesis +
           theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9)),
         width = 9, height = 5, dpi = 400)
}


# ==============================================================================
# GUARDADO
# ==============================================================================

saveRDS(
  list(
    tabla_cv       = tabla_cv,
    cv_largo       = cv_largo,
    desplazamiento = desplazamiento,
    coeficientes   = coeficientes,
    sensibilidad_k = sensibilidad_k,
    cv_umbrales    = CV_UMBRALES,
    k_principal    = rob$k_principal,
    session        = sessionInfo()
  ),
  file.path(ruta_out, "robusto_diagnostico.rds"))

cat("\n============================================================\n")
cat("DIAGNOSTICO COMPLETADO\n")
cat("============================================================\n")
