# ==============================================================================
# 01_1_verificacion_imputacion.R
# ------------------------------------------------------------------------------
# Revisión, puntos 14.1 y 14.2: faltantes e imputación.
#
# La matriz tiene seis valores faltantes en variables distintas de la distancia
# por carretera (que no se trata como faltante, ver 01):
#
#   analfabetismo_2018 y analfabetismo_rural_2018 en Popayán y Cali. Se
#     recuperan de los planes de desarrollo municipal; como ninguno reporta el
#     dato rural, la tasa total se usa como aproximación (proxy).
#   sisben_grupoA_2024 en Jambaló. Se imputa por KNN (VIM::kNN, k = 5).
#
# El script documenta lo que la revisión pide y no modifica ninguna matriz.
#
#   A. KNN de Jambaló. Distancia y estandarización de los predictores,
#      donantes, valor imputado, sensibilidad a k y validación dejando un
#      municipio fuera (leave-one-out).
#   B. Proxy de analfabetismo rural en Popayán y Cali: relación entre la tasa
#      rural y la total en los demás municipios.
#   C. Presencia de estas variables en las especificaciones seleccionadas.
#
# DISTANCIA. VIM::kNN usa la distancia de Gower: promedio, sobre los
# predictores, de la diferencia absoluta dividida por el rango de cada uno.
# Esa división estandariza los predictores (todos pesan igual y ninguno domina
# por su escala), por lo que no se estandariza aparte. Aquí se programa con R
# base para mostrar los donantes; el script comprueba que el valor reproducido
# coincida con el de la matriz que dejó VIM::kNN.
#
# Debe correrse después de 01. Lee:
#   output/matriz_sae_final_v2.rds
#   output/seleccion_stepwise_completo_both.rds   (opcional, parte C)
#
# SALIDAS (output/)
#   imputacion_knn_donantes.csv
#   imputacion_knn_sensibilidad_k.csv
#   imputacion_knn_validacion.csv
#   proxy_analfabetismo_rural.csv
# ==============================================================================

library(here)

ruta_out <- here("output")

matriz <- readRDS(file.path(ruta_out, "matriz_sae_final_v2.rds"))
matriz <- as.data.frame(matriz)


# ==============================================================================
# A. KNN DE JAMBALÓ
# ==============================================================================

COD_JAMBALO <- "19364"

# Las mismas predictoras que usa 01 en la imputación.
PREDICTORAS <- c("pct_indigena_2024", "pct_rural", "pob_total",
                 "tasa_dependencia", "acueducto_2018", "internet_2018")

# La variable imputada es un conteo (personas en Sisbén grupo A). La matriz
# final guarda su proporción respecto a la población; el conteo se recupera
# de ella.
matriz$conteo_sisben_A <- matriz$pct_sisben_A / 100 * matriz$pob_total

# Distancia de Gower entre todos los pares de municipios.
X   <- matriz[, PREDICTORAS]
rng <- sapply(X, function(z) diff(range(z)))

D <- Reduce(`+`, lapply(PREDICTORAS, function(v) {
  abs(outer(X[[v]], X[[v]], "-")) / rng[[v]]
})) / length(PREDICTORAS)

rownames(D) <- colnames(D) <- matriz$cod_mun

jambalo   <- which(matriz$cod_mun == COD_JAMBALO)
donantes  <- setdiff(seq_len(nrow(matriz)), jambalo)

# Imputación como la hace VIM::kNN: mediana de los k donantes más cercanos.
# modo "conteo"     : mediana de los conteos (lo que hace 01).
# modo "proporcion" : mediana de las proporciones, llevada a la población del
#                     municipio (alternativa de sensibilidad).
imputar <- function(i, k, modo, pool) {
  vec <- pool[order(D[i, pool])][seq_len(k)]
  if (modo == "conteo") {
    median(matriz$conteo_sisben_A[vec]) / matriz$pob_total[i] * 100
  } else {
    median(matriz$pct_sisben_A[vec])
  }
}

# ---- A1. Donantes y valor imputado -------------------------------------------
vec5 <- donantes[order(D[jambalo, donantes])][1:5]

tabla_donantes <- data.frame(
  cod_mun      = matriz$cod_mun[vec5],
  municipio    = matriz$Municipio[vec5],
  distancia    = round(D[jambalo, vec5], 3),
  poblacion    = matriz$pob_total[vec5],
  conteo       = round(matriz$conteo_sisben_A[vec5]),
  pct_sisben_A = round(matriz$pct_sisben_A[vec5], 1)
)

pct_replica <- imputar(jambalo, 5, "conteo", donantes)
pct_matriz  <- matriz$pct_sisben_A[jambalo]

cat("== KNN DE JAMBALÓ (k = 5, distancia de Gower) ==\n")
print(tabla_donantes, row.names = FALSE)
cat("\nConteo imputado (mediana de los donantes):",
    median(matriz$conteo_sisben_A[vec5]), "\n")
cat("Población de Jambaló:", matriz$pob_total[jambalo], "\n")
cat("Proporción imputada (réplica):", round(pct_replica, 2), "%\n")
cat("Proporción en la matriz (VIM::kNN):", round(pct_matriz, 2), "%\n\n")

if (abs(pct_replica - pct_matriz) > 0.05) {
  warning("La réplica no reproduce el valor de la matriz: revisar la ",
          "versión de VIM o las predictoras.")
}

cat("Proporción de Sisbén A en los demás municipios: mediana",
    round(median(matriz$pct_sisben_A[donantes]), 1), "%, rango",
    round(min(matriz$pct_sisben_A[donantes]), 1), "a",
    round(max(matriz$pct_sisben_A[donantes]), 1), "%\n\n")

write.csv(tabla_donantes, file.path(ruta_out, "imputacion_knn_donantes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- A2. Sensibilidad a k ----------------------------------------------------
sens_k <- do.call(rbind, lapply(c(3, 5, 7, 10), function(k) {
  data.frame(k = k,
             pct_conteo     = round(imputar(jambalo, k, "conteo", donantes), 2),
             pct_proporcion = round(imputar(jambalo, k, "proporcion", donantes), 2))
}))

cat("== SENSIBILIDAD DEL VALOR IMPUTADO A k ==\n")
print(sens_k, row.names = FALSE)
cat("\n")

write.csv(sens_k, file.path(ruta_out, "imputacion_knn_sensibilidad_k.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- A3. Validación dejando un municipio fuera -------------------------------
# Se oculta, uno a uno, el valor de cada municipio observado, se imputa con los
# demás donantes y se compara con el valor real. La referencia es imputar con la
# mediana de todos los municipios.
real <- matriz$pct_sisben_A[donantes]

validar <- function(k, modo) {
  pred <- sapply(donantes, function(i) imputar(i, k, modo, setdiff(donantes, i)))
  data.frame(k = k, metodo = modo,
             MAE  = mean(abs(pred - real)),
             RMSE = sqrt(mean((pred - real)^2)))
}

validacion <- do.call(rbind, c(
  lapply(c(1, 3, 5, 7, 10, 15), function(k) validar(k, "conteo")),
  lapply(c(1, 3, 5, 7, 10, 15), function(k) validar(k, "proporcion"))
))
validacion <- rbind(
  validacion,
  data.frame(k = NA, metodo = "mediana global",
             MAE  = mean(abs(median(real) - real)),
             RMSE = sqrt(mean((median(real) - real)^2)))
)
validacion[, c("MAE", "RMSE")] <- round(validacion[, c("MAE", "RMSE")], 2)

cat("== VALIDACIÓN LEAVE-ONE-OUT (puntos porcentuales) ==\n")
print(validacion, row.names = FALSE)
cat("\n")

write.csv(validacion, file.path(ruta_out, "imputacion_knn_validacion.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# B. PROXY DE ANALFABETISMO RURAL EN POPAYÁN Y CALI
# ==============================================================================
# En Popayán y Cali la tasa rural se aproxima con la total. Para ver cuánto
# puede alejarse de la tasa rural efectiva, se examina la relación entre ambas
# tasas en los demás municipios.

COD_CAPITALES <- c("19001", "76001")
resto <- matriz[!matriz$cod_mun %in% COD_CAPITALES, ]

cociente <- resto$analfabetismo_rural_2018 / resto$analfabetismo_2018
ajuste   <- lm(analfabetismo_rural_2018 ~ analfabetismo_2018, data = resto)

capitales <- matriz[matriz$cod_mun %in% COD_CAPITALES, ]

proxy <- data.frame(
  cod_mun          = capitales$cod_mun,
  municipio        = capitales$Municipio,
  tasa_total       = capitales$analfabetismo_2018,
  proxy_usado      = capitales$analfabetismo_rural_2018,
  rural_esperada   = round(predict(ajuste,
                                   data.frame(analfabetismo_2018 = capitales$analfabetismo_2018)), 2),
  rural_cociente   = round(median(cociente) * capitales$analfabetismo_2018, 2)
)

cat("== PROXY DE ANALFABETISMO RURAL ==\n")
cat("Municipios de referencia:", nrow(resto), "\n")
cat("Tasa rural >= tasa total en", sum(cociente >= 1), "de", nrow(resto), "\n")
cat("Cociente rural/total: mediana", round(median(cociente), 2),
    ", rango intercuartílico", round(quantile(cociente, .25), 2), "a",
    round(quantile(cociente, .75), 2), "\n")
cat("Correlación rural-total:",
    round(cor(resto$analfabetismo_2018, resto$analfabetismo_rural_2018), 2),
    ", R2 de la recta:", round(summary(ajuste)$r.squared, 2), "\n")
cat("Rango de la tasa total en los municipios de referencia:",
    range(resto$analfabetismo_2018), "\n\n")
print(proxy, row.names = FALSE)
cat("\n")

write.csv(proxy, file.path(ruta_out, "proxy_analfabetismo_rural.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ==============================================================================
# C. ¿ENTRAN ESTAS VARIABLES EN LAS ESPECIFICACIONES SELECCIONADAS?
# ==============================================================================
# En las matrices de selección la proporción de Sisbén A se llama pct_sisben_A.

ruta_sel <- file.path(ruta_out, "seleccion_stepwise_completo_both.rds")

if (file.exists(ruta_sel)) {
  
  sel <- readRDS(ruta_sel)$tabla
  vars_imputadas <- c("pct_sisben_A", "analfabetismo_2018",
                      "analfabetismo_rural_2018")
  
  cat("== PRESENCIA EN LAS ESPECIFICACIONES SELECCIONADAS ==\n")
  for (i in seq_len(nrow(sel))) {
    covs <- strsplit(sel$covariables[i], " \\+ ")[[1]]
    esta <- intersect(vars_imputadas, covs)
    cat(sel$modelo[i], "(", sel$criterio[i], "):",
        if (length(esta) == 0) "ninguna de las tres variables"
        else paste(esta, collapse = ", "), "\n")
  }
  cat("\n")
  
} else {
  cat("No se encontró", ruta_sel, "; parte C omitida.\n\n")
}

cat("Verificación de la imputación completada.\n")