# ==============================================================================
# 02_4e_comparacion_mse_pseudo_boot.R
# ------------------------------------------------------------------------------
# Estimación de la pobreza monetaria municipal mediante modelos de áreas
# pequeñas. Cauca y Valle del Cauca, 2024.
#
# ETAPA II.4e. Comparación del error cuadrático medio por pseudolinealización
# y por bootstrap en REBLUP y REBLUP-BC.
#
# Requiere output/robusto.rds, producido por 02_4_robusto.R, que contiene los
# ajustes principales con bootstrap (ajustes_robustos) y los de
# pseudolinealización (ajustes_robustos_pseudo).
#
# ------------------------------------------------------------------------------
# OBJETIVO
# ------------------------------------------------------------------------------
# Describir la sensibilidad de la precisión estimada al estimador del error
# cuadrático medio, manteniendo fija la especificación del modelo y la
# estimación puntual.
#
# La comparación no interpreta ninguno de los dos estimadores del ECM como
# "verdadero". Se estudia cuánto cambia la precisión estimada al sustituir la
# pseudolinealización por bootstrap.
#
# Se verifica y reporta:
#
#   1. misma especificación, constante de sintonización y estimaciones
#      puntuales bajo ambos métodos;
#
#   2. razón ECM_boot / ECM_pseudo por municipio;
#
#   3. distribución del coeficiente de variación bajo cada estimador;
#
#   4. identidad:
#
#          CV_boot / CV_pseudo
#          =
#          sqrt(ECM_boot / ECM_pseudo),
#
#      dado que las estimaciones puntuales son iguales;
#
#   5. comportamiento de REBLUP-BC frente a REBLUP en los municipios
#      identificados mediante la regla de recorte;
#
#   6. correspondencia entre municipios recortados y municipios en los que
#      cambia el ECM pseudolineal de REBLUP-BC.
#
# ------------------------------------------------------------------------------
# INTERPRETACION
# ------------------------------------------------------------------------------
# Las diferencias entre los estimadores del ECM se interpretan como sensibilidad
# de la medida de precisión al procedimiento de estimación del ECM.
#
# El bootstrap no se toma como referencia de verdad ni como criterio de
# superioridad. La comparación es descriptiva.
#
# ------------------------------------------------------------------------------
# MUNICIPIOS RECORTADOS
# ------------------------------------------------------------------------------
# El recorte de REBLUP-BC se determina a partir de la regla del predictor:
#
#     REBLUP_i < y_i - c * ee_i
#     o
#     REBLUP_i > y_i + c * ee_i
#
# donde:
#
#     c  = mult_constant
#     ee = sqrt(psi_i)
#
# La identificación del recorte es independiente del ECM.
#
# ------------------------------------------------------------------------------
# FUENTES METODOLOGICAS
# ------------------------------------------------------------------------------
# Chambers, Chandra y Tzavidis (2011).
# Warnholz (2016).
#
# ------------------------------------------------------------------------------
# SALIDAS
# ------------------------------------------------------------------------------
# output/mse_pseudo_boot_verificacion.csv
# output/mse_pseudo_boot_detalle.csv
# output/mse_pseudo_boot_resumen.csv
# output/mse_pseudo_boot_recorte_resumen.csv
# output/mse_pseudo_boot_recorte_municipios.csv
# ==============================================================================


# ==============================================================================
# 0. LIBRERIAS Y OBJETOS
# ==============================================================================

library(dplyr)
library(tidyr)
library(here)

select <- dplyr::select
filter <- dplyr::filter

ruta_out <- here("output")

rob <- readRDS(file.path(ruta_out, "robusto.rds"))

datos            <- rob$datos
ESPECIFICACIONES <- rob$especificaciones
K                <- rob$k_principal
C                <- rob$mult_constant
METODOS          <- c("reblup", "reblupbc")


# ==============================================================================
# 1. CONFIGURACION
# ==============================================================================

# Umbral descriptivo para identificar diferencias superiores al 20 % entre
# los estimadores del ECM.
#
# NO corresponde a un contraste estadistico.
UMBRAL_DIF <- 0.20

# Tolerancia relativa para considerar iguales dos ECM.
TOL_REL <- 1e-6

# Tolerancia absoluta para considerar identicas las estimaciones puntuales.
TOL_EST <- 1e-10


# ==============================================================================
# 2. DATOS DIRECTOS
# ==============================================================================

y  <- as.numeric(datos$pobreza_monetaria)
psi <- as.numeric(datos$varianza_pobreza)
ee <- sqrt(psi)

CV_directo <- ee / y

n_dom <- length(y)


# ==============================================================================
# 3. FUNCIONES AUXILIARES
# ==============================================================================

clave_pseudo <- function(nm, metodo) {
  paste(nm, metodo, format(K), sep = "_")
}

# En robusto.rds ambas versiones comparten la clave; difieren en la lista.
clave_boot <- function(nm, metodo) {
  paste(nm, metodo, format(K), sep = "_")
}


extraer <- function(aj, etiqueta) {
  
  if (is.null(aj)) {
    stop("No existe el ajuste: ", etiqueta)
  }
  
  dominio <- as.character(aj$ind$Domain)
  
  if (!identical(dominio, as.character(datos$cod_mun))) {
    stop(
      "El orden de dominios de ",
      etiqueta,
      " no coincide con datos$cod_mun."
    )
  }
  
  list(
    pred = as.numeric(aj$ind$FH),
    mse  = as.numeric(aj$MSE$FH)
  )
}


cv_de <- function(pred, mse) {
  sqrt(mse) / pred
}


# ==============================================================================
# 4. VERIFICACION DE QUE SOLO CAMBIA EL ESTIMADOR DEL ECM
# ==============================================================================

verificacion <- list()

for (nm in names(ESPECIFICACIONES)) {
  
  for (metodo in METODOS) {
    
    clave_p <- clave_pseudo(nm, metodo)
    clave_b <- clave_boot(nm, metodo)
    
    aj_p <- rob$ajustes_robustos_pseudo[[clave_p]]
    aj_b <- rob$ajustes_robustos[[clave_b]]
    
    ep <- extraer(aj_p, clave_p)
    eb <- extraer(aj_b, clave_b)
    
    # --------------------------------------------------------------------------
    # Misma especificacion
    # --------------------------------------------------------------------------
    
    k_igual <- isTRUE(
      all.equal(
        aj_p$model$k,
        aj_b$model$k
      )
    )
    
    mult_igual <- isTRUE(
      all.equal(
        aj_p$model$mult_constant,
        aj_b$model$mult_constant
      )
    )
    
    # --------------------------------------------------------------------------
    # Mismas estimaciones puntuales
    # --------------------------------------------------------------------------
    
    dif_estimacion <- max(
      abs(ep$pred - eb$pred)
    )
    
    # --------------------------------------------------------------------------
    # Metodo declarado para el ECM
    # --------------------------------------------------------------------------
    
    mse_pseudo_metodo <- aj_p$method$MSE_method
    mse_boot_metodo   <- aj_b$method$MSE_method
    
    verificacion[[paste(nm, metodo, sep = "_")]] <- data.frame(
      
      especificacion = nm,
      metodo = metodo,
      
      mse_pseudo = mse_pseudo_metodo,
      mse_boot   = mse_boot_metodo,
      
      k_igual = k_igual,
      mult_igual = mult_igual,
      
      dif_max_estimacion =
        dif_estimacion,
      
      stringsAsFactors = FALSE
    )
  }
}

verificacion <- do.call(rbind, verificacion)


# ------------------------------------------------------------------------------
# Comprobaciones
# ------------------------------------------------------------------------------

stopifnot(
  "La pseudolinealizacion no esta identificada correctamente" =
    all(
      verificacion$mse_pseudo ==
        "pseudo linearization"
    ),
  
  "El bootstrap no esta identificado correctamente" =
    all(
      verificacion$mse_boot ==
        "bootstrap"
    ),
  
  "La constante k difiere entre ajustes" =
    all(verificacion$k_igual),
  
  "La constante multiplicativa difiere entre ajustes" =
    all(verificacion$mult_igual),
  
  "Las estimaciones puntuales difieren entre ajustes" =
    all(
      verificacion$dif_max_estimacion <
        TOL_EST
    )
)


cat(
  "\nOK  Mismas especificaciones, misma configuracion y mismas ",
  "estimaciones puntuales.\n",
  "    La unica diferencia es el estimador del ECM.\n",
  sep = ""
)


# ==============================================================================
# 5. DETALLE POR MUNICIPIO
# ==============================================================================

detalle <- list()

for (nm in names(ESPECIFICACIONES)) {
  
  for (metodo in METODOS) {
    
    clave_p <- clave_pseudo(nm, metodo)
    clave_b <- clave_boot(nm, metodo)
    
    ep <- extraer(
      rob$ajustes_robustos_pseudo[[clave_p]],
      clave_p
    )
    
    eb <- extraer(
      rob$ajustes_robustos[[clave_b]],
      clave_b
    )
    
    cv_p <- cv_de(
      ep$pred,
      ep$mse
    )
    
    cv_b <- cv_de(
      eb$pred,
      eb$mse
    )
    
    # --------------------------------------------------------------------------
    # Razones
    # --------------------------------------------------------------------------
    
    razon_mse <- eb$mse / ep$mse
    
    razon_cv <- cv_b / cv_p
    
    detalle[[paste(nm, metodo, sep = "_")]] <-
      data.frame(
        
        especificacion = nm,
        metodo = metodo,
        
        cod_mun = datos$cod_mun,
        Municipio = datos$Municipio,
        
        estimacion = ep$pred,
        
        mse_pseudo = ep$mse,
        mse_boot   = eb$mse,
        
        cv_pseudo = cv_p,
        cv_boot   = cv_b,
        
        razon_mse = razon_mse,
        razon_cv  = razon_cv,
        
        cv_directo = CV_directo,
        
        pseudo_menor_directo =
          cv_p < CV_directo,
        
        boot_menor_directo =
          cv_b < CV_directo,
        
        stringsAsFactors = FALSE
      )
  }
}

detalle <- do.call(rbind, detalle)


# ==============================================================================
# 6. COMPROBACION DE LA IDENTIDAD ENTRE ECM Y CV
# ==============================================================================

# Dado que:
#
#       CV = sqrt(ECM) / theta_hat
#
# y theta_hat es igual bajo ambos metodos:
#
#       CV_boot / CV_pseudo
#       =
#       sqrt(ECM_boot / ECM_pseudo)

detalle$dif_identidad <-
  abs(
    detalle$razon_cv^2 -
      detalle$razon_mse
  )

dif_identidad_max <-
  max(
    detalle$dif_identidad,
    na.rm = TRUE
  )

cat(
  "\nMaxima diferencia numerica en la identidad ",
  "razones ECM-CV: ",
  signif(dif_identidad_max, 6),
  "\n",
  sep = ""
)

stopifnot(
  "No se cumple la identidad entre razon de ECM y razon de CV" =
    dif_identidad_max < 1e-10
)


# ==============================================================================
# 7. RESUMEN DEL ECM Y DEL CV
# ==============================================================================

resumen <- detalle %>%
  
  group_by(
    especificacion,
    metodo
  ) %>%
  
  summarise(
    
    n = dplyr::n(),
    
    # --------------------------------------------------------------------------
    # Razon ECM_boot / ECM_pseudo
    # --------------------------------------------------------------------------
    
    razon_mse_p10 =
      quantile(
        razon_mse,
        0.10,
        names = FALSE,
        na.rm = TRUE
      ),
    
    razon_mse_mediana =
      median(
        razon_mse,
        na.rm = TRUE
      ),
    
    razon_mse_p90 =
      quantile(
        razon_mse,
        0.90,
        names = FALSE,
        na.rm = TRUE
      ),
    
    razon_mse_max =
      max(
        razon_mse,
        na.rm = TRUE
      ),
    
    n_boot_mayor_20 =
      sum(
        razon_mse > 1 + UMBRAL_DIF,
        na.rm = TRUE
      ),
    
    n_boot_menor_20 =
      sum(
        razon_mse < 1 - UMBRAL_DIF,
        na.rm = TRUE
      ),
    
    # --------------------------------------------------------------------------
    # CV pseudolinealizacion
    # --------------------------------------------------------------------------
    
    cv_pseudo_media =
      mean(
        cv_pseudo,
        na.rm = TRUE
      ) * 100,
    
    cv_pseudo_mediana =
      median(
        cv_pseudo,
        na.rm = TRUE
      ) * 100,
    
    cv_pseudo_p90 =
      quantile(
        cv_pseudo,
        0.90,
        names = FALSE,
        na.rm = TRUE
      ) * 100,
    
    cv_pseudo_max =
      max(
        cv_pseudo,
        na.rm = TRUE
      ) * 100,
    
    # --------------------------------------------------------------------------
    # CV bootstrap
    # --------------------------------------------------------------------------
    
    cv_boot_media =
      mean(
        cv_boot,
        na.rm = TRUE
      ) * 100,
    
    cv_boot_mediana =
      median(
        cv_boot,
        na.rm = TRUE
      ) * 100,
    
    cv_boot_p90 =
      quantile(
        cv_boot,
        0.90,
        names = FALSE,
        na.rm = TRUE
      ) * 100,
    
    cv_boot_max =
      max(
        cv_boot,
        na.rm = TRUE
      ) * 100,
    
    # --------------------------------------------------------------------------
    # Comparacion con el directo
    # --------------------------------------------------------------------------
    
    n_pseudo_menor_directo =
      sum(
        pseudo_menor_directo,
        na.rm = TRUE
      ),
    
    n_boot_menor_directo =
      sum(
        boot_menor_directo,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    across(
      where(is.double),
      ~ round(.x, 3)
    )
  ) %>%
  
  as.data.frame()


# ==============================================================================
# 8. IDENTIFICACION DE MUNICIPIOS RECORTADOS
# ==============================================================================
#
# IMPORTANTE:
#
# El recorte no se define a partir del MSE.
# Se define a partir de la regla del REBLUP:
#
#       reblup < y - c*ee
#       o
#       reblup > y + c*ee
#
# ------------------------------------------------------------------------------

datos_recorte <- data.frame(
  
  cod_mun = datos$cod_mun,
  Municipio = datos$Municipio,
  directo = y,
  ee_directo = ee,
  limite_inferior = y - C * ee,
  limite_superior = y + C * ee,
  
  stringsAsFactors = FALSE
)


# ==============================================================================
# 9. REBLUP-BC FRENTE A REBLUP
# ==============================================================================

ancho <- detalle %>%
  
  select(
    especificacion,
    metodo,
    cod_mun,
    Municipio,
    estimacion,
    mse_pseudo,
    cv_pseudo,
    cv_boot
  ) %>%
  
  pivot_wider(
    names_from = metodo,
    
    values_from = c(
      estimacion,
      mse_pseudo,
      cv_pseudo,
      cv_boot
    )
  ) %>%
  
  left_join(
    datos_recorte %>%
      select(
        cod_mun,
        directo,
        ee_directo,
        limite_inferior,
        limite_superior
      ),
    by = "cod_mun"
  ) %>%
  
  mutate(
    
    # --------------------------------------------------------------------------
    # Regla de recorte
    # --------------------------------------------------------------------------
    
    recortado =
      
      estimacion_reblup < limite_inferior |
      estimacion_reblup > limite_superior,
    
    # --------------------------------------------------------------------------
    # Cambio relativo del ECM pseudolineal de BC frente a REBLUP
    # --------------------------------------------------------------------------
    
    razon_mse_bc_rb_pseudo =
      
      mse_pseudo_reblupbc /
      mse_pseudo_reblup,
    
    # --------------------------------------------------------------------------
    # Cambio relativo del CV de BC frente a REBLUP
    # --------------------------------------------------------------------------
    
    razon_bc_rb_pseudo =
      
      cv_pseudo_reblupbc /
      cv_pseudo_reblup,
    
    razon_bc_rb_boot =
      
      cv_boot_reblupbc /
      cv_boot_reblup,
    
    # --------------------------------------------------------------------------
    # Verificacion independiente
    # --------------------------------------------------------------------------
    #
    # Fuera del recorte, el MSE pseudolineal de BC debería coincidir con
    # REBLUP. Dentro del recorte debe diferir.
    # --------------------------------------------------------------------------
    
    cambio_mse_pseudo =
      
      abs(
        mse_pseudo_reblupbc /
          mse_pseudo_reblup -
          1
      ) > TOL_REL
  )


# ------------------------------------------------------------------------------
# Correspondencia recorte vs cambio en ECM
# ------------------------------------------------------------------------------

coincidencias_recorte <- sum(
  ancho$recortado ==
    ancho$cambio_mse_pseudo
)

cat(
  "\nCorrespondencia entre municipios recortados y municipios ",
  "con cambio en el ECM pseudolineal:\n",
  coincidencias_recorte,
  " de ",
  nrow(ancho),
  "\n",
  sep = ""
)

stopifnot(
  "El recorte no coincide con el cambio en el ECM pseudolineal" =
    coincidencias_recorte == nrow(ancho)
)


# ==============================================================================
# 10. RESUMEN DEL COMPORTAMIENTO POR RECORTE
# ==============================================================================

resumen_recorte <- ancho %>%
  
  group_by(
    especificacion,
    recortado
  ) %>%
  
  summarise(
    
    n = dplyr::n(),
    
    # --------------------------------------------------------------------------
    # REBLUP-BC / REBLUP bajo pseudolinealizacion
    # --------------------------------------------------------------------------
    
    razon_pseudo_med =
      median(
        razon_bc_rb_pseudo,
        na.rm = TRUE
      ),
    
    razon_pseudo_p10 =
      quantile(
        razon_bc_rb_pseudo,
        0.10,
        names = FALSE,
        na.rm = TRUE
      ),
    
    razon_pseudo_p90 =
      quantile(
        razon_bc_rb_pseudo,
        0.90,
        names = FALSE,
        na.rm = TRUE
      ),
    
    razon_pseudo_max =
      max(
        razon_bc_rb_pseudo,
        na.rm = TRUE
      ),
    
    # --------------------------------------------------------------------------
    # REBLUP-BC / REBLUP bajo bootstrap
    # --------------------------------------------------------------------------
    
    razon_boot_med =
      median(
        razon_bc_rb_boot,
        na.rm = TRUE
      ),
    
    razon_boot_p10 =
      quantile(
        razon_bc_rb_boot,
        0.10,
        names = FALSE,
        na.rm = TRUE
      ),
    
    razon_boot_p90 =
      quantile(
        razon_bc_rb_boot,
        0.90,
        names = FALSE,
        na.rm = TRUE
      ),
    
    razon_boot_max =
      max(
        razon_bc_rb_boot,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    across(
      where(is.double),
      ~ round(.x, 3)
    )
  ) %>%
  
  as.data.frame()


# ==============================================================================
# 11. DETALLE DE LOS MUNICIPIOS RECORTADOS
# ==============================================================================

municipios_recortados <- ancho %>%
  
  filter(recortado) %>%
  
  transmute(
    
    especificacion,
    
    cod_mun,
    
    Municipio,
    
    distancia_ee =
      round(
        (estimacion_reblup - directo) /
          ee_directo,
        3
      ),
    
    cv_reblup_pseudo =
      round(
        cv_pseudo_reblup * 100,
        2
      ),
    
    cv_reblupbc_pseudo =
      round(
        cv_pseudo_reblupbc * 100,
        2
      ),
    
    cv_reblup_boot =
      round(
        cv_boot_reblup * 100,
        2
      ),
    
    cv_reblupbc_boot =
      round(
        cv_boot_reblupbc * 100,
        2
      ),
    
    razon_bc_rb_pseudo =
      round(
        razon_bc_rb_pseudo,
        3
      ),
    
    razon_bc_rb_boot =
      round(
        razon_bc_rb_boot,
        3
      ),
    
    cambio_mse_pseudo
  ) %>%
  
  arrange(
    especificacion,
    desc(razon_bc_rb_pseudo)
  )


# ==============================================================================
# 12. RESULTADOS
# ==============================================================================

cat(
  "\n============================================================\n",
  "COMPARACION PSEUDOLINEALIZACION FRENTE A BOOTSTRAP\n",
  "============================================================\n\n",
  sep = ""
)


cat(
  "1. VERIFICACION DE LA CONFIGURACION\n\n"
)

print(
  verificacion,
  row.names = FALSE
)


cat(
  "\n------------------------------------------------------------\n",
  "2. RAZON ECM_boot / ECM_pseudo\n",
  "------------------------------------------------------------\n\n",
  sep = ""
)

print(
  resumen %>%
    select(
      especificacion,
      metodo,
      razon_mse_p10,
      razon_mse_mediana,
      razon_mse_p90,
      razon_mse_max,
      n_boot_mayor_20,
      n_boot_menor_20
    ),
  row.names = FALSE
)


cat(
  "\n------------------------------------------------------------\n",
  "3. COEFICIENTE DE VARIACION (%)\n",
  "------------------------------------------------------------\n\n",
  sep = ""
)

print(
  resumen %>%
    select(
      especificacion,
      metodo,
      cv_pseudo_media,
      cv_boot_media,
      cv_pseudo_mediana,
      cv_boot_mediana,
      cv_pseudo_p90,
      cv_boot_p90,
      cv_pseudo_max,
      cv_boot_max,
      n_pseudo_menor_directo,
      n_boot_menor_directo
    ),
  row.names = FALSE
)


cat(
  "\n------------------------------------------------------------\n",
  "4. REBLUP-BC FRENTE A REBLUP SEGUN RECORTE\n",
  "------------------------------------------------------------\n\n",
  sep = ""
)

print(
  resumen_recorte,
  row.names = FALSE
)


cat(
  "\n------------------------------------------------------------\n",
  "5. MUNICIPIOS RECORTADOS\n",
  "------------------------------------------------------------\n\n",
  sep = ""
)

print(
  municipios_recortados,
  row.names = FALSE
)


# ==============================================================================
# 13. LECTURA DESCRIPTIVA
# ==============================================================================

cat(
  "\n------------------------------------------------------------\n",
  "LECTURA\n",
  "------------------------------------------------------------\n",
  sep = ""
)

cat(
  "\nLa razon ECM_boot / ECM_pseudo se interpreta como una medida descriptiva ",
  "de sensibilidad del ECM al procedimiento de estimacion.\n",
  sep = ""
)

cat(
  "\nUna razon superior a 1 indica un ECM estimado mayor mediante bootstrap; ",
  "una razon inferior a 1 indica un ECM estimado menor.\n",
  sep = ""
)

cat(
  "\nEl umbral del 20 % se utiliza unicamente con fines descriptivos y no ",
  "corresponde a un contraste estadistico.\n",
  sep = ""
)

cat(
  "\nEn REBLUP-BC, la comparacion entre municipios recortados y no recortados ",
  "permite determinar si las diferencias observadas entre pseudolinealizacion ",
  "y bootstrap se concentran en las areas en las que actua la correccion.\n",
  sep = ""
)


# ==============================================================================
# 14. GUARDADO
# ==============================================================================

write.csv(
  verificacion,
  file.path(
    ruta_out,
    "mse_pseudo_boot_verificacion.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  detalle,
  file.path(
    ruta_out,
    "mse_pseudo_boot_detalle.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  resumen,
  file.path(
    ruta_out,
    "mse_pseudo_boot_resumen.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  resumen_recorte,
  file.path(
    ruta_out,
    "mse_pseudo_boot_recorte_resumen.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  municipios_recortados,
  file.path(
    ruta_out,
    "mse_pseudo_boot_recorte_municipios.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ==============================================================================
# 15. FINAL
# ==============================================================================

cat(
  "\n============================================================\n",
  "COMPARACION COMPLETADA\n",
  "============================================================\n",
  sep = ""
)