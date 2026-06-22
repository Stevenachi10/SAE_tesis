# Estimación de Áreas Pequeñas (SAE) para Pobreza Monetaria en Cauca y Valle del Cauca

Este repositorio contiene el *pipeline* metodológico y el código en R para la estimación de la pobreza monetaria a nivel municipal en los departamentos del Cauca y Valle del Cauca (Colombia), utilizando modelos de Estimación de Áreas Pequeñas (SAE).

El proyecto implementa un enfoque de modelamiento a nivel de área (Fay-Herriot) para obtener mejores predictores lineales insesgados empíricos (EBLUP), combinando estimaciones directas con covariables demográficas, económicas, espaciales y de conflicto armado.

## 🗂️ Estructura del Repositorio

El proyecto está diseñado para ser reproducible mediante rutas relativas con el paquete `here`.

* `SAEtesis.Rproj`: Archivo principal del proyecto en R.
* `scripts/`: Código fuente, organizado secuencialmente.
* `data/`: Directorio de datos (parcialmente ignorado en GitHub por tamaño).
  * `raw/`: Bases originales y matriz de varianza del estimador directo.
  * `pivoteadas/`: Covariables procesadas en formato `.rds` y cartografía municipal.
* `output/`: Resultados exportados: figuras, tablas de comparación y estimaciones EBLUP finales.

*(Nota: los microdatos y algunas bases espaciales pesadas no se incluyen en el repositorio público por su tamaño. El flujo de trabajo parte de las matrices consolidadas.)*

## ⚙️ Pipeline Metodológico (orden de ejecución)

**Preparación de datos**

* **`01_consolidacion_matriz_sae.R`** — Ensamblaje de la base municipal. Cruza la variable objetivo (estimación directa de pobreza y su varianza, calculada con el método delta sobre `cvlog`) con las covariables (DNP, DANE, Terridata, distancias espaciales).
* **`02_transformacion_variables.R`** — *Feature engineering* y manejo de asimetría: evaluación del sesgo con el estadístico *Medcouple* (`robustbase`), transformaciones logarítmicas y creación de variables indicadoras para covariables con exceso de ceros.
* **`03_seleccion_variables.R`** — Selección de variables por consenso de tres métodos: LASSO (validación cruzada), *forward stepwise* (AIC) y *best subset* (BIC).

**Estimación**

* **`04_estimacion_fay_herriot.R`** — Ajuste del modelo Fay-Herriot clásico. Compara los modelos candidatos por coeficiente de variación (CV), significancia de los coeficientes y criterios de ajuste; incluye pruebas de robustez (codificación de Ley 617 y transformaciones de la respuesta) y extrae el EBLUP final (modelo M09).

**Análisis exploratorio e influencia**

* **`05_descriptivas.R`** — Estadística descriptiva y exploración: mapas de conflicto y presencia de coca, y mapa de calor de correlaciones de Spearman entre covariables.
* **`06_atipicos.R`** — Detección de observaciones atípicas por distancia de Mahalanobis y análisis de influencia (reajuste del modelo excluyendo capitales y municipios atípicos).

**Visualización**

* **`07_mapasdepobreza.R`** — Mapas de coropletas de la pobreza estimada (EBLUP) y del coeficiente de variación, con `ggplot2`, `sf` y `ggspatial`.

**Diagnósticos y validación**

* **`08_supuestos.R`** — Diagnóstico de supuestos: normalidad de residuos y efectos aleatorios (Shapiro-Wilk y gráficos cuantil-cuantil), homocedasticidad y autocorrelación espacial (índice de Moran). Incluye la comparación con el Fay-Herriot espacial (SFH).
* **`08.1_multicolinealidad.R`** — Análisis de multicolinealidad de las covariables mediante el factor de inflación de varianza (VIF).
* **`09_validacion.R`** — Validación de la precisión: comparación del MSE analítico (Prasad-Rao) frente al bootstrap paramétrico, y cuantificación de la ganancia de eficiencia respecto al estimador directo.

## 🔗 Correspondencia con el documento

| Script | Sección del documento |
|---|---|
| 01 – 03 | Metodología — Tratamiento de datos y selección de variables |
| 04 | Resultados — Modelo seleccionado y estimación EBLUP |
| 05 – 06 | Resultados — Análisis descriptivo y de influencia |
| 07 | Resultados — Distribución territorial de la pobreza (mapas) |
| 08, 08.1 | Resultados — Validación de supuestos |
| 09 | Resultados — Precisión y validación |

## 🛠️ Herramientas y Paquetes Principales

* **Manejo espacial y datos:** `sf`, `geodata`, `ggspatial`, `dplyr`, `tidyr`, `here`. Los límites municipales provienen del Marco Geoestadístico Nacional (DANE).
* **Selección y modelado:** `glmnet`, `leaps`, `sae`, `emdi`, `car`.
* **Diagnóstico espacial:** `spdep`.
* **Estadística robusta:** `robustbase`.
* **Visualización:** `ggplot2`, `patchwork`.
