# Estimación de Áreas Pequeñas (SAE) para Pobreza Monetaria en Cauca y Valle del Cauca

Este repositorio contiene el *pipeline* metodológico y el código en R para la estimación de la pobreza monetaria a nivel municipal en los departamentos del Cauca y Valle del Cauca (Colombia), utilizando modelos de Estimación de Áreas Pequeñas (SAE). 

El proyecto implementa un enfoque de modelamiento a nivel de área (Fay-Herriot) para obtener predictores empíricos insesgados (EBLUP), combinando estimaciones directas con covariables demográficas, económicas, espaciales y de conflicto armado.

## 🗂️ Estructura del Repositorio

El proyecto está diseñado para ser 100% reproducible utilizando rutas relativas mediante el paquete `here`.

*   `SAEtesis.Rproj`: Archivo principal del proyecto en R.
*   `scripts/`: Contiene el código fuente organizado secuencialmente.
*   `data/`: Directorio de datos (ignorado en GitHub por tamaño).
    *   `raw/`: Bases de datos originales y matriz de varianza del estimador directo.
    *   `pivoteadas/`: Covariables procesadas en formato `.rds`.
*   `output/`: Resultados exportados, tablas de comparación y estimaciones finales EBLUP.

*(Nota: Los microdatos y bases de datos espaciales no se incluyen en este repositorio público por restricciones de peso El flujo de trabajo parte de las matrices consolidadas).*

## ⚙️ Pipeline Metodológico (Orden de Ejecución)

El análisis está segmentado en cuatro etapas claras y reproducibles:

1.  **`01_consolidacion_matriz_sae.R`**
    Ensamblaje secuencial de la base de datos municipal. Cruza la variable objetivo (estimaciones directas de pobreza y su varianza real calculada mediante el método delta / `cvlog`) con múltiples covariables (DNP, DANE, Terridata, distancias espaciales).
2.  **`02_transformacion_variables.R`**
    *Feature engineering* y manejo de asimetría. Incluye evaluación del sesgo continuo real mediante el estadístico *Medcouple* (`robustbase`), transformaciones logarítmicas y creación de variables indicadoras (dummies) para covariables con exceso de ceros.
3.  **`03_seleccion_variables.R`**
    Filtrado . Implementa un consenso entre tres métodos estadísticos:
    *   Regularización LASSO (Validación cruzada).
    *   *Forward Stepwise* (Criterio de Información de Akaike - AIC).
    *   *Best Subset Selection* (Criterio de Información Bayesiano - BIC).
4.  **`04_estimacion_fay_herriot.R`**
    Ajuste y validación del modelo espacial base. Compara modelos candidatos evaluando el Coeficiente de Variación (CV), la significancia individual y los Factores de Inflación de Varianza (VIF). Extrae las predicciones finales (EBLUP) para el mapeo y las conclusiones.

## 🛠️ Herramientas y Paquetes Principales

El desarrollo de este *pipeline* se basa en los siguientes paquetes del entorno de R:
*   **Manejo Espacial y Datos:** `sf`, `geodata`, `dplyr`, `tidyr`, `here`.
*   **Selección y Modelado:** `glmnet`, `leaps`, `sae`, `car`.
*   **Estadística Robusta:** `robustbase`.
*   **Visualización:** `ggplot2`.
