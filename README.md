Estimación en áreas pequeñas de la pobreza monetaria municipal en Cauca y Valle del Cauca

Código en R del trabajo de grado de la Maestría en Ciencia de Datos de la Pontificia Universidad Javeriana Cali. El proyecto estima la tasa de pobreza monetaria de 2024 en los 84 municipios de Cauca y Valle del Cauca mediante modelos de estimación en áreas pequeñas a nivel de área (Fay–Herriot). Las estimaciones directas del DANE se combinan con covariables municipales de capital humano, demografía, territorio, capacidad fiscal, seguridad social y conflicto armado.

Se evalúan tres familias de modelos sobre las mismas especificaciones:

Fay–Herriot convencional, en escala original y bajo las transformaciones logarítmica y logit.
Fay–Herriot espacial, con efectos de área autorregresivos y contigüidad de primer orden.
Fay–Herriot robusto, con los predictores REBLUP y REBLUP-BC.
Estructura
SAEtesis.Rproj
scripts/            Código del análisis, en orden de ejecución
data/
  raw/              Estimaciones directas de pobreza y su varianza (DANE, 2024)
  pivoteadas/       Covariables municipales ya depuradas y cartografía (ver su README)
output/             Salidas de los scripts: matrices, tablas, objetos de ajuste y figuras
  figuras/          Diagnósticos de normalidad

Las rutas se resuelven con el paquete here, por lo que el proyecto debe abrirse desde SAEtesis.Rproj.

La reproducción parte de las covariables ya depuradas de data/pivoteadas/. El procesamiento previo de las descargas originales no se incluye en el repositorio; el contenido y la fuente de cada insumo se documentan en data/pivoteadas/README.md. El trabajo usa las estimaciones y varianzas publicadas por el DANE, no los microdatos de la encuesta.

Orden de ejecución
Script	Etapa	Contenido	Salida principal
00_ejecutar_pipeline.R	—	Ejecuta los scripts siguientes en orden, cada uno en un proceso de R separado	—
01_etapa1_construccion_base.R	I	Integración de fuentes por código Divipola, medidas relativas, tratamiento de faltantes, diagnóstico de asimetría (medcouple), transformaciones y análisis descriptivo	matriz_sae_transformada_v2.rds
02_1c_seleccion_stepwise.R	II	Selección de covariables sobre el modelo Fay–Herriot, por búsqueda en ambas direcciones con los criterios KIC, KICc y KICb2 (Marhuenda et al., 2014), bajo máxima verosimilitud y con VIF < 10	seleccion_stepwise_completo_both.rds
02_2a_modelado_variantes.R	II	Ajuste convencional por REML en escala original, logarítmica y logit, y sensibilidad al método de estimación de la varianza	modelado_variantes.rds
02_2a_modelado_variantes_diagnostico.R	II	Coeficientes de variación y factores de contracción de las variantes anteriores	modelado_diagnostico.rds
02_3_espacial.R	II	Moran y Geary sobre el estimador directo y las covariables, y ajuste del modelo espacial	espacial.rds
02_4_robusto.R	II	Ajuste robusto con k = 1,345 y sensibilidad a k y a la constante de corrección	robusto.rds
02_4_robusto_diagnostico.R	II	Desplazamiento de las estimaciones, estabilidad de coeficientes, comparación entre REBLUP y REBLUP-BC y sensibilidad	robusto_diagnostico.rds
03_1_normalidad.R	III	Diagnóstico descriptivo de normalidad de los efectos de área y los residuos con las herramientas de emdi	sup_normalidad.csv

La selección produce tres especificaciones anidadas, que se usan en todos los scripts posteriores:

Especificación	Criterio	Covariables
M1	KICb2	6
M2	KICc	7
M3	KIC	9
Cómo ejecutar
Abrir SAEtesis.Rproj en RStudio.
Ejecutar completo scripts/00_ejecutar_pipeline.R. En su vector PASOS se elige qué scripts correr.

02_1c_seleccion_stepwise.R tarda cerca de tres horas por el bootstrap de KICb2 (B = 500). Su salida está incluida en output/ y el script maestro lo deja desactivado por defecto.

Los scripts que usan plot() de emdi no deben ejecutarse línea por línea: el paquete pide pulsar Enter entre figuras, y esas pausas consumen las líneas siguientes del script. El script maestro evita el problema.

Paquetes
Modelado: emdi
Datos espaciales: sf, spdep
Manejo de datos: dplyr, tidyr, readxl, here
Preparación de covariables: VIM (imputación KNN), robustbase (medcouple), car (VIF)
Salidas: ggplot2, knitr
Ejecución: callr

Desarrollado en R 4.5.

Estado

Etapas I y II completas; Etapa III en curso (diagnóstico de supuestos completado; precisión, consistencia agregada y selección del modelo final pendientes). La Etapa IV (determinantes, mapas y tablero interactivo) está pendiente.
