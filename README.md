# Pobreza monetaria municipal en Cauca y Valle del Cauca, 2024

Código de la tesis de Maestría en Ciencia de Datos (Pontificia Universidad
Javeriana Cali): estimación de la pobreza monetaria en los 84 municipios de
Cauca y Valle del Cauca mediante modelos Fay-Herriot.

## Versión de la tesis

Los resultados del documento corresponden a la versión **`v1.0-tesis`** de este
repositorio (pestaña *Releases*). La rama principal puede cambiar después de la
entrega; para reproducir la tesis, descargue esa versión y no la rama.

## Estructura

```
scripts/            Scripts del análisis, en orden de ejecución
data/raw/           Referencia municipal DANE y cifras departamentales
data/pivoteadas/    Covariables auxiliares ya procesadas y cartografía MGN 2025
output/             Salidas (se regeneran al ejecutar los scripts)
```

## Requisitos

- R 4.5.1. Las versiones de los paquetes con que se obtuvieron los resultados
  están en `versiones_paquetes.txt`, generado con
  `scripts/00_versiones_paquetes.R`.

## Cómo ejecutarlo

1. Abrir `SAEtesis.Rproj` en RStudio.
2. Ejecutar `scripts/00_ejecutar_pipeline.R`.

## Entradas

- `data/raw/varianza estimador.xlsx`: estimaciones municipales de pobreza
  monetaria publicadas por el DANE para 2024 y sus errores estándar (referencia
  municipal). Variable de interés del estudio.
- `data/raw/Agregados.xlsx`: cifras departamentales oficiales, usadas en la
  consistencia agregada y el benchmarking.
- `data/pivoteadas/`: covariables auxiliares por fuente y cartografía MGN 2025.

Las salidas de `output/` no se incluyen en el repositorio porque se regeneran
con los scripts. Se conservan las siguientes:

- `seleccion_stepwise_completo_both.rds`, cuya obtención tarda cerca de tres
  horas (el paso `02_1c` está desactivado por defecto en el pipeline).
- `final_estimaciones_municipales_m1_log.csv`, con las estimaciones finales de
  pobreza por municipio, su CV e intervalo al 95 %.
- `figuras/`, con los gráficos de diagnóstico del anexo de normalidad.
- `seleccion_kicb2_log.rds`, la selección en escala logarítmica del anexo de
  sensibilidad (el paso `02_1d` está desactivado por defecto).
- `estabilidad_kicb2/`, con los resúmenes de la estabilidad de la selección
  (57 corridas ejecutadas en una máquina virtual; el paso `02_1e` está
  desactivado por defecto).
- `dashboard_pobreza.html`, el tablero interactivo. Se abre con doble clic en
  cualquier navegador, sin conexión ni instalación.

## Etapas

| Scripts | Contenido |
|---|---|
| `01` y `01_1` | Construcción de la base integrada y verificación de la imputación |
| `02_1c` a `02_4g` | Selección de covariables (con su estabilidad y la selección en escala logarítmica) y ajuste de los modelos Fay-Herriot estándar, espacial y robusto |
| `03_0` a `04_1` | Supuestos, precisión, incertidumbre de la referencia, consistencia con las cifras departamentales, sensibilidad a la escala de selección y benchmarking |
| `05_1` | Contribución de las variables auxiliares mediante valores de Shapley (KICb2, MSE y CV) |
| `06_0` y `06_1` | Estimaciones finales, coeficientes y mapas |
| `07_2` | Análisis exploratorio de las especificaciones |
| `08` y `09` | Tablero interactivo (`output/dashboard_pobreza.html`) |

## Correspondencia con la tesis

Script que produce cada grupo de resultados y archivo de `output/` donde queda.

| Resultados de la tesis | Script | Salida principal |
|---|---|---|
| Base integrada, valores faltantes e imputación | `01_etapa1_construccion_base.R`, `01_1_verificacion_imputacion.R` | `matriz_sae_transformada_v2.rds`, `diagnostico_faltantes_pre.csv`, `imputacion_knn_*.csv` |
| Selección de covariables (M1, M2, M3) | `02_1c_seleccion_stepwise.R` | `seleccion_stepwise_completo_both.rds` |
| Estabilidad de la selección | `02_1e_estabilidad_seleccion_kicb2.R` | `estabilidad_kicb2/estabilidad_*.csv` |
| Selección en escala logarítmica (anexo) | `02_1d_seleccion_kicb2_log.R`, `03_9_evaluacion_seleccion_log.R` | `seleccion_kicb2_log.rds`, `sens_escala_*.csv` |
| FH bajo cada escala y diagnóstico distribucional | `02_2a_modelado_variantes.R`, `02_2a_modelado_variantes_diagnostico.R`, `03_0_normalidad.R` | `modelado_variantes.rds`, `mod_cv.csv`, `mod_normalidad.csv`, `sup_normalidad.csv` |
| Modelo espacial | `02_3_espacial.R`, `02_3b_ganancia_espacial.R` | `espacial.rds`, `esp_ajuste_queen.csv`, `esp_ganancia_mse.csv` |
| Modelo robusto y sus estimadores del MSE | `02_4_robusto.R` a `02_4g_figura_coeficientes.R` | `robusto.rds`, `rob_*.csv`, `mse_pseudo_boot_*.csv`, `warnholz_verificacion_*.csv` |
| Precisión, concordancia e intervalos (Etapa III) | `03_1_precision_diagnosticos.R` | `prec_cv_resumen.csv`, `prec_bias_*.csv`, `prec_brown_gof.csv`, `prec_coverage.csv`, `prec_intervalos.csv` |
| Consistencia con las cifras departamentales | `03_2_agregado_departamental.R` | `agr_comparacion.csv`, `agr_resumen_compatibilidad.csv` |
| Benchmarking (sensibilidad) | `04_1_benchmarking_m1log.R` | `bench_*.csv` |
| Estimaciones finales, coeficientes y mapas | `06_0_estimaciones_coef_cv_m1_log.R`, `06_1_mapas_m1_log.R` | `final_estimaciones_municipales_m1_log.csv`, `final_coeficientes_m1_log.csv`, `fig_mapas_m1_log.png` |
| Valores de Shapley (Etapa IV) | `05_1_shapley_etapa4_completo.R` | `shapley_kicb2_m1*_ranking.csv`, `shapley_mse_m1_log_ranking.csv`, `shapley_cv_m1_log_ranking.csv`, `auditoria_shapley_*.csv` |
| Tablero interactivo | `08_exportar_datos_dashboard.R`, `09_dashboard_construir.R` | `dashboard_pobreza.html` |
