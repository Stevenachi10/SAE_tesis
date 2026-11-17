# Pobreza monetaria municipal en Cauca y Valle del Cauca, 2024

Código de la tesis de Maestría en Ciencia de Datos (Pontificia Universidad
Javeriana Cali): estimación de la pobreza monetaria en los 84 municipios de
Cauca y Valle del Cauca mediante modelos Fay-Herriot.

## Estructura

```
scripts/            Scripts del análisis, en orden de ejecución
data/               Estimaciones directas del DANE y covariables auxiliares
output/             Salidas (se regeneran al ejecutar los scripts)
```

## Cómo ejecutarlo

1. Abrir `SAEtesis.Rproj` en RStudio.
2. Ejecutar `scripts/00_ejecutar_pipeline.R`.

Las salidas de `output/` no se incluyen en el repositorio porque se regeneran
con los scripts. Se conservan tres:

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
