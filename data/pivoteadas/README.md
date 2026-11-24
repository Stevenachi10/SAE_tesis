# Covariables auxiliares (`data/pivoteadas/`)

Archivos de entrada de `scripts/01_etapa1_construccion_base.R`. Cada archivo
tiene una fila por municipio (84 municipios de Cauca y Valle del Cauca) y
contiene las descargas de las fuentes oficiales ya depuradas y pivoteadas a
formato ancho. Los scripts de ese procesamiento previo no forman parte del
repositorio. Las fuentes y los años de cada variable se describen en el
capítulo de Datos de la tesis.

| Archivo | Columnas | Contenido |
|---|---|---|
| `afiliados_regimen_2022.rds` | 8 | Afiliados por régimen de salud, 2022 |
| `cov_carretera.rds` | 4 | Distancia por carretera, distancia de enganche a la red vial e indicador de aislamiento vial |
| `cov_demografia_etnia.rds` | 9 | Población indígena y afrodescendiente, menores de 15 años, índice de envejecimiento y Sisbén grupo A, 2024 |
| `cov_demografica_cauca_valle.rds` | 8 | Población total, proporción rural y urbana y tasa de dependencia |
| `cov_descripciongeneral.rds` | 9 | Categoría Ley 617, categoría de ruralidad, extensión, población, densidad y variación intercensal |
| `cov_desempeñofiscal.rds` | 5 | Índice de Desempeño Fiscal, 2024 |
| `cov_distancias.rds` | 5 | Distancias geodésicas a Popayán, a Cali y a la capital departamental |
| `cov_educacion_cauca_valle.rds` | 10 | Cobertura neta, deserción, analfabetismo y puntajes Saber 11 |
| `cov_ejecucion_ingresos.rds` | 10 | Ingresos totales, corrientes, tributarios y no tributarios, transferencias, SGP y recursos de capital |
| `cov_fiscal.rds` | 6 | Pagos administrativos, en educación, en salud, en servicios públicos y totales, 2024 |
| `cov_laboral.rds` | 7 | Cotizantes, empresas formales y ocupados formales, 2016 |
| `cov_pib.rds` | 9 | Composición sectorial y valor agregado municipal, 2023 |
| `cov_salud.rds` | 18 | Mortalidad, vacunación, controles prenatales, partos atendidos por personal calificado y fecundidad adolescente, 2020 |
| `cov_seguridad.rds` | 8 | Tasas de hurto, homicidio, extorsión y violencia intrafamiliar, 2019 |
| `cov_territorial_cauca_valle.rds` | 10 | Impuesto predial, predios y avalúo catastral |
| `cov_vivienda.rds` | 8 | Coberturas de acueducto, alcantarillado, energía e internet, 2018 |
| `terridata_conflicto_2023.rds` | 8 | Cultivos de coca, incidencia del conflicto armado, desplazamiento y víctimas, 2023 |
| `MGN2025_MPIO_GRAFICO/` | | Cartografía municipal MGN 2025 (DANE), recortada a Cauca y Valle del Cauca |
