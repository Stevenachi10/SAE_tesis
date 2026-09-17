# Insumos pivoteados

Archivos de entrada del script `01_etapa1_construccion_base.R`, salvo la
cartografía, que usa `02_3_espacial.R`.
Corresponden a las descargas de las fuentes oficiales ya depuradas
y organizadas en formato ancho, con una fila por municipio para los
84 dominios de estudio (Cauca y Valle del Cauca).

La variable dependiente no está en esta carpeta: se lee de
`data/raw/varianza estimador.xlsx` (DANE, estimaciones SAE de pobreza
monetaria municipal, 2024).

## afiliados_regimen_2022.rds

- Fuente: TerriData (ADRES)
- Dimensiones: 84 filas x 8 columnas
- Columnas: CodDepto, Departamento, CodMunicipio, Municipio, C_2022, E_2022, S_2022, I_2022

## cov_carretera.rds

- Fuente: elaboración propia. Distancia por vía terrestre a la capital departamental calculada con OSRM desde los centroides municipales de GADM
- Dimensiones: 84 filas x 4 columnas
- Columnas: municipio, dist_carretera_km, snap_km, aislado_vial
- Nota: no incluye la columna de departamento; ver la limitación descrita en el script 01

## cov_demografia_etnia.rds

- Fuente: TerriData (DANE); Sisbén IV (DNP) para `sisben_grupoA_2024`
- Dimensiones: 84 filas x 9 columnas
- Columnas: Código Departamento, Departamento, Código Entidad, Entidad, pct_indigena_2024, pct_afro_2024, pct_menores15_2024, indice_envejecimiento_2024, sisben_grupoA_2024

## cov_demografica_cauca_valle.rds

- Fuente: DANE, proyecciones de población municipal
- Dimensiones: 84 filas x 8 columnas
- Columnas: DP, DPNOM, MPIO, DPMP, pob_total, pct_rural, pct_urbano, tasa_dependencia

## cov_descripciongeneral.rds

- Fuente: TerriData (DNP, IGAC)
- Dimensiones: 84 filas x 9 columnas
- Columnas: Departamento, Código Entidad, Entidad, Categoría ley 617 de 2000_2024, Categoría de ruralidad_2000, Extensión_2017, Población total_2024, Densidad poblacional_2024, Variación porcentual intercensal 2005 - 2018_2018

## cov_desempeñofiscal.rds

- Fuente: DNP, Índice de Desempeño Fiscal
- Dimensiones: 84 filas x 5 columnas
- Columnas: Código.Departamento, Departamento, Código.Entidad, Entidad, Indice de Desempeño Fiscal 2024

## cov_distancias.rds

- Fuente: elaboración propia. Distancias geodésicas calculadas desde los centroides municipales de GADM
- Dimensiones: 84 filas x 5 columnas
- Columnas: departamento, municipio, dist_popayan_km, dist_cali_km, dist_capital_km

## cov_educacion_cauca_valle.rds

- Fuente: TerriData (MEN, ICFES, DANE - Censo 2018)
- Dimensiones: 84 filas x 10 columnas
- Columnas: Código.Departamento, Departamento, Código.Entidad, Entidad, Cobertura neta en educación - Total_2024, Tasa de deserción intra-anual del sector oficial en educación básica y media (Desde transición hasta once)_2024, Tasa de Analfabetismo (Censo)_2018, Tasa de Analfabetismo Rural (Censo)_2018, Puntaje promedio Pruebas Saber 11 - Matemáticas_2022, Puntaje promedio Pruebas Saber 11 - Lectura crítica_2022

## cov_ejecucion_ingresos.rds

- Fuente: CHIP (CUIPO), ejecución presupuestal de ingresos
- Dimensiones: 84 filas x 10 columnas
- Columnas: Departamento, Código DANE, Entidad, INGRESOS, INGRESOS CORRIENTES, INGRESOS TRIBUTARIOS, INGRESOS NO TRIBUTARIOS, TRANSFERENCIAS CORRIENTES, SISTEMA GENERAL DE PARTICIPACIONES, RECURSOS DE CAPITAL

## cov_fiscal.rds

- Fuente: CHIP (CUIPO), ejecución presupuestal de gastos
- Dimensiones: 84 filas x 6 columnas
- Columnas: Código.DANE, pagos_admin_2024, pagos_educacion_2024, pagos_salud_2024, pagos_serv_publicos_2024, pagos_total_2024

## cov_laboral.rds

- Fuente: TerriData (COMPLETAR entidad de origen)
- Dimensiones: 84 filas x 7 columnas
- Columnas: Código Departamento, Departamento, Código Entidad, Entidad, cotizantes_2016, empresas_formales_2016, pct_ocupados_form_2016

## cov_pib.rds

- Fuente: DANE, valor agregado municipal
- Dimensiones: 84 filas x 9 columnas
- Columnas: cod_dep, departamento, cod_mun, municipio, pct_primario_2023, pct_secundario_2023, pct_terciario_2023, valor_agregado_2023, peso_relativo_2023

## cov_salud.rds

- Fuente: TerriData (Ministerio de Salud, DANE - estadísticas vitales)
- Dimensiones: 84 filas x 18 columnas
- Columnas: Código.Departamento, Departamento, Código.Entidad, Entidad, tasa_mortalidad_2020, mort_materna_2020, mort_infantil_1a_2020, cobertura_pentavalente_2020, bajo_peso_2020, control_prenatal_4_2020, promedio_prenatales_2020, mort_infantil_5a_2020, mort_eda_2020, mort_ira_2020, mort_desnutricion_2020, fecundidad_15_19_2020, mort_neonatal_2020, partos_calificados_2020

## cov_seguridad.rds

- Fuente: TerriData (Policía Nacional)
- Dimensiones: 84 filas x 8 columnas
- Columnas: Código Departamento, Departamento, Código Entidad, Entidad, tasa_hurto_2019, tasa_homicidios_2019, tasa_extorsion_2019, tasa_violencia_intra_2019

## cov_territorial_cauca_valle.rds

- Fuente: TerriData (IGAC)
- Dimensiones: 84 filas x 10 columnas
- Columnas: Código Departamento, Departamento, Código Entidad, Entidad, Predial_2020, Recaudo efectivo por impuesto predial_2020, Total de predios_2022, Avalúo catastral total_2020, Porcentaje de predios urbanos_2022, Porcentaje de predios rurales_2022

## cov_vivienda.rds

- Fuente: TerriData (DANE - Censo 2018)
- Dimensiones: 84 filas x 8 columnas
- Columnas: Código Entidad, Entidad, Código Departamento, Departamento, acueducto_2018, alcantarillado_2018, energia_2018, internet_2018

## terridata_conflicto_2023.rds

- Fuente: TerriData (UNODC - SIMCI, UARIV, DNP)
- Dimensiones: 84 filas x 8 columnas
- Columnas: Código Entidad, Entidad, Código Departamento, Departamento, coca_2023, iica_2023, desplaz_2023, victimas_2023

## MGN2025_MPIO_GRAFICO/

- Fuente: DANE, Marco Geoestadístico Nacional 2025, nivel municipal (versión gráfica)
- Contenido: shapefile `MGN_ADM_MPIO_GRAFICO.shp` con sus archivos asociados (.dbf, .shx, .prj, .cpg), recortado a los municipios de Cauca (19) y Valle del Cauca (76) para reducir su tamaño. Las geometrías de esos municipios no se modificaron.
- Uso: construcción de la matriz de contigüidad en `02_3_espacial.R`

