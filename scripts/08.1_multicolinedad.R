library(car); library(dplyr); library(here)

matriz <- readRDS(here("output", "matriz_sae_transformada_v2.rds"))

mod_vif <- lm(pobreza_monetaria ~ extension + saber11_lectura_2022 +
                partos_calificados_2020 + dist_popayan_km + dist_cali_km +
                pct_contributivo + rural_disperso + tiene_coca +
                desplaz_tasa_1000 + categoria_617, data = matriz)

vif_tab <- data.frame(Variable = names(vif(mod_vif)),
                      VIF = round(as.numeric(vif(mod_vif)), 2))
print(vif_tab, row.names = FALSE)
cat(sprintf("\nVIF máximo: %.2f\n", max(vif_tab$VIF)))