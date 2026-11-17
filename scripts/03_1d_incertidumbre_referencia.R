# ==============================================================================
# 03_1d_incertidumbre_referencia.R
# ------------------------------------------------------------------------------
# La incertidumbre publicada de la referencia municipal DANE frente al nivel de
# pobreza y la ruralidad. Respalda el parrafo de precision sobre el agrupamiento
# por terciles de la referencia.
#
# SALIDA (por consola): correlaciones de Spearman en los 84 municipios.
# ==============================================================================
library(here)

m <- as.data.frame(readRDS(here("output", "matriz_sae_transformada_v2.rds")))

p      <- m$pobreza_monetaria
rmse   <- sqrt(m$varianza_pobreza) * 100          # error estandar publicado, pp
cv     <- sqrt(m$varianza_pobreza) / p * 100      # CV de la referencia, %
rural  <- m$pct_rural
resto  <- m$varianza_pobreza / (p * (1 - p))      # varianza sin el componente binomial

sp <- function(a, b) round(cor(a, b, method = "spearman"), 2)

res <- data.frame(
  par = c("RMSE referencia ~ pobreza",
          "RMSE referencia ~ ruralidad",
          "pobreza ~ ruralidad",
          "CV referencia ~ pobreza",
          "CV referencia ~ ruralidad",
          "varianza / [p(1-p)] ~ pobreza",
          "varianza / [p(1-p)] ~ ruralidad"),
  spearman = c(sp(rmse, p), sp(rmse, rural), sp(p, rural),
               sp(cv, p), sp(cv, rural), sp(resto, p), sp(resto, rural))
)
print(res, row.names = FALSE)

cat("\nRango de pobreza:", round(range(p), 3), "\n")
cat("Municipios con pobreza > 0,5:", sum(p > 0.5), "de", length(p), "\n")