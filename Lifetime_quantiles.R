#   ____________________________________________________________________________
#   Compute groupwise longevity quantiles for quantile plots                ####

sepkoski00 = readRDS("sepkoski00.rds")

group_lab = c("Trilobita",
              "Brachiopoda",
              "Gastropoda",
              "Cephalopoda",
              "Anthozoa",
              "Bivalvia")

sepkoski00filtered = subset(sepkoski00, sepkoski00$grupoid < 7 &
                              sepkoski00$longevidad > 0)
sepkoski00filtered$FA = ceiling(sepkoski00filtered$FA)

lifetime_quantiles = matrix(NA, nrow = 6, ncol = n_quantiles_grp)
row.names(lifetime_quantiles) = group_lab
colnames(lifetime_quantiles) = quantiles_grp

for (j in 1:6) {
  subset_longevities = sepkoski00filtered[sepkoski00filtered$grupoid == j,
                                          "longevidad"]
  lifetime_quantiles[j,] = quantile(subset_longevities, quantiles_grp,
                                    na.rm = TRUE)
}

saveRDS(lifetime_quantiles, "lifetime_quantiles.rds")
