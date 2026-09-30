# 52_fix_dii_rebuild.R -- post-hoc fix: recompute DII (skip NA params, fix vitE name),
#                         rebuild the 5-cycle analysis sample WITHOUT re-downloading.
suppressPackageStartupMessages({ library(dplyr) })

df <- readRDS("/root/autodl-tmp/nhanes/nhanes_merged_full.rds")
cat("full rds nrow =", nrow(df), "\n")

# --- DII reference (17 params), with vit E corrected to DR1TATOC ---
dii_ref <- data.frame(
  param  = c("energy","protein","total_fat","sat_fat","trans_fat","cholesterol",
             "carbohydrate","fiber","vit_b6","vit_b12","vit_c","vit_d","vit_e",
             "iron","magnesium","zinc","alcohol"),
  var    = c("DR1TKCAL","DR1TPROT","DR1TTFAT","DR1TSFAT","DR1TTRF","DR1TCHOL",
             "DR1TCARB","DR1TFIBE","DR1TVB6","DR1TVB12","DR1TVC","DR1TVD","DR1TATOC",
             "DR1TIRON","DR1TMAGN","DR1TZINC","DR1TALCO"),
  effect = c( 0.180, 0.021, 0.298, 0.373, 0.229, 0.110,
              0.097,-0.663,-0.365, 0.106,-0.424,-0.446,-0.419,
              0.032,-0.484,-0.313,-0.278),
  mean   = c(2056, 79.4, 71.4, 28.6, 3.15, 279.4,
             272.2, 18.8, 1.47, 5.15, 118.2, 6.26, 8.73,
             13.35, 310.1, 9.84, 13.98),
  sd     = c(338, 13.9, 19.4, 8.0, 3.75, 51.2,
             40.0, 4.9, 0.74, 2.70, 43.46, 2.21, 1.49,
             3.71, 139.4, 2.19, 3.72),
  stringsAsFactors = FALSE
)

calc_dii_fixed <- function(df) {
  dii <- rep(0, nrow(df))
  any_valid <- rep(FALSE, nrow(df))
  n_used <- rep(0L, nrow(df))
  for (i in seq_len(nrow(dii_ref))) {
    v <- dii_ref$var[i]
    if (!v %in% names(df)) next
    x <- as.numeric(df[[v]])
    z  <- (x - dii_ref$mean[i]) / dii_ref$sd[i]
    cc <- 2 * pnorm(z) - 1
    contrib <- cc * dii_ref$effect[i]
    valid <- !is.na(contrib)
    any_valid <- any_valid | valid
    n_used <- n_used + as.integer(valid)
    dii <- dii + ifelse(valid, contrib, 0)
  }
  dii[!any_valid] <- NA_real_
  attr(dii, "n_used") <- n_used
  dii
}

# which nutrient params are actually available (columns present)?
avail <- intersect(dii_ref$var, names(df))
cat("available DII params (columns):", length(avail), "of 17\n")
cat("  ", paste(avail, collapse=", "), "\n")

# recompute DII
df$DII <- calc_dii_fixed(df)

# --- rebuild analysis sample (same filter as 20_nhanes_prep.R) ---
df2 <- df %>% filter(
  !is.na(DII), !is.na(crp_mgL), !is.na(phq9_total),
  !is.na(weight), !is.na(SDMVPSU), !is.na(SDMVSTRA)
)
cat("\n=== rebuilt analysis sample ===\n")
cat("n =", nrow(df2), "\n")
cat("by cycle:\n"); print(table(df2$cycle, useNA="ifany"))

# prevalence (unweighted)
cat("\nsomatic_high unweighted prevalence =", round(mean(df2$somatic_high, na.rm=TRUE), 4), "\n")
cat("cognitive_high unweighted prevalence =", round(mean(df2$cognitive_high, na.rm=TRUE), 4), "\n")
cat("DII mean =", round(mean(df2$DII, na.rm=TRUE), 3), " sd =", round(sd(df2$DII, na.rm=TRUE), 3), "\n")
cat("CRP geometric mean =", round(exp(mean(log(df2$crp_mgL), na.rm=TRUE)), 3), "\n")

# save
saveRDS(df2, "/root/autodl-tmp/nhanes/nhanes_merged_analysis_v2.rds")
cat("\nSaved nhanes_merged_analysis_v2.rds (n =", nrow(df2), ")\n")
cat("DONE\n")
