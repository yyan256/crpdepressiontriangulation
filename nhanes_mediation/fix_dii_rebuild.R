# 52_fix_dii_rebuild.R -- post-hoc fix: recompute DII (skip NA params, fix vitE name),
#                         fix smoker/education/race coding, rebuild the 5-cycle
#                         analysis sample WITHOUT re-downloading the whole pipeline.
suppressPackageStartupMessages({ library(dplyr); library(nhanesA) })

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

# ---- 修复 smoker / education / race 编码（与 nhanes_prep.R 保持完全一致）----
#   smoker：SMQ020(一生吸过≥100支) + SMQ040(现在是否吸烟) 正确构造"当前吸烟者"
#   education：DMDEDUC2 正确标签（3=高中毕业，非"some college"）
#   race：合并 Mexican American + Other Hispanic → Hispanic（论文 4 类）

# --- 下载 SMQ040（full.rds 仅含 SMQ020，需补当前吸烟状态）---
smq_cycles <- data.frame(
  cyc = c("2005-2006","2007-2008","2009-2010","2015-2016","2017-2018"),
  pfx = c("D","E","F","I","J"), stringsAsFactors = FALSE)
safe_smq <- function(tbl, max_try = 5) {
  out <- NULL
  for (k in seq_len(max_try)) {
    out <- tryCatch(as.data.frame(nhanes(tbl, translated = FALSE)), error = function(e) NULL)
    if (!is.null(out) && nrow(out) > 0) return(out)
    Sys.sleep(3)
  }
  NULL
}
smq040_all <- lapply(seq_len(nrow(smq_cycles)), function(i) {
  pfx <- smq_cycles$pfx[i]
  tb <- safe_smq(paste0("SMQ_", pfx))
  if (is.null(tb) || !"SMQ040" %in% names(tb)) return(NULL)
  data.frame(SEQN = tb$SEQN, SMQ040 = suppressWarnings(as.numeric(as.character(tb$SMQ040))))
})
smq040_all <- bind_rows(Filter(Negate(is.null), smq040_all))
if (nrow(smq040_all) > 0) {
  df <- merge(df, smq040_all, by = "SEQN", all.x = TRUE)
  cat("SMQ040 downloaded & merged:", nrow(smq040_all), "rows across",
      length(unique(smq040_all$SMQ040)), "categories\n")
} else {
  cat("WARNING: could not download SMQ040 (network?); smoker will be NA\n")
}

if ("SMQ020" %in% names(df) && "SMQ040" %in% names(df)) {
  df$smoker <- ifelse(df$SMQ020 == 2, 0,
               ifelse(df$SMQ040 %in% c(1, 2), 1,
               ifelse(df$SMQ040 == 3, 0, NA)))
} else if ("SMQ020" %in% names(df)) {
  # 兜底：若 SMQ040 缺失（网络失败），明确标注 smoker 不可靠，避免静默产出 99.97%
  warning("SMQ040 not present; smoker left as NA. Re-run to download SMQ040.")
  df$smoker <- NA_real_
}
if ("DMDEDUC2" %in% names(df)) {
  df$education <- factor(dplyr::case_when(
    df$DMDEDUC2 %in% c(1, 2) ~ "< High school",
    df$DMDEDUC2 == 3 ~ "High school graduate",
    df$DMDEDUC2 %in% c(4, 5) ~ "> High school",
    TRUE ~ "Other"
  ))
}
if ("RIDRETH1" %in% names(df)) {
  df$race <- factor(dplyr::case_when(
    df$RIDRETH1 == 3 ~ "Non-Hispanic White",
    df$RIDRETH1 == 4 ~ "Non-Hispanic Black",
    df$RIDRETH1 %in% c(1, 2) ~ "Hispanic",
    TRUE ~ "Other"
  ))
}

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
cat("current smoker (unweighted) =", round(mean(df2$smoker, na.rm=TRUE), 4),
    " (should be ~0.20; historical bug gave ~0.9997)\n")
cat("education levels:", paste(levels(df2$education), collapse=", "), "\n")

# save
saveRDS(df2, "/root/autodl-tmp/nhanes/nhanes_merged_analysis_v2.rds")
cat("\nSaved nhanes_merged_analysis_v2.rds (n =", nrow(df2), ")\n")
cat("DONE\n")
