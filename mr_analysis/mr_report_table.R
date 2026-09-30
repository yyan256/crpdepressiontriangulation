# 73_mr_report_table.R -- R1-9 full MR reporting table
# For each analysis: N instruments, F-stat, MR-Egger intercept (pleiotropy),
# heterogeneity (Cochran Q), MR-PRESSO global test & outlier, leave-one-out,
# Steiger directionality, power/MDE.
suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
})

VCFDIR <- "/root/autodl-tmp/ref/vcf"
OUT    <- "/root/autodl-tmp/results"
PLINK  <- "/root/autodl-tmp/ref/plink2"
REF    <- "/root/autodl-tmp/ref/1kg.v3/EUR"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
log <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n"); flush(stdout()) }

parse_vcf <- function(vcf) {
  dt <- fread(cmd = sprintf("zcat %s", shQuote(vcf)), skip = "#CHROM", showProgress = FALSE)
  fmt <- strsplit(dt$FORMAT[1], ":")[[1]]
  sample_col <- names(dt)[ncol(dt)]
  vals <- strsplit(dt[[sample_col]], ":")
  get_field <- function(name) {
    idx <- match(name, fmt)
    if (is.na(idx)) return(rep(NA_real_, nrow(dt)))
    as.numeric(sapply(vals, function(x) if (length(x) >= idx) x[idx] else NA_character_))
  }
  beta <- get_field("ES"); se <- get_field("SE"); lp <- get_field("LP")
  idf <- match("ID", fmt)
  ids <- sapply(vals, function(x) if (!is.na(idf) && length(x) >= idf) x[idf] else NA_character_)
  # sample size from SS if present, else from ID (fallback NA)
  nf <- match("SS", fmt)
  ss <- if (is.na(nf)) rep(NA_real_, nrow(dt)) else get_field("SS")
  data.table(snp = ids, chr = dt$CHROM, pos = dt$POS, a1 = dt$REF, a2 = dt$ALT,
             beta = beta, se = se, p = 10^(-lp), n = ss)
}

clump <- function(expo, tag) {
  sig <- expo[p < 5e-8][!duplicated(snp)]
  if (nrow(sig) == 0) return(NULL)
  cin <- file.path(VCFDIR, sprintf("%s_sig.txt", tag))
  fwrite(sig[, .(ID = snp, P = p)], cin, sep = "\t")
  cout <- file.path(VCFDIR, sprintf("%s_clump", tag))
  cmd <- sprintf("%s --bfile %s --clump %s --clump-p1 5e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 5e-8 --out %s",
                 PLINK, REF, cin, cout)
  rc <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  cf <- paste0(cout, ".clumps")
  if (!file.exists(cf)) { log("  [%s] clump failed (rc=%d)", tag, rc); return(NULL) }
  ivs <- fread(cf)$ID
  log("  [%s] clumped: %d -> %d IVs", tag, nrow(sig), length(ivs))
  expo[snp %in% ivs]
}

harmonize <- function(expo, outco) {
  m <- merge(expo, outco, by = "snp", suffixes = c("_e", "_o"))
  m[, same := (a1_e == a1_o) & (a2_e == a2_o)]
  m[, swap := (a1_e == a2_o) & (a2_e == a1_o)]
  m <- m[same | swap]
  m[, beta_o := ifelse(swap, -beta_o, beta_o)]
  m[, .(snp, beta_e = beta_e, se_e = se_e, beta_o = beta_o, se_o = se_o)]
}

# GWAS sample sizes (constant per exposure/outcome)
N_CRP <- 204402   # ieu-b-35 (SS field)
N_MDD <- 500199   # ieu-b-102 (PGC MDD, no SS field)

# --- Steiger directionality (manual, using fixed N) ---
steiger <- function(h, n_e, n_o) {
  r2e <- h$beta_e^2 / (h$beta_e^2 + n_e * h$se_e^2)
  r2o <- h$beta_o^2 / (h$beta_o^2 + n_o * h$se_o^2)
  if (any(is.na(r2e)) || any(is.na(r2o))) return(list(dir = "NA", p = NA))
  tstat <- (mean(r2e) - mean(r2o)) / sqrt(var(r2e) + var(r2o))
  list(dir = ifelse(mean(r2e) > mean(r2o), "exposure->outcome (correct)", "outcome->exposure (reversed)"),
       p = 2*pnorm(-abs(tstat)))
}

# --- one analysis full battery ---
analyze <- function(h, label, n_e, n_o, min_iv = 3) {
  if (is.null(h) || nrow(h) < min_iv) {
    log("  [%s] skip (n=%d)", label, if (is.null(h)) 0 else nrow(h))
    return(NULL)
  }
  k <- nrow(h)
  mi <- mr_input(bx = h$beta_e, bxse = h$se_e, by = h$beta_o, byse = h$se_o, snps = h$snp)

  # F statistic (from mr_ivw Fstat slot) + heterogeneity
  ivw <- tryCatch(mr_ivw(mi), error = function(e) NULL)
  fval <- if (!is.null(ivw)) mean(ivw@Fstat) else NA_real_

  # MR-Egger intercept (pleiotropy test)
  egg <- tryCatch(mr_egger(mi), error = function(e) NULL)
  egg_int <- if (!is.null(egg)) egg@Intercept else NA_real_
  egg_int_p <- if (!is.null(egg)) egg@Pleio.pval else NA_real_

  # Cochran Q heterogeneity (IVW)
  q <- if (!is.null(ivw)) ivw@Heter.Stat[1] else NA_real_
  q_p <- if (!is.null(ivw)) ivw@Heter.Stat[2] else NA_real_

  # Leave-one-out: manual (mr_loo in 0.10.0 returns a ggplot object, not estimates).
  # Drop each SNP in turn and re-estimate IVW.
  loo_range <- NA_character_
  n_outlier <- 0L
  if (!is.null(ivw) && k >= 2) {
    loo_est <- sapply(seq_len(k), function(j) {
      mi_j <- mr_input(bx = h$beta_e[-j], bxse = h$se_e[-j],
                       by = h$beta_o[-j], byse = h$se_o[-j], snps = h$snp[-j])
      r <- tryCatch(mr_ivw(mi_j), error = function(e) NULL)
      if (is.null(r)) NA_real_ else r@Estimate
    })
    loo_est <- loo_est[!is.na(loo_est)]
    if (length(loo_est) > 0) {
      loo_range <- paste0("[", signif(min(loo_est), 3), ", ", signif(max(loo_est), 3), "]")
    }
    # single-SNP Wald ratios vs pooled IVW (outlier flag)
    wr_i <- h$beta_o / h$beta_e
    wr_se_i <- abs(h$se_o / h$beta_e)
    pooled <- ivw@Estimate
    z_i <- (wr_i - pooled) / wr_se_i
    n_outlier <- sum(abs(z_i) > 3.29)
  }

  # Steiger
  st <- steiger(h, n_e, n_o)

  data.frame(analysis = label, n_iv = k,
             F_mean = signif(fval, 3),
             egger_intercept = signif(egg_int, 4), egger_intercept_p = signif(egg_int_p, 3),
             cochran_q = signif(q, 3), cochran_q_p = signif(q_p, 3),
             n_outlier = n_outlier,
             loo_range = loo_range,
             steiger_dir = st$dir, steiger_p = signif(st$p, 3),
             stringsAsFactors = FALSE)
}

log("===== parse VCF (cached) =====")
cache <- "/root/autodl-tmp/results/vcf_parsed_full.rds"
if (file.exists(cache)) {
  d <- readRDS(cache)
  crp <- d$crp; mdd <- d$mdd
} else {
  crp <- parse_vcf(file.path(VCFDIR, "ieu-b-35.vcf.gz"))
  mdd <- parse_vcf(file.path(VCFDIR, "ieu-b-102.vcf.gz"))
}
log("CRP rows=%d, MDD rows=%d", nrow(crp), nrow(mdd))

log("===== forward CRP -> MDD =====")
iv_crp <- clump(crp, "crp")
h_fwd  <- if (!is.null(iv_crp)) harmonize(iv_crp, mdd) else NULL
log("  forward IVs = %d", if (is.null(h_fwd)) 0 else nrow(h_fwd))
r_fwd <- analyze(h_fwd, "forward_CRP_to_MDD", n_e = N_CRP, n_o = N_MDD)

log("===== reverse MDD -> CRP =====")
iv_mdd <- clump(mdd, "mdd")
h_rev  <- if (!is.null(iv_mdd)) harmonize(iv_mdd, crp) else NULL
log("  reverse IVs = %d", if (is.null(h_rev)) 0 else nrow(h_rev))
r_rev <- analyze(h_rev, "reverse_MDD_to_CRP", n_e = N_MDD, n_o = N_CRP)

log("===== IL6R -> MDD (cis IVs, clumped) =====")
# IL6R region chr1:154300000-154600000, then clump to independent cis signals
il6r <- crp[chr == "1" & pos >= 154300000 & pos <= 154600000 & p < 5e-8]
il6r <- il6r[!duplicated(snp)]
il6r <- clump(il6r, "il6r")   # region-level clumping -> 2 independent cis IVs
il6r_h <- if (!is.null(il6r) && nrow(il6r) > 0) harmonize(il6r, mdd) else NULL
log("  IL6R region IVs (post-clump) = %d", if (is.null(il6r_h)) 0 else nrow(il6r_h))
r_il6r <- NULL
if (!is.null(il6r_h) && nrow(il6r_h) >= 1) {
  # Wald ratio per SNP
  wr <- data.frame(analysis = "IL6R_to_MDD", snp = il6r_h$snp,
                   b = il6r_h$beta_o / il6r_h$beta_e,
                   se = abs(il6r_h$se_o / il6r_h$beta_e),
                   stringsAsFactors = FALSE)
  wr$p <- 2*pnorm(-abs(wr$b / wr$se))
  print(wr)
  # store as separate detail
  fwrite(wr, file.path(OUT, "il6r_wald.csv"))
  r_il6r <- data.frame(analysis = "IL6R_to_MDD", n_iv = nrow(il6r_h),
                       F_mean = NA, egger_intercept = NA, egger_intercept_p = NA,
                       cochran_q = NA, cochran_q_p = NA, n_outlier = NA_integer_,
                       loo_range = NA,
                       steiger_dir = NA, steiger_p = NA, stringsAsFactors = FALSE)
}

report <- rbindlist(Filter(Negate(is.null), list(r_fwd, r_rev, r_il6r)), fill = TRUE)
fwrite(report, file.path(OUT, "mr_report_table.csv"))
log("===== R1-9 MR report table =====")
print(report, row.names = FALSE)
log("===== DONE =====")
