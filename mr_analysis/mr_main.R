# 60_mr_main_v2.R -- full two-sample MR (VCF direct parse + plink2 clumping + harmonize + MR)
# Forward: CRP(ieu-b-35) -> MDD(ieu-b-102);  Reverse: MDD(ieu-b-102) -> CRP(ieu-b-35)
suppressPackageStartupMessages({ library(data.table); library(MendelianRandomization) })

OUT    <- "/root/autodl-tmp/results"
VCFDIR <- "/root/autodl-tmp/ref/vcf"
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
  data.table(snp = ids, chr = dt$CHROM, pos = dt$POS, a1 = dt$REF, a2 = dt$ALT,
             beta = beta, se = se, p = 10^(-lp))
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
  cl <- fread(cf)
  ivs <- cl$ID
  log("  [%s] clumped: %d significant -> %d independent IVs", tag, nrow(sig), length(ivs))
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

run_mr <- function(h, label) {
  if (is.null(h) || nrow(h) < 3) { log("  [%s] IVs < 3, skip", label); return(NULL) }
  mi <- mr_input(bx = h$beta_e, bxse = h$se_e, by = h$beta_o, byse = h$se_o, snps = h$snp)
  res <- tryCatch(mr_allmethods(mi, method = "main"), error = function(e) NULL)
  if (is.null(res)) { log("  [%s] MR estimation failed", label); return(NULL) }
  v <- res@Values
  v <- v[!grepl("intercept", v$Method, ignore.case = TRUE), ]
  out <- data.frame(analysis = label, method = v$Method, b = v$Estimate,
             se = v[["Std Error"]], lo = v$Estimate - 1.96*v[["Std Error"]],
             hi = v$Estimate + 1.96*v[["Std Error"]], p = v[["P-value"]],
             n_snp = nrow(h), stringsAsFactors = FALSE)
  # 显式随机效应 IVW（论文报告的主估计；mr_allmethods 的 IVW 默认为固定效应）
  ivw_re <- tryCatch(mr_ivw(mi, model = "random"), error = function(e) NULL)
  if (!is.null(ivw_re)) {
    out <- rbind(out, data.frame(analysis = label, method = "IVW (random effects)",
                 b = ivw_re@Estimate, se = ivw_re@StdError,
                 lo = ivw_re@Estimate - 1.96*ivw_re@StdError,
                 hi = ivw_re@Estimate + 1.96*ivw_re@StdError,
                 p = ivw_re@Pvalue, n_snp = nrow(h), stringsAsFactors = FALSE))
  }
  out
}

log("===== parse VCF =====")
crp <- parse_vcf(file.path(VCFDIR, "ieu-b-35.vcf.gz"))
mdd <- parse_vcf(file.path(VCFDIR, "ieu-b-102.vcf.gz"))
log("CRP SNPs = %d ; MDD SNPs = %d", nrow(crp), nrow(mdd))

log("===== forward: CRP -> MDD =====")
iv_crp <- clump(crp, "crp")
h_fwd  <- if (!is.null(iv_crp)) harmonize(iv_crp, mdd) else NULL
log("  forward harmonized IVs = %d", if (is.null(h_fwd)) 0 else nrow(h_fwd))
mr_fwd <- run_mr(h_fwd, "forward_CRP_to_MDD")

log("===== reverse: MDD -> CRP =====")
iv_mdd <- clump(mdd, "mdd")
h_rev  <- if (!is.null(iv_mdd)) harmonize(iv_mdd, crp) else NULL
log("  reverse harmonized IVs = %d", if (is.null(h_rev)) 0 else nrow(h_rev))
mr_rev <- run_mr(h_rev, "reverse_MDD_to_CRP")

mr_all <- rbindlist(Filter(Negate(is.null), list(mr_fwd, mr_rev)), fill = TRUE)
fwrite(mr_all, file.path(OUT, "mr_vcf_results.csv"))
log("========== MR summary ==========")
print(mr_all)
log("========== DONE ==========")
