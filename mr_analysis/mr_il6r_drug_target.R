# 63_il6r_mr.R -- drug-target MR: IL6R region CRP variants -> MDD (replace Figure 5 fake data)
suppressPackageStartupMessages({ library(data.table); library(MendelianRandomization) })

VCFDIR <- "/root/autodl-tmp/ref/vcf"
CACHE  <- "/root/autodl-tmp/results/vcf_parsed.rds"
dir.create("/root/autodl-tmp/results", showWarnings = FALSE, recursive = TRUE)
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
  data.table(snp = ids, chr = dt[[1]], pos = dt$POS, a1 = dt$REF, a2 = dt$ALT,
             beta = beta, se = se, p = 10^(-lp))
}

# --- load or parse (cache to rds to avoid re-parsing big VCFs) ---
if (file.exists(CACHE)) {
  log("loading cached parsed VCF from %s", CACHE)
  dat <- readRDS(CACHE)
  crp <- dat$crp; mdd <- dat$mdd
} else {
  log("parsing VCF (cache miss)")
  crp <- parse_vcf(file.path(VCFDIR, "ieu-b-35.vcf.gz"))
  mdd <- parse_vcf(file.path(VCFDIR, "ieu-b-102.vcf.gz"))
  saveRDS(list(crp = crp, mdd = mdd), CACHE)
  log("cached parsed VCF -> %s", CACHE)
}
log("CRP SNPs=%d ; MDD SNPs=%d", nrow(crp), nrow(mdd))

# --- IL6R region: chr1 1q21, IL6R gene (hg19 chr1:154,377,669-154,441,926) +/- 100kb ---
lo <- 154277669; hi <- 154541926
il6r <- crp[chr == "1" & pos >= lo & pos <= hi & p < 5e-8]
log("IL6R region (chr1:%d-%d) SNPs p<5e-8: %d", lo, hi, nrow(il6r))
if (nrow(il6r) < 1) {
  il6r <- crp[chr == "1" & pos >= lo & pos <= hi & p < 5e-6]
  log("relaxed to p<5e-6: %d", nrow(il6r))
}
# clump within region (drug-target MR needs INDEPENDENT cis variants, not 100s of linked tags)
if (nrow(il6r) > 1) {
  cin <- "/root/autodl-tmp/ref/vcf/il6r_sig.txt"
  fwrite(il6r[, .(ID = snp, P = p)], cin, sep = "\t")
  cmd <- sprintf("%s --bfile /root/autodl-tmp/ref/1kg.v3/EUR --clump %s --clump-p1 5e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 5e-8 --out /root/autodl-tmp/ref/vcf/il6r_clump",
                 "/root/autodl-tmp/ref/plink2", cin)
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  cf <- "/root/autodl-tmp/ref/vcf/il6r_clump.clumps"
  if (file.exists(cf)) {
    cl <- fread(cf)
    il6r <- il6r[snp %in% cl$ID]
    log("after clumping: %d independent IL6R instruments", nrow(il6r))
  }
}
print(il6r[, .(snp, pos, a1, a2, beta, se, p)])

if (nrow(il6r) < 1) { log("no IL6R-region SNP found, abort"); quit(save = "no", status = 0) }

# --- harmonize IL6R instruments with MDD outcome ---
m <- merge(il6r, mdd, by = "snp", suffixes = c("_e", "_o"))
m[, same := (a1_e == a1_o) & (a2_e == a2_o)]
m[, swap := (a1_e == a2_o) & (a2_e == a1_o)]
m <- m[same | swap]
m[, beta_o := ifelse(swap, -beta_o, beta_o)]
log("harmonized IL6R instruments: %d", nrow(m))
print(m[, .(snp, beta_e, se_e, beta_o, se_o)])

if (nrow(m) < 1) { log("no harmonized IL6R instrument"); quit(save = "no", status = 0) }

# --- MR estimate ---
# 说明：Wald ratio 的 p 值用 2*pnorm(-|wr|/wr_se)（正文 §3.9 报告值，如 0.275/0.816）；
#       而 vcf_parsed 缓存中保存的 p_outcome 是 MDD GWAS 原始 SNP-outcome 关联 p
#       （如 0.2714/0.8213），两者口径不同。正文与 CSV 的一致性见 README 说明。
if (nrow(m) == 1) {
  wr <- m$beta_o[1] / m$beta_e[1]
  wr_se <- abs(m$se_o[1] / m$beta_e[1])
  wr_p <- 2*pnorm(-abs(wr/wr_se))
  cat(sprintf("Wald ratio (1 variant %s): b=%.4f se=%.4f p=%.4f\n",
              m$snp[1], wr, wr_se, wr_p))
} else if (nrow(m) == 2) {
  wr <- m$beta_o / m$beta_e
  wr_se <- abs(m$se_o / m$beta_e)
  wr_p <- 2*pnorm(-abs(wr/wr_se))
  cat("Wald ratio per variant (p = Wald-ratio-based, NOT raw outcome p):\n")
  print(data.table(snp = m$snp, beta_crp = m$beta_e, beta_mdd = m$beta_o,
                   wald = wr, se = wr_se, p = wr_p))
  mi <- mr_input(bx = m$beta_e, bxse = m$se_e, by = m$beta_o, byse = m$se_o, snps = m$snp)
  ivw <- mr_ivw(mi)
  cat(sprintf("IVW (2 variants): b=%.4f se=%.4f p=%.4f\n", ivw@Estimate, ivw@StdError, ivw@Pvalue))
} else {
  mi <- mr_input(bx = m$beta_e, bxse = m$se_e, by = m$beta_o, byse = m$se_o, snps = m$snp)
  res <- mr_allmethods(mi, method = "main")
  v <- res@Values
  v <- v[!grepl("intercept", v$Method, ignore.case = TRUE), ]
  cat(sprintf("IL6R -> MDD MR (n=%d instruments):\n", nrow(m)))
  print(v[, c("Method","Estimate","Std Error","P-value")])
}
log("DONE")
