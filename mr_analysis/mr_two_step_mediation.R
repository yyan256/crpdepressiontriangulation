# 64_two_step_mr.R -- two-step MR mediation: nutrient -> CRP -> MDD (replace Figure 3 placeholder)
suppressPackageStartupMessages({ library(data.table); library(MendelianRandomization) })

VCFDIR <- "/root/autodl-tmp/ref/vcf"
PLINK  <- "/root/autodl-tmp/ref/plink2"
REF    <- "/root/autodl-tmp/ref/1kg.v3/EUR"
crp_mdd_cache <- "/root/autodl-tmp/results/vcf_parsed.rds"
full_cache    <- "/root/autodl-tmp/results/vcf_parsed_full.rds"
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

clump <- function(expo, tag) {
  sig <- expo[p < 5e-8][!duplicated(snp)]
  if (nrow(sig) == 0) return(NULL)
  cin <- file.path(VCFDIR, sprintf("%s_sig.txt", tag))
  fwrite(sig[, .(ID = snp, P = p)], cin, sep = "\t")
  cout <- file.path(VCFDIR, sprintf("%s_clump", tag))
  cmd <- sprintf("%s --bfile %s --clump %s --clump-p1 5e-8 --clump-r2 0.001 --clump-kb 10000 --clump-p2 5e-8 --out %s",
                 PLINK, REF, cin, cout)
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  cf <- paste0(cout, ".clumps")
  if (!file.exists(cf)) return(NULL)
  cl <- fread(cf)
  log("  [%s] clumped %d -> %d IVs", tag, nrow(sig), length(cl$ID))
  expo[snp %in% cl$ID]
}

harmonize <- function(expo, outco) {
  m <- merge(expo, outco, by = "snp", suffixes = c("_e", "_o"))
  m[, same := (a1_e == a1_o) & (a2_e == a2_o)]
  m[, swap := (a1_e == a2_o) & (a2_e == a1_o)]
  m <- m[same | swap]
  m[, beta_o := ifelse(swap, -beta_o, beta_o)]
  m[, .(snp, beta_e = beta_e, se_e = se_e, beta_o = beta_o, se_o = se_o)]
}

run_ivw <- function(expo, outco, tag) {
  iv <- clump(expo, tag)
  if (is.null(iv) || nrow(iv) < 1) { log("  [%s] no instruments", tag); return(NULL) }
  h <- harmonize(iv, outco)
  if (nrow(h) < 1) { log("  [%s] no harmonized IVs", tag); return(NULL) }
  mi <- mr_input(bx = h$beta_e, bxse = h$se_e, by = h$beta_o, byse = h$se_o, snps = h$snp)
  r <- mr_ivw(mi)
  data.frame(tag = tag, b = r@Estimate, se = r@StdError, p = r@Pvalue, n = nrow(h))
}

# --- load / parse (reuse crp/mdd cache; parse fiber/vitd once) ---
if (file.exists(full_cache)) {
  dat <- readRDS(full_cache)
} else {
  if (file.exists(crp_mdd_cache)) {
    d <- readRDS(crp_mdd_cache); crp <- d$crp; mdd <- d$mdd
  } else {
    crp <- parse_vcf(file.path(VCFDIR, "ieu-b-35.vcf.gz"))
    mdd <- parse_vcf(file.path(VCFDIR, "ieu-b-102.vcf.gz"))
  }
  log("parsing fiber + vitd VCF (large)")
  fiber <- parse_vcf(file.path(VCFDIR, "ukb-b-19085.vcf.gz"))
  vitd  <- parse_vcf(file.path(VCFDIR, "ieu-b-4812.vcf.gz"))
  dat <- list(fiber = fiber, vitd = vitd, crp = crp, mdd = mdd)
  saveRDS(dat, full_cache)
  log("cached -> %s", full_cache)
}
fiber <- dat$fiber; vitd <- dat$vitd; crp <- dat$crp; mdd <- dat$mdd
log("fiber=%d vitd=%d crp=%d mdd=%d SNPs", nrow(fiber), nrow(vitd), nrow(crp), nrow(mdd))

# --- Step 1: nutrient -> CRP (alpha) ---
log("Step 1: nutrient -> CRP")
a_fiber <- run_ivw(fiber, crp, "fiber_to_CRP")
a_vitd  <- run_ivw(vitd,  crp, "vitd_to_CRP")

# --- Step 2: CRP -> MDD (beta) ---
log("Step 2: CRP -> MDD")
b_crp <- run_ivw(crp, mdd, "CRP_to_MDD")

# --- indirect effect = alpha * beta (delta-method SE) ---
calc <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NULL)
  ind <- a$b * b$b
  se  <- sqrt(a$b^2 * b$se^2 + b$b^2 * a$se^2)
  data.frame(path = paste0(a$tag, " x ", b$tag), alpha = a$b, alpha_se = a$se,
             beta = b$b, beta_se = b$se, indirect = ind, se = se,
             p = 2*pnorm(-abs(ind/se)))
}

log("===== results =====")
cat("-- alpha (Step 1) --\n"); print(a_fiber); print(a_vitd)
cat("-- beta (Step 2) --\n"); print(b_crp)
cat("-- indirect effect (alpha*beta) --\n")
print(rbind(calc(a_fiber, b_crp), calc(a_vitd, b_crp)))
log("DONE")
