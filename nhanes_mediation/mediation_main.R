# 53_med_v2.R -- mediation on corrected 5-cycle sample (n=24019)
suppressPackageStartupMessages({ library(survey); library(glmnet); library(dplyr) })
options(survey.lonely.psu = "adjust", warn = -1)

df <- readRDS("/root/autodl-tmp/nhanes/nhanes_merged_analysis_v2.rds")
cat(sprintf("analysis nrow=%d\n", nrow(df)))
for (v in c("race","education","cycle")) df[[v]] <- factor(df[[v]])
design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight, data = df, nest = TRUE)

COV_BASE <- c("RIDAGEYR","sex","race","education","smoker","drinker","cycle")
COV_BMI  <- c(COV_BASE, "BMXBMI")

fit_glmnet_coef <- function(outcome, mediator, exposure, covs, lambda=1e-4) {
  ff <- as.formula(sprintf("%s ~ %s + %s + %s", outcome, mediator, exposure, paste(covs, collapse=" + ")))
  allv <- all.vars(ff); cc <- complete.cases(df[, allv, drop=FALSE]); dd <- df[cc, ]
  x <- model.matrix(ff, data=dd); y <- dd[[outcome]]; w <- dd$weight
  fit <- glmnet(x, y, family="binomial", alpha=0, lambda=lambda, weights=w, intercept=TRUE, standardize=FALSE)
  cf <- as.numeric(coef(fit)); names(cf) <- rownames(coef(fit))
  list(b=cf[mediator], cp=cf[exposure], n=sum(cc))
}

run_one <- function(outcome, covs, label) {
  fm_med <- as.formula(sprintf("%s ~ %s + %s", "log_crp", "DII", paste(covs, collapse=" + ")))
  m_med <- svyglm(fm_med, design=design, family=gaussian())
  a <- coef(m_med)["DII"]
  o <- fit_glmnet_coef(outcome, "log_crp", "DII", covs)
  b <- o$b; cp <- o$cp
  indirect <- a*b; denom <- indirect + cp
  prop <- if (abs(denom)>1e-9) indirect/denom else NA
  cat(sprintf("  %s: a=%.4f b=%.4f c'=%.4f indirect=%.5f prop=%.2f%% (n=%d)\n",
              label, a, b, cp, indirect, prop*100, o$n))
  data.frame(label=label, a=a, b=b, cp=cp, indirect=indirect, prop=prop, n=o$n, stringsAsFactors=FALSE)
}

cat("===== 主分析 (corrected 5-cycle, glmnet ridge) =====\n")
res <- rbind(
  run_one("somatic_high",  COV_BASE, "somatic_noBMI"),
  run_one("somatic_high",  COV_BMI,  "somatic_BMI"),
  run_one("cognitive_high",COV_BASE, "cognitive_noBMI"),
  run_one("cognitive_high",COV_BMI,  "cognitive_BMI")
)
cat("\n=== 对照论文 ===\n")
cat("  论文 somatic: noBMI 13.6% / BMI 6.0%\n")
cat("  论文 cognitive: 1.2% (CI -0.8~3.5)\n")
cat("DONE\n")
