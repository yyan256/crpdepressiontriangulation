# 72_sensitivity_continuous.R -- R2-6 连续结局敏感性分析
# 结局：PHQ-9 连续躯体维度、连续认知维度、总 PHQ-9、probable depression (PHQ-9>=10)
# 中介：DII -> log_crp -> 结局（连续结局用 svyglm gaussian，二分结局用 glmnet 岭 logistic）
suppressPackageStartupMessages({ library(survey); library(glmnet); library(dplyr) })
options(survey.lonely.psu = "adjust", warn = -1)

df <- readRDS("/root/autodl-tmp/nhanes/nhanes_merged_analysis_v2.rds")
cat(sprintf("analysis nrow=%d\n", nrow(df)))
for (v in c("race","education","cycle")) df[[v]] <- factor(df[[v]])
design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight, data = df, nest = TRUE)

COV_BASE <- c("RIDAGEYR","sex","race","education","smoker","drinker","cycle")
COV_BMI  <- c(COV_BASE, "BMXBMI")

# ---- 连续结局的 svyglm gaussian 系数提取 ----
fit_svy_coef <- function(outcome, mediator, exposure, covs) {
  ff <- as.formula(sprintf("%s ~ %s + %s + %s", outcome, mediator, exposure, paste(covs, collapse=" + ")))
  m <- svyglm(ff, design=design, family=gaussian())
  cc <- coef(m)
  n <- nrow(df[complete.cases(df[, all.vars(ff), drop=FALSE]), ])
  list(b=as.numeric(cc[mediator]), cp=as.numeric(cc[exposure]), n=n)
}

# ---- 二分结局的 glmnet 岭 logistic 系数提取 ----
fit_glmnet_coef <- function(outcome, mediator, exposure, covs, lambda=1e-4) {
  ff <- as.formula(sprintf("%s ~ %s + %s + %s", outcome, mediator, exposure, paste(covs, collapse=" + ")))
  allv <- all.vars(ff); cc <- complete.cases(df[, allv, drop=FALSE]); dd <- df[cc, ]
  x <- model.matrix(ff, data=dd); y <- dd[[outcome]]; w <- dd$weight
  fit <- glmnet(x, y, family="binomial", alpha=0, lambda=lambda, weights=w, intercept=TRUE, standardize=FALSE)
  cf <- as.numeric(coef(fit)); names(cf) <- rownames(coef(fit))
  list(b=cf[mediator], cp=cf[exposure], n=sum(cc))
}

# ---- a 路径（暴露->中介，gaussian，与结局无关）----
fit_a <- function(covs) {
  ff <- as.formula(sprintf("log_crp ~ DII + %s", paste(covs, collapse=" + ")))
  m <- svyglm(ff, design=design, family=gaussian())
  as.numeric(coef(m)["DII"])
}

run_cont <- function(outcome, covs, label, binary=FALSE) {
  a <- fit_a(covs)
  if (binary) {
    o <- fit_glmnet_coef(outcome, "log_crp", "DII", covs)
  } else {
    o <- fit_svy_coef(outcome, "log_crp", "DII", covs)
  }
  b <- o$b; cp <- o$cp
  indirect <- a*b; denom <- indirect + cp
  prop <- if (abs(denom)>1e-9) indirect/denom else NA
  cat(sprintf("  %s: a=%.4f b=%.4f c'=%.4f indirect=%.5f prop=%.2f%% (n=%d)\n",
              label, a, b, cp, indirect, prop*100, o$n))
  data.frame(label=label, a=a, b=b, cp=cp, indirect=indirect, prop=prop, n=o$n, stringsAsFactors=FALSE)
}

cat("===== R2-6 敏感性分析：连续结局 =====\n")
cat("--- (1) 连续躯体维度 (continuous somatic) ---\n")
r1 <- rbind(
  run_cont("phq9_somatic", COV_BASE, "somatic_cont_noBMI"),
  run_cont("phq9_somatic", COV_BMI,  "somatic_cont_BMI")
)
cat("--- (2) 连续认知维度 (continuous cognitive) ---\n")
r2 <- rbind(
  run_cont("phq9_cognitive", COV_BASE, "cognitive_cont_noBMI"),
  run_cont("phq9_cognitive", COV_BMI,  "cognitive_cont_BMI")
)
cat("--- (3) 总 PHQ-9 (total PHQ-9) ---\n")
r3 <- rbind(
  run_cont("phq9_total", COV_BASE, "total_noBMI"),
  run_cont("phq9_total", COV_BMI,  "total_BMI")
)
cat("--- (4) probable depression PHQ-9>=10 (binary) ---\n")
r4 <- rbind(
  run_cont("phq9_ge10", COV_BASE, "ge10_noBMI", binary=TRUE),
  run_cont("phq9_ge10", COV_BMI,  "ge10_BMI",  binary=TRUE)
)
cat("\nDONE\n")
