# 68_delta_ci_v2.R -- delta-method 95% CI using glmnet ridge + sandwich vcov
suppressPackageStartupMessages({ library(survey); library(glmnet); library(dplyr) })
options(survey.lonely.psu = "adjust", warn = -1)

df <- readRDS("/root/autodl-tmp/nhanes/nhanes_merged_analysis_v2.rds")
for (v in c("race","education","cycle")) df[[v]] <- factor(df[[v]])

COV_BASE <- c("RIDAGEYR","sex","race","education","smoker","drinker","cycle")
COV_BMI  <- c(COV_BASE, "BMXBMI")
LAMBDA <- 1e-4

# weighted sandwich vcov for ridge-penalized logistic (glmnet)
# V = H^{-1} B H^{-1}, H = X'WX + 2*lambda*I, B = X' diag(w^2 * mu) X
sandwich_vcov <- function(x, y, cf, lambda, w) {
  eta <- as.numeric(x %*% cf)
  p <- 1/(1 + exp(-eta))
  mu <- pmax(p * (1 - p), 1e-8)
  Wv <- mu * w
  H <- crossprod(x, x * Wv) + 2*lambda*diag(ncol(x))
  Wv2 <- mu * w^2
  B <- crossprod(x, x * Wv2)
  Hi <- solve(H)
  V <- Hi %*% B %*% Hi
  colnames(V) <- rownames(V) <- colnames(x)
  V
}

# multivariate delta method for prop = a*b/(a*b+cp)
delta_ci <- function(a, b, cp, Va, Vb, Vc, Cov_bc) {
  # strip names so c(prop=...) does not become "prop.DII"
  a <- as.numeric(a); b <- as.numeric(b); cp <- as.numeric(cp)
  Va <- as.numeric(Va); Vb <- as.numeric(Vb); Vc <- as.numeric(Vc); Cov_bc <- as.numeric(Cov_bc)
  ind <- a*b
  tot <- ind + cp
  prop <- ind/tot
  ga <- b*cp/tot^2
  gb <- a*cp/tot^2
  gc <- -ind/tot^2
  Vp <- ga^2*Va + gb^2*Vb + gc^2*Vc + 2*gb*gc*Cov_bc
  se <- sqrt(max(Vp, 0))
  Vi <- b^2*Va + a^2*Vb
  se_i <- sqrt(max(Vi, 0))
  c(prop=prop, prop_lo=prop-1.96*se, prop_hi=prop+1.96*se,
    indirect=ind, ind_lo=ind-1.96*se_i, ind_hi=ind+1.96*se_i)
}

run_one <- function(outcome, covs, label) {
  allv <- unique(c(outcome, "log_crp", "DII", covs))
  cc <- complete.cases(df[, allv, drop=FALSE])
  dd <- df[cc, ]
  # a path: svyglm gaussian (survey-corrected SE)
  dsub <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~weight, data=dd, nest=TRUE)
  fm_a <- as.formula(sprintf("log_crp ~ DII + %s", paste(covs, collapse=" + ")))
  m_a <- svyglm(fm_a, design=dsub, family=gaussian())
  a <- coef(m_a)["DII"]; Va <- vcov(m_a)["DII","DII"]
  # b/c' path: glmnet ridge (no separation) + sandwich vcov
  fm_b <- as.formula(sprintf("%s ~ log_crp + DII + %s", outcome, paste(covs, collapse=" + ")))
  x <- model.matrix(fm_b, data=dd)
  y <- dd[[outcome]]; w <- dd$weight
  gr <- glmnet(x, y, family="binomial", alpha=0, lambda=LAMBDA, weights=w, intercept=FALSE, standardize=FALSE)
  cfm <- coef(gr)
  # glmnet emits a zero pseudo-intercept row even with intercept=FALSE -> drop it
  cf <- as.numeric(cfm)[-1]; names(cf) <- rownames(cfm)[-1]
  b <- cf["log_crp"]; cp <- cf["DII"]
  V <- sandwich_vcov(x, y, cf, LAMBDA, w)
  Vb <- V["log_crp","log_crp"]; Vc <- V["DII","DII"]; Cov_bc <- V["log_crp","DII"]
  dc <- delta_ci(a, b, cp, Va, Vb, Vc, Cov_bc)
  cat(sprintf("== %s (n=%d) ==\n", label, nrow(dd)))
  cat(sprintf("  a=%.4f (se=%.4f)  b=%.4f (se=%.4f)  c'=%.4f (se=%.4f)\n",
      a, sqrt(Va), b, sqrt(Vb), cp, sqrt(Vc)))
  cat(sprintf("  prop=%.2f%% [%.2f, %.2f]  indirect=%.5f [%.5f, %.5f]\n",
      dc["prop"]*100, dc["prop_lo"]*100, dc["prop_hi"]*100,
      dc["indirect"], dc["ind_lo"], dc["ind_hi"]))
  cat("\n")
}

cat("===== delta-method 95% CI (glmnet ridge + sandwich vcov, 5-cycle n=24019) =====\n")
run_one("somatic_high",  COV_BASE, "somatic_noBMI")
run_one("somatic_high",  COV_BMI,  "somatic_BMI")
run_one("cognitive_high",COV_BASE, "cognitive_noBMI")
run_one("cognitive_high",COV_BMI,  "cognitive_BMI")
cat("DONE\n")
