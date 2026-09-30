# mediation_main.R -- primary mediation analysis (DII -> hs-CRP -> depressive symptom dimensions)
#
# =============================================================================
# Methodological notes (responding to R1-6 "product-of-coefficients validity on
# the reported scale; whether strata/PSUs were incorporated")
# -----------------------------------------------------------------------------
# PRIMARY analysis uses the survey package so that the NHANES complex-sampling
# design is correctly incorporated into the standard errors:
#   path a : log_crp ~ DII + covs            (Gaussian link, svyglm)
#   path b : outcome ~ log_crp + DII + covs  (probit link, svyglm)
#   path c': the DII coefficient in the path-b model (direct effect)
#   indirect := a*b ; total := indirect + c' ; prop := indirect/total
# The design is specified explicitly as strata (SDMVSTRA), PSUs (SDMVPSU), and
# survey weights, so clustering and stratification enter the covariance matrix.
# The indirect effect, its proportion, and their 95% CIs are obtained by the
# delta method from the design-based (cluster-robust) vcov matrices.
#
# SENSITIVITY analyses:
#   (i)  lavaan WLSMV (probit) with sampling.weights -- weights observations
#        but does NOT additionally adjust for clustering; used only as a
#        robustness cross-check.
#   (ii) glmnet ridge logistic (logit scale; intercept=FALSE) -- a robustness
#        check against complete separation; note the b path is on the logit
#        scale here, so the proportion is not the primary report.
# =============================================================================
suppressPackageStartupMessages({ library(survey); library(glmnet); library(lavaan) })
options(survey.lonely.psu = "adjust", warn = -1)

df <- readRDS("/root/autodl-tmp/nhanes/nhanes_merged_analysis_v2.rds")
cat(sprintf("analysis nrow=%d\n", nrow(df)))
for (v in c("race","education","cycle")) df[[v]] <- factor(df[[v]])

design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight,
                    data = df, nest = TRUE)
cat(sprintf("design degf=%s  (PSU clusters=%d, strata=%d)\n",
    format(degf(design)),
    length(unique(interaction(df$SDMVSTRA, df$SDMVPSU))),
    length(unique(df$SDMVSTRA))))

COV_BASE <- c("RIDAGEYR","sex","race","education","smoker","drinker","cycle")
COV_BMI  <- c(COV_BASE, "BMXBMI")

# ---- delta method (multivariate; a independent of (b, c')) ----
delta_ci <- function(a, b, cp, Va, Vb, Vc, Cov_bc) {
  a <- as.numeric(a); b <- as.numeric(b); cp <- as.numeric(cp)
  Va <- as.numeric(Va); Vb <- as.numeric(Vb); Vc <- as.numeric(Vc); Cov_bc <- as.numeric(Cov_bc)
  ind <- a * b
  tot <- ind + cp
  prop <- ind / tot
  Vi <- b^2*Va + a^2*Vb
  se_i <- sqrt(max(Vi, 0))
  ga <- b*cp / tot^2
  gb <- a*cp / tot^2
  gc <- -ind / tot^2
  Vp <- ga^2*Va + gb^2*Vb + gc^2*Vc + 2*gb*gc*Cov_bc
  se_p <- sqrt(max(Vp, 0))
  c(ind = ind, ind_lo = ind - 1.96*se_i, ind_hi = ind + 1.96*se_i,
    prop = prop, prop_lo = prop - 1.96*se_p, prop_hi = prop + 1.96*se_p)
}

# -----------------------------------------------------------------------------
# PRIMARY: survey svyglm (Gaussian path a + probit path b/c')
# -----------------------------------------------------------------------------
run_svy <- function(outcome, covs, label) {
  allv <- unique(c(outcome, "log_crp", "DII", covs))
  cc <- complete.cases(df[, allv, drop = FALSE])
  dd <- df[cc, ]
  dsub <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight,
                    data = dd, nest = TRUE)
  fm_a <- as.formula(sprintf("log_crp ~ DII + %s", paste(covs, collapse = " + ")))
  m_a <- svyglm(fm_a, design = dsub, family = gaussian())
  a <- coef(m_a)["DII"]; Va <- vcov(m_a)["DII", "DII"]
  fm_b <- as.formula(sprintf("%s ~ log_crp + DII + %s", outcome, paste(covs, collapse = " + ")))
  m_b <- svyglm(fm_b, design = dsub, family = binomial(link = "probit"))
  b  <- coef(m_b)["log_crp"]; Vb  <- vcov(m_b)["log_crp", "log_crp"]
  cp <- coef(m_b)["DII"];      Vc  <- vcov(m_b)["DII", "DII"]
  Cov_bc <- vcov(m_b)["log_crp", "DII"]
  dc <- delta_ci(a, b, cp, Va, Vb, Vc, Cov_bc)
  cat(sprintf("  %-16s a=%.4f b=%.4f c'=%.4f  prop=%.2f%% [%.2f, %.2f]  indirect=%.5f [%.5f, %.5f]  (n=%d)\n",
      label, a, b, cp, dc["prop"]*100, dc["prop_lo"]*100, dc["prop_hi"]*100,
      dc["ind"], dc["ind_lo"], dc["ind_hi"], nrow(dd)))
  data.frame(label = label, method = "svyglm_probit", n = nrow(dd),
             a = a, a_se = sqrt(Va), b = b, b_se = sqrt(Vb), cp = cp, cp_se = sqrt(Vc),
             indirect = dc["ind"], ind_lo = dc["ind_lo"], ind_hi = dc["ind_hi"],
             prop = dc["prop"], prop_lo = dc["prop_lo"], prop_hi = dc["prop_hi"],
             stringsAsFactors = FALSE)
}

cat("===== PRIMARY: survey svyglm (probit; strata/PSU/weights incorporated) =====\n")
res_primary <- rbind(
  run_svy("somatic_high",   COV_BASE, "somatic_noBMI"),
  run_svy("somatic_high",   COV_BMI,  "somatic_BMI"),
  run_svy("cognitive_high", COV_BASE, "cognitive_noBMI"),
  run_svy("cognitive_high", COV_BMI,  "cognitive_BMI")
)

# -----------------------------------------------------------------------------
# SENSITIVITY (i): lavaan WLSMV (probit) + sampling.weights (no clustering)
# -----------------------------------------------------------------------------
enc <- function(var) {
  x <- as.character(df[[var]])
  lv <- sort(unique(x[!is.na(x)]))
  rest <- setdiff(lv, lv[1])
  out <- character(0)
  for (lev in rest) {
    nm <- paste0(var, "_d", match(lev, rest))
    df[[nm]] <<- as.integer(x == lev)
    out <- c(out, nm)
  }
  out
}
rc <- enc("race"); ec <- enc("education"); cc <- enc("cycle")
replace_covs <- function(covs) {
  for (pair in list(c("race", rc), c("education", ec), c("cycle", cc))) {
    if (pair[1] %in% covs) covs <- c(setdiff(covs, pair[1]), pair[2])
  }
  covs
}
COV_BASE_LAV <- replace_covs(COV_BASE)
COV_BMI_LAV  <- replace_covs(COV_BMI)

run_lavaan <- function(outcome, covs, label) {
  covstr <- paste(covs, collapse = " + ")
  model <- sprintf(
    "log_crp ~ a*DII + %s\n%s ~ b*log_crp + cp*DII + %s\nindirect := a*b\ntotal := indirect + cp\nprop := indirect/total",
    covstr, outcome, covstr)
  fit <- tryCatch(
    sem(model = model, data = df, ordered = outcome, estimator = "WLSMV",
        sampling.weights = "weight"),
    error = function(e) { cat("  lavaan ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(invisible(NULL))
  pe <- parameterEstimates(fit)
  a <- pe$est[pe$label == "a"]; b <- pe$est[pe$label == "b"]
  cp <- pe$est[pe$label == "cp"]; prop <- pe$est[pe$label == "prop"]
  ci_prop <- pe[pe$label == "prop", c("ci.lower","ci.upper")]
  cat(sprintf("  %-16s a=%.4f b=%.4f c'=%.4f  prop=%.2f%% [%.2f, %.2f]  (lavaan)\n",
      label, a, b, cp, prop*100, ci_prop[[1]]*100, ci_prop[[2]]*100))
  data.frame(label = label, method = "lavaan_wlsmv_w", a = a, b = b, cp = cp,
             prop = prop, prop_lo = ci_prop[[1]], prop_hi = ci_prop[[2]],
             stringsAsFactors = FALSE)
}

cat("\n===== SENSITIVITY (i): lavaan WLSMV + sampling weights =====\n")
res_lavaan <- rbind(
  run_lavaan("somatic_high",   COV_BASE_LAV, "somatic_noBMI"),
  run_lavaan("somatic_high",   COV_BMI_LAV,  "somatic_BMI"),
  run_lavaan("cognitive_high", COV_BASE_LAV, "cognitive_noBMI"),
  run_lavaan("cognitive_high", COV_BMI_LAV,  "cognitive_BMI")
)

# -----------------------------------------------------------------------------
# SENSITIVITY (ii): glmnet ridge logistic (logit scale; intercept=FALSE)
# -----------------------------------------------------------------------------
fit_glmnet_coef <- function(outcome, mediator, exposure, covs, lambda = 1e-4) {
  ff <- as.formula(sprintf("%s ~ %s + %s + %s", outcome, mediator, exposure, paste(covs, collapse = " + ")))
  allv <- all.vars(ff); cc <- complete.cases(df[, allv, drop = FALSE]); dd <- df[cc, ]
  x <- model.matrix(ff, data = dd); y <- dd[[outcome]]; w <- dd$weight
  fit <- glmnet(x, y, family = "binomial", alpha = 0, lambda = lambda, weights = w,
                intercept = FALSE, standardize = FALSE)
  cf <- as.numeric(coef(fit))[-1]; names(cf) <- rownames(coef(fit))[-1]
  list(b = cf[mediator], cp = cf[exposure], n = sum(cc))
}

run_sens <- function(outcome, covs, label) {
  fm_med <- as.formula(sprintf("log_crp ~ DII + %s", paste(covs, collapse = " + ")))
  m_med <- svyglm(fm_med, design = design, family = gaussian())
  a <- coef(m_med)["DII"]
  o <- fit_glmnet_coef(outcome, "log_crp", "DII", covs)
  b <- o$b; cp <- o$cp
  indirect <- a * b; denom <- indirect + cp
  prop <- if (abs(denom) > 1e-9) indirect / denom else NA
  cat(sprintf("  %-16s a=%.4f b=%.4f c'=%.4f  prop=%.2f%%  (logit sensitivity, n=%d)\n",
      label, a, b, cp, prop*100, o$n))
  data.frame(label = label, method = "glmnet_logit", a = a, b = b, cp = cp,
             indirect = indirect, prop = prop, n = o$n, stringsAsFactors = FALSE)
}

cat("\n===== SENSITIVITY (ii): glmnet ridge logistic (logit scale) =====\n")
res_sens <- rbind(
  run_sens("somatic_high",   COV_BASE, "somatic_noBMI"),
  run_sens("somatic_high",   COV_BMI,  "somatic_BMI"),
  run_sens("cognitive_high", COV_BASE, "cognitive_noBMI"),
  run_sens("cognitive_high", COV_BMI,  "cognitive_BMI")
)

dir.create("/root/autodl-tmp/results", showWarnings = FALSE)
write.csv(res_primary, "/root/autodl-tmp/results/mediation_primary_svyglm_probit.csv", row.names = FALSE)
cat("\nDONE\n")
