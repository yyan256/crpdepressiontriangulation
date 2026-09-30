# 57_lavaan_repro2.R -- lavaan WLSMV reproduction (dummy-encoded covariates)
# WLSMV = probit link, the only estimator lavaan supports for binary outcomes.
suppressPackageStartupMessages({ library(lavaan) })
df <- readRDS("/root/autodl-tmp/nhanes/nhanes_merged_analysis_v2.rds")

# --- inspect coding ---
cat("race:");      print(table(df$race, useNA = "ifany"))
cat("education:"); print(table(df$education, useNA = "ifany"))
cat("cycle:");     print(table(df$cycle, useNA = "ifany"))
cat("sex:");       print(table(df$sex, useNA = "ifany"))
cat("smoker:");    print(table(df$smoker, useNA = "ifany"))
cat("drinker:");   print(table(df$drinker, useNA = "ifany"))

# --- dummy-encode multi-level unordered factors (WLSMV cannot take them directly) ---
enc <- function(var) {
  x <- as.character(df[[var]])
  lv <- sort(unique(x[!is.na(x)]))
  rest <- setdiff(lv, lv[1])           # first level = reference
  out <- character(0)
  for (lev in rest) {
    nm <- paste0(var, "_d", match(lev, rest))
    df[[nm]] <<- as.integer(x == lev)
    out <- c(out, nm)
  }
  out
}
rc <- enc("race"); ec <- enc("education"); cc <- enc("cycle")

COV_BASE <- c("RIDAGEYR","sex","smoker","drinker", rc, ec, cc)
COV_BMI  <- c(COV_BASE, "BMXBMI")

run_one <- function(outcome, covs, label, use_w = FALSE) {
  covstr <- paste(covs, collapse = " + ")
  model <- sprintf("log_crp ~ a*DII + %s\n%s ~ b*log_crp + cp*DII + %s\nindirect := a*b\ntotal := indirect + cp\nprop := indirect/total",
                   covstr, outcome, covstr)
  args <- list(model = model, data = df, ordered = outcome, estimator = "WLSMV")
  if (use_w) args$sampling.weights <- "weight"
  t0 <- Sys.time()
  fit <- tryCatch(do.call(sem, args),
                  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(fit)) return(invisible(NULL))
  pe <- parameterEstimates(fit)
  a   <- pe$est[pe$label == "a"]
  b   <- pe$est[pe$label == "b"]
  cp  <- pe$est[pe$label == "cp"]
  ind <- pe$est[pe$label == "indirect"]
  prop <- pe$est[pe$label == "prop"]
  cat(sprintf("  %-14s %s: a=%.4f b=%.4f cp=%.4f ind=%.5f prop=%.2f%% (%.1fs)\n",
              label, if (use_w) "WLSMV+w" else "WLSMV", a, b, cp, ind, prop*100, dt))
  flush.console()
  invisible(NULL)
}

cat("\n===== WLSMV (probit, no weights) =====\n")
run_one("somatic_high",  COV_BASE, "somatic_noBMI")
run_one("somatic_high",  COV_BMI,  "somatic_BMI")
run_one("cognitive_high",COV_BASE, "cognitive_noBMI")
run_one("cognitive_high",COV_BMI,  "cognitive_BMI")

cat("\n===== WLSMV + sampling.weights =====\n")
run_one("somatic_high",  COV_BASE, "somatic_noBMI",  TRUE)
run_one("somatic_high",  COV_BMI,  "somatic_BMI",    TRUE)
run_one("cognitive_high",COV_BASE, "cognitive_noBMI",TRUE)
run_one("cognitive_high",COV_BMI,  "cognitive_BMI",  TRUE)

cat("\npaper: somatic 13.6% (noBMI) / 6.0% (BMI); cognitive 1.2% (CI -0.8~3.5)\n")
cat("glmnet ridge (de novo): somatic 14.09%/6.46%; cognitive 11.64%/6.16%\n")
cat("DONE\n")
