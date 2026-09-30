# 53_med_v2.R -- mediation on corrected 5-cycle sample (n=24019)
#
# =============================================================================
# 方法学说明（回应 R1-6 "product-of-coefficients 是否在报告尺度上有效"）
# -----------------------------------------------------------------------------
# 主结果统一到 **probit 链接**：用 lavaan WLSMV + sampling weights 同时估计
#     log_crp ~ a*DII + covs            （path a：线性，log-CRP 连续中介）
#     outcome ~ b*log_crp + cp*DII + covs（path b / c'：probit 链接）
#     indirect := a*b ; total := indirect + cp ; prop := indirect/total
#   a 与 b 在同一模型内、同一数据、同一协变量集、同一 survey 权重下联合估计，
#   链接函数与尺度明确且一致（log-CRP 为高斯尺度、结局为 probit 标准正态尺度），
#   间接效应与占比由 delta 法（lavaan 内置）直接给出 CI —— 不再出现
#   "线性系数 × logit 系数" 的尺度混用问题。
#
# 敏感性分析：glmnet ridge logistic（解决完全分离）。注意这里 b 路径是
#   logit 尺度，与主结果的 probit 尺度不同，仅作为稳健性交叉验证，
#   不用于主报告占比。
# =============================================================================
suppressPackageStartupMessages({ library(survey); library(glmnet); library(dplyr); library(lavaan) })
options(survey.lonely.psu = "adjust", warn = -1)

df <- readRDS("/root/autodl-tmp/nhanes/nhanes_merged_analysis_v2.rds")
cat(sprintf("analysis nrow=%d\n", nrow(df)))
for (v in c("race","education","cycle")) df[[v]] <- factor(df[[v]])

COV_BASE <- c("RIDAGEYR","sex","race","education","smoker","drinker","cycle")
COV_BMI  <- c(COV_BASE, "BMXBMI")

# -----------------------------------------------------------------------------
# 主结果：lavaan WLSMV（probit），dummy 编码无序多分类协变量
# -----------------------------------------------------------------------------
# dummy-encode multi-level unordered factors (WLSMV cannot take factors directly)
enc <- function(var) {
  x <- as.character(df[[var]])
  lv <- sort(unique(x[!is.na(x)]))
  rest <- setdiff(lv, lv[1])            # first level = reference
  out <- character(0)
  for (lev in rest) {
    nm <- paste0(var, "_d", match(lev, rest))
    df[[nm]] <<- as.integer(x == lev)
    out <- c(out, nm)
  }
  out
}
rc <- enc("race"); ec <- enc("education"); cc <- enc("cycle")
# 将 covs 中的 factor 列替换为 dummy 列（sex 是 2 水平可保留，但统一转 numeric 更稳妥）
dummy_map <- c(race = list(rc), education = list(ec), cycle = list(cc))
replace_covs <- function(covs) {
  out <- covs
  for (nm in names(dummy_map)) {
    if (nm %in% out) out <- c(setdiff(out, nm), dummy_map[[nm]])
  }
  out
}
COV_BASE_LAV <- replace_covs(COV_BASE)
COV_BMI_LAV  <- replace_covs(COV_BMI)

run_lavaan <- function(outcome, covs, label, use_w = FALSE) {
  covstr <- paste(covs, collapse = " + ")
  model <- sprintf(
    "log_crp ~ a*DII + %s\n%s ~ b*log_crp + cp*DII + %s\nindirect := a*b\ntotal := indirect + cp\nprop := indirect/total",
    covstr, outcome, covstr)
  args <- list(model = model, data = df, ordered = outcome, estimator = "WLSMV")
  # 主结果使用 NHANES 复杂抽样的 sampling weights（与 Table 1 / svyglm 一致）
  if (use_w) args$sampling.weights <- "weight"
  fit <- tryCatch(
    do.call(sem, args),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(invisible(NULL))
  pe <- parameterEstimates(fit)
  a <- pe$est[pe$label == "a"];  b  <- pe$est[pe$label == "b"]
  cp <- pe$est[pe$label == "cp"]; ind <- pe$est[pe$label == "indirect"]
  prop <- pe$est[pe$label == "prop"]
  ci_ind <- pe[pe$label == "indirect", c("ci.lower","ci.upper")]
  ci_prop <- pe[pe$label == "prop", c("ci.lower","ci.upper")]
  cat(sprintf("  %-14s %s: a=%.4f b=%.4f cp=%.4f indirect=%.5f [%.5f, %.5f] prop=%.2f%% [%.2f, %.2f]\n",
              label, if (use_w) "WLSMV+w" else "WLSMV", a, b, cp, ind,
              ci_ind[[1]], ci_ind[[2]], prop*100, ci_prop[[1]]*100, ci_prop[[2]]*100))
  data.frame(label=label, weighted=use_w, a=a, b=b, cp=cp, indirect=ind,
             ind_lo=ci_ind[[1]], ind_hi=ci_ind[[2]],
             prop=prop, prop_lo=ci_prop[[1]], prop_hi=ci_prop[[2]],
             stringsAsFactors=FALSE)
}

cat("===== 主结果：lavaan WLSMV (probit, sampling weights) =====\n")
res_main <- rbind(
  run_lavaan("somatic_high",  COV_BASE_LAV, "somatic_noBMI",  TRUE),
  run_lavaan("somatic_high",  COV_BMI_LAV,  "somatic_BMI",    TRUE),
  run_lavaan("cognitive_high",COV_BASE_LAV, "cognitive_noBMI", TRUE),
  run_lavaan("cognitive_high",COV_BMI_LAV,  "cognitive_BMI",   TRUE)
)
cat("\n===== 稳健性：lavaan WLSMV (probit, 未加权) =====\n")
res_main_unw <- rbind(
  run_lavaan("somatic_high",  COV_BASE_LAV, "somatic_noBMI"),
  run_lavaan("somatic_high",  COV_BMI_LAV,  "somatic_BMI"),
  run_lavaan("cognitive_high",COV_BASE_LAV, "cognitive_noBMI"),
  run_lavaan("cognitive_high",COV_BMI_LAV,  "cognitive_BMI")
)

# -----------------------------------------------------------------------------
# 敏感性：glmnet ridge logistic（logit 尺度；intercept=FALSE 避免双重截距）
# -----------------------------------------------------------------------------
fit_glmnet_coef <- function(outcome, mediator, exposure, covs, lambda=1e-4) {
  ff <- as.formula(sprintf("%s ~ %s + %s + %s", outcome, mediator, exposure, paste(covs, collapse=" + ")))
  allv <- all.vars(ff); cc <- complete.cases(df[, allv, drop=FALSE]); dd <- df[cc, ]
  x <- model.matrix(ff, data=dd); y <- dd[[outcome]]; w <- dd$weight
  # intercept=FALSE：model.matrix 已含 intercept 列，避免双重截距（历史 bug）
  fit <- glmnet(x, y, family="binomial", alpha=0, lambda=lambda, weights=w, intercept=FALSE, standardize=FALSE)
  cf <- as.numeric(coef(fit))[-1]; names(cf) <- rownames(coef(fit))[-1]  # 去掉 glmnet 伪 intercept 行
  list(b=cf[mediator], cp=cf[exposure], n=sum(cc))
}

run_sens <- function(outcome, covs, label) {
  fm_med <- as.formula(sprintf("%s ~ %s + %s", "log_crp", "DII", paste(covs, collapse=" + ")))
  m_med <- svyglm(fm_med, design=svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~weight, data=df, nest=TRUE), family=gaussian())
  a <- coef(m_med)["DII"]
  o <- fit_glmnet_coef(outcome, "log_crp", "DII", covs)
  b <- o$b; cp <- o$cp
  indirect <- a*b; denom <- indirect + cp
  prop <- if (abs(denom)>1e-9) indirect/denom else NA
  cat(sprintf("  %-14s (logit敏感): a=%.4f b=%.4f c'=%.4f indirect=%.5f prop=%.2f%% (n=%d)\n",
              label, a, b, cp, indirect, prop*100, o$n))
  data.frame(label=paste0(label,"_sens"), a=a, b=b, cp=cp, indirect=indirect, prop=prop, n=o$n, stringsAsFactors=FALSE)
}

cat("\n===== 敏感性：glmnet ridge (logit 尺度，仅交叉验证) =====\n")
res_sens <- rbind(
  run_sens("somatic_high",  COV_BASE, "somatic_noBMI"),
  run_sens("somatic_high",  COV_BMI,  "somatic_BMI"),
  run_sens("cognitive_high",COV_BASE, "cognitive_noBMI"),
  run_sens("cognitive_high",COV_BMI,  "cognitive_BMI")
)

cat("\n=== 对照论文 ===\n")
cat("  论文 somatic: noBMI 13.6% / BMI 6.0%\n")
cat("  论文 cognitive: 1.2% (CI -0.8~3.5)\n")
cat("DONE\n")
