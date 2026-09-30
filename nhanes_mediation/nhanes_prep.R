# =============================================================================
# 20_nhanes_prep.R — Phase 3 数据准备（de novo）
# 论文: Dietary Anti-Inflammatory Benefits Bypass CRP→Depression Pathway
# 目标: 从 NHANES 原始数据重建分析样本，计算 DII / PHQ-9 / hs-CRP /
#       BMI / 协变量 / survey 权重，保存为 rds，并输出各阶段排除人数（回应 R1-19 流程）。
# 依赖: nhanesA, data.table, dplyr, tidyr
# 环境: AutoDL Linux（已装 nhanesA 1.4.1 / data.table / dplyr）
# 用法: Rscript 20_nhanes_prep.R
# 预计: 1–2 小时（5 周期 × 多文件下载）
# 重要：NHANES 在 2011-2012(G)、2013-2014(H) 未测 CRP（官方数据缺口），
#       故仅纳入 CRP 可用的 5 个周期（与同类文献一致，见 README_Phase3.md）。
# =============================================================================

suppressPackageStartupMessages({
  library(nhanesA)
  library(data.table)
  library(dplyr)
})

OUT <- "/root/autodl-tmp/nhanes"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
                         sprintf(...), "\n")

# -----------------------------------------------------------------------------
# 0. 5 个 CRP 可用周期 + 文件名前缀
#    NHANES 未在 2011-2012(G)、2013-2014(H) 测定 CRP，故排除这两个周期
# -----------------------------------------------------------------------------
cycles <- data.frame(
  cyc = c("2005-2006","2007-2008","2009-2010","2015-2016","2017-2018"),
  pfx = c("D","E","F","I","J"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 1. DII 17 参数参考表（Shivappa et al. 2014, Public Health Nutrition）
#    列: 参数名 / DR1TOT 变量名 / 炎症效应分数 / 全球均值 / 全球 SD
# -----------------------------------------------------------------------------
dii_ref <- data.frame(
  param  = c("energy","protein","total_fat","sat_fat","trans_fat","cholesterol",
             "carbohydrate","fiber","vit_b6","vit_b12","vit_c","vit_d","vit_e",
             "iron","magnesium","zinc","alcohol"),
  var    = c("DR1TKCAL","DR1TPROT","DR1TTFAT","DR1TSFAT","DR1TTRF","DR1TCHOL",
             "DR1TCARB","DR1TFIBE","DR1TVB6","DR1TVB12","DR1TVC","DR1TVD","DR1TATOC",
             "DR1TIRON","DR1TMAGN","DR1TZINC","DR1TALCO"),
  effect = c( 0.180, 0.021, 0.298, 0.373, 0.229, 0.110,
              0.097,-0.663,-0.365, 0.106,-0.424,-0.446,-0.419,
              0.032,-0.484,-0.313,-0.278),
  mean   = c(2056, 79.4, 71.4, 28.6, 3.15, 279.4,
             272.2, 18.8, 1.47, 5.15, 118.2, 6.26, 8.73,
             13.35, 310.1, 9.84, 13.98),
  sd     = c(338, 13.9, 19.4, 8.0, 3.75, 51.2,
             40.0, 4.9, 0.74, 2.70, 43.46, 2.21, 1.49,
             3.71, 139.4, 2.19, 3.72),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 2. 计算 DII 分数（对一行营养素向量）
#    算法（Shivappa 2014）：
#      z = (摄入 - 全球均值) / 全球SD
#      percentile = pnorm(z)
#      centered = 2 * percentile - 1          # 范围 [-1, 1]
#      score = centered * 炎症效应分数
#      DII = sum(score)  （对可用的参数求和；缺参数自动跳过）
# -----------------------------------------------------------------------------
calc_dii_vec <- function(df) {
  # 向量化计算：对整列逐参数累加（缺参数自动跳过）
  dii <- rep(0, nrow(df))
  for (i in seq_len(nrow(dii_ref))) {
    v <- dii_ref$var[i]
    if (!v %in% names(df)) next
    x <- as.numeric(df[[v]])
    z  <- (x - dii_ref$mean[i]) / dii_ref$sd[i]
    cc <- 2 * pnorm(z) - 1
    dii <- dii + cc * dii_ref$effect[i]
  }
  dii
}

# -----------------------------------------------------------------------------
# 3. 稳健下载（nhanes 失败返回 NULL，不中断整个循环）
# -----------------------------------------------------------------------------
safe_nhanes <- function(tbl, max_try = 5, sleep_s = 3) {
  # translated = FALSE：返回原始整数编码（而非 "Male"/"Female" 等文字），
  # 与下游 RIAGENDR==1、RIDRETH1==3 等数值比较严格对齐
  # 关键：NHANES CDC 服务器对单个表偶发 404/超时/空返回（已实测 ALQ_E/ALQ_F 间歇
  #      "not available" 或返回 0 行），必须重试并把 0 行判为失败重试，否则整周期协变量丢失。
  out <- NULL
  for (k in seq_len(max_try)) {
    out <- tryCatch(as.data.frame(nhanes(tbl, translated = FALSE)),
                    error = function(e) NULL)
    if (!is.null(out) && nrow(out) > 0) return(out)
    Sys.sleep(sleep_s)
  }
  NULL
}

# 稳健转数值：无论返回 numeric / character / factor，都取底层编码
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# -----------------------------------------------------------------------------
# 4. 主循环：下载 + 构建每周期数据
# -----------------------------------------------------------------------------
all_rows <- vector("list", nrow(cycles))

for (i in seq_len(nrow(cycles))) {
  pfx <- cycles$pfx[i]; cyc <- cycles$cyc[i]
  log("===== 周期 %s (%s) =====", cyc, pfx)

  # ---- 4.1 饮食 Day1 总量（DII 17 参数来源）----
  dr1 <- safe_nhanes(paste0("DR1TOT_", pfx))

  # ---- 4.2 抑郁 PHQ-9 ----
  dpq <- safe_nhanes(paste0("DPQ_", pfx))

  # ---- 4.3 CRP：表名随周期变化，靠变量名判断单位 ----
  #   2005-2010 = CRP_<p>（LBXCRP, mg/dL）
  #   2015-2018 = HSCRP_<p>（LBXHSCRP, mg/L）
  crp_candidates <- if (pfx %in% c("D","E","F")) {
    paste0("CRP_", pfx)
  } else {
    c(paste0("HSCRP_", pfx), paste0("CRP_", pfx))
  }
  crp <- NULL; crp_unit <- NULL
  for (ct in crp_candidates) {
    tmp <- safe_nhanes(ct)
    if (is.null(tmp)) next
    if ("LBXCRP" %in% names(tmp)) {
      crp <- tmp[, c("SEQN","LBXCRP")]; crp_unit <- "mg/dL"; break
    }
    if ("LBXHSCRP" %in% names(tmp)) {
      crp <- tmp[, c("SEQN","LBXHSCRP")]; crp_unit <- "mg/L"; break
    }
  }
  if (is.null(crp)) log("  [警告] %s 未找到 CRP 表，该周期 CRP 将全缺失", cyc)

  # ---- 4.4 BMI ----
  bmx <- safe_nhanes(paste0("BMX_", pfx))

  # ---- 4.5 人口学 + survey 变量 ----
  demo <- safe_nhanes(paste0("DEMO_", pfx))

  # ---- 4.6 吸烟 ----
  smq <- safe_nhanes(paste0("SMQ_", pfx))

  # ---- 4.7 饮酒 ----
  alq <- safe_nhanes(paste0("ALQ_", pfx))

  # ---- 4.8 体力活动（扩充协变量，可选）----
  paq <- safe_nhanes(paste0("PAQ_", pfx))

  # ---- 4.9 睡眠（扩充协变量，可选）----
  slq <- safe_nhanes(paste0("SLQ_", pfx))

  if (is.null(dr1) || is.null(dpq) || is.null(demo)) {
    log("  [跳过] %s 关键文件缺失", cyc); next
  }

  # ---- 提取并重命名变量（统一内部名）----
  dat <- data.frame(SEQN = dr1$SEQN, stringsAsFactors = FALSE)

  # DII 营养素：逐列提取，缺失则不纳入
  for (k in seq_len(nrow(dii_ref))) {
    v <- dii_ref$var[k]
    if (v %in% names(dr1)) dat[[v]] <- to_num(dr1[[v]])
  }

  # PHQ-9 9 项（DPQ010..DPQ090，取值 0-3；7=拒绝/9=不知 → NA）
  # 关键：DPQ 行数少于 DR1TOT，必须先按 SEQN 合并再取值，不能按行直接赋值
  phq_vars <- sprintf("DPQ%03d", 1:9 * 10)
  phq_vars <- intersect(phq_vars, names(dpq))
  if (length(phq_vars) > 0) {
    dat <- merge(dat, dpq[, c("SEQN", phq_vars)], by = "SEQN", all.x = TRUE)
  }
  for (q in 1:9) {
    vn <- sprintf("DPQ%03d", q * 10)
    nm <- paste0("phq", q)
    if (vn %in% names(dat)) {
      x <- to_num(dat[[vn]])
      x[x %in% c(7, 9)] <- NA
      dat[[nm]] <- x
    } else {
      dat[[nm]] <- NA_real_   # 保险：缺项补 NA，避免后续 phq9_total 求和报错
    }
  }

  # hs-CRP（统一到 mg/L：LBXCRP 是 mg/dL → ×10；LBXHSCRP 已是 mg/L → 不变）
  if (!is.null(crp)) {
    if (crp_unit == "mg/dL") {
      crp$crp_mgL <- to_num(crp[[setdiff(names(crp), "SEQN")[1]]]) * 10
    } else {
      crp$crp_mgL <- to_num(crp[[setdiff(names(crp), "SEQN")[1]]])
    }
    dat <- merge(dat, crp[, c("SEQN","crp_mgL")], by = "SEQN", all.x = TRUE)
  } else {
    dat$crp_mgL <- NA_real_
  }

  # BMI
  if (!is.null(bmx) && "BMXBMI" %in% names(bmx)) {
    dat <- merge(dat, bmx[, c("SEQN","BMXBMI")], by = "SEQN", all.x = TRUE)
  } else {
    dat$BMXBMI <- NA_real_
  }

  # 人口学 + survey 变量（用 WTMEC2YR 体检权重，因 hs-CRP/BMI 为 MEC 测量）
  demo_vars <- intersect(c("RIAGENDR","RIDAGEYR","RIDRETH1","DMDEDUC2","INDFMPIR",
                           "SDMVPSU","SDMVSTRA","WTMEC2YR"), names(demo))
  if (length(demo_vars) > 0) {
    dat <- merge(dat, demo[, c("SEQN", demo_vars)], by = "SEQN", all.x = TRUE)
    for (v in demo_vars) dat[[v]] <- to_num(dat[[v]])
  }

  # 吸烟（当前吸烟者）—— 正确构造，修复历史编码 bug：
  #   SMQ020 = "一生是否吸过 ≥100 支烟" (1=是, 2=否, 7/9=拒绝/不知)
  #   SMQ040 = "现在是否吸烟" (1=每天, 2=有些天, 3=不吸烟, 7/9=拒绝/不知; 仅 SMQ020==1 者填答)
  #   错误旧逻辑曾把 SMQ020 误当"当前吸烟频率"，导致 smoker=1 占 ~99.97%。
  #   正确构造：SMQ020==2(从未吸过100支) → 非当前吸烟者=0；
  #             SMQ020==1 且 SMQ040∈{1,2}(每天/有些天) → 当前吸烟者=1；
  #             SMQ020==1 且 SMQ040==3(已戒烟) → 非当前吸烟者=0；其余 → NA。
  if (!is.null(smq) && "SMQ020" %in% names(smq)) {
    smq_cols <- intersect(c("SEQN","SMQ020","SMQ040"), names(smq))
    dat <- merge(dat, smq[, smq_cols, drop = FALSE], by = "SEQN", all.x = TRUE)
    if ("SMQ020" %in% names(dat)) dat$SMQ020 <- to_num(dat$SMQ020)
    if ("SMQ040" %in% names(dat)) dat$SMQ040 <- to_num(dat$SMQ040)
  }

  # 饮酒：合并终生/过去12月酒精变量，跨周期一致构造 drinker（见下方派生）
  #   D/E/F/I: ALQ101=过去12月≥12饮(1=是,2=否); ALQ110=终生≥12饮(1=是,2=否, 仅 ALQ101=2 时填答)
  #   J      : ALQ111=终生≥12饮(1=是,2=否)  [2017-2018 问卷重构, 无 ALQ101/ALQ110]
  # 注意：ALQ110 存在严重跳问缺失（仅 ALQ101=2 时填答，整体 ~71% NA），不能直接当 drinker；
  #       正确做法是先用 ALQ101(过去12月, ~91% 完整) 作主体，再用 ALQ110 把"终生曾饮但过去12月未饮"
  #       的人补回为饮酒=1；J 周期用 ALQ111(终生筛选器, ~93% 完整)。统一为"终生是否饮酒"。
  alq_vars <- intersect(c("ALQ101","ALQ110","ALQ111"), names(alq))
  if (length(alq_vars) > 0) {
    dat <- merge(dat, alq[, c("SEQN", alq_vars)], by = "SEQN", all.x = TRUE)
    for (v in alq_vars) dat[[v]] <- to_num(dat[[v]])
  }

  # 体力活动 PAQ650 (vigorous work) / PAQ665 (moderate work) — 作为扩充协变量
  if (!is.null(paq)) {
    pa_vars <- intersect(c("PAQ650","PAQ665","PAQ610","PAQ635"), names(paq))
    if (length(pa_vars) > 0) dat <- merge(dat, paq[, c("SEQN", pa_vars)], by = "SEQN", all.x = TRUE)
  }

  # 睡眠 SLD010H（每晚睡眠小时）
  if (!is.null(slq) && "SLD010H" %in% names(slq)) {
    dat <- merge(dat, slq[, c("SEQN","SLD010H")], by = "SEQN", all.x = TRUE)
    dat$SLD010H <- to_num(dat$SLD010H)
  }

  dat$cycle <- cyc
  all_rows[[i]] <- dat
  log("  %s 合并完成: n=%d", cyc, nrow(dat))
}

# -----------------------------------------------------------------------------
# 5. 合并 5 周期 + 派生变量
# -----------------------------------------------------------------------------
df <- bind_rows(all_rows)
log("合并后总行数: %d", nrow(df))

# 记录各阶段排除人数（回应 R1-19 流程表）
flow <- data.frame(step = character(), n = integer(), stringsAsFactors = FALSE)
flow <- rbind(flow, data.frame(step = "合并后原始行数", n = nrow(df)))

# 成人 ≥18（排除 NA 年龄）
df <- df %>% filter(!is.na(RIDAGEYR), RIDAGEYR >= 18)
flow <- rbind(flow, data.frame(step = "年龄>=18", n = nrow(df)))

# 计算 DII（向量化）
df$DII <- calc_dii_vec(df)
flow <- rbind(flow, data.frame(step = "DII 非缺失", n = sum(!is.na(df$DII))))

# PHQ-9 维度（连续）+ 二分结局 + 总 PHQ-9
df <- df %>% mutate(
  phq9_somatic  = phq3 + phq4 + phq5 + phq8,                 # 睡眠+疲劳+食欲+精神运动
  phq9_cognitive = phq1 + phq2 + phq6 + phq7 + phq9,         # 快感缺失+低落+无价值+注意+自杀
  phq9_total    = phq1 + phq2 + phq3 + phq4 + phq5 + phq6 + phq7 + phq8 + phq9,
  somatic_high  = ifelse(phq9_somatic  >= 5, 1, 0),
  cognitive_high = ifelse(phq9_cognitive >= 5, 1, 0),
  phq9_ge10     = ifelse(phq9_total >= 10, 1, 0)
)

# log 变换 hs-CRP（用 base::log，避免与上方自定义的 log 打印函数重名冲突）
df$log_crp <- base::log(df$crp_mgL + 0.01)

# 协变量编码
df <- df %>% mutate(
  sex   = factor(ifelse(RIAGENDR == 1, "Male", ifelse(RIAGENDR == 2, "Female", NA))),
  race  = factor(case_when(
    RIDRETH1 == 3 ~ "Non-Hispanic White",
    RIDRETH1 == 4 ~ "Non-Hispanic Black",
    RIDRETH1 %in% c(1, 2) ~ "Hispanic",
    TRUE ~ "Other"
  )),
  education = factor(case_when(
    DMDEDUC2 %in% c(1,2) ~ "< High school",
    DMDEDUC2 == 3 ~ "High school graduate",
    DMDEDUC2 %in% c(4,5) ~ "> High school",
    TRUE ~ "Other"
  )),
  smoker = ifelse(SMQ020 == 2, 0,
           ifelse(SMQ040 %in% c(1, 2), 1,
           ifelse(SMQ040 == 3, 0, NA))),
  # 终生是否饮酒(ever drinker)二分：跨 5 周期一致构造
  #   D/E/F/I: ALQ101==1(过去12月饮) -> 1; ALQ101==2 且 ALQ110==1(终生曾饮) -> 1; ALQ101==2 且 ALQ110==2 -> 0
  #   J      : ALQ111==1(终生曾饮) -> 1; ALQ111==2 -> 0
  #   拒绝/不知(7,9)或缺失 -> NA
  drinker = case_when(
    ALQ101 == 1 ~ 1,
    ALQ101 == 2 & ALQ110 == 1 ~ 1,
    ALQ101 == 2 & ALQ110 == 2 ~ 0,
    ALQ111 == 1 ~ 1,
    ALQ111 == 2 ~ 0,
    TRUE ~ NA_real_
  ),
  bmi_cat = case_when(
    BMXBMI < 25 ~ "Normal",
    BMXBMI >= 25 & BMXBMI < 30 ~ "Overweight",
    BMXBMI >= 30 ~ "Obese",
    TRUE ~ NA_character_
  )
)

# survey 权重：5 周期合并 = WTMEC2YR / 5（NHANES 多周期合并标准做法）
df$weight <- ifelse(is.na(df$WTMEC2YR), NA_real_, df$WTMEC2YR / 5)

# ---- 验证 drinker 修复（回应 R1-8 数值矛盾排查）----
log("drinker 总体 NA 比例 = %.3f (修复目标 <0.10)", mean(is.na(df$drinker), na.rm = FALSE))
if ("cycle" %in% names(df)) {
  tab <- table(df$cycle, is.na(df$drinker), useNA = "ifany")
  log("drinker NA 按周期:\n%s", paste(capture.output(print(tab)), collapse = "\n"))
  # 数值分布(0/1)按周期，验证 J 周期 ALQ111 编码方向未被反转
  dval <- table(df$cycle, df$drinker, useNA = "ifany")
  log("drinker 取值分布(0=未饮,1=曾饮)按周期:\n%s", paste(capture.output(print(dval)), collapse = "\n"))
}


# 完整分析样本：valid DII + hs-CRP + PHQ-9（9 项全非缺失）+ 协变量
df_complete <- df %>% filter(
  !is.na(DII),
  !is.na(crp_mgL),
  !is.na(phq9_total),
  !is.na(weight), !is.na(SDMVPSU), !is.na(SDMVSTRA)
)
flow <- rbind(flow, data.frame(step = "完整分析样本(DII+CRP+PHQ9+权重)", n = nrow(df_complete)))

# -----------------------------------------------------------------------------
# 6. 保存
# -----------------------------------------------------------------------------
saveRDS(df, file.path(OUT, "nhanes_merged_full.rds"))
saveRDS(df_complete, file.path(OUT, "nhanes_merged_analysis.rds"))
write.csv(flow, file.path(OUT, "sample_flow.csv"), row.names = FALSE)

log("========== 完成 ==========")
log("完整合并样本(含缺失): %d", nrow(df))
log("分析样本: %d", nrow(df_complete))
print(flow)
cat("\n关键文件:\n  ", file.path(OUT, "nhanes_merged_analysis.rds"), "\n")
cat("  ", file.path(OUT, "sample_flow.csv"), "\n")
