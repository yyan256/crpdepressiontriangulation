#!/usr/bin/env Rscript
# ============================================================
# setup_R_packages.R — 稳健安装 R 包（多镜像自动探测 + 可重复执行）
# 解决两个问题：
#   1) AutoDL「学术加速」的代理会拦断国内 CRAN 镜像 → 安装前先清代理
#   2) 单一清华源偶尔不可达 → 逐个探测镜像，用第一个能下的
# 用法：Rscript setup_R_packages.R
# ============================================================

# ---- 1) 清理代理环境变量（关键：学术加速代理不转发 CRAN 国内镜像）----
for (v in c("http_proxy","https_proxy","HTTP_PROXY","HTTPS_PROXY",
            "all_proxy","ALL_PROXY","no_proxy","NO_PROXY")) {
  Sys.unsetenv(v)
}
cat("[proxy] 已清理代理环境变量\n")

options(Ncpus = parallel::detectCores(), timeout = 300)

# ---- 2) 候选镜像（国内优先，官方兜底）----
mirrors <- c(
  "https://mirrors.tuna.tsinghua.edu.cn/CRAN",
  "https://mirrors.ustc.edu.cn/CRAN",
  "https://mirrors.aliyun.com/CRAN",
  "https://mirrors.cloud.tencent.com/CRAN",
  "https://cloud.r-project.org"
)

# ---- 3) 探测：找一个能成功拉取 PACKAGES 索引的镜像 ----
chosen <- NULL
for (m in mirrors) {
  cat("[probe] 探测", m, "... ")
  ok <- tryCatch({
    ap <- available.packages(repos = m, type = "source", quiet = TRUE)
    is.matrix(ap) && nrow(ap) > 10000
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (isTRUE(ok)) { chosen <- m; cat("OK\n"); break } else cat("FAIL\n")
}
if (is.null(chosen)) stop("所有 CRAN 镜像均不可达，请检查网络后重跑本脚本")
options(repos = c(CRAN = chosen))
cat("== 选用镜像:", chosen, "==\n")

# ---- 4) 待装包 ----
need <- c(
  "TwoSampleMR",            # 主 MR 框架
  "MendelianRandomization", # IVW / Egger / PRESSO 等估计量
  "ieugwasr",               # OpenGWAS API 客户端（新版 api.opengwas.io）
  "survey",                 # NHANES 加权分析
  "lavaan",                 # SEM 中介
  "nhanesA",                # NHANES 数据直接下载
  "data.table", "dplyr", "tidyr", "ggplot2", "Matrix"
)
missing <- need[!need %in% rownames(installed.packages())]
if (length(missing) > 0) {
  cat("\n[install] 开始安装", length(missing), "个包:\n")
  print(missing)
  install.packages(missing)
} else {
  cat("\n[install] 所有目标包已安装，跳过\n")
}

# ---- 5) 打印最终结果（只打印已装上的）----
cat("\n--- 已安装版本 ---\n")
inst <- need[need %in% rownames(installed.packages())]
for (p in inst) cat(sprintf("%-24s %s\n", p, as.character(packageVersion(p))))
miss2 <- setdiff(need, inst)
if (length(miss2)) cat("\n[WARN] 仍未装上（需排查）:", paste(miss2, collapse = ", "), "\n")
