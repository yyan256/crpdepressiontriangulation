#!/usr/bin/env Rscript
# ============================================================================
# 09_download_vcf.R — 用签名 URL 下载 5 个 GWAS 的完整 VCF（绕过 Pipeline 故障）
#
# 背景：直接匿名访问 opengwas.io/files/... 返回 502（文件服务与 Pipeline 同故障），
#       但 ieugwasr 的 gwasinfo_files() 走 Metadata 服务（gwasinfo 已确认可用），
#       返回带签名的 S3 下载 URL（2 小时有效）。用签名 URL + curl 即可下到真数据。
#
# 用法：Rscript 09_download_vcf.R
# 产物：/root/autodl-tmp/ref/vcf/{id}.vcf.gz + .tbi（每个几十 MB~数 GB）
# ============================================================================

suppressPackageStartupMessages({ library(ieugwasr) })

ids  <- c("ieu-b-35", "ieu-b-102", "ieu-a-1187", "ukb-b-19085", "ieu-b-4812")
dest <- "/root/autodl-tmp/ref/vcf"
dir.create(dest, showWarnings = FALSE, recursive = TRUE)

MIN_SIZE <- 1e6   # 1MB，小于此视为错误页，重下

for (id in ids) {
  cat("\n========== ", id, " ==========\n", sep = "")

  # 1) 拿签名下载 URL
  f <- tryCatch(gwasinfo_files(id),
                error = function(e) { cat("  [ERROR] gwasinfo_files:", conditionMessage(e), "\n"); NULL })
  if (is.null(f) || (is.data.frame(f) && nrow(f) == 0)) {
    cat("  [失败] gwasinfo_files 无返回（可能该服务也受影响）\n")
    next
  }
  cat("  返回列名:", paste(names(f), collapse = ", "), "\n")
  print(f)

  # 2) 通用提取所有 http(s) URL
  urls <- character()
  for (col in names(f)) {
    vals <- as.character(f[[col]])
    urls <- c(urls, vals[grepl("^https?://", vals)])
  }
  urls <- unique(urls)
  if (length(urls) == 0) { cat("  [无 URL]\n"); next }

  # 3) 逐个下载（只下 .vcf.gz 与 .vcf.gz.tbi）
  for (u in urls) {
    fn <- basename(sub("\\?.*$", "", u))
    if (!grepl("\\.(vcf\\.gz|vcf\\.gz\\.tbi)$", fn)) {
      cat("  [跳过非VCF]", fn, "\n"); next
    }
    out <- file.path(dest, fn)
    if (file.exists(out) && file.info(out)$size >= MIN_SIZE) {
      cat("  [已存在] ", fn, " (", round(file.info(out)$size/1e6, 1), "MB)\n", sep = ""); next
    }
    cat("  [下载] ", fn, "\n", sep = "")
    # 优先 aria2c 多线程（Oracle 对象存储支持 range 分段，速度可达 curl 的 5-10 倍），
    # 否则退回 curl 断点续传。签名 URL 2h 过期，必须尽量快下完大文件。
    if (nzchar(Sys.which("aria2c"))) {
      rc <- system(sprintf(
        "aria2c -x 16 -s 16 -k 1M --file-allocation=none --summary-interval=10 -d '%s' -o '%s' '%s'",
        dest, fn, u))
    } else {
      rc <- system(sprintf("curl -L -C - --retry 3 --retry-delay 5 -o '%s' '%s'", out, u))
    }
    if (rc == 0 && file.exists(out) && file.info(out)$size >= MIN_SIZE) {
      cat("    OK ", round(file.info(out)$size/1e6, 1), "MB\n", sep = "")
    } else {
      sz <- if (file.exists(out)) file.info(out)$size else 0
      cat("    FAIL rc=", rc, " size=", sz, "（若 size 很小=又是错误页）\n", sep = "")
      if (file.exists(out) && sz < MIN_SIZE) file.remove(out)
    }
  }
}

cat("\n===== 最终清单 =====\n")
if (dir.exists(dest)) {
  fl <- list.files(dest, full.names = TRUE)
  for (p in fl) cat(sprintf("  %8.1f MB  %s\n", file.info(p)$size/1e6, basename(p)))
}

cat("\n===== 校验表头（前 2 行数据）=====\n")
for (id in ids) {
  vf <- file.path(dest, paste0(id, ".vcf.gz"))
  if (file.exists(vf) && file.info(vf)$size >= MIN_SIZE) {
    cat("\n----- ", id, " -----\n", sep = "")
    system(sprintf("zcat '%s' 2>/dev/null | grep -v '^##' | head -3", vf))
  } else {
    cat("\n----- ", id, " : [缺失或过小] -----\n", sep = "")
  }
}
cat("\n===== 结束 =====\n")
