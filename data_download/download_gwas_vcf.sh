#!/bin/bash
# 08_download_vcf.sh — 直接下载 5 个 GWAS 的完整 VCF（绕过 API 的 tophits/associations）
#
# 背景：OpenGWAS 的 Pipeline 服务（tophits/associations）处于计划内维护，已挂一整夜。
#       但完整 GWAS summary statistics 以 VCF 文件托管在文件服务器，可 curl 直下。
#       探测（07_probe_vcf.sh）确认 gwas.mrcieu.ac.uk/files/{id}/{id}.vcf.gz 返回 301
#       重定向（有效地址），需 curl -L 跟随。
#
# 用法：bash 08_download_vcf.sh
# 产物：/root/autodl-tmp/ref/vcf/{id}.vcf.gz + .tbi

set -u

BASE="https://gwas.mrcieu.ac.uk/files"
DEST="/root/autodl-tmp/ref/vcf"
mkdir -p "$DEST"

IDS=("ieu-b-35" "ieu-b-102" "ieu-a-1187" "ukb-b-19085" "ieu-b-4812")

echo "===== 开始下载 5 个 GWAS 的 VCF 文件 ====="
echo "目标目录: $DEST"
echo ""

for id in "${IDS[@]}"; do
  for ext in "vcf.gz" "vcf.gz.tbi"; do
    url="$BASE/$id/$id.$ext"
    out="$DEST/$id.$ext"

    if [ -f "$out" ] && [ -s "$out" ]; then
      echo "[skip] 已存在且非空: $id.$ext"
      continue
    fi

    echo "[down] $id.$ext"
    # -L 跟随 301 重定向；-C - 断点续传；--retry 自动重试；-o 落盘
    curl -L -C - --retry 3 --retry-delay 5 -o "$out" "$url" 2>/dev/null
    rc=$?
    if [ $rc -eq 0 ] && [ -s "$out" ]; then
      sz=$(du -h "$out" | cut -f1)
      echo "       OK  size=$sz"
    else
      echo "       FAIL (rc=$rc) — 删除残留，稍后单独重试"
      rm -f "$out"
    fi
  done
done

echo ""
echo "===== 下载完成，文件清单 ====="
ls -lh "$DEST"

echo ""
echo "===== 校验：各 VCF 能否正常解压读取（表头 + 前 2 行数据）====="
for id in "${IDS[@]}"; do
  f="$DEST/$id.vcf.gz"
  if [ -s "$f" ]; then
    echo ""
    echo "----- $id -----"
    zcat "$f" 2>/dev/null | grep -v '^##' | head -3
  else
    echo "----- $id : [文件缺失！] -----"
  fi
done

echo ""
echo "===== 结束 ====="
