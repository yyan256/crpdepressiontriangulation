#!/usr/bin/env bash
# ============================================================
# 00_download_ref_panel.sh
# 下载真实 LD clumping 所需的两样东西：
#   1. plink2 二进制（本地 clumping 工具）
#   2. 1000G EUR 参考面板（1kg.v3，bed/bim/fam，用于计算 LD r²）
#
# 用途：修复 B1 缺陷 —— 原 PR_v4.py 的 ld_clump 是"按 bp 窗口去重"的假 clump，
#       本脚本下载真实参考面板后，10_mr_main.R 会用 plink2 做真 LD clumping
#       （r²<0.001，window 10Mb），满足审稿人 R1-9 对"clump 前后工具数"的报告要求。
#
# 用法：bash 00_download_ref_panel.sh
# 产出：
#   /root/autodl-tmp/ref/plink2          (可执行文件)
#   /root/autodl-tmp/ref/1kg.v3/EUR.bed  (参考面板)
#   /root/autodl-tmp/ref/1kg.v3/EUR.bim
#   /root/autodl-tmp/ref/1kg.v3/EUR.fam
# ============================================================

set -uo pipefail   # 注意：不用 -e，下载失败要给 fallback 机会

REFDIR="/root/autodl-tmp/ref"
mkdir -p "$REFDIR"

echo "===== [1/3] 下载 plink2 二进制 ====="
PLINK2="$REFDIR/plink2"
if [ -x "$PLINK2" ]; then
  echo "  plink2 已存在，跳过"
else
  # 官方最新版（AWS 上的稳定资源）
  cd "$REFDIR" || exit 1
  wget -q --show-progress -O plink2.zip \
    "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip" \
    || wget -q --show-progress -O plink2.zip \
    "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250428.zip"
  if [ -f plink2.zip ]; then
    unzip -o plink2.zip plink2 2>/dev/null || unzip -o plink2.zip
    chmod +x plink2 2>/dev/null
    rm -f plink2.zip
  fi
fi

if [ -x "$PLINK2" ]; then
  echo "  [PASS] plink2 就绪：$("$PLINK2" --version 2>&1 | head -1)"
else
  echo "  [WARN] plink2 未就绪，稍后可用 conda 安装：conda install -y -c bioconda plink2"
fi

echo ""
echo "===== [2/3] 下载 1000G EUR 参考面板（1kg.v3，约 1.5GB）====="
LDDIR="$REFDIR/1kg.v3"
mkdir -p "$LDDIR"

if [ -f "$LDDIR/EUR.bed" ] && [ -f "$LDDIR/EUR.bim" ] && [ -f "$LDDIR/EUR.fam" ]; then
  echo "  EUR 参考面板已存在，跳过"
else
  cd "$REFDIR" || exit 1

  # 关键：先开学术加速，否则 fileserve/S3 会卡在几十 KB/s
  if [ -f /etc/network_turbo ]; then
    source /etc/network_turbo 2>/dev/null && echo "  [加速] 已开启学术加速"
  else
    echo "  [提示] 未检测到学术加速脚本，若速度 <500KB/s 请手动 source /etc/network_turbo"
  fi

  # 优先用 aria2c 多线程（16 连接），没有就退回 wget 断点续传
  if command -v aria2c >/dev/null 2>&1; then
    echo "  使用 aria2c 多线程下载（16 连接，断点续传）..."
    aria2c -c -x 16 -s 16 -k 1M \
      "https://mrcieu.s3.amazonaws.com/ld/1kg.v3.tgz" \
      "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz" \
      -d "$REFDIR" -o 1kg.v3.tgz
  else
    # wget -c 断点续传：中断后重跑不会从头下
    echo "  使用 wget 断点续传（优先 https S3 源）..."
    wget -c -q --show-progress -T 60 -t 3 -O 1kg.v3.tgz \
      "https://mrcieu.s3.amazonaws.com/ld/1kg.v3.tgz" \
      || wget -c -q --show-progress -T 60 -t 3 -O 1kg.v3.tgz \
      "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
  fi

  # 校验完整性：1.5GB 的 tgz 若没下完会 tar 报错，据此判断
  if [ -f 1kg.v3.tgz ]; then
    sz=$(stat -c%s 1kg.v3.tgz 2>/dev/null || stat -f%z 1kg.v3.tgz 2>/dev/null || echo 0)
    echo "  已下载 1kg.v3.tgz：$((sz/1024/1024)) MB"
    if [ "$sz" -gt 1400000000 ]; then
      echo "  解压中（约 1-3 分钟）..."
      tar -xzf 1kg.v3.tgz && rm -f 1kg.v3.tgz
      # 自动归位：MRCIEU 的 1kg.v3.tgz 解压后是"平铺"的（EUR.bed 直接落在 ref/ 下），
      # 若 EUR.bed 不在 1kg.v3/ 子目录，则统一归位过去（10_mr_main.R 引用的是 1kg.v3/EUR）
      if [ ! -f "$LDDIR/EUR.bed" ]; then
        EUR_BED=$(find "$REFDIR" -maxdepth 3 -name "EUR.bed" 2>/dev/null | head -1)
        if [ -n "$EUR_BED" ]; then
          echo "  检测到 EUR.bed 在 $(dirname "$EUR_BED")，归位到 $LDDIR/ ..."
          mkdir -p "$LDDIR"
          for pop in EUR EAS AFR AMR SAS; do
            for ext in bed bim fam; do
              [ -f "$(dirname "$EUR_BED")/$pop.$ext" ] && mv "$(dirname "$EUR_BED")/$pop.$ext" "$LDDIR/"
            done
          done
        fi
      fi
    else
      echo "  [WARN] 文件不完整（<1.4GB），未解压；重跑本脚本会断点续传"
    fi
  fi
fi

if [ -f "$LDDIR/EUR.bed" ] && [ -f "$LDDIR/EUR.bim" ] && [ -f "$LDDIR/EUR.fam" ]; then
  echo "  [PASS] EUR 参考面板就绪："
  ls -lh "$LDDIR/EUR.bed" "$LDDIR/EUR.bim" "$LDDIR/EUR.fam"
else
  echo "  [WARN] EUR 参考面板未就绪。备选方案："
  echo "    1) 用学术加速后重跑本脚本（source /etc/network_turbo）"
  echo "    2) 手动下载解压到 $LDDIR/（需 EUR.bed/bim/fam 三件套）"
fi

echo ""
echo "===== [3/3] 校验汇总 ====="
if [ -x "$PLINK2" ] && [ -f "$LDDIR/EUR.bed" ]; then
  echo "  [PASS] clumping 环境全部就绪，可运行 10_mr_main.R"
else
  echo "  [FAIL] 尚有缺失，10_mr_main.R 在 clump 阶段会清晰报错并提示补装"
fi
