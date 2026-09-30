#!/bin/bash
# ============================================================
# 00_env_setup.sh — AutoDL 实例环境一键安装
# 项目：PSY-D-26-02386 大修 de novo 重跑（MR + NHANES 中介）
# 适配：Ubuntu + Miniconda 镜像（RTX 4080 SUPER / 12 vCPU / 62GB 实测可用）
#
# 用法：建议在 screen 里跑（SSH 断线不影响）：
#   screen -S setup
#   bash 00_env_setup.sh
#   （Ctrl+A 再按 D 退出 screen；重连用 screen -r setup）
# 脚本可重复执行：已装过的包会自动跳过，中断后直接重跑即可。
# 预计耗时：20–40 分钟（取决于源速度）
#
# ⚠ 注意：必须在「正常开机」状态执行！
#   无卡模式只有 0.5 核 / 2GB 内存，装 R 包会 OOM 失败。
# ============================================================
set -euo pipefail

# 开启 AutoDL 学术加速（对 conda/pip/wget 的外网下载都有效；重复执行无副作用）
source /etc/network_turbo 2>/dev/null || echo "(未找到学术加速脚本，继续用默认网络)"

echo "===== [1/6] 系统包 ====="
# 换 AutoDL 内部 apt 源通常已配好，直接装
apt-get update -qq || sudo apt-get update -qq
apt-get install -y -qq tabix wget curl git bzip2 libcurl4-openssl-dev \
  libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev \
  libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  libgmp-dev libmpfr-dev || \
  sudo apt-get install -y -qq tabix wget curl git bzip2 libcurl4-openssl-dev \
  libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev \
  libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  libgmp-dev libmpfr-dev
echo "tabix 版本: $(tabix --version 2>&1 | head -1)"

echo "===== [2/6] 确认 conda 可用 ====="
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source /root/miniconda3/etc/profile.d/conda.sh
conda --version
# 关闭 libmamba solver 的报错（可消除 "conda-libmamba-solver ... QueryFormat" 提示，纯装饰性）
conda config --set solver classic 2>/dev/null || true

echo "===== [3/6] 安装 R 4.3+（conda-forge）====="
# 若已有 r-base 可跳过
if ! command -v R >/dev/null 2>&1; then
  conda install -y -c conda-forge "r-base>=4.3" r-essentials
fi
R --version | head -1

echo "===== [4/6] 安装 R 包（多镜像自动探测，可重复执行）====="
# 独立脚本会先清代理、再逐个镜像探测，解决"清华源连不上 / not available"问题
if [ -f setup_R_packages.R ]; then
  Rscript setup_R_packages.R
else
  echo "未找到 setup_R_packages.R，请与 00_env_setup.sh 一起上传"
  exit 1
fi

echo "===== [5/6] 写入 OpenGWAS token 到 ~/.Renviron ====="
# Token 有效期至 2026-10-11 左右（14 天有效，过期需到 api.opengwas.io 重新生成）
if ! grep -q "OPENGWAS_JWT" ~/.Renviron 2>/dev/null; then
  echo "请把 jwt token 粘贴到下面（一行，以 eyJ 开头）：" >&2
  read -r JWT
  echo "OPENGWAS_JWT=${JWT}" >> ~/.Renviron
fi
echo "Renviron 已配置: $(grep -c OPENGWAS_JWT ~/.Renviron) 条"

echo "===== [6/6] Python 包（pip 走 AutoDL 学术加速）====="
source /etc/network_turbo 2>/dev/null || true
pip install -q --upgrade pip
pip install -q pandas numpy scipy matplotlib statsmodels lifelines openpyxl
python -c "import pandas, numpy, scipy, matplotlib; print('python 包 OK:', pandas.__version__, numpy.__version__)"

echo ""
echo "=============================================="
echo " 环境安装完成。下一步：bash 01_smoke_test.sh"
echo "=============================================="
