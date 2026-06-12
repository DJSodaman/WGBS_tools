#!/bin/bash
# WGBS_tools 一键环境安装脚本
# 使用 mamba 构建环境并安装本软件包
set -e

MAMBA=$(which mamba 2>/dev/null || echo "")
if [ -z "$MAMBA" ]; then
    echo "[错误] 未找到 mamba，请先安装 mamba 或 conda"
    echo "  安装方法: conda install -c conda-forge mamba -y"
    exit 1
fi

echo "============================================"
echo "  WGBS_tools 环境安装"
echo "  mamba: $MAMBA"
echo "============================================"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_YAML="${SCRIPT_DIR}/env.yaml"

echo "[1/2] 用 mamba 创建环境 WGBS_tools ..."
mamba env create -f "$ENV_YAML" -y 2>&1 | tail -20

echo ""
echo "[2/2] 安装 WGBS_tools 本身 (pip install -e) ..."
# 获取 mamba/conda 安装路径
CONDA_PREFIX=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
ENV_PATH="${CONDA_PREFIX}/envs/WGBS_tools"
"${ENV_PATH}/bin/pip" install -e "$SCRIPT_DIR" 2>&1 | tail -5

echo ""
echo "============================================"
echo "  WGBS_tools 环境安装完成！"
echo ""
echo "  激活环境: conda activate WGBS_tools"
echo "  运行分析: python -m wgbs_tools --help"
echo "============================================"
