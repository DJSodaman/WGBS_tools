"""
WGBS_tools 集中配置
====================
所有可配置的路径、参数、并行策略在此定义。
"""

import os
from dataclasses import dataclass, field
from typing import Optional, List


# ============================================================
# 分析参数 (匹配公司 WGBS pipeline)
# ============================================================

@dataclass
class AnalysisConfig:
    """分析配置，可通过 CLI 或样本 manifest 覆盖。"""

    # --- Bismark/Bowtie2 比对参数 ---
    bowtie2_score_min: str = "L,0,-0.2"
    max_insert: int = 700
    dovetail: bool = True

    # --- 甲基化提取参数 ---
    ignore_r1: int = 5
    ignore_r2: int = 5  # 仅 PE 使用

    # --- 统计参数 ---
    bisulfite_non_conversion_rate: float = 0.01  # H0 非转化率
    bh_correction_grouping: str = "by_context"   # "by_context" | "global"

    # --- fastp 修剪参数 (公司参数) ---
    fastp_params: dict = field(default_factory=lambda: {
        "cut_front_window_size": 1,
        "cut_front_mean_quality": 3,
        "cut_tail_window_size": 1,
        "cut_tail_mean_quality": 3,
        "cut_right_window_size": 4,
        "cut_right_mean_quality": 15,
        "length_required": 36,
    })

    # --- 并行策略 ---
    bismark_parallel_instances: int = 8
    bismark_bowtie2_threads: int = 13
    methyl_extractor_multicore: int = 16
    max_parallel_samples: int = 3
    total_cores: int = 112


# ============================================================
# 默认路径 (运行时由 CLI --ref_dir / --output_dir 覆盖)
# ============================================================

DEFAULT_REF_DIR = "/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/ref"
DEFAULT_OUTPUT_DIR = "/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/WGBS_tools_new_result"

# 参考基因组文件名 (在 ref_dir 内)
GENOME_FNA = "GCF_000182925.2_NC12_genomic.fna"
GENOME_GFF = "GCF_000182925.2_NC12_genomic.gff"
