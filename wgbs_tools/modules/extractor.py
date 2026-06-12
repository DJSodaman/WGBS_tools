"""
甲基化提取模块 — bismark_methylation_extractor wrapper
支持单端和双端数据。
"""

import os
from ..utils.shell import run_cmd


def run_methylation_extractor(
    dedup_bam: str,
    output_dir: str,
    sample_name: str,
    paired_end: bool = True,
    multicore: int = 16,
    ignore_r1: int = 5,
    ignore_r2: int = 5,
    dry_run: bool = False,
) -> str:
    """
    运行 bismark_methylation_extractor --comprehensive。

    参数:
        dedup_bam: 去重后的 BAM 文件
        output_dir: 输出目录
        sample_name: 样本名
        paired_end: 是否为配对端
        multicore: --multicore 参数
        ignore_r1: 忽略 read1 前 N bp
        ignore_r2: 忽略 read2 前 N bp (仅 PE)

    返回:
        甲基化提取输出目录
    """
    os.makedirs(output_dir, exist_ok=True)

    mode_str = "PE" if paired_end else "SE"

    cmd_parts = [
        "bismark_methylation_extractor",
        f"--multicore {multicore}",
        "--comprehensive",
        "--gzip",
        f"--ignore {ignore_r1}",
        f"--output_dir '{output_dir}'",
    ]

    if paired_end:
        cmd_parts.append("--paired-end")
        cmd_parts.append("--no_overlap")
        cmd_parts.append(f"--ignore_r2 {ignore_r2}")
    else:
        cmd_parts.append("--single-end")

    cmd_parts.append(f"'{dedup_bam}'")
    cmd = " ".join(cmd_parts)

    run_cmd(cmd, log_prefix=f"extract:{sample_name}({mode_str})", dry_run=dry_run)

    return output_dir
