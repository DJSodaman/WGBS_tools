"""
Read 修剪模块 — fastp wrapper
支持单端和双端数据。
"""

import os
from ..utils.shell import run_cmd


def run_fastp(
    read1: str,
    output_dir: str,
    sample_name: str,
    read2: str = "",
    threads: int = 8,
    fastp_params: dict = None,
    dry_run: bool = False,
) -> tuple:
    """
    运行 fastp 修剪。

    参数:
        read1: 输入 read1 FASTQ 路径
        read2: 输入 read2 FASTQ 路径 (空字符串 = 单端)
        output_dir: 输出目录
        sample_name: 样本名
        threads: 线程数
        fastp_params: fastp 参数字典

    返回:
        (trimmed_read1_path, trimmed_read2_path_or_empty)
    """
    if fastp_params is None:
        fastp_params = {}

    os.makedirs(output_dir, exist_ok=True)

    paired_end = bool(read2)
    out1 = os.path.join(output_dir, f"{sample_name}_trimmed_1.fq.gz")
    out2 = os.path.join(output_dir, f"{sample_name}_trimmed_2.fq.gz") if paired_end else ""

    html_report = os.path.join(output_dir, f"{sample_name}_fastp.html")
    json_report = os.path.join(output_dir, f"{sample_name}_fastp.json")

    # 构建 fastp 命令
    cmd_parts = [
        "fastp",
        f"--in1 '{read1}'",
    ]
    if paired_end:
        cmd_parts.append(f"--in2 '{read2}'")
        cmd_parts.append(f"--out1 '{out1}'")
        cmd_parts.append(f"--out2 '{out2}'")
    else:
        cmd_parts.append(f"--out1 '{out1}'")

    # 修剪参数
    for key, val in fastp_params.items():
        if val is True:
            cmd_parts.append(f"--{key}")
        else:
            cmd_parts.append(f"--{key}={val}")

    cmd_parts.extend([
        f"--thread {threads}",
        f"--html '{html_report}'",
        f"--json '{json_report}'",
    ])

    cmd = " ".join(cmd_parts)
    run_cmd(cmd, log_prefix=f"trim:{sample_name}", dry_run=dry_run)

    return (out1, out2)
