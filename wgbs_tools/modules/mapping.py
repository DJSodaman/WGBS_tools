"""
Bismark 比对模块 — 单端和双端数据
"""

import os
from ..utils.shell import run_cmd


def run_bismark_alignment(
    read1: str,
    ref_dir: str,
    output_dir: str,
    sample_name: str,
    read2: str = "",
    parallel_instances: int = 8,
    bowtie2_threads: int = 13,
    score_min: str = "L,0,-0.2",
    max_insert: int = 700,
    dovetail: bool = True,
    dry_run: bool = False,
) -> str:
    """
    运行 Bismark 比对。

    参数:
        read1: 输入 read1 FASTQ (或单端唯一输入)
        read2: 输入 read2 FASTQ (空 = 单端)
        ref_dir: 参考基因组目录 (含 Bisulfite_Genome/)
        output_dir: BAM 输出目录
        sample_name: 样本名
        parallel_instances: Bismark --parallel 参数
        bowtie2_threads: bowtie2 -p 参数
        score_min: bowtie2 --score_min 参数
        max_insert: 最大插入长度 (仅 PE)
        dovetail: 是否允许 dovetail 比对 (仅 PE)

    返回:
        生成的 BAM 文件路径
    """
    os.makedirs(output_dir, exist_ok=True)

    paired_end = bool(read2)
    mode_str = "PE" if paired_end else "SE"

    # 注意: --basename 与 --parallel 不兼容 (Bismark 0.24.x)
    # 不指定 --basename, Bismark 根据输入文件名自动命名
    cmd_parts = [
        "bismark",
        "--bowtie2",
        f"--genome_folder '{ref_dir}'",
        f"--score_min {score_min}",
        f"--parallel {parallel_instances}",
        f"-p {bowtie2_threads}",
        "--bam",
        "--quiet",
        f"-o '{output_dir}'",
    ]

    if paired_end:
        cmd_parts.append(f"-1 '{read1}' -2 '{read2}'")
        cmd_parts.append(f"-X {max_insert}")
        if dovetail:
            cmd_parts.append("--dovetail")
    else:
        cmd_parts.append(f"'{read1}'")

    cmd = " ".join(cmd_parts)
    run_cmd(cmd, log_prefix=f"bismark:{sample_name}({mode_str})", dry_run=dry_run)

    # Bismark 输出 BAM 查找 (无 --basename 时自动命名)
    # 期望命名: {read1_basename}_bismark_bt2.bam
    read1_base = os.path.splitext(os.path.basename(read1))[0]
    # 处理 .fq.gz -> .fq -> 基名
    for ext in [".fq.gz", ".fastq.gz", ".fq", ".fastq"]:
        if read1_base.endswith(ext):
            read1_base = read1_base[:-len(ext)]
            break
    expected_bam = os.path.join(output_dir, f"{read1_base}_bismark_bt2.bam")

    if os.path.exists(expected_bam):
        return expected_bam

    # 备用: 搜索包含样本名的 BAM
    for f in sorted(os.listdir(output_dir), reverse=True):
        if f.endswith("_bismark_bt2.bam"):
            path = os.path.join(output_dir, f)
            return path

    raise FileNotFoundError(
        f"Bismark BAM 输出未找到: 期望 {expected_bam}\n"
        f"  输出目录内容: {os.listdir(output_dir)}"
    )
