"""
去重模块 — deduplicate_bismark wrapper
支持单端和双端数据。
"""

import os
from ..utils.shell import run_cmd


def deduplicate(
    bam_file: str,
    output_dir: str,
    sample_name: str,
    paired_end: bool = True,
    dry_run: bool = False,
) -> str:
    """
    运行 deduplicate_bismark。

    参数:
        bam_file: 输入 BAM 文件
        output_dir: 输出目录
        sample_name: 样本名
        paired_end: 是否为配对端

    返回:
        去重后的 BAM 文件路径
    """
    os.makedirs(output_dir, exist_ok=True)

    mode_flag = "--paired" if paired_end else "--single"
    mode_str = "PE" if paired_end else "SE"

    cmd = (
        f"deduplicate_bismark "
        f"{mode_flag} "
        f"--output_dir '{output_dir}' "
        f"--bam '{bam_file}'"
    )

    run_cmd(cmd, log_prefix=f"dedup:{sample_name}({mode_str})", dry_run=dry_run)

    # deduplicate_bismark 输出: {bam_base}.deduplicated.bam
    base = os.path.splitext(os.path.basename(bam_file))[0]
    # 去掉 _bismark_bt2 后缀 (由 deduplicate_bismark 自动处理)
    for suffix in ["_bismark_bt2", "_bismark"]:
        if base.endswith(suffix):
            base = base[:-len(suffix)]
            break

    expected = os.path.join(output_dir, f"{base}.deduplicated.bam")

    if os.path.exists(expected):
        return expected

    # 搜索可能的输出
    for f in os.listdir(output_dir):
        if f.endswith(".deduplicated.bam"):
            return os.path.join(output_dir, f)

    raise FileNotFoundError(f"去重 BAM 未找到: {expected}")
