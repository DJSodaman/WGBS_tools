"""
参考基因组处理模块
- Bismark Bowtie2 索引构建
- 基因组 C 统计计算与缓存
"""

import os
import pickle
from ..utils.shell import run_cmd
from ..utils.genome_utils import compute_genome_c_stats, get_genome_totals


def build_bismark_index(ref_dir: str, dry_run: bool = False) -> bool:
    """
    构建 Bismark + Bowtie2 基因组索引。
    在 ref_dir 中检测 Bisulfite_Genome/ 目录。
    若已存在则跳过，否则运行 bismark_genome_preparation。

    返回 True 表示索引已就绪。
    """
    bisulfite_dir = os.path.join(ref_dir, "Bisulfite_Genome")
    ct_dir = os.path.join(bisulfite_dir, "CT_conversion")
    ga_dir = os.path.join(bisulfite_dir, "GA_conversion")

    if os.path.isdir(ct_dir) and os.path.isdir(ga_dir):
        print(f"[genome] Bismark 索引已存在: {bisulfite_dir}/")
        return True

    print(f"[genome] 构建 Bismark Bowtie2 索引 (首次)...")
    print(f"[genome] 参考基因组目录: {ref_dir}")

    cmd = (
        f"bismark_genome_preparation --bowtie2 --verbose "
        f"'{ref_dir}'"
    )

    run_cmd(cmd, log_prefix="genome:index", dry_run=dry_run)

    if os.path.isdir(ct_dir) and os.path.isdir(ga_dir):
        print(f"[genome] Bismark 索引构建完成: {bisulfite_dir}/")
        return True
    else:
        print(f"[genome] 警告: 未检测到 Bisulfite_Genome/CT_conversion，索引可能构建失败")

    return False


def compute_and_cache_c_stats(ref_dir: str, genome_fna: str,
                              cache_dir: str = "") -> dict:
    """
    计算基因组 C 统计。结果缓存在 pickle 中。

    返回:
        {"per_chrom": {...}, "totals": {"total_C": N, "C_CG": N, ...}}
    """
    fasta_path = os.path.join(ref_dir, genome_fna)
    cache_path = os.path.join(
        cache_dir or ref_dir,
        f"{genome_fna}.c_stats.pkl"
    )

    # 检查缓存
    if os.path.exists(cache_path):
        print(f"[genome] 从缓存加载 C 统计: {cache_path}")
        with open(cache_path, "rb") as f:
            return pickle.load(f)

    print(f"[genome] 计算基因组 C 统计 (双链)...")
    per_chrom = compute_genome_c_stats(fasta_path)
    totals = get_genome_totals(per_chrom)

    result = {"per_chrom": per_chrom, "totals": totals}
    print(f"[genome] 基因组 C 统计: "
          f"CG={totals['C_CG']:,} "
          f"CHG={totals['C_CHG']:,} "
          f"CHH={totals['C_CHH']:,} "
          f"Total={totals['total_C']:,}")

    # 写缓存
    with open(cache_path, "wb") as f:
        pickle.dump(result, f)
    print(f"[genome] C 统计已缓存: {cache_path}")

    return result
